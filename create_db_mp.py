#!/usr/bin/env python3
# author          : Pavel
# date            : 13.07.16
# version         : 0.1
# python_version  : 3
# copyright       : Pavel 2016
# license         : GPL3
#==============================================================================

import os
import sys
import argparse
import sqlite3 as lite
import marshal

from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from read_input import read_input
from multiprocessing import cpu_count, Pool
from pharmacophore import Pharmacophore, read_smarts_feature_file, load_multi_conf_mol


def create_tables(cursor, bin_step, smarts, store_coords, fp, nohash):
    cursor.execute("CREATE TABLE conformers(conf_id INTEGER PRIMARY KEY AUTOINCREMENT, mol_name TEXT, "
                   "mol_stereo_id TEXT)")
    if not nohash:
        cursor.execute("ALTER TABLE conformers ADD COLUMN ph_hash BLOB")
    if fp:
        cursor.execute("ALTER TABLE conformers ADD COLUMN fp BLOB")
    cursor.execute("CREATE INDEX mol_name_idx ON conformers(mol_name)")
    cursor.execute("CREATE TABLE settings(bin_step NUMERIC)")
    cursor.execute("INSERT INTO settings VALUES(?)", (bin_step, ))
    cursor.execute("CREATE TABLE smarts_features(label TEXT, value TEXT)")
    if smarts:
        for k, values in smarts.items():
            for v in values:
                cursor.execute("INSERT INTO smarts_features VALUES(?, ?)", (k, Chem.MolToSmarts(v)))
    if store_coords:
        cursor.execute("CREATE TABLE feature_coords(conf_id INTEGER, feature_label TEXT, "
                       "x NUMERIC, y NUMERIC, z NUMERIC, "
                       "FOREIGN KEY(conf_id) REFERENCES conformers(conf_id))")


def insert_res_db(cur, res, store_coords, stereo_id):
    for item in res:   # mol_name, hash, coords, fp
        if stereo_id:
            mol_name, mol_stereo_id = item[0].rsplit("_", 1)
        else:
            mol_name = item[0]
            mol_stereo_id = '0'
        data = [mol_name, mol_stereo_id]
        if item[1] is not None:
            data.append(item[1])
        if item[3] is not None:
            data.append(item[3])
        cur.execute("INSERT INTO conformers VALUES(NULL, %s)" % ','.join(['?'] * len(data)), data)
        if store_coords and item[2] is not None:
            cur.execute("SELECT MAX(conf_id) FROM conformers")
            conf_id = cur.fetchone()[0]
            cur.executemany("INSERT INTO feature_coords VALUES (?, ?, ?, ?, ?)", ((conf_id, i[0], *i[1]) for i in item[2]))


def compress_db(cursor, store_coords):
    # remove duplicated hashes for the same mol
    cursor.execute("DELETE FROM conformers WHERE rowid NOT IN (SELECT MIN(rowid) "
                   "FROM conformers GROUP BY mol_name, mol_stereo_id, ph_hash)")
    if store_coords:
        cursor.execute("DELETE FROM feature_coords WHERE conf_id NOT IN (SELECT conf_id "
                       "FROM conformers)")


def insert_res_txt(f, res, lines_set, stereo_id):
    # function is called if hash is not None, no need for additional check
    for item in res:   # mol_name, hash, coords, fp
        if stereo_id:
            mol_name, mol_stereo_id = item[0].rsplit("_", 1)
        else:
            mol_name = item[0]
            mol_stereo_id = '0'
        record = (mol_name, mol_stereo_id, item[1])
        if record not in lines_set:
            f.write("%s\t%s\t%s\n" % record)
            lines_set.add(record)


def prep_input(fname, id_field_name, smarts, bin_step, store_coords, multiconf, tolerance, fp, nohash):
    if fname is None:
        read_iterator = read_input(fname, 'sdf', id_field_name, removeHs=False)
    else:
        read_iterator = read_input(fname, id_field_name=id_field_name, removeHs=False)
    for mol, mol_name in read_iterator:
        yield mol, mol_name, smarts, bin_step, store_coords, multiconf, tolerance, fp, nohash


def map_process_mol(args):
    return process_mol(*args)


def process_mol(mol, mol_name, smarts, bin_step, store_coords, multiconf, tolerance, fp, nohash):
    # process_factory is within process scope only
    if multiconf:
        if 'process_factory' in globals():
            ps = load_multi_conf_mol(mol, factory=process_factory, bin_step=bin_step)
        else:
            ps = load_multi_conf_mol(mol, smarts_features=smarts, bin_step=bin_step)
        output = []
        for p in ps:
            coords = p.get_feature_coords() if store_coords else None
            hash = p.get_signature_md5(tol=tolerance) if not nohash else None
            fp_bin = marshal.dumps(p.get_fp_on_bits()) if fp else None
            output.append((mol_name, hash, coords, fp_bin))
        return output
    else:
        p = Pharmacophore(bin_step)
        if 'process_factory' in globals():
            p.load_from_feature_factory(mol, process_factory)
        elif smarts:
            p.load_from_smarts(mol, smarts)
        coords = p.get_feature_coords() if store_coords else None
        hash = p.get_signature_md5(tol=tolerance) if not nohash else None
        fp_bin = marshal.dumps(p.get_fp_on_bits()) if fp else None
        return [(mol_name, hash, coords, fp_bin)]


def pool_init(fdef_fname):
    global process_factory
    process_factory = ChemicalFeatures.BuildFeatureFactory(fdef_fname) if fdef_fname else None


def main_params(conformers_fname, out_fname, dbout_fname, bin_step, rewrite_db, store_coords, id_field_name, stereo_id,
                smarts_features_fname, rdkit_factory, fp, nohash, ncpu, tolerance, verbose):

    if out_fname is None and dbout_fname is None:
        raise ValueError("You should specify at least one of possible outputs: text file or SQLite DB.")

    # check DB existence
    if dbout_fname is not None and os.path.isfile(dbout_fname):
        if rewrite_db:
            os.remove(dbout_fname)
        else:
            raise FileExistsError("DB exists. To rewrite it add -r key to the command line call.")

    # load smarts features
    smarts = None
    if smarts_features_fname and rdkit_factory is None:
        smarts = read_smarts_feature_file(smarts_features_fname)

    multiconf = conformers_fname.lower().endswith('.pkl')

    # open db
    if dbout_fname is not None:
        conn = lite.connect(dbout_fname)
        cur = conn.cursor()
        create_tables(cur, bin_step, smarts, store_coords, fp, nohash)

    # open text file
    if out_fname is not None and not nohash:
        ftxt = open(out_fname, 'wt')
        ftxt.write("mol_id\tstereo_id\thash\n")

    nprocess = max(min(ncpu, cpu_count()), 1)

    if rdkit_factory:
        p = Pool(nprocess, initializer=pool_init, initargs=[rdkit_factory])
    else:
        p = Pool(nprocess)

    lines_set = set()

    counter = 0

    try:
        for i, res in enumerate(p.imap_unordered(map_process_mol,
                                                 prep_input(conformers_fname, id_field_name, smarts, bin_step, store_coords, multiconf, tolerance, fp, nohash),
                                                 chunksize=10), 1):
            counter += len(res)
            if dbout_fname is not None:
                insert_res_db(cur, res, store_coords, stereo_id)
                if i % 1000 == 0:
                    conn.commit()
            if out_fname is not None and not nohash:
                insert_res_txt(ftxt, res, lines_set, stereo_id)
            if verbose and i % 100 == 0:
                sys.stderr.write('\rprocessed %i molecules and %i conformers/pharmacophores' % (i, counter))
                sys.stderr.flush()

        if dbout_fname is not None:
            # compress_db(cur, store_coords)
            conn.commit()

    finally:
        p.close()
        if dbout_fname is not None:
            conn.close()
        if out_fname is not None:
            ftxt.close()

    if verbose:
        sys.stderr.write("\n")


def main():
    parser = argparse.ArgumentParser(description='Create DB with pharmacophore representation of input compound '
                                                 'library.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', metavar='conformers.sdf', required=False, default=None,
                        help='SDF or SDF.GZ file with conformers structures. '
                             'Conformers of the same molecule should have identical names. '
                             'Alternatively PKL files can be used as input, which are stored pickled tuples of '
                             'multi-conformer molecules and their names. '
                             'If omitted STDIN will be parsed as SDF file.')
    parser.add_argument('-o', '--out', metavar='output.txt', required=False, default=None,
                        help='output text file, which will contain: mol_id, stereo_id, pharm_hash, stereo_sign. '
                             'Coordinates will not be stored in text file.')
    parser.add_argument('-d', '--dbout', metavar='output.db', required=False, default=None,
                        help='output DB SQLite file.')
    parser.add_argument('-b', '--bin_step', default=1,
                        help='binning step. Default: 1.')
    parser.add_argument('-t', '--tolerance', default=0,
                        help='tolerance volume for the calculation of the stereo sign. If the volume of the '
                             'tetrahedron created by four points less than tolerance then those points are considered '
                             'lying on the same plane (flat; stereo sign is 0). Default: 0.')
    parser.add_argument('-s', '--smarts_features', metavar='smarts.txt', default=None, nargs='?',
                        const=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'smarts_features.txt'),
                        help='text file with definition of pharmacophore features. Tab-separated two columns: '
                             'first - SMARTS, second - feature name. If file name is not specified the default file '
                             'from the script dir will be used.')
    parser.add_argument('--rdkit_factory', metavar='features.fdef', default=None, nargs='?',
                        const=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'smarts_features.fdef'),
                        help='text file with definition of pharmacophore features in RDKit format. If file name is not '
                             'specified the default file from the script dir will be used. This option has '
                             'a priority over smarts_features.')
    parser.add_argument('-f', '--id_field_name', metavar='field_name', default=None,
                        help='field name of compound ID (sdf). If omitted for sdf molecule titles will be used or '
                             'auto-generated names. Please note if you use SDTIN as input, molecule names should be '
                             'stored in a title field of a MOL block, property fields will not be read.')
    parser.add_argument('--stereo_id', action='store_true', default=False,
                        help='set this option if mol names contain stereo_id after last underscore character '
                             '(e.g. MolName_1, Mol_name_2, etc). Then different stereoisomers will be considered '
                             'together when duplicated pharmacophores will be removed. Otherwise each stereoisomer '
                             'will be treated as an individual compound.')
    parser.add_argument('--save_coords', action='store_true', default=False,
                        help='store coordinates of features for each conformer to retrieve full pharmacophore if '
                             'needed. Coordinates can be stored only in DB output. This will enlarge your database.')
    parser.add_argument('--fp', action='store_true', default=False,
                        help='calculate pseudo-3D pharmacophore RDKit fingerprints to speed up screening. '
                             'Works only if rdkit_factory is supplied.')
    parser.add_argument('--nohash', action='store_true', default=False,
                        help='do not generate pharmacophore hashes.')
    parser.add_argument('-r', '--rewrite', action='store_true', default=False,
                        help='rewrite existed DB.')
    parser.add_argument('-c', '--ncpu', default=1,
                        help='number of cpu used. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": conformers_fname = v
        if o == "out": out_fname = v
        if o == "dbout": dbout_fname = v
        if o == "bin_step": bin_step = float(v)
        if o == "smarts_features": smarts_features = v
        if o == "id_field_name": id_field_name = v
        if o == "rewrite": rewrite_db = v
        if o == "verbose": verbose = v
        if o == "ncpu": ncpu = int(v)
        if o == "save_coords": store_coords = v
        if o == "stereo_id": stereo_id = v
        if o == "tolerance": tolerance = float(v)
        if o == "rdkit_factory": rdkit_factory = v
        if o == "fp": fp = v
        if o == "nohash": nohash = v

    if rdkit_factory is None and smarts_features is None:
        rdkit_factory = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'smarts_features.fdef')

    main_params(out_fname=out_fname,
                dbout_fname=dbout_fname,
                smarts_features_fname=smarts_features,
                rdkit_factory=rdkit_factory,
                conformers_fname=conformers_fname,
                bin_step=bin_step,
                rewrite_db=rewrite_db,
                store_coords=store_coords,
                id_field_name=id_field_name,
                stereo_id=stereo_id,
                fp=fp,
                nohash=nohash,
                verbose=verbose,
                ncpu=ncpu,
                tolerance=tolerance)


if __name__ == '__main__':
    main()
