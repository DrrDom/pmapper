#!/usr/bin/env python3
# author          : Pavel
# date            : 12.05.16
# version         : 0.1
# python_version  : 3
# copyright       : Pavel 2016
# license         : GPL3
#==============================================================================

import os
import sys
import time
import argparse
import marshal
import sqlite3 as lite
from collections import defaultdict
from multiprocessing import Process

from pharmacophore import Pharmacophore


# def load_queries(f):
#     q = dict()
#     p = Pharmacophore()
#     # if os.path.splitext(f)[1].lower() == '.pml':
#     #     p.load_ls_model(f)
#     if os.path.splitext(f)[1].lower() == '.pma':
#         p.load_from_pma(f)
#     q[os.path.basename(f)] = p
#     return q


def compare_fp(query_fp, fp):
    return (query_fp & fp) == query_fp


def save_result_of_screen(database, queries_path, file_ts, bin_step, out_path, n_omit):
    start_time = time.time()

    q_list_of_dicts = []
    list_outfiles = []
    for i, filename in enumerate(sorted(os.listdir(queries_path))):
        q_list_of_dicts.append(load_queries(os.path.join(queries_path, filename), bin_step))
        list_outfiles.append(open(os.path.join(out_path, 'screen_{}_{}.txt'.format(os.path.basename(database).split('.')[0],
                                                                                   filename.split('.')[0])), 'w'))
        list_outfiles[i].write('compound\tconf_id\tmodel\tfeature_ids\n')

    with open(file_ts) as fats:
        mol_ts = [el.strip() for el in fats]

    with lite.connect(database) as con:
        cur = con.cursor()
        cur.execute("SELECT DISTINCT(mol_name) FROM conformers")
        mol_names = set(v[0] for v in cur.fetchall()).difference(mol_ts)

        for mol_name in mol_names:
            # load pharmacophores for all conformers
            cur.execute("SELECT conf_id, feature_label, x, y, z FROM feature_coords WHERE conf_id IN "
                        "(SELECT conf_id FROM conformers WHERE mol_name = ?)", (mol_name,))
            res = cur.fetchall()
            confs = defaultdict(list)
            for r in res:
                confs[r[0]].append((r[1], tuple(r[2:])))
            confs_pharm = {}
            for k, v in confs.items():
                confs_pharm[k] = Pharmacophore()
                confs_pharm[k].load_from_feature_coords(v)

            # count features
            features = confs_pharm[list(confs_pharm.keys())[0]].get_features_count()

            for i, q in enumerate(q_list_of_dicts):
                for p_name, p_obj in q.items():
                    # if number of features is less than in query model then pharmacophore cannot fit model and should be omitted
                    if all(features[k] >= v for k, v in p_obj['feature_count'].items()):
                        for conf_id, p in confs_pharm.items():
                            p.update(bin_step=p_obj['model'].get_bin_step())
                            res = p.fit_model(p_obj['model'], n_omitted=n_omit)
                            if res:
                                list_outfiles[i].write('%s\t%s\t%s\t%s\n' % (mol_name, conf_id, p_name, ','.join(map(str, res))))
                                break

    return 'screening {}: ({}s)\n'.format(os.path.basename(database), round(time.time() - start_time,3))


def main(db_fname, query_fname, out_fname, verbose):

    if out_fname is not None:
        out_f = open(out_fname, 'wt')

    # load query
    q = Pharmacophore()
    q.load_from_pma(query_fname)
    q_fp = q.get_fp()

    conn = lite.connect(db_fname)
    cur = conn.cursor()

    cur.execute("SELECT bin_step FROM settings")
    db_bin_step = cur.fetchone()[0]

    # check bin steps
    if q.get_bin_step() != db_bin_step:
        sys.stderr.write('Model has a different bin step from compounds in database. '
                         'It would be skipped.\n')
        raise Exception('Incompatible bin step')

    cur.execute("SELECT DISTINCT(mol_name) FROM conformers")
    mol_names = [i[0] for i in cur.fetchall()]

    for mol_name in mol_names:

        # load fp for all conformers
        cur.execute("SELECT conf_id, fp FROM conformers WHERE mol_name = ?", (mol_name, ))
        data = cur.fetchall()

        conf_ids = []
        for conf_id, fp in data:
            if compare_fp(q_fp, marshal.loads(fp)):
                conf_ids.append(conf_id)

        if conf_ids:
            # load pharmacophores for selected conformers
            sql = "SELECT conf_id, feature_label, x, y, z FROM feature_coords WHERE conf_id IN (%s)" % \
                  ','.join(['?'] * len(conf_ids))
            cur.execute(sql, conf_ids)
            res = cur.fetchall()
            confs = defaultdict(list)
            for r in res:
                confs[r[0]].append((r[1], tuple(r[2:])))
            for conf_id, feature_coords in confs.items():
                p = Pharmacophore()
                p.load_from_feature_coords(feature_coords)
                res = p.fit_model(q)
                if res:
                    res_string = '\t'. join(map(str, (mol_name, conf_id, ','.join(map(str, res))))) + '\n'
                    if out_fname is not None:
                        out_f.write(res_string)
                    else:
                        sys.stdout.write(res_string)
                    break  # match the first conformer


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Screen SQLite DB with compounds against pharmacophore queries. '
                                     'Output is printed out in STDOUT.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--database', metavar='input.db', required=True,
                        help='SQLite input file with compounds to screen.')
    parser.add_argument('-q', '--query', metavar='model.pma', required=True,
                        help='pharmacophore models. Several models can be specified. Model format is recognized by '
                             'file extension. Currently pma (pmapper) and pml (LigandScout) formats are supported.')
    parser.add_argument('-o', '--output', default=None,
                        help='path to output text file with screening results.')
    # parser.add_argument('-n', '--n_omit', metavar='NUMBER', default=0,
    #                     help='number of omitted features.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "database": db_fname = v
        if o == "query": query_fname = v
        if o == "output": out_fname = v
        # if o == "n_omit": n_omit = int(v)
        if o == "verbose": verbose = v

    main(db_fname=db_fname,
         query_fname=query_fname,
         out_fname=out_fname,
         # n_omit=n_omit,
         verbose=verbose)
