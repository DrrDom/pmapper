#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 14-02-2019
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2019
# license         : 
#==============================================================================

import os
import sys
import argparse
from multiprocessing import Pool, cpu_count
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures
from .pharmacophore import Pharmacophore as P


def get_count(smi, name):
    m = Chem.MolFromSmiles(smi)
    if m:
        p = P()
        q = p._get_features_atom_ids_factory(m, process_factory)
        feat = {k: len(set(v)) for k, v in q.items()}
        for k in process_factory.GetFeatureFamilies():
            if k not in feat:
                feat[k] = 0
        return name, feat
    else:
        return name, dict(zip(process_factory.GetFeatureFamilies(), ['NA'] * len(process_factory.GetFeatureFamilies())))


def get_count_mp(items):
    return get_count(*items)


def read_smi(fname, sep="\t"):
    with open(fname) as f:
        for line in f:
            items = line.strip().split(sep)
            if len(items) == 1:
                yield items[0], items[0]
            else:
                yield items[0], items[1]


def pool_init(fdef_fname):
    global process_factory
    process_factory = ChemicalFeatures.BuildFeatureFactory(fdef_fname)


def entry_point():
    parser = argparse.ArgumentParser(description='Count the number of pharmacophore features of each type.')
    parser.add_argument('-i', '--in', metavar='input.smi', required=True,
                        help='input SMILES file. Should contain mol title as a second field. '
                             'Fields are tab-separated. No header.')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output text file with calculated number of features. '
                             'If input molecules cannot be parsed NA values will be inserted in output.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, default=1,
                        help='Number of CPU cores to use. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in": in_fname = v
        if o == "out": out_fname = v
        if o == "ncpu": ncpu = int(v)
        if o == "verbose": verbose = v

    p = Pool(min(ncpu, cpu_count()),
             initializer=pool_init,
             initargs=[os.path.join(os.path.dirname(os.path.abspath(__file__)), 'smarts_features.fdef')])

    with open(out_fname, 'wt') as f:
        for i, (name, res) in enumerate(p.imap(get_count_mp, read_smi(in_fname), chunksize=100)):
            if i == 0:
                f.write('Name\t' + '\t'.join(['n' + k for k in sorted(res)]) + '\n')  # header
            f.write(name + '\t' + '\t'.join(map(str, [v for _, v in sorted(res.items())])) + '\n')
            if verbose and (i + 1) % 100 == 0:
                sys.stderr.write('\r%i molecules passed' % (i + 1))
                sys.stderr.flush()


if __name__ == '__main__':
    entry_point()
