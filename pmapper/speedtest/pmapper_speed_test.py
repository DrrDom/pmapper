#!/usr/bin/env python3

import argparse
import pickle
import os
from collections import defaultdict
from timeit import default_timer as timer
from pmapper.utils import load_multi_conf_mol


def __read_pkl(fname):
    with open(fname, 'rb') as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break


def entry_point():
    parser = argparse.ArgumentParser(description='Measure performance of calculation of 3D pharmacophore hashes for '
                                                 'DrugBank compounds. Enabled caching increase speed for repeated '
                                                 'calculation of hashes for the same pharmacophores.')

    print('========== Reading of conformers of molecules ==========')
    mols = dict()
    start = timer()
    for i, (mol, mol_name) in enumerate(__read_pkl(os.path.join(os.path.dirname(__file__), 'confs.pkl'))):
        mols[mol_name] = mol
    print(f'{len(mols)} molecules were read in {round(timer() - start, 5)} s')

    print('\n========== Creation of pharmacophores (with enabled caching) ==========')
    d = defaultdict(list)
    start = timer()
    for i, (mol_name, mol) in enumerate(mols.items()):
        p = load_multi_conf_mol(mol, bin_step=1, cached=True)
        nfeat = sum(p[0].get_features_count().values())
        d[nfeat].extend(p)
        sum(len(i) for i in d.values())
    print(f'{sum(len(i) for i in d.values())} pharmacophores were created in {round(timer() - start, 5)} s')

    print('\n========== First calculation of hashes ==========')
    for k in sorted(d.keys()):
        start = timer()
        for p in d[k]:
            p.get_signature_md5()
        t = timer() - start
        print(f'{len(d[k])} pharmacophores with {k} features - {round(t, 5)}s or {round(t / len(d[k]), 5)}s per pharmacophore')

    print('\n========== Second calculation of hashes of the same pharmacophores ==========')
    for k in sorted(d.keys()):
        start = timer()
        for p in d[k]:
            p.get_signature_md5()
        t = timer() - start
        print(f'{len(d[k])} pharmacophores with {k} features - {round(t, 5)}s or {round(t / len(d[k]), 5)}s per pharmacophore')


if __name__ == '__main__':
    entry_point()

