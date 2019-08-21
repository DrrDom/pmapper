#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 21-08-2019
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2019
# license         : 
#==============================================================================

from os import path
from rdkit import Chem
from rdkit.Chem import ChemicalFeatures


def load_smarts(filename=None):
    """
    Loads custom SMARTS patterns of features from a file.

    :param filename: name of a text containing SMARTS patterns of pharmacophore features.
                     If None the default patterns will be loaded. Default: None.
    :type filename: str
    :return: dictionary where keys are feature labels and values are tuples of corresponding SMARTS patterns in
             RDKit Mol format
    :rtype: dict

    """
    output = dict()
    if filename is None:
        filename = path.join(path.abspath(path.dirname(__file__)), 'smarts_features.txt')
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line and line[0] != '#':
                tmp = line.split()
                q = Chem.MolFromSmarts(tmp[0])
                if tmp[1] not in output.keys():
                    output[tmp[1]] = [q]
                else:
                    output[tmp[1]].append(q)
    output = {k: tuple(v) for k, v in output.items()}
    return output


def load_factory(filename=None):
    """
    Loads RDKit factory with custom feature patterns from a file.

    :param filename: file name of fdef format file. If None the default patterns will be loaded. Default: None.
    :type filename: str
    :return: object of MolChemicalFeatureFactory class

    """
    if filename is None:
        filename = path.join(path.abspath(path.dirname(__file__)), 'smarts_features.fdef')
    return ChemicalFeatures.BuildFeatureFactory(filename)
