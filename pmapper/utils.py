#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 17-07-2019
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2019
# license         : 
#==============================================================================

from os import path
from rdkit import Chem
from rdkit.Chem import AllChem, ChemicalFeatures
from .pharmacophore import Pharmacophore


def load_smarts(filename=path.join(path.abspath(path.dirname(__file__)), 'smarts_features.txt')):
    """
    Load feature SMARTS patterns from a file.

    :param filename: name of a text containing SMARTS patterns of pharmacophore features
    :type filename: str
    :return: dictionary where keys are feature labels and values are tuples of corresponding SMARTS patterns in
             RDKit Mol format

    """
    output = dict()
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


def load_factory(filename=path.join(path.abspath(path.dirname(__file__)), 'smarts_features.fdef')):
    """
    Load RDKit factory with feature patterns from a file.

    :param filename: file name of fdef format file
    :type filename: str
    :return: object of MolChemicalFeatureFactory class
    """
    return ChemicalFeatures.BuildFeatureFactory(filename)


def load_multi_conf_mol(mol, smarts_features=None, factory=None, bin_step=1, cached=False):
    """
    Convenience function which loads all conformers of a molecule into a list of pharmacophore objects.

    :param mol: RDKit Mol
    :param smarts_features: dictionary of SMARTS of features obtained with `load_smarts` function from `pmapper.util`
                            module
    :param factory: RDKit MolChemicalFeatureFactory loaded with `load_factory` function from `pmapper.util` module
    :param bin_step: binning step
    :param cached: whether or not to cache intermediate computation results. This substantially increases speed
                   of repeated computation of a hash or fingerprints.
    :return: list of pharmacophore objects

    """
    # factory or smarts_features should be None to select only one procedure
    if smarts_features is not None and factory is not None:
        raise ValueError("Only one options should be not None (smarts_features or factory)")
    output = []
    p = Pharmacophore(bin_step, cached)
    if smarts_features is not None:
        ids = p._get_features_atom_ids(mol, smarts_features)
    elif factory is not None:
        ids = p._get_features_atom_ids_factory(mol, factory)
    else:
        return output
    for conf in mol.GetConformers():
        p = Pharmacophore(bin_step, cached)
        p.load_from_atom_ids(mol, ids, conf.GetId())
        output.append(p)
    return output


def get_rms(p1, p2):
    """
    Calculates RMS based om RDKit Mol representations of pharmacophores.

    :param p1: the first Pharmacophore class object
    :param p2: the second Pharmacophore class object
    :return: best rms value between two pharmacophores. Value -1 will be returned if two
             pharmacophores have different sets of features, e.g. aaHD and aaHDD
    :rtype: float

    """

    def get_pharm_str(p):
        return tuple(sorted(p.get_features_count().items()))

    if get_pharm_str(p1) == get_pharm_str(p2):
        res = AllChem.GetBestRMS(p1.get_mol(), p2.get_mol(), 0, 0)
    else:
        res = -1
    return res

