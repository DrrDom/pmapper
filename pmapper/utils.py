#!/usr/bin/env python3
#==============================================================================
# author          : Pavel Polishchuk
# date            : 17-07-2019
# version         : 
# python_version  : 
# copyright       : Pavel Polishchuk 2019
# license         : 
#==============================================================================

from rdkit.Chem import AllChem
from .pharmacophore import Pharmacophore
from .customize import load_smarts


def load_multi_conf_mol(mol, smarts_features=None, factory=None, bin_step=1, cached=False):
    """
    Convenience function which loads all conformers of a molecule into a list of pharmacophore objects.

    :param mol: RDKit Mol
    :param smarts_features: dictionary of SMARTS of features obtained with `load_smarts` function from `pmapper.util`
                            module. Default: None.
    :param factory: RDKit MolChemicalFeatureFactory loaded with `load_factory` function from `pmapper.util` module.
                    Default: None.
    :param bin_step: binning step
    :param cached: whether or not to cache intermediate computation results. This substantially increases speed
                   of repeated computation of a hash or fingerprints.
    :return: list of pharmacophore objects

    Note: if both arguments `smarts_features` and `factory` are None the default patterns will be used.

    """
    # factory or smarts_features should be None to select only one procedure
    if smarts_features is not None and factory is not None:
        raise ValueError("Only one options should be not None (smarts_features or factory)")
    if smarts_features is None and factory is None:
        smarts_features = __smarts_patterns
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
    Calculates best RMSD between two pharmacophores.

    :param p1: the first Pharmacophore class object
    :param p2: the second Pharmacophore class object
    :return: best RMSD value between two pharmacophores. Value -1 will be returned if two
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


__smarts_patterns = load_smarts()
