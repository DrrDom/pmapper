#!/usr/bin/env python
# author          : Pavel
# date            : 30.05.16
# license         : BSD-3
#==============================================================================

import networkx as nx
import json
import numpy as np
import random
import sys
import warnings

from rdkit import Chem
from rdkit.Chem import Conformer, rdMolAlign
from rdkit.Geometry import Point3D
from collections import Counter, defaultdict, OrderedDict
from itertools import combinations, product
from hashlib import md5
from xml.dom import minidom
from networkx.algorithms import isomorphism as iso
from math import sqrt, asin, pi

from .customize import load_smarts

_smarts_patterns = load_smarts()


class __PharmacophoreBase():

    """
    Basic class which stores and manages pharmacophore object

    """

    def __init__(self, bin_step=1, cached=False):
        """
        Initializes Pharmacophore instance

        :param bin_step: binning step
        :type bin_step: float
        :param cached: whether or not to cache intermediate computation results. This substantially increases speed
                       of repeated computation of a hash or fingerprints.
        :type cached: bool

        """
        self.__g = nx.Graph()
        self.__bin_step = bin_step
        self.__cached = cached
        self.__cache = dict()

    def __drop_cache(self):
        self.__cache = dict()

    def _get_cached(self):
        return self.__cached

    @staticmethod
    def __remove_dupl(ls):
        # remove duplicates preserving the order of items
        seen = set()
        seen_add = seen.add
        output = [x for x in ls if not (x in seen or seen_add(x))]
        if len(output) < len(ls):
            warnings.warn('Features with identical labels and coordinates were found and duplicates were removed')
        return output

    def load_from_feature_coords(self, feature_coords):
        # feature_coords: [('A', (1.23, 2.34, 3.45)), ('A', (4.56, 5.67, 6.78)), ...]
        # remove full duplicates from features
        feature_coords = self.__remove_dupl(feature_coords)
        self.__g.clear()
        for i, (label, coords) in enumerate(feature_coords):
            self.__g.add_node(i, label=label, xyz=coords)
        self.__update_dists()

    def add_feature(self, label, xyz):
        """
        Add a feature with its coordinates to the pharmacophore
        :param label: text label of the feature
        :param xyz: tuple of xyz coordinates of the feature
        :return: id of the added feature
        """
        if self.__g.nodes:
            feature_id = max(list(self.__g.nodes)) + 1
        else:
            feature_id = 0
        self.__g.add_node(feature_id, label=label, xyz=tuple(xyz))
        self.__drop_cache()
        self.__update_dists()
        return feature_id

    def __update_dists(self, bin_step=None):
        for i, j in combinations(self.__g.nodes(), 2):
            dist = self.__dist(self.__g.nodes[i]['xyz'], self.__g.nodes[j]['xyz'], bin_step)
            self.__g.add_edge(i, j, dist=dist)

    def __dist(self, coord1, coord2, bin_step=None):
        # coord1, coord2 - tuples of (x, y, z)
        # bin_step: None means default, 0 means no transformation, other number will use as regular
        tmp = sum([(i - j) ** 2 for i, j in zip(coord1, coord2)]) ** 0.5
        if bin_step == 0 or self.__bin_step == 0:
            return tmp
        elif bin_step is None:
            return int(tmp // self.__bin_step)
        else:
            return int(tmp // bin_step)

    def __get_canon_feature_signatures2(self, ids):

        feature_labels = dict(zip(ids, (self.__g.nodes[i]['label'] for i in ids)))
        feature_signatures = []
        for i in ids:
            sign = []
            for j in ids:
                if i != j:
                    sign.append('%s%i' % (feature_labels[j], self.__g.edges[i, j]['dist']))
            feature_signatures.append(feature_labels[i] + ''.join(sorted(sign)))
        return tuple(feature_signatures)

    def _get_ids(self, ids=None):
        if ids is None:
            ids = self.__g.nodes()
        return tuple(sorted(set(ids)))

    def __get_signature_dict(self, ids, tol, ncomb=4):
        d = defaultdict(int)
        for qudruplet_ids in combinations(ids, min(len(ids), ncomb)):
            if self.__cached:
                try:
                    res = self.__cache[qudruplet_ids + (tol,)]
                except KeyError:
                    res = str(self.__gen_quadruplet_canon_name_stereo(qudruplet_ids, tol) + (tol,))
                    self.__cache[qudruplet_ids + (tol,)] = res
            else:
                res = str(self.__gen_quadruplet_canon_name_stereo(qudruplet_ids, tol) + (tol,))
            d[res] += 1
        return d

    def __get_full_hash(self, ids=None, tol=0, hex=True):
        d = self.__get_signature_dict(ids, tol)
        h = md5(str(tuple(sorted(d.items()))).encode())
        if hex:
            return h.hexdigest()
        else:
            return h.digest()

    def __gen_quadruplet_canon_name_stereo(self, feature_ids, tol=0):
        # return canon quadruplet signature and stereo

        def sign_dihedral_angle(coords):
            # Alina code
            b1 = [i1 - i2 for i1, i2 in zip(coords[0], coords[1])]
            b2 = [i1 - i2 for i1, i2 in zip(coords[1], coords[2])]
            b3 = [i1 - i2 for i1, i2 in zip(coords[2], coords[3])]

            n1 = ((b1[1]*b2[2] - b1[2]*b2[1]), -(b1[0]*b2[2] - b1[2]*b2[0]), (b1[0]*b2[1] - b1[1]*b2[0])) # normal to plane 1
            n2 = ((b2[1]*b3[2] - b2[2]*b3[1]), -(b2[0]*b3[2] - b2[2]*b3[0]), (b2[0]*b3[1] - b2[1]*b3[0])) # normal to plane 2

            cos_angle = (n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]) / (sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2]) + sqrt(n2[0]*n2[0] + n2[1]*n2[1] + n2[2]*n2[2]))

            if cos_angle > 0:
                res = 1
            elif cos_angle < 0:
                res = -1
            else:
                res = 0
            return res

        # feature_names = self.__get_canon_feature_signatures(ids=feature_ids, short_version=True)
        feature_names = self.__get_canon_feature_signatures2(feature_ids)

        # less than 4 unique feature coordinates
        if len(set(tuple(coords for (label, coords) in self.get_feature_coords(feature_ids)))) < 4:

            stereo = 0

        else:

            c = Counter(feature_names)

            # system AAAA or AAAB or AABC is achiral
            if len(c) == 1 or any(v == 3 for k, v in c.items()) or len(c) == 3:
                stereo = 0

            else:

                names, ids = self.__sort_two_lists(feature_names, feature_ids)

                if len(c) == len(feature_names):  # system ABCD
                    stereo = self.__get_quadruplet_stereo(coord=tuple(self.__g.nodes[i]['xyz'] for i in ids), tol=tol)

                else:  # system AABB

                    # if A1-B1 == A1-B2 and A2-B1 == A2-B2 distances or A1-B1 == A2-B1 and A1-B2 == A2-B2 then simplex is achiral
                    if (self.__g.edges[ids[0], ids[2]]['dist'] == self.__g.edges[ids[0], ids[3]]['dist'] and
                        self.__g.edges[ids[1], ids[2]]['dist'] == self.__g.edges[ids[1], ids[3]]['dist']) or \
                       (self.__g.edges[ids[0], ids[2]]['dist'] == self.__g.edges[ids[1], ids[2]]['dist'] and
                        self.__g.edges[ids[0], ids[3]]['dist'] - self.__g.edges[ids[1], ids[3]]['dist']):
                        stereo = 0
                    else:  # swap B vertices to put on the higher position B vertex with a shorter distance to the first A vertex
                        if self.__g.edges[ids[0], ids[2]]['dist'] > self.__g.edges[ids[0], ids[3]]['dist']:
                            ids[2], ids[3] = ids[3], ids[2]
                        stereo = self.__get_quadruplet_stereo(coord=tuple(self.__g.nodes[i]['xyz'] for i in ids), tol=tol)
                        # modifies the sign to distinguish trapeze and parallelogram-like quadruplets
                        stereo += 10 * sign_dihedral_angle(tuple(self.__g.nodes[ids[i]]['xyz'] for i in [0, 2, 3, 1]))

        return '|'.join(sorted(feature_names)), stereo

    @staticmethod
    def __sort_two_lists(primary, secondary):
        # sort two lists by order of elements of the primary list
        paired_sorted = sorted(zip(primary, secondary), key=lambda x: x[0])
        return map(list, zip(*paired_sorted))  # two lists

    @staticmethod
    def __get_quadruplet_stereo(coord, tol=0):
        # coord - tuple of tuples containing coordinates of four points in a specific order
        # ((x1, y1, z1), (x2, y2, z2), (x3, y3, z3), (x4, y4, z4))
        # triple product is calculated for 1-2, 1-3 and 1-4 vectors

        def get_angles(p, k):
            # get the angles of the lines which ends in k-point and a plane
            # p is (4 x 3) numpy array
            # http://kitchingroup.cheme.cmu.edu/blog/2015/01/18/Equation-of-a-plane-through-three-points/
            pp = np.delete(p, k, 0)
            v1 = pp[1] - pp[0]
            v2 = pp[2] - pp[0]
            cp = np.cross(v1, v2)
            d = np.dot(cp, pp[2])
            dist = (np.dot(cp, p[k]) - d) / sqrt(sum(cp ** 2))  # signed distance to plane
            return tuple(180 * asin(dist / np.linalg.norm(p[k] - pp[i])) / pi for i in range(3))

        def check_tolerance(p, tol=0):
            # return True if quadruplet within the tolerance range
            # look for a minimal angle regardless sign
            for i in range(4):
                res = get_angles(np.array(p), i)
                if any(-tol <= value <= tol for value in res):
                    return True
            return False
            # res = np.array([get_angles(np.array(p), i) for i in range(4)])
            # return ((-tol <= res) & (res <= tol)).any()

        def det(a):
            return (a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2])
                    - a[1][0] * (a[0][1] * a[2][2] - a[2][1] * a[0][2])
                    + a[2][0] * (a[0][1] * a[1][2] - a[1][1] * a[0][2]))

        if len(set(coord)) < 4 or tol and check_tolerance(coord, tol):
            return 0
        else:
            # calc difference between coord of the first point and all other points
            b = [[j - coord[0][i] for i, j in enumerate(elm)] for elm in coord[1:]]
            # calc the signed volume of parallelepiped
            d = det(b)
            if d > 0:
                return 1
            if d < 0:
                return -1
            else:
                return 0

    def get_bin_step(self):
        """
        Returns current binning step value

        :return: value of a binning step of a pharmacophore
        :rtype: float

        """
        return self.__bin_step

    def get_graph(self):
        """
        Returns a copy of a pharmacophore graph

        :return: a copy of a NetworkX graph object of a pharmacophore

        """
        return self.__g.copy()

    def get_features_count(self):
        """
        Returns number of each feature type

        :return: Counter (dict-like) object with feature labels and their feature count of a pharmacophore

        """
        data = self.__g.nodes(data=True)
        return Counter([item[1]['label'] for item in data])

    def get_signature_md5(self, ids=None, tol=0, hex=True):
        """
        Returns pharmacophore hash.

        :param ids: iterable with feature ids to be used to compute pharmacophore hash
        :type ids: iterable (int)
        :param tol: tolerance value to ignore small deviation of quadruplets of features from planarity. Minimal angle
                    between an edge and a plane formed by other three features. Quadruplets having at least one angle
                    less than tolerance value are assigned 0 chirality.
        :type tol: float
        :param hex: return hex representation (True) or binary (False)
        :type hex: bool
        :return: md5 hash of a pharmacophore
        :rtype: str

        """
        return self.__get_full_hash(self._get_ids(ids), tol, hex)

    def get_feature_coords(self, ids=None):
        """
        Returns coordinates of features.

        :param ids: iterable with feature ids to be used
        :type ids: iterable (int)
        :return: list of 2-tuples where each tuple consists of a label and a 3-tuple with coordinates

        """
        if ids is None:
            return [(v['label'], v['xyz']) for k, v in self.__g.nodes(data=True)]
        else:
            return [(v['label'], v['xyz']) for k, v in self.__g.nodes(data=True) if k in set(ids)]

    def get_feature_ids(self):
        """
        Returns a dict of feature labels and a tuple of corresponding features ids

        :return: dist {label: (id1, id2, ...), ...}
        """
        output = defaultdict(list)
        for i, data in self.__g.nodes(data=True):
            output[data['label']].append(i)
        return {k: tuple(v) for k, v in output.items()}

    def update(self, bin_step=None, cached=None):
        """
        Changes parameters of the pharmacophore instance.

        :param bin_step: binning step.
        :type bin_step: float
        :param cached: whether or not to cache intermediate computation results. This substantially increases speed
                       of repeated computation of a hash or fingerprints.
        :type cached: bool

        """
        if bin_step is not None and bin_step != self.__bin_step:
            self.__bin_step = bin_step
            self.__update_dists(bin_step)
        if cached is not None:
            self.__cached = cached

    def iterate_pharm(self, min_features=1, max_features=None, tol=0, return_feature_ids=True):
        """
        Iterates over subsets of features to get their hashes.

        :param min_features: minimum number of features in a subset
        :type min_features: int
        :param max_features: maximum number of features in a subset
        :type max_features: int
        :param tol: tolerance
        :type tol: float
        :param return_feature_ids: whether or not return feature ids
        :type return_feature_ids: bool
        :return: generator over hashes of feature subsets or over 2-tuples with a hash and a tuple with feature ids.

        """
        ids = self._get_ids()
        if max_features is None:
            max_features = len(self.__g.nodes())
        else:
            max_features = min(max_features, len(self.__g.nodes()))
        for n in range(min_features, max_features + 1):
            for comb in combinations(ids, n):
                if return_feature_ids:
                    yield self.__get_full_hash(ids=comb, tol=tol), comb
                else:
                    yield self.__get_full_hash(ids=comb, tol=tol)

    def iterate_pharm1(self, fix_ids, tol=0, return_feature_ids=True):
        """
        Iterates over subsets of features created by addition of a single feature to the input list of features.

        :param fix_ids: iterable with feature ids which will be used as a constant part of enumerated feature subsets
        :type fix_ids: iterable (int)
        :param tol: tolerance
        :type tol: float
        :param return_feature_ids: whether or not return feature ids
        :type return_feature_ids: bool
        :return: generator over hashes of feature subsets or over 2-tuples with a hash and a tuple with feature ids.

        """
        visited_id_comb = set()
        for ids in fix_ids:
            add_ids = set(self.__g.nodes()).difference(ids)
            for add_id in add_ids:
                i = tuple(sorted(tuple(ids) + (add_id, )))
                if i not in visited_id_comb:
                    visited_id_comb.add(i)
                    if return_feature_ids:
                        yield self.__get_full_hash(ids=i, tol=tol), i
                    else:
                        yield self.__get_full_hash(ids=i, tol=tol)

    def get_fp(self, min_features=3, max_features=3, tol=0, nbits=2048, activate_bits=1):
        """
        Returns a bitstring fingerprint of a pharmacophore encoded by subsets of features

        :param min_features: minimum number of features in a subset
        :type min_features: int
        :param max_features: maximum number of features in a subset
        :type max_features: int
        :param tol: tolerance
        :type tol: float
        :param nbits: length of a bit string
        :type nbits: int
        :param activate_bits: number of activated bits per feature subset
        :type activate_bits: int
        :return: set of numbers of activated bits (bitstring)
        :rtype: set

        """
        output = set()
        for h in self.iterate_pharm(min_features, max_features, tol, False):
            random.seed(int(h, 16))
            for i in range(activate_bits):
                output.add(random.randrange(nbits))
        return output

    def get_fp2(self, min_features=3, max_features=3, tol=(0, ), nbits=(2048, ), activate_bits=(1, )):
        """
        Returns a bitstring fingerprint of a pharmacophore encoded by subsets of features obtained with different setups

        :param min_features: minimum number of features in a subset
        :type min_features: int
        :param max_features: maximum number of features in a subset
        :type max_features: int
        :param tol: iterable with tolerance values
        :type tol: iterable (float)
        :param nbits: iterable with lengths of a bit string
        :type nbits: iterable (int)
        :param activate_bits: iterable with numbers of activated bits per feature subset
        :type activate_bits: iterable (int)
        :return: dictionary where key is a tuple of (nbits, activate_bits, tol) and value is a set of
                 numbers of activated bits (bitstring)
        :rtype: dict

        """
        output = defaultdict(set)
        for tol_ in tol:
            for h in self.iterate_pharm(min_features, max_features, tol_, False):
                seed = int(h, 16)
                for nbits_, act_bits_ in product(nbits, activate_bits):
                    random.seed(seed)
                    for i in range(act_bits_):
                        output[(nbits_, act_bits_, tol_)].add(random.randrange(nbits_))
        return output

    def get_descriptors(self, tol=0, ncomb=4, ids=None):
        """
        Returns count-based descriptor string of a pharmacophore.

        :param tol: tolerance
        :type tol: float
        :param ncomb: number of feature combinations in descriptors, can be an integer from 1 to 4
        :type ncomb: int
        :param ids: feature ids used to calculate descriptors
        :type ids: list
        :return: dictionary where keys are signatures of quadruplets and values are counts of quadruples with
                 identical signatures.
                 Signature of a quadruplet is for example 'AA2A2D2|AA2A3D0|AA2A3D3|DA0A2A3|0|0|1';
                 first blocks are signatures of separate features, then configuration sign, tolerance and
                 binning step value.
        :rtype: dict

        """
        if not isinstance(ncomb, int) or ncomb < 1 or ncomb > 4:
            return tuple()
        ids = self._get_ids(ids)
        d = self.__get_signature_dict(ids, tol, ncomb)
        return {k[2:-1].replace("', ", '|').replace(", ", '|') + f'|{self.__bin_step}': v for k, v in d.items()}


class __PharmacophoreMol(__PharmacophoreBase):

    """
    Class which represents pharmacophore as RDKit Mol object

    """

    __feat_dict_mol = {'A': 89, 'P': 15, 'N': 7, 'H': 1, 'D': 66, 'a': 10, 'T': 11}

    def __init__(self, bin_step=1, cached=False):
        super().__init__(bin_step, cached)

    def get_mol(self, ids=None):
        """
        Returns RDKit Mol object of a pharmacophore where features are replaced with atoms

        :param ids: iterable with feature ids to be used
        :type ids: iterable (int)
        :return: RDKit RWMol

        """
        pmol = Chem.RWMol()
        all_coords = self.get_feature_coords(ids=ids)
        for item in all_coords:
            a = Chem.Atom(self.__feat_dict_mol[item[0]])
            pmol.AddAtom(a)
        c = Conformer(len(all_coords))
        for i, coords in enumerate(all_coords):
            c.SetAtomPosition(i, Point3D(*coords[1]))
        pmol.AddConformer(c, True)
        return pmol


class __PharmacophoreMatch(__PharmacophoreMol):

    """
    Class which implements matching of pharmacopores

    """

    def __init__(self, bin_step=1, cached=False):
        super().__init__(bin_step, cached)
        self.__nm = iso.categorical_node_match('label', '_')
        self.__em = iso.numerical_edge_match('dist', 0)

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['_PharmacophoreMatch__nm']
        del state['_PharmacophoreMatch__em']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.__nm = iso.categorical_node_match('label', '_')
        self.__em = iso.numerical_edge_match('dist', 0)

    def __get_transformation_matrix_and_rms(self, model, mapping):
        return rdMolAlign.GetAlignmentTransform(self.get_mol(), model.get_mol(), atomMap=tuple(mapping.items()))

    def __fit_graph(self, model):
        if self.get_bin_step() != 0:
            gm = iso.GraphMatcher(self._PharmacophoreBase__g, model, node_match=self.__nm, edge_match=self.__em)
        else:
            gm = iso.GraphMatcher(self._PharmacophoreBase__g, model, node_match=self.__nm, edge_match=iso.numerical_edge_match('dist', 0, atol=0.75))
        return gm

    def fit_model(self, model, n_omitted=0, essential_features=None, tol=0, get_transform_matrix=False, get_rms=False):
        """
        Matches the supplied pharmacophore model.

        :param model: a pharmacophore model which is used for matching (it should be a subgraph of the current
                      pharmacophore graph).
        :type model: Pharmacophore
        :param n_omitted: a number of simultaneously omitted features in a model pharmacophore
        :type n_omitted: int
        :param essential_features: a list of ids of features which will not be omitted in a model pharmacophore.
        :type essential_features: list
        :param tol: tolerance
        :type tol: float
        :param get_transform_matrix: if set, the function will return a transformation matrix as an additional output
                                     to align the pharmacophore to a model
        :type get_transform_matrix: bool
        :param get_rms: if set True, the function will return rms value of matched pharmacophore
        :return: None if a model does not match; if matched the output depends on set arguments.
                 By default output is a dictionary of current pharmacophore feature ids (keys) and mapped model feature
                 ids (values). If get_transform_matrix or get_rms arguments were set to True the output will be
                 2 or 3-tuple, where the first item if a tuple of matched features ids, the next items are
                 a transformation matrix and/or rms value.

        """

        if self.get_bin_step() != model.get_bin_step():
            raise ValueError('bin_steps in both pharmacophores are different, %f and %f' %
                             (self.get_bin_step(), model.get_bin_step()))

        ids = set(model._get_ids())
        if essential_features:
            if not ids.issuperset(essential_features):
                raise Exception('bin_steps in both pharmacophores are different, %f and %f' %
                                (self.get_bin_step(), model.get_bin_step()))
            optional_features = ids.difference(essential_features)
        else:
            optional_features = ids

        # check first for graph isomorphism, if yes iteratively evaluate hashes and compare with reference
        if n_omitted:
            for n in range(1, n_omitted + 1):
                for i in combinations(optional_features, n):
                    gm = self.__fit_graph(model._PharmacophoreBase__g.subgraph(ids.difference(i)))
                    for j, mapping in enumerate(gm.subgraph_isomorphisms_iter()):
                        if j == 0:
                            ref = model.get_signature_md5(ids=tuple(mapping.values()), tol=tol)
                        if self.get_signature_md5(ids=tuple(mapping.keys()), tol=tol) == ref:
                            output = mapping
                            if get_transform_matrix or get_rms:
                                output = (output,)
                                rms, matrix = self.__get_transformation_matrix_and_rms(model, mapping)
                                if get_transform_matrix:
                                    output += (matrix, )
                                if get_rms:
                                    output += (rms, )
                                return output
                            else:
                                return output
        else:
            gm = self.__fit_graph(model._PharmacophoreBase__g)
            for j, mapping in enumerate(gm.subgraph_isomorphisms_iter()):
                if j == 0:
                    ref = model.get_signature_md5(tol=tol)
                if self.get_signature_md5(ids=tuple(mapping.keys()), tol=tol) == ref:
                    output = mapping
                    if get_transform_matrix or get_rms:
                        output = (output, )
                        rms, matrix = self.__get_transformation_matrix_and_rms(model, mapping)
                        if get_transform_matrix:
                            output += (matrix, )
                        if get_rms:
                            output += (rms, )
                        return output
                    else:
                        return output
        return None


class __PharmacophoreLoadedMol(__PharmacophoreMatch):

    """
    This class manages a loaded molecule
    """

    def __init__(self, bin_step=1, cached=False):
        super().__init__(bin_step, cached)
        self.atom_ids = None    # will store list of tuples (('A', (1,)), ('P', (5,6,7)), ...),
                                # the order should be preserved while reading to keep the correspondence with
                                # feature ids which will be computed in load_from_feature_coords

    def load_from_mol(self, mol):
        """
        Creates pharmacophore from RDKit Mol. Uses default definition of feature SMARTS.

        :param mol: RDKit Mol object
        :return: nothing

        """
        self.load_from_smarts(mol, _smarts_patterns)

    def load_from_smarts(self, mol, smarts):
        """
        Creates pharmacophore from RDKit Mol and features encoded by custom SMARTS.

        :param mol: RDKit Mol object
        :param smarts: dictionary of SMARTS of features obtained with `load_smarts` function from `pmapper.util` module
        :return: nothing

        """
        features_atom_ids = self._get_features_atom_ids(mol, smarts)
        self.load_from_atom_ids(mol, features_atom_ids)

    def load_from_feature_factory(self, mol, factory):
        """
        Creates pharmacophore from RDKit Mol and features encoded by custom RDKit feature factory.

        :param mol: RDKit Mol object
        :param factory: object of MolChemicalFeatureFactory class
        :return: nothing

        """
        # factory is ChemicalFeatures.BuildFeatureFactory(fdefName)
        features_atom_ids = self._get_features_atom_ids_factory(mol, factory)
        self.load_from_atom_ids(mol, features_atom_ids)

    @staticmethod
    def __filter_features(atom_ids):
        # function exclude tuples if they are a complete subset of a bigger tuples to avoid superfluous features
        # for example the hydrophobic features has ids (1,2) and (1,2,3,4,5), only the second will be kept
        # another example the hydrophobic features has ids (1,2) and (2,3,4,5), both will be kept
        # atom ids is a list of tuples with atom ids [(1,), (1,2), (2,3,4)]
        # OUTPUT: tuple of tuples ((1,2), (2,3,4))
        tmp = sorted(tuple(set(atom_ids)), key=len, reverse=True)
        res = [tmp[0]]
        for item in tmp[1:]:
            item = set(item)
            if any(item.issubset(bigger) for bigger in res):
                continue
            else:
                res.append(tuple(item))
        return tuple(res)

    def _get_features_atom_ids_factory(self, mol, factory):
        """
        function works with RDKit feature factory
        output is dict
        {'N': ((1,), (4,5,6)), 'P': ((2,)), 'a': ((7,8,9,10,11)), ...}
        """
        features = factory.GetFeaturesForMol(Chem.AddHs(mol))
        output = OrderedDict()
        for f in features:
            name = f.GetFamily()
            if name in output:
                output[name].append(f.GetAtomIds())
            else:
                output[name] = [f.GetAtomIds()]
        output = {k: self.__filter_features(v) for k, v in output.items()}
        return output

    def _get_features_atom_ids(self, mol, smarts_features):
        """
        function works with explicit smarts definitions of features
        output is dict
        {'N': ((1,), (4,5,6)), 'P': ((2,)), 'a': ((7,8,9,10,11)), ...}
        """

        def get_map_atom_ids_list(m, query):
            """
            Returns list of tuples with atom ids if match and None otherwise
            ids started from 0, 1, 2, 3, ...
            """
            ids = m.GetSubstructMatches(query)
            ids = [tuple(sorted(i)) for i in ids]
            return ids

        output = OrderedDict()

        m = Chem.AddHs(mol)

        for name, smarts_tuple in smarts_features.items():
            tmp = []
            for smarts in smarts_tuple:
                ids = get_map_atom_ids_list(m, smarts)
                if ids:
                    tmp.extend(ids)
            # exclude features if their atom ids is full subset of another feature of the same type
            if tmp:
                res = self.__filter_features(tmp)
                output[name] = tuple(res)  # tuple of tuples with atom ids
            else:
                output[name] = tuple(tmp)  # empty set of ids
        return output

    def load_from_atom_ids(self, mol, atom_features_ids, confId=-1):
        """
        Creates pharmacophore from RDKit Mol and atom ids subsets associated with particular features

        :param mol: RDKit Mol object
        :param atom_features_ids: dictionary where keys are feature labels and values are lists of tuples with atom
                                  ids of individual features, e.g.
                                  {'A': [(12,), (14,)], 'H': [(11,12,13,14,15,16)], ...}
        :param confId: id of a conformer in a molecule
        :return: nothing

        """
        # atom_features_ids - dict of feature labels and list of ids tuples
        # {'A': [(12,), (14,)], 'H': [(11,12,13,14,15,16)], ...}
        # return list of tuples containing feature labels with their coordinates
        # [('A', (1.23, 2.34, 3.45)), ('A', (4.56, 5.67, 6.78)), ...]

        def get_feature_coord(molecule, ids):
            if len(ids) == 1:
                return tuple(molecule.GetConformer(confId).GetAtomPosition(ids[0]))
            else:
                x = 0
                y = 0
                z = 0
                for i in ids:
                    coord = molecule.GetConformer(confId).GetAtomPosition(i)
                    x += coord.x
                    y += coord.y
                    z += coord.z
                return x / len(ids), y / len(ids), z / len(ids)

        self.atom_ids = []
        feature_coords = []
        for name, features in atom_features_ids.items():
            for feature in features:
                feature_coords.append((name, get_feature_coord(mol, feature)))
                self.atom_ids.append((name, feature))
        self.load_from_feature_coords(feature_coords)
        self.atom_ids = tuple(self.atom_ids)
        return 0

    def get_descriptors_atom_exclusion(self, natoms, tol=0, ncomb=4):
        """
        Return descriptors for subsets of features excluding every atom in an input molecule

        :param natoms: total number of heavy atoms in a molecule
        :type natoms: int
        :param tol: tolerance
        :type tol: float
        :param ncomb: number of feature combinations in descriptors, can be an integer from 1 to 4
        :type ncomb: int
        :return: list of dictionaries where keys are signatures of quadruplets and values are counts of quadruples with
                 identical signatures.
                 Signature of a quadruplet is for example 'AA2A2D2|AA2A3D0|AA2A3D3|DA0A2A3|0|0|1';
                 first blocks are signatures of separate features, then configuration sign, tolerance and
                 binning step value.
                 Every entry of the list is the descriptors calculated with removal of an atom of
                 the corresponding index, the first entry - the first atom was removed.
        :rtype: dict
        """
        output = []
        for i in range(natoms):
            ids = []
            for j, item in enumerate(self.atom_ids):
                if i not in item[1]:
                    ids.append(j)
            output.append(self.get_descriptors(tol=tol, ncomb=ncomb, ids=ids))
        return output


class __PharmacophoreFiles(__PharmacophoreLoadedMol):

    """
    The class to support reading and writing pharmacophores in different formats

    """

    __feat_dict_ls = {"A": "HBA", "H": "H", "D": "HBD", "P": "PI", "N": "NI", "a": "AR", "e": "exclusion"}

    def load_ls_model(self, pml_fname):
        """
        Reads pharmacophore from LigandScout pml-file

        :param pml_fname: file name of a LigandScout pml-file
        :return: nothing

        """

        feature_names_dict = {'HBD': 'D', 'HBA': 'A', 'NI': 'N', 'PI': 'P', 'H': 'H', 'AR': 'a', 'exclusion': 'e'}
        coord = []

        d = minidom.parse(pml_fname)
        for p in d.getElementsByTagName('point'):
            feature_name = feature_names_dict[p.getAttribute('name')]
            optional = p.getAttribute('optional')
            disabled = p.getAttribute('disabled')
            if optional == disabled == 'false':
                for pos in p.getElementsByTagName('position'):
                    c = tuple(float(pos.getAttribute(i)) for i in ('x3', 'y3', 'z3'))
                    coord.append((feature_name, c))
        for p in d.getElementsByTagName('plane'):
            feature_name = feature_names_dict[p.getAttribute('name')]
            optional = p.getAttribute('optional')
            disabled = p.getAttribute('disabled')
            if optional == disabled == 'false':
                for pos in p.getElementsByTagName('position'):
                    c = tuple(float(pos.getAttribute(i)) for i in ('x3', 'y3', 'z3'))
                    coord.append((feature_name, c))
        for p in d.getElementsByTagName('vector'):
            feature_name = feature_names_dict[p.getAttribute('name')]
            optional = p.getAttribute('optional')
            disabled = p.getAttribute('disabled')
            pointsToLigand = p.getAttribute('pointsToLigand')
            if optional == disabled == 'false':
                if pointsToLigand == 'false':
                    poses = p.getElementsByTagName('origin')
                elif pointsToLigand == 'true':
                    poses = p.getElementsByTagName('target')
                for pos in poses:
                    c = tuple(float(pos.getAttribute(i)) for i in ('x3', 'y3', 'z3'))
                    coord.append((feature_name, c))

        import pprint
        pprint.pprint(coord)

        self.load_from_feature_coords(coord)

    def save_ls_model(self, fname, name="pmapper_pharmcophore"):
        """
        Saves pharmacophore to LigandScout pml-file.

        :param fname: pml-file name
        :param name: name of a pharmacophore which would be stored in a file and will be displayed in LigandScout
        :return: nothing

        """
        coords = self.get_feature_coords()
        doc = minidom.Document()
        root = doc.createElement('pharmacophore')
        root.setAttribute('name', name)
        root.setAttribute('pharmacophoreType', 'LIGAND_SCOUT')
        doc.appendChild(root)
        for i, feature in enumerate(coords):
            if feature[0] in self.__feat_dict_ls:
                if feature[0] == "e":
                    point = doc.createElement('volume')
                elif feature[0] != "a":
                    point = doc.createElement('point')
                else:
                    point = doc.createElement('plane')
                point.setAttribute('name', self.__feat_dict_ls[feature[0]])
                point.setAttribute('featureId', str(i))
                point.setAttribute('optional', 'false')
                point.setAttribute('disabled', 'false')
                point.setAttribute('weight', '1.0')
                point.setAttribute('id', 'feature' + str(i))
                root.appendChild(point)
                position = doc.createElement('position')
                for k, j in zip(['x3', 'y3', 'z3'], feature[1]):
                    position.setAttribute(k, str(j))
                position.setAttribute('tolerance', '1')
                point.appendChild(position)
                if feature[0] == "a":
                    normal = doc.createElement('normal')
                    normal.setAttribute('x3', '1')
                    normal.setAttribute('y3', '1')
                    normal.setAttribute('z3', '1')
                    normal.setAttribute('tolerance', '0.5')
                    point.appendChild(normal)

        with open(fname, 'w') as f:
            f.write(doc.toprettyxml(indent="  "))

    def save_to_pma(self, fname, feature_ids=None):
        """
        Saves pharmacophore in json format. This is a native way to store pharmacophore objects in a readable format.

        :param fname: pma-file name
        :param feature_ids: ids of features which should be stored. Default: None (all features).
        :return: nothing

        """
        coords = self.get_feature_coords(feature_ids)
        obj = {'bin_step': self.get_bin_step(), 'feature_coords': coords}
        with open(fname, 'wt') as f:
            f.write(json.dumps(obj))

    def load_from_pma(self, fname):
        """
        Reads pharmacophore from a pma-file.

        :param fname: pma-file name
        :return: nothing

        """
        with open(fname) as f:
            d = json.loads(f.readline().strip())
            feature_coords = tuple((feature[0], tuple(feature[1])) for feature in d['feature_coords'])
            self.load_from_feature_coords(feature_coords)
            self.update(d['bin_step'])

    def load_from_xyz(self, fname):
        """
        Reads pharmacophore from xyz-file.

        :param fname: xyz-file name
        :return: nothing

        """
        with open(fname) as f:
            feature_coords = []
            f.readline()
            line = f.readline().strip()
            if line:
                opts = dict(item.split('=') for item in line.split(';'))
                if 'bin_step' in opts:
                    self.update(bin_step=float(opts['bin_step']))
            for line in f:
                label, *coords = line.strip().split()
                coords = tuple(map(float, coords))
                feature_coords.append((label, coords))
            self.load_from_feature_coords(tuple(feature_coords))

    def save_to_xyz(self, fname, feature_ids=None, ndigits=2):
        """
        Save pharmacophoe to xyz format

        :param fname: xyz-file name
        :param feature_ids: ids of features which should be stored. Default: None (all features).
        :param ndigits: number of digits after the dot to save coordinates
        :return: nothing
        """
        with open(fname, 'wt') as f:
            f.write('\n')
            f.write(f'bin_step={self.get_bin_step()}\n')
            for i, (x, y, z) in self.get_feature_coords(feature_ids):
                f.write(f'{i} {round(x, ndigits)} {round(y, ndigits)} {round(z, ndigits)}\n')

    def load_from_file(self, fname):
        """
        Reads pharmacophore from file. File format will be recognized by extension.

        :param fname: file name
        :return: nothing
        """
        ext = fname.rsplit('.', 1)[1]
        if ext == 'pml':
            self.load_ls_model(fname)
        elif ext == 'pma':
            self.load_from_pma(fname)
        elif ext == 'xyz':
            self.load_from_xyz(fname)
        else:
            sys.stderr.write('Unknown extension %s. Pharmacophore was not loaded. '
                             'Only pml, pma and xyz formats are supported.' % fname)

    def load_from_pharmit(self, fname):
        """
        Reads pharmacophore from pharmit json-file. Only enabled features will be read. Directed features will be
        reduced to undirected ones (arrows to spheres). Radius will be ignored at read.

        :param fname: json-file name with a pharmit pharmacophore model
        :return: nothing
        """
        feature_names_dict = {'HydrogenDonor': 'D', 'HydrogenAcceptor': 'A', 'NegativeIon': 'N', 'PositiveIon': 'P',
                              'Hydrophobic': 'H', 'Aromatic': 'a'}

        coords = []
        p = json.load(open(fname))['points']
        for i in p:
            if i['enabled']:
                # use pmapper feature labels if they are in the list of pre-defined (standard) features
                label = feature_names_dict[i['name']] if i['name'] in feature_names_dict.keys() else i['name']
                coords.append((label, (i['x'], i['y'], i['z'])))
        self.load_from_feature_coords(tuple(coords))

    def save_to_pharmit(self, fname, feature_ids=None, ndigits=2):
        """
        Save pharmacophore in pharmit format. All radius will be set to bin_step of the pharmacophore. Labels absent
        in the standard feature dictionary will be saved as is.

        :param fname: json-file name
        :param feature_ids: ids of features which should be stored. Default: None (all features).
        :param ndigits: number of digits after the dot to save coordinates
        :return: nothing
        """
        feature_names_dict = {'D': 'HydrogenDonor', 'A': 'HydrogenAcceptor', 'N': 'NegativeIon', 'P': 'PositiveIon',
                              'H': 'Hydrophobic', 'a': 'Aromatic'}

        data = {'points': []}
        bin_step = self.get_bin_step()
        for label, (x, y, z) in self.get_feature_coords(ids=feature_ids):
            label = feature_names_dict[label] if label in feature_names_dict.keys() else label
            data['points'].append({'name': label,
                                   'x': round(x, ndigits),
                                   'y': round(y, ndigits),
                                   'z': round(z, ndigits),
                                   'radius': bin_step,
                                   'enabled': True})
        with open(fname, 'wt') as f:
            json.dump(data, f, indent=4)


class Pharmacophore(__PharmacophoreFiles):

    """
    Main class

    """

    def get_mirror_pharmacophore(self):
        """
        Returns a new mirrored Pharmacophore instance. Bin step and cached arguments will be the same as for
        the source instance.

        :return: a new instance of a Pharmacophore class with all features mirrored in yz-plane.

        """
        p = Pharmacophore(bin_step=self.get_bin_step(), cached=self._get_cached())
        coords = tuple((label, (-x, y, z)) for label, (x, y, z) in self.get_feature_coords())
        p.load_from_feature_coords(coords)
        return p

    def get_subpharmacophore(self, ids=None):
        """
        Returns a new Pharmacophore instance with chosen features. Bin step and cached arguments will be the same as
        for the source instance.

        :param ids: a list of chosen feature ids or None (all features)
        :return: a new instance containing only selected feature.

        """
        p = Pharmacophore(bin_step=self.get_bin_step(), cached=self._get_cached())
        p.load_from_feature_coords(self.get_feature_coords(ids=ids))
        return p
