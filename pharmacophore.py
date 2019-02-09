#!/usr/bin/env python
# author          : Pavel
# date            : 30.05.16
# license         : BSD-3
#==============================================================================

import networkx as nx
import pickle
import json
import numpy as np
import random

from rdkit import Chem
from rdkit.Chem import Conformer, rdMolAlign
from rdkit.Geometry import Point3D
from collections import Counter, defaultdict, OrderedDict
from itertools import combinations, product, permutations
from hashlib import md5
from xml.dom import minidom
from networkx.algorithms import isomorphism as iso
from math import sqrt, asin, pi


def read_smarts_feature_file(file_name):
    output = dict()
    with open(file_name) as f:
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


def load_multi_conf_mol(mol, smarts_features=None, factory=None, bin_step=1, cached=False):
    # factory or smarts_featurs should be None to select only one procedure
    if smarts_features is not None and factory is not None:
        raise Exception("Only one options should be not None (smarts_features or factory)")
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


class PharmacophoreBase():

    __primes_vertex = {'a': 2, 'H': 3, 'A': 5, 'D': 7, 'P': 11,'N': 13}
    __primes_edge = (31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137,
                     139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229)

    def __init__(self, bin_step=1, cached=False):
        self.__g = nx.Graph()
        self.__bin_step = bin_step
        self.__cached = cached
        self.__cache = dict()
        self.__nx_version = int(nx.__version__.split('.')[0])

    @staticmethod
    def __remove_dupl(ls):
        # remove duplicates preserving the order of items
        seen = set()
        seen_add = seen.add
        return [x for x in ls if not (x in seen or seen_add(x))]

    def load_from_feature_coords(self, feature_coords):
        # feature_coords: [('A', (1.23, 2.34, 3.45)), ('A', (4.56, 5.67, 6.78)), ...]
        # remove full duplicates from features
        feature_coords = self.__remove_dupl(feature_coords)
        self.__g.clear()
        for i, (label, coords) in enumerate(feature_coords):
            self.__g.add_node(i, label=label, xyz=coords, plabel=PharmacophoreBase.__primes_vertex[label])
        self.__update_dists()

    def __update_dists(self, bin_step=None):
        if self.__nx_version == 2:
            for i, j in combinations(self.__g.nodes(), 2):
                dist = self.__dist(self.__g.nodes[i]['xyz'], self.__g.nodes[j]['xyz'], bin_step)
                self.__g.add_edge(i, j, dist=dist, pdist=PharmacophoreBase.__primes_edge[dist])
        else:
            for i, j in combinations(self.__g.nodes(), 2):
                dist = self.__dist(self.__g.node[i]['xyz'], self.__g.node[j]['xyz'], bin_step)
                self.__g.add_edge(i, j, dist=dist, pdist=PharmacophoreBase.__primes_edge[dist])

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

    # def __get_feature_signatures(self, ids=None, feature_labels=None):
    #     """
    #     :param ids: tuple/list of feature ids
    #     :param feature_labels: tuple/list of feature labels
    #       if both parameters are supplied they must have the same order and length
    #     :return: tuple of feature signatures in the same order as supplied to the function
    #     """
    #     if ids is not None and feature_labels is not None and len(ids) != len(feature_labels):
    #         raise Exception('The number of feature ids and labels must be the same')
    #     ids = self._get_ids(ids)
    #     if feature_labels is None:
    #         feature_labels = tuple(self.__g.node[i]['label'] for i in ids)
    #     feature_signatures = []
    #     for num_i, i in enumerate(ids):
    #         sign = []
    #         for num_j, j in enumerate(ids):
    #             if i != j:
    #                 if self.__nx_version == 2:
    #                     sign.append((feature_labels[num_j], self.__g.edges[i, j]['dist']))
    #                 else:
    #                     sign.append((feature_labels[num_j], self.__g.edge[i][j]['dist']))
    #         feature_signatures.append((feature_labels[num_i],) + tuple(sorted(sign)))
    #     return tuple(feature_signatures)

    # def __get_canon_feature_signatures(self, ids=None, feature_labels=None, short_version=False):
    #     ids = self._get_ids(ids)
    #     f = self.__get_feature_signatures(ids=ids, feature_labels=feature_labels)
    #     # the second iteration is reasonable to use for 5 and more point graphs
    #     # for 4-point graph 1 iteration is enough
    #     if not short_version:
    #         f = self.__get_feature_signatures(ids=ids, feature_labels=f)
    #     return f

    def __get_canon_feature_signatures2(self, ids):

        feature_labels = dict(zip(ids, (self.__g.node[i]['label'] for i in ids)))
        feature_signatures = []
        for i in ids:
            sign = []
            for j in ids:
                if i != j:
                    if self.__nx_version == 2:
                        sign.append('%s%i' % (feature_labels[j], self.__g.edges[i, j]['dist']))
                    else:
                        sign.append('%s%i' % (feature_labels[j], self.__g.edge[i][j]['dist']))
            feature_signatures.append(feature_labels[i] + ''.join(sorted(sign)))
        return tuple(feature_signatures)

    def _get_ids(self, ids=None):
        if ids is None:
            ids = self.__g.nodes()
        return tuple(sorted(set(ids)))

    # def __get_graph_signature(self, ids=None):
    #     ids = self._get_ids(ids)
    #     return tuple(sorted(self.__get_canon_feature_signatures(ids=ids)))

    # def __get_graph_signature_md5(self, ids=None):
    #     s = self.__get_graph_signature(ids=ids)
    #     return md5(pickle.dumps(repr(s)))

    def __get_signature_dict(self, ids, tol):
        d = defaultdict(int)
        for qudruplet_ids in combinations(ids, min(len(ids), 4)):
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

    def __get_full_hash(self, ids=None, tol=0):
        d = self.__get_signature_dict(ids, tol)
        return md5(pickle.dumps(str(tuple(sorted(d.items()))))).hexdigest()

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

                if self.__nx_version == 2:

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

                else:

                    if len(c) == len(feature_names):   # system ABCD
                        stereo = self.__get_quadruplet_stereo(coord=tuple(self.__g.node[i]['xyz'] for i in ids), tol=tol)

                    else:   # system AABB

                        # if A1-B1 == A1-B2 and A2-B1 == A2-B2 distances or A1-B1 == A2-B1 and A1-B2 == A2-B2 then simplex is achiral
                        if (self.__g.edge[ids[0]][ids[2]]['dist'] == self.__g.edge[ids[0]][ids[3]]['dist'] and
                            self.__g.edge[ids[1]][ids[2]]['dist'] == self.__g.edge[ids[1]][ids[3]]['dist']) or \
                           (self.__g.edge[ids[0]][ids[2]]['dist'] == self.__g.edge[ids[1]][ids[2]]['dist'] and
                            self.__g.edge[ids[0]][ids[3]]['dist'] - self.__g.edge[ids[1]][ids[3]]['dist']):
                            stereo = 0
                        else: # swap B vertices to put on the higher position B vertex with a shorter distance to the first A vertex
                            if self.__g.edge[ids[0]][ids[2]]['dist'] > self.__g.edge[ids[0]][ids[3]]['dist']:
                                ids[2], ids[3] = ids[3], ids[2]
                            stereo = self.__get_quadruplet_stereo(coord=tuple(self.__g.node[i]['xyz'] for i in ids), tol=tol)
                            # modifies the sign to distinguish trapeze and parallelogram-like quadruplets
                            stereo += 10 * sign_dihedral_angle(tuple(self.__g.node[ids[i]]['xyz'] for i in [0, 2, 3, 1]))

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
        return self.__bin_step

    def get_graph(self):
        return self.__g.copy()

    def get_features_count(self):
        data = self.__g.nodes(data=True)
        return Counter([item[1]['label'] for item in data])

    def get_signature_md5(self, ids=None, tol=0):
        return self.__get_full_hash(self._get_ids(ids), tol)

    def get_feature_coords(self, ids=None):
        if ids is None:
            return [(v['label'], v['xyz']) for k, v in self.__g.nodes(data=True)]
        else:
            return [(v['label'], v['xyz']) for k, v in self.__g.nodes(data=True) if k in set(ids)]

    def get_mirror_pharmacophore(self):
        p = Pharmacophore()
        coords = tuple((label, (-x, y, z)) for label, (x, y, z) in self.get_feature_coords())
        p.load_from_feature_coords(coords)
        return p

    def update(self, bin_step=None, cached=None):
        if bin_step is not None and bin_step != self.__bin_step:
            self.__bin_step = bin_step
            self.__update_dists(bin_step)
        if cached is not None:
            self.__cached = cached

    def iterate_pharm(self, min_features=1, max_features=None, tol=0, return_feature_ids=True):
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
        take list of feature ids and add one more feature for them to enumerate all possible combinations

        :param fix_ids: list of tuples/lists with feature ids
        :param tol:
        :param return_feature_ids:
        :return:
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
        output = set()
        for h in self.iterate_pharm(min_features, max_features, tol, False):
            random.seed(int(h, 16))
            for i in range(activate_bits):
                output.add(random.randrange(nbits))
        return output

    def get_fp2(self, min_features=3, max_features=3, tol=(0, ), nbits=(2048, ), activate_bits=(1, )):
        # return dict of {(nbits, activate_bits): {bit set}, ...}
        output = defaultdict(set)
        for tol_ in tol:
            for h in self.iterate_pharm(min_features, max_features, tol_, False):
                seed = int(h, 16)
                for nbits_, act_bits_ in product(nbits, activate_bits):
                    random.seed(seed)
                    for i in range(act_bits_):
                        output[(nbits_, act_bits_, tol_)].add(random.randrange(nbits_))
        return output

    def get_descriptors(self, tol=0):
        """

        :param tol: tolerance
        :return: dict of subpharmacophore names (keys) and counts (values)
        """
        ids = self._get_ids(None)
        d = self.__get_signature_dict(ids, tol)
        return {k[2:-1].replace("', ", '|').replace(", ", '|'): v for k, v in d.items()}


class PharmacophoreMol(PharmacophoreBase):

    __feat_dict_mol = {'A': 89, 'P': 15, 'N': 7, 'H': 1, 'D': 66, 'a': 10}

    def __init__(self, bin_step=1, cached=False):
        super().__init__(bin_step, cached)

    def get_mol(self, ids=None):
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


class PharmacophoreMatch(PharmacophoreMol):

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

    def __get_transformation_matrix(self, model, mapping):
        return rdMolAlign.GetAlignmentTransform(self.get_mol(), self.get_mol(model), atomMap=tuple(mapping.items()))[1]

    def __fit_graph(self, model):
        if self.get_bin_step() != 0:
            gm = iso.GraphMatcher(self._PharmacophoreBase__g, model, node_match=self.__nm, edge_match=self.__em)
        else:
            gm = iso.GraphMatcher(self._PharmacophoreBase__g, model, node_match=self.__nm, edge_match=iso.numerical_edge_match('dist', 0, atol=0.75))
        return gm

    def fit_model(self, model, n_omitted=0, essential_features=None, tol=0, get_transform_matrix=False):
        """
        model: a pharmacophore model which is used for matching (it should be a subgraph of the current
            pharmacophore graph).
        n_omitted: a number of simultaneously omitted features in a model pharmacophore
        essential_features: a list of ids of features which will not be omitted in a model pharmacophore,
            not mentioned features can be omitted iteratively (optional features).
            Default: None - means all features are optional.
        tol: tolerance when define stereoconfiguration
        get_transform_matrix: if set, the function will return a transformation matrix as an additional output

        return: tuple of feature ids of a target (query) model fit to the current pharmacophore, additionally a
            transformation matrix can be returned
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
                            if get_transform_matrix:
                                return tuple(mapping.values()), self.__get_transformation_matrix(model, mapping)
                            else:
                                return tuple(mapping.values())
        else:
            gm = self.__fit_graph(model._PharmacophoreBase__g)
            for j, mapping in enumerate(gm.subgraph_isomorphisms_iter()):
                if j == 0:
                    ref = model.get_signature_md5(tol=tol)
                if self.get_signature_md5(ids=tuple(mapping.keys()), tol=tol) == ref:
                    if get_transform_matrix:
                        return tuple(ids), self.__get_transformation_matrix(model, mapping)
                    else:
                        return tuple(ids)
        return None


class Pharmacophore(PharmacophoreMatch):

    __feat_dict_ls = {"A": "HBA", "H": "H", "D": "HBD", "P": "PI", "N": "NI", "a": "AR"}

    def load_from_smarts(self, mol, smarts):
        features_atom_ids = self._get_features_atom_ids(mol, smarts)
        self.load_from_atom_ids(mol, features_atom_ids)

    def load_from_feature_factory(self, mol, factory):
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

    @staticmethod
    def _get_features_atom_ids_factory(mol, factory):
        """
        function works with RDKit feature factory
        output is dict
        {'N': ((1,), (4,5,6)), 'P': ((2,)), 'a': ((7,8,9,10,11)), ...}
        """
        features = factory.GetFeaturesForMol(Chem.AddHs(mol))
        output = dict()
        for f in features:
            name = f.GetFamily()
            if name in output:
                output[name].append(f.GetAtomIds())
            else:
                output[name] = [f.GetAtomIds()]
        output = {k: Pharmacophore.__filter_features(v) for k, v in output.items()}
        return output

    @staticmethod
    def _get_features_atom_ids(mol, smarts_features):
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

        output = dict()

        m = Chem.AddHs(mol)

        for name, smarts_tuple in smarts_features.items():
            tmp = []
            for smarts in smarts_tuple:
                ids = get_map_atom_ids_list(m, smarts)
                if ids:
                    tmp.extend(ids)
            # exclude features if their atom ids is full subset of another feature of the same type
            if tmp:
                res = Pharmacophore.__filter_features(tmp)
                output[name] = tuple(res)  # tuple of tuples with atom ids
            else:
                output[name] = tuple(tmp)  # empty set of ids
        return output

    def load_from_atom_ids(self, mol, atom_features_ids, confId=-1):
        # atom_features_ids - dict of feature labels and list of ids tuples
        # {'A': [(12,), (14)], 'H': [(11,12,13,14,15,16)], ...}
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

        feature_coords = []
        for name, features in atom_features_ids.items():
            for feature in features:
                feature_coords.append((name, get_feature_coord(mol, feature)))
        self.load_from_feature_coords(feature_coords)
        return 0

    def load_ls_model(self, pml_fname):

        feature_names_dict = {'HBD': 'D', 'HBA': 'A', 'NI': 'N', 'PI': 'P', 'H': 'H', 'AR': 'a'}
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
        coords = self.get_feature_coords()
        doc = minidom.Document()
        root = doc.createElement('pharmacophore')
        root.setAttribute('name', name)
        root.setAttribute('pharmacophoreType', 'LIGAND_SCOUT')
        doc.appendChild(root)
        for i, feature in enumerate(coords):
            if feature[0] in self.__feat_dict_ls:
                if feature[0] != "a":
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
        coords = self.get_feature_coords(feature_ids)
        obj = {'bin_step': self.get_bin_step(), 'feature_coords': coords}
        with open(fname, 'wt') as f:
            f.write(json.dumps(obj))

    def load_from_pma(self, fname):
        with open(fname) as f:
            d = json.loads(f.readline().strip())
            feature_coords = tuple((feature[0], tuple(feature[1])) for feature in d['feature_coords'])
            self.load_from_feature_coords(feature_coords)
            self.update(d['bin_step'])

