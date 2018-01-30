#!/usr/bin/env python
# author          : Pavel
# date            : 30.05.16
# version         : 0.1
# python_version  : 3
# copyright       : Pavel 2016
# license         : GPL3
#==============================================================================

import networkx as nx
import pickle
import json
import numpy as np

from rdkit import Chem
from collections import Counter, defaultdict
from itertools import combinations
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


def load_multi_conf_mol(mol, smarts_features, factory=None, bin_step=2):
    # factory or smarts_featurs should be None to select only one procedure
    if smarts_features is not None and factory is not None:
        raise Exception("Only one options should be not None (smarts_features or factory)")
    output = []
    p = Pharmacophore(bin_step)
    if smarts_features is not None:
        ids = p._get_features_atom_ids(mol, smarts_features)
    elif factory is not None:
        ids = p._get_features_atom_ids_factory(mol, factory)
    else:
        return output
    for conf in mol.GetConformers():
        p = Pharmacophore(bin_step)
        p.load_from_atom_ids(mol, ids, conf.GetId())
        output.append(p)
    return output


class PharmacophoreBase():

    def __init__(self, bin_step=2):
        self.__g = nx.Graph()
        self.__bin_step = bin_step
        self.__nx_version = int(nx.__version__.split('.')[0])
        self.__cached = False
        self.__cached_ids = None
        self.__cached_canon_feature_signatures = None

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
            self.__g.add_node(i, label=label, xyz=coords)
        self.__update_dists()

    def __update_dists(self):
        if self.__nx_version == 2:
            for i, j in combinations(self.__g.nodes(), 2):
                self.__g.add_edge(i, j, dist=self.__dist(self.__g.nodes[i]['xyz'], self.__g.nodes[j]['xyz']))
        else:
            for i, j in combinations(self.__g.nodes(), 2):
                self.__g.add_edge(i, j, dist=self.__dist(self.__g.node[i]['xyz'], self.__g.node[j]['xyz']))

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

    def __get_feature_signatures(self, ids=None, feature_labels=None):
        """
        :param ids: tuple/list of feature ids
        :param feature_labels: tuple/list of feature labels
          if both parameters are supplied they must have the same order and length
        :return: tuple of feature signatures in the same order as supplied to the function
        """
        if ids is not None and feature_labels is not None and len(ids) != len(feature_labels):
            raise Exception('The number of feature ids and labels must be the same')
        if ids is None:
            ids = self.__g.nodes()
        if feature_labels is None:
            feature_labels = tuple(self.__g.node[i]['label'] for i in ids)
        feature_signatures = []
        for num_i, i in enumerate(ids):
            sign = []
            for num_j, j in enumerate(ids):
                if i != j:
                    if self.__nx_version == 2:
                        sign.append((feature_labels[num_j], self.__g.edges[i, j]['dist']))
                    else:
                        sign.append((feature_labels[num_j], self.__g.edge[i][j]['dist']))
            feature_signatures.append((feature_labels[num_i],) + tuple(sorted(sign)))
        return tuple(feature_signatures)

    def __get_canon_feature_signatures(self, ids=None, feature_labels=None, cache_results=False):
        ids = self._get_ids(ids)
        if self.__cached and self.__cached_ids == ids:
            return self.__cached_canon_feature_signatures
        else:
            f = self.__get_feature_signatures(ids=ids, feature_labels=feature_labels)
            # f = tuple(self.__get_dense_feature_sign(s) for s in f)
            f = self.__get_feature_signatures(ids=ids, feature_labels=f)
            if cache_results:
                self.__cached = True
                self.__cached_ids = ids
                self.__cached_canon_feature_signatures = f
            return f

    # @staticmethod
    # def __get_dense_feature_sign(feature_signature):
    #     return ''.join(str(j) for i in feature_signature for j in i)

    def _get_ids(self, ids=None):
        if ids is None:
            ids = self.__g.nodes()
        return tuple(sorted(set(ids)))

    def __get_signature(self, ids=None):
        ids = self._get_ids(ids)
        if self.__cached and self.__cached_ids == ids:
            return tuple(sorted(self.__cached_canon_feature_signatures))
        else:
            return tuple(sorted(self.__get_canon_feature_signatures(ids=ids, cache_results=True)))

    def __get_signature_md5(self, ids=None):
        s = self.__get_signature(ids=ids)
        return md5(pickle.dumps(repr(s)))

    def _get_stereo(self, ids=None, tol=0):
        ids = self._get_ids(ids)
        if len(ids) > 3:
            # sort to get sorted ids
            canon_names = self.__get_canon_feature_signatures(ids=ids, cache_results=True)
            canon_names, ids = self.__sort_two_lists(canon_names, ids)

            if len(set(canon_names)) == len(canon_names):
                # if all points are different sign can be determined by the first four points which are not in one plane
                stereo = 0
                for comb in combinations(ids, 4):
                    if self.__nx_version == 2:
                        s = self.__get_stereo_sign(coord=tuple(self.__g.nodes[i]['xyz'] for i in comb), tol=tol)
                    else:
                        s = self.__get_stereo_sign(coord=tuple(self.__g.node[i]['xyz'] for i in comb), tol=tol)
                    if s != 0:
                        stereo = s
                        break
            else:
                stereo = self.__calc_full_stereo(ids, tol)

        else:
            stereo = 0

        return stereo

    def __calc_full_stereo(self, ids, tol=0):
        # input features and ids are already sorted by feature names (canon_names)

        # sum stereo of individual simplexes
        # achiral objects should have all 0 (right and left simplexes should compensate each other)
        # for chiral objects the stereo is defined by the first non-zero simplex (simplexes are sorted by priority)
        d = defaultdict(int)
        # labels = self.__get_canon_feature_signatures(ids, cache_results=True)
        for comb in combinations(range(len(ids)), 4):
            simplex_ids = tuple(ids[i] for i in comb)
            name, stereo = self.__gen_canon_simplex_name(simplex_ids,
                                                         self.__get_canon_feature_signatures(ids=simplex_ids, cache_results=False),
                                                         # self.__get_canon_feature_signatures(ids=simplex_ids, feature_labels=tuple(labels[i] for i in comb), cache_results=False),
                                                         tol)
            d[name] += stereo

        for k, v in sorted(d.items()):
            if v > 0:
                return 1
            elif v < 0:
                return -1
        return 0

    def __gen_canon_simplex_name(self, feature_ids, feature_names, tol=0):
        # return canon simplex signature and stereo

        c = Counter(feature_names)

        # system AAAA or AAAB or AABC is achiral
        if len(c) == 1 or any(v == 3 for k, v in c.items()) or len(c) == 3:
            stereo = 0

        else:

            names, ids = self.__sort_two_lists(feature_names, feature_ids)

            if self.__nx_version == 2:

                if len(c) == len(names):  # system ABCD
                    stereo = self.__get_stereo_sign(coord=tuple(self.__g.nodes[i]['xyz'] for i in ids), tol=tol)

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
                            names[2], names[3] = names[3], names[2]
                        stereo = self.__get_stereo_sign(coord=tuple(self.__g.nodes[i]['xyz'] for i in ids), tol=tol)

            else:

                if len(c) == len(names):   # system ABCD
                    stereo = self.__get_stereo_sign(coord=tuple(self.__g.node[i]['xyz'] for i in ids), tol=tol)

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
                            names[2], names[3] = names[3], names[2]
                        stereo = self.__get_stereo_sign(coord=tuple(self.__g.node[i]['xyz'] for i in ids), tol=tol)

        return tuple(sorted(feature_names)), stereo

    @staticmethod
    def __sort_two_lists(primary, secondary):
        # sort two lists by order of elements of the primary list
        paired_sorted = sorted(zip(primary, secondary), key=lambda x: x[0])
        return map(list, zip(*paired_sorted))  # two lists

    @staticmethod
    def __get_stereo_sign(coord, tol=0):
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

    def get_signature(self):
        return self.__get_signature()

    def get_signature_md5(self):
        return self.__get_signature_md5().hexdigest()

    def get_signature_md5bin(self):
        return self.__get_signature_md5().digest()

    def get_stereo(self, tol=0):
        return self._get_stereo(tol=tol)

    def get_feature_coords(self):
        return [(v['label'], v['xyz']) for k, v in self.__g.nodes(data=True)]

    def update(self, bin_step=None):
        if bin_step is not None and bin_step != self.__bin_step:
            self.__bin_step = bin_step
            self.__update_dists()
            self.__cached = False

    def iterate_pharm(self, min_features=1, max_features=None, tol=0, return_feature_ids=True):
        ids = self._get_ids()
        if max_features is None:
            max_features = len(self.__g.nodes())
        for n in range(min_features, max_features + 1):
            for comb in combinations(ids, n):
                if return_feature_ids:
                    yield self.__get_signature_md5(ids=comb).hexdigest(), \
                          self._get_stereo(ids=comb, tol=tol), \
                          comb
                else:
                    yield self.__get_signature_md5(ids=comb).hexdigest(), \
                          self._get_stereo(ids=comb, tol=tol)


class PharmacophoreMatch(PharmacophoreBase):

    def __init__(self, bin_step=2):
        super().__init__(bin_step=bin_step)
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

    def __fit_graph(self, g, mapping=False):
        if self.get_bin_step() != 0:
            gm = iso.GraphMatcher(self._Pharmacophore__g, g, node_match=self.__nm, edge_match=self.__em)
        else:
            gm = iso.GraphMatcher(self._Pharmacophore__g, g, node_match=self.__nm, edge_match=iso.numerical_edge_match('dist', 0, atol=0.75))
        if not mapping:
            return gm.subgraph_is_isomorphic()
        else:
            return gm.subgraph_is_isomorphic(), gm.mapping

    def fit_model(self, target, n_omitted=1, essential_features=None):
        """
        target is a target pharmacophore model which is used for matching (it should be a subgraph of the current
            pharmacophore graph).
        n_omitted is a number of possible simultaneously omitted features in target pharmacophore.
        essential_features is a list of ids of features which will not be omitted in target pharmacophore,
            not mentioned features will be omitted iteratively (optional features).
            Default: None - means all features are optional.
        polar_only if True then only polar features of target pharmacophore will be used for matching.

        return: tuple of feature ids of a target (query) model fit to the current pharmacophore
        """

        if self.get_bin_step() != target.get_bin_step():
            raise ValueError('bin_steps in both pharmacophores are different, %f and %f' %
                             (self.get_bin_step(), target.get_bin_step()))

        ids = set(target._get_ids())
        if essential_features:
            ids = ids.union(essential_features)
            optional_features = ids.difference(essential_features)
        else:
            optional_features = ids

        if self.__fit_graph(target.get_graph().subgraph(ids)):
            return tuple(ids)

        for n in range(1, n_omitted + 1):
            for i in combinations(optional_features, n):
                res, mapping = self.__fit_graph(target.get_graph().subgraph(ids.difference(i)), mapping=True)
                if res and self._get_stereo(ids=tuple(mapping.keys())) == target._get_stereo(ids=tuple(ids.difference(i))):
                    return tuple(ids.difference(i))
        return None


class Pharmacophore(PharmacophoreMatch):

    def __init__(self, bin_step=2):
        super().__init__(bin_step=bin_step)

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

    def save_to_pma(self, fname, feature_ids=None):
        coords = self.get_feature_coords()
        if feature_ids:
            coords = [v for i, v in coords if i in feature_ids]
        obj = {'bin_step': self.get_bin_step(), 'feature_coords': coords}
        with open(fname, 'wt') as f:
            f.write(json.dumps(obj))

    def load_from_pma(self, fname):
        with open(fname)as f:
            d = json.loads(f.readline().strip())
            feature_coords = tuple((feature[0], tuple(feature[1])) for feature in d['feature_coords'])
            self.load_from_feature_coords(feature_coords)
            self.update(d['bin_step'])
