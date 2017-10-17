#!/usr/bin/env python3
# author          : Pavel Polishchuk
# date            : 14.11.16
# version         : 0.1
# python_version  : 3
# copyright       : Pavel Polishchuk, 2016
# license         : GPL3
#==============================================================================

__author__ = 'Pavel Polishchuk'

import sys
import argparse
import operator
from collections import defaultdict, Counter, OrderedDict

# constants
inf = float("inf")


def read_txt(fname, stereo_id):
    d = defaultdict(set)
    with open(fname) as f:
        f.readline()
        for line in f:
            tmp = line.strip().split()
            if stereo_id:
                d["%s_%s" % (tmp[0], tmp[1])].add((tmp[2], tmp[3]))
            else:
                d[tmp[0]].add((tmp[2], tmp[3]))
    return d


def rank_naive_split(d, ratio_only=False, verbose=False, binned_names=None, pharm_number=-1):
    # binned_names is a dict where values are sets of mol_names in each bin
    if binned_names is None:
        output = rank_naive(d, ratio_only, verbose, pharm_number)
    else:
        output = dict()
        for v in binned_names.values():
            output.update(rank_naive(get_bin(d, v), ratio_only, verbose, pharm_number))
    return output


def rank_naive(d, ratio_only=False, verbose=False, pharm_number=-1):
    # pharm_number should be called with non -1 value only for naive ranking!!!
    output = dict()
    pharm = Counter(v for values in d.values() for v in values)
    uniq_pharm = set(k for k, v in pharm.items() if v == 1)
    for i, (k, v) in enumerate(d.items()):
        u_num = len(uniq_pharm.intersection(v))
        ratio = u_num/len(v) if v else 0  # for empty set of patterns ratio will always be 0
        # set ratio of compound with low number of pharmacopohres to 1, for higher priority selection
        if len(v) <= pharm_number:
            ratio = 1
        if ratio_only:
            output[k] = ratio
        else:
            output[k] = {'uniq': u_num, 'all': len(v), 'ratio': ratio}
        if verbose and i % 1000 == 0:
            sys.stderr.write('\r%i' % i)
            sys.stderr.flush()
    return output


def split_on_bins(d, split):
    splits = list(split) + [float("inf")]
    output = defaultdict(set)
    for k, v in d.items():
        for s in splits:
            if len(v) < s:
                output[s].add(k)
                break
    return output


def get_bin(d, names):
    output = dict()
    for n in names:
        if n in d:
            output[n] = d[n]
    return output


def set_cover_solver(d, min_step, verbose, split=None, pharm_number=-1):

    output = OrderedDict()
    eliminated_mol_ids = []
    uncovered_patterns = set(v for values in d.values() for v in values)
    total_patterns_count = len(uncovered_patterns)
    total_mol_count = len(d)
    step = 0

    binned_names = split_on_bins(d, split) if split is not None else None

    while uncovered_patterns:

        step += 1
        r = rank_naive_split(d, True, False, binned_names)
        max_r = max(r.values())
        if sum(tuple(v == max_r for v in r.values())) < min_step:
            # select top min_step mol_ids according to their R value
            mol_ids = tuple(k for k, v in sorted(r.items(), key=operator.itemgetter(1), reverse=True)[:min_step])
        else:
            mol_ids = tuple(k for k, v in r.items() if v == max_r)

        # add molecules with low number of pharmacophores directly on the first step
        if step == 1 and pharm_number > 0:
            mol_ids += tuple(k for k, v in d.items() if len(v) <= pharm_number)
            mol_ids = tuple(set(mol_ids))

        # remove patterns covered by selected mols
        covered_patterns = set(v for mol_id in mol_ids for v in d[mol_id])
        uncovered_patterns.difference_update(covered_patterns)
        chem_coverage = (len(output) + len(mol_ids)) / total_mol_count
        pharm_coverage = (total_patterns_count - len(uncovered_patterns)) / total_patterns_count
        for mol_id in mol_ids:
            output[mol_id] = {'step': step, 'chemspace_coverage': chem_coverage, 'pharmspace_coverage': pharm_coverage}

        # remove selected mols and covered patterns from the remaining mols
        for mol_id in mol_ids:
            del(d[mol_id])
        for k, v in list(d.items()):
            tmp = v.difference(covered_patterns)
            if tmp:
                d[k] = tmp
            else:
                del d[k]
                eliminated_mol_ids.append(k)

        if verbose:
            sys.stderr.write('step: %i, chem space coverage: %f, pharm space coverage: %f\n' %
                             (step, chem_coverage, pharm_coverage))
            sys.stderr.flush()

    # add redundant mols to the end of output
    step += 1
    for k in d:
        output[k] = {'step': step, 'chemspace_coverage': 1, 'pharmspace_coverage': 1}
    for mol_id in eliminated_mol_ids:
        output[mol_id] = {'step': step, 'chemspace_coverage': 1, 'pharmspace_coverage': 1}

    return output


def get_scores(d, algorithm, covered_patterns):
    output = dict()
    for k, v in d.items():
        new_patterns = len(v.difference(covered_patterns))
        if new_patterns:
            if algorithm == 1:
                output[k] = len(v) ** 0.95 / new_patterns
            elif algorithm == 2:
                output[k] = len(v) ** 1.05 / new_patterns
        else:
            output[k] = inf
    return output


def get_unique_mols(d):
    output = list()
    pharm = Counter(v for values in d.values() for v in values)
    uniq_pharm = set(k for k, v in pharm.items() if v == 1)
    for k, v in d.items():
        u_num = len(uniq_pharm.intersection(v))
        ratio = u_num/len(v) if v else 0  # for empty set of patterns ratio will always be 0
        if ratio == 1:
            output.append(k)
    return output


def set_cover_solver2(d, algorithm, unique, step_size, verbose):

    output = OrderedDict()
    covered_patterns = set()
    discarded_mols = list()
    uncovered_patterns = set(v for values in d.values() for v in values)
    total_patterns_count = len(uncovered_patterns)
    total_mol_count = len(d)
    step = 0

    while uncovered_patterns:

        step += 1

        # select compounds with all unique pharmacophores on the first step
        if step == 1 and unique:

            mol_ids = get_unique_mols(d)
            if not mol_ids:
                max_len = max(map(len, d.values()))
                mol_ids = tuple(k for k, v in d.items() if len(v) == max_len)

        else:

            scores = get_scores(d, algorithm, covered_patterns)

            # select redundant mols and remove them
            new_discarded_mols = tuple(k for k, v in scores.items() if v == inf)
            for m in new_discarded_mols:
                del scores[m]
                del d[m]
            discarded_mols.extend(new_discarded_mols)

            # select mol(s) with the lowest score
            mol_ids = [ k for k, v in sorted(scores.items(), key=operator.itemgetter(1))[:step_size] ]
            # best_score = inf
            # best_mol = None
            # for k, s in scores.items():
            #     if s < best_score:
            #         best_mol = k
            #         best_score = s
            # mol_ids = [best_mol]

        new_covered_patterns = set(v for mol_id in mol_ids for v in d[mol_id])
        covered_patterns.update(new_covered_patterns)
        uncovered_patterns.difference_update(new_covered_patterns)
        chem_coverage = (len(output) + len(mol_ids)) / total_mol_count
        pharm_coverage = (total_patterns_count - len(uncovered_patterns)) / total_patterns_count
        for mol_id in mol_ids:
            output[mol_id] = {'step': step, 'chemspace_coverage': chem_coverage, 'pharmspace_coverage': pharm_coverage}

        # remove selected mols
        for mol_id in mol_ids:
            del d[mol_id]

        if verbose:
            sys.stderr.write('step: %i, selected mols: %i, discarded mols: %i, chem space coverage: %f, pharm space coverage: %f\n' %
                             (step, len(output), len(discarded_mols), chem_coverage, pharm_coverage))
            sys.stderr.flush()

    # add redundant mols to the end of output
    step += 1
    for k in d:
        output[k] = {'step': step, 'chemspace_coverage': 1, 'pharmspace_coverage': 1}
    for mol_id in discarded_mols:
        output[mol_id] = {'step': step, 'chemspace_coverage': 1, 'pharmspace_coverage': 1}

    return output


def save_txt_output(fname, compound_dict):
    with open(fname, 'wt') as f:
        f.write('Compound\tuniq_patterns\tall_patterns\tratio\n')
        for k, v in compound_dict.items():
            f.write('%s\t%i\t%i\t%f\n' % (k, v['uniq'], v['all'], v['ratio']))


def save_txt_output2(fname, compound_dict):
    with open(fname, 'wt') as f:
        f.write('Compound\tstep\tchem_space_coverage\tpharm_space_coverage\n')
        for k, v in compound_dict.items():
            f.write('%s\t%i\t%f\t%f\n' % (k, v['step'], v['chemspace_coverage'], v['pharmspace_coverage']))


def main_params(in_fname, out_fname, algorithm, unique, step_size, stereo_id, verbose):

    d = read_txt(in_fname, stereo_id)
    res = set_cover_solver2(d, algorithm, unique, step_size, verbose)
    save_txt_output2(out_fname, res)


def main():
    parser = argparse.ArgumentParser(description='Rank compounds based on their pharmacophore signatures.')
    parser.add_argument('-i', '--input', metavar='pharmacophores.txt', required=True,
                        help='Input text file (or SQLite DB) with calculated pharmacophore signatures.')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='output text file with ranked compounds.')
    parser.add_argument('-a', '--algorithm', metavar='NUMBER', required=False, default=2,
                        help='currently compounds are ranked based on greedy solution of set cover problem using '
                             'two cost functions. Specify 1 if you would like to prioritize selection of compounds '
                             'covering larger number of pharmacophores; 2 if you would like to prioritize compounds '
                             'with smaller number of covered pharmacophores. Default: 2.')
    parser.add_argument('-u', '--unique', action='store_true', default=False,
                        help='set this option if you want to remove all compounds covering only unique pharmacophores '
                             'directly to the output on the first step. Otherwise all input compounds will be ranked '
                             'that may substantially increase the time of execution.')
    parser.add_argument('-s', '--step_size', metavar='NUMBER', required=False, default=1,
                        help='the number of compounds selected on each iteration (step). Default: 1.')
    parser.add_argument('--stereo_id', action='store_true', default=False,
                        help='set this option if you want to consider stereoisomers (records having the same mol_name '
                             'but different stereo_id) as different molecules.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress.')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "input": in_fname = v
        if o == "out": out_fname = v
        if o == "verbose": verbose = v
        if o == "stereo_id": stereo_id = v
        if o == "algorithm": algorithm = int(v)
        if o == "unique": unique = v
        if o == "step_size": step_size = int(v)

    main_params(in_fname=in_fname,
                out_fname=out_fname,
                algorithm=algorithm,
                unique=unique,
                step_size=step_size,
                stereo_id=stereo_id,
                verbose=verbose)


if __name__ == '__main__':
    main()
