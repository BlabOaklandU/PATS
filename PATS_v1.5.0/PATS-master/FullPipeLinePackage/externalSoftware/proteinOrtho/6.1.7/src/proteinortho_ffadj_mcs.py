#!/usr/bin/env python3

from sys import stderr, exit, argv, maxsize
from copy import deepcopy
from bisect import bisect
from itertools import product
from random import randint
from math import ceil
import logging as log

class BothStrands:

    def __eq__(self, x):
        return x == '+' or x == '-' or isinstance(x, BothStrands)

    def __str__(self):
        return '+/-'


class Run:
    direction = None
    startG1 = None
    startG2 = None
    endG1 = None
    endG2 = None
    weight = None

    def __init__(self, startG1, startG2, weight, direction):
        self.direction = direction
        self.startG1 = startG1
        self.startG2 = startG2
        self.endG1 = startG1
        self.endG2 = startG2
        self.weight = list()
        self.weight.append(weight)

    def getWeight(self, alpha):
        adjTerm = 0
        if len(self.weight) > 1:
            adjTerm = sum([self.weight[i] * self.weight[i+1] for i in
                           range(len(self.weight)-1)])
        edgeTerm = sum([w**2 for w in self.weight])
        return alpha * adjTerm + (1-alpha) * edgeTerm

    def extendRun(self, nextG1, nextG2, weight):
        if self.direction == '+':
            self.endG1 = nextG1
            self.endG2 = nextG2
            self.weight.append(weight)
        else:
            self.endG1 = nextG1
            self.startG2 = nextG2
            self.weight.append(weight)

    def __len__(self):
        return len(self.weight)

    def __str__(self):
        return 'G1:%s-%s G2:%s-%s %s' % (self.startG1, self.endG1,
                                                self.startG2, self.endG2,
                                                self.direction)


def readDistsAndOrder(dist_file, edgeThreshold):
    res = dict()
    hasMultipleChromosomes = False

    g1_chrom_pos, g2_chrom_pos = {}, {}
    chr1, chr2 = '0', '0'
    with open(dist_file) as data:
        for line in data:
            cells = line.rstrip().split('\t')
            if not res:
                hasMultipleChromosomes = len(cells) == 6

            if hasMultipleChromosomes:
                chr1 = cells[0]
                g1 = int(cells[1])
                chr2 = cells[2]
                g2 = int(cells[3])
                direction = cells[4]
                edgeWeight = float(cells[5])
            else:
                g1 = int(cells[0])
                g2 = int(cells[1])
                direction = cells[2]
                edgeWeight = float(cells[3])

            if edgeWeight < edgeThreshold:
                continue

            if chr1 not in g1_chrom_pos:
                g1_chrom_pos[chr1] = set()
            if chr2 not in g2_chrom_pos:
                g2_chrom_pos[chr2] = set()

            g1_chrom_pos[chr1].add(g1)
            g2_chrom_pos[chr2].add(g2)

            l0 = (chr1, g1)
            l1 = (chr2, g2)

            if l0 not in res:
                res[l0] = dict()
            # append mapping pos in mappedGenome and the weight of the corresponding edge
            res[l0][l1] = (direction == '1' and '+' or '-', edgeWeight)

    # construct genome order
    tel1, g1 = sort_genome(g1_chrom_pos)
    tel2, g2 = sort_genome(g2_chrom_pos)

    # add telomeres
    for t1, t2 in product(tel1, tel2):
        if t1 not in res:
            res[t1] = dict()
        res[t1][t2] = (BothStrands(), 1)

    return hasMultipleChromosomes, g1, g2, res


def sort_genome(chrom_pos):
    g = list()
    telomeres = set()
    for k in sorted(chrom_pos.keys()):
        g.append((k, -1))
        telomeres.add((k, -1))
        g.extend([(k, i) for i in sorted(chrom_pos[k])])
        g.append((k, maxsize))
        telomeres.add((k, maxsize))
    return telomeres, g


def insertIntoRunList(runs, runList, alpha):
    keys = [x.getWeight(alpha) for x in runList]
    for run in runs:
        i = bisect(keys, run.getWeight(alpha))
        keys.insert(i, run.getWeight(alpha))
        runList.insert(i, run)


def checkMatching(g1, g2, g1_runs, g2_runs, runs, dist):
    g1pos = dict(zip(g1, range(len(g1))))
    g2pos = dict(zip(g2, range(len(g2))))

    if len(g1) != len(g2):
        log.error('unequal length: |G1| = %d, |G2| = %d' % (len(g1), len(g2)))
    if len(g1) != len(g1_runs) or len(g2) != len(g2_runs):
        log.error(('Annotation vector length does not match with genome ' +
                   'length: len(G1) = %s, len(g1_runs) = %s, len(G2) = %s, len(' +
                   'g2_runs) = %s') % (len(g1), len(g1_runs), len(g2), len(g2_runs)))

    all_included = set()
    r_counter = 0
    prev_run = None
    c_adj = 0
    for i in range(len(g1)):
        if not g1_runs[i]:
            log.error('Gene %s is not included in any run' % g1[i])
            continue
        if len(g1_runs[i]) > 1:
            log.error('Gene %s is included in more than one run: %s' % (g1[i],
                      ', '.join(map(str, g1_runs[i]))))
            continue

        r = list(g1_runs[i])[0]

        if prev_run != r:
            c_adj += len(r.weight)-1
            if r not in runs:
                log.error('Run %s not included in run list.' % r)
            if r in all_included:
                log.error(('Run %s occurs twice in G1. Current gene position: %s') % (r, g1[i]))
            r_counter += len(r.weight)
            prev_run = r

        all_included.add(r)
        k = i-g1pos[r.startG1]
        if r.direction == '+':
            g2j = g2[g2pos[r.startG2] + k]
        else:
            g2j = g2[g2pos[r.endG2] - k]
        eWgt = dist[g1[i]][g2j][1]

        if r.weight[k] != eWgt:
            log.error('Edge weight of %s-%s differs in run %s, should be %.6f but is %.6f' %
                      (g1[i], g2j, r, eWgt, r.weight[k]))

    missing_runs = all_included.symmetric_difference(runs)
    if missing_runs:
        log.error('Additional runs in runslist that are not part in the matching: %s' %
                  (list(map(str, missing_runs))))

    log.info('Number of adjacencies is %s in matching of size %s.' % (c_adj, len(g1)))

    if r_counter != len(g1):
        log.error(('Sum of run lengths does not equal matching size! Sum ' +
                   'of run lengths: %s, matching size: %s') % (r_counter, len(g1)))

    for j in range(len(g2)):
        if not g2_runs[j]:
            log.error('Gene %s is not included in any run' % g2[j])
        if len(g2_runs[j]) > 1:
            log.error('Gene %s is included in more than one run: %s' % (g2[j],
                      ', '.join(map(str, g2_runs[j]))))
        if g2_runs[j].difference(all_included):
            log.error('G2 differs in runs from G1 on position %s: %s' % (g2[j],
                      ', '.join(map(str, g2_runs[j].difference(all_included)))))

    for r in runs:
        if r.startG1 not in g1pos or r.endG1 not in g1pos or r.startG2 not in \
                g2pos or r.endG2 not in g2pos:
            log.error('Positions of run %s can not be mapped back to the genomes.' % r)
            continue
        if len(g1) <= g1pos[r.startG1] or len(g1) <= g1pos[r.endG1] or \
                len(g2) <= g2pos[r.startG2] or len(g2) <= g2pos[r.endG2]:
            log.error('Positions of run %s exceed borders of the genomes' % r)
            continue
        if g1[g1pos[r.startG1]] != r.startG1 or g2[g2pos[r.startG2]] != \
                r.startG2:
            log.error(('Start of run %s is not coherent with genome ' + \
                    'position on %s (G1) or %s (G2)') %(r, g1[g1pos[r.startG1]],
                        g2[g2pos[r.startG2]]))
        if g1[g1pos[r.endG1]] != r.endG1 or g2[g2pos[r.endG2]] != r.endG2:
            log.error(('End of run %s is not coherent with genome ' + \
                    'position on %s (G1) or %s (G2)') %(r, g1[g1pos[r.endG1]],
                        g2[g2pos[r.endG2]]))
        if g1pos[r.endG1] - g1pos[r.startG1] != g2pos[r.endG2] - \
                g2pos[r.startG2] or g1pos[r.endG1] - g1pos[r.startG1] < 0:
            log.error(('Length of run %s is erroneous: %s (on G1), %s ' + \
                    '(on G2)') %(r, g1pos[r.endG1] - g1pos[r.startG1],
                        g2pos[r.endG2] - g2pos[r.startG2]))
        if len(r.weight) != g1pos[r.endG1] - g1pos[r.startG1] + 1:
            log.error(('Number of weights does not comply with run length. ' + \
                    'Weights: %s, run length: %s, run: %s') %(len(r.weight),
                        g1pos[r.endG1] - g1pos[r.startG1], r))

        g1_chromosomes = set([x[0] for x in g1[g1pos[r.startG1]:g1pos[r.endG1]+1]])
        g2_chromosomes = set([x[0] for x in g2[g2pos[r.startG2]:g2pos[r.endG2]+1]])
        if len(g1_chromosomes) != 1 and len(g2_chromosomes) != 1:
            log.error(('Number of chromosomes on G1 (#chrs: %s) or G2 ' + \
                    '(#chrs: %s) in run %s is not 1 (Meaning that possibly' + \
                    ' the run extends over two or more chromosomes, which ' + \
                    'shouldn\'t be allowed).') %(len(g1_chromosomes),
                        len(g2_chromosomes), r))

    # are all runs merged that can be merged?
    run_ends = dict()
    for r in runs:
        if r.direction == '+':
            run_ends[r.startG1] = (r.direction, r.startG2)
            run_ends[r.endG1] = (r.direction, r.endG2)
        else:
            run_ends[r.startG1] = (r.direction, r.endG2)
            run_ends[r.endG1] = (r.direction, r.startG2)

    for i in range(len(g1)-1):
        g1i = g1[i]
        g1i2 = g1[i+1]
        if g1i in run_ends and g1i2 in run_ends and run_ends[g1i][0] == \
                run_ends[g1i2][0] and g1_runs[i] != g1_runs[i+1]:
            direction = run_ends[g1i][0]
            g2i = run_ends[g1i][1]
            g2i2 = run_ends[g1i2][1]
            if direction == '+' and g2pos[g2i] == g2pos[g2i2]-1:
                log.error('Runs %s and %s could be merged, but are not!' %
                          (list(map(str, g1_runs[i]))[0], list(map(str, g1_runs[i+1]))[0]))
            elif direction == '-' and g2pos[g2i] == g2pos[g2i2]+1:
                log.error('Runs %s and %s could be merged, but are not!' %
                          (list(map(str, g1_runs[i]))[0], list(map(str, g1_runs[i+1]))[0]))


def getAllRuns(g1, g2, d):
    g2pos = dict(zip(g2, range(len(g2))))
    g1_runs = [set() for _ in g1]
    g2_runs = [set() for _ in g2]
    activeRuns = list()
    reportedRuns = list()

    for i in range(len(g1)):
        curPos = g1[i]
        newRunList = list()
        forbiddenRunStarts = list()

        # check if link exists, otherwise terminate all runs
        e = curPos in d
        # iterate over all runs
        for r in activeRuns:
            jEnd = g2pos[r.endG2]
            jStart = g2pos[r.startG2]
            if r.startG1[0] != curPos[0]:
                # run could not be extended
                log.info('terminating run %s, continue with %s.' % (r, curPos))
                reportedRuns.append(r)
                continue
            # extend to the right
            if e and r.direction == '+' and len(g2) > jEnd + 1 and \
                    g2[jEnd+1] in d[curPos] and d[curPos][g2[jEnd+1]][0] == \
                    '+' and g2[jEnd+1][0] == r.endG2[0]:
                g2_gene_r = g2[jEnd+1]
                r.extendRun(curPos, g2_gene_r, d[curPos][g2_gene_r][1])
                newRunList.append(r)
                forbiddenRunStarts.append(('+', g2_gene_r))
                g1_runs[i].add(r)
                g2_runs[jEnd+1].add(r)
                log.debug('Extended run %s to the right' % r)

            # extend to the left
            elif e and r.direction == '-' and jStart > 0 and \
                    g2[jStart-1] in d[curPos] and d[curPos][g2[jStart-1]][0] == \
                    '-' and g2[jStart-1][0] == r.startG2[0]:
                g2_gene_l = g2[jStart-1]
                r.extendRun(curPos, g2_gene_l, d[curPos][g2_gene_l][1])
                newRunList.append(r)
                g1_runs[i].add(r)
                g2_runs[jStart-1].add(r)
                forbiddenRunStarts.append(('-', g2_gene_l))
                log.debug('Extended run %s to the left' % r)
            else:
                # run could not be extended
                log.info(('Terminate and report run %s, because %s has '
                          'no further consecutive edge.') % (r, curPos))
                reportedRuns.append(r)

        # if no edge exists, nothing has to be done...
        if e:
            for (g2_gene, (direction, weight)) in list(d[curPos].items()):
                if (direction, g2_gene) not in forbiddenRunStarts:
                    j = g2pos[g2_gene]
                    if isinstance(direction, BothStrands):
                        r = Run(curPos, g2_gene, weight, '+')
                        newRunList.append(r)
                        g1_runs[i].add(r)
                        g2_runs[j].add(r)
                        log.debug(('Start new (%s) run %s') % (direction, r))
                        r = Run(curPos, g2_gene, weight, '-')
                        newRunList.append(r)
                        g1_runs[i].add(r)
                        g2_runs[j].add(r)
                        log.debug(('Start new (%s) run %s') % (direction, r))
                    else:
                        r = Run(curPos, g2_gene, weight, direction)
                        newRunList.append(r)
                        g1_runs[i].add(r)
                        g2_runs[j].add(r)
                        log.debug(('Start new (%s) run %s') % (direction, r))
        activeRuns = newRunList
    reportedRuns.extend(activeRuns)
    return (g1_runs, g2_runs, reportedRuns)


def replaceByNew(g1_runs, g2_runs, i, j, r_old, r_new):
    while r_old in g1_runs[i]:
        g1_runs[i].remove(r_old)
        g1_runs[i].add(r_new)
        g2_runs[j].remove(r_old)
        g2_runs[j].add(r_new)
        i += 1
        j += 1
        if len(g1_runs) <= i or len(g2_runs) <= j:
            break


def doMatching(g1, g2, g1_runs, g2_runs, m, runList, alpha):
    g1pos = dict(zip(g1, range(len(g1))))
    g2pos = dict(zip(g2, range(len(g2))))
    newRuns = set()

    for k in range(g1pos[m.endG1] - g1pos[m.startG1] + 1):
        i = g1pos[m.startG1] + k
        j = g2pos[m.startG2] + k

        for r in set(g1_runs[i]):
            if r == m:
                continue
            g1_runs[i].remove(r)

            if r in runList:
                runList.remove(r)

            if g1pos[r.startG1] < i:
                overlap = g1pos[r.endG1] - i
                log.info(('Run %s overlaps with selected run %s by %s ' +
                          'at position G1:%s.') % (r, m, overlap+1, g1[i]))
                r_new = deepcopy(r)
                r_new.endG1 = g1[i-1]
                if r.direction == '+':
                    # check weight
                    r_new.endG2 = g2[g2pos[r.endG2] - overlap - 1]
                    r_new.weight = r.weight[:-overlap-1]
                    r.weight = r.weight[-overlap-1:]
                    r.startG2 = g2[g2pos[r.endG2]-overlap]
                    g2_runs[g2pos[r.startG2]].remove(r)
                else:
                    r_new.startG2 = g2[g2pos[r.startG2] + overlap + 1]
                    r_new.weight = r.weight[:-overlap-1]
                    r.weight = r.weight[-overlap-1:]
                    r.endG2 = g2[g2pos[r.startG2] + overlap]
                    g2_runs[g2pos[r.endG2]].remove(r)
                r.startG1 = g1[i]
                log.info('Divided overlapping run in %s and %s' % (r_new, r))
                # do you see that r.startG2 is already at the right position?
                replaceByNew(g1_runs, g2_runs, g1pos[r_new.startG1],
                             g2pos[r_new.startG2], r, r_new)
                newRuns.add(r_new)

            elif g1pos[r.startG1] == i:
                if r.direction == '+':
                    g2_runs[g2pos[r.startG2]].remove(r)
                else:
                    g2_runs[g2pos[r.endG2]].remove(r)
            if len(g1) > i+1 and i < g1pos[r.endG1]:
                # run start cannot be larger than i
                log.info(('Run %s interfers with current run %s at ' +
                          'position G1:%s. Shifting.') % (r, m, g1[i]))
                r.startG1 = g1[i+1]
                del r.weight[0]
                if r.direction == '+':
                    r.startG2 = g2[g2pos[r.startG2]+1]
                else:
                    r.endG2 = g2[g2pos[r.endG2]-1]

                log.info('Shifted run is now located at %s' % r)
                newRuns.add(r)
            elif r in newRuns:
                newRuns.remove(r)

        for r in set(g2_runs[j]):
            if r == m:
                continue
            g2_runs[j].remove(r)

            if r in runList:
                runList.remove(r)

            if g2pos[r.startG2] < j:
                overlap = g2pos[r.endG2] - j
                log.info(('Run %s overlaps with selected run %s by %s ' +
                          'at position G2:%s.') % (r, m, overlap+1, g2[j]))
                r_new = deepcopy(r)
                r_new.endG2 = g2[j-1]
                if r.direction == '+':
                    r_new.endG1 = g1[g1pos[r.endG1]-overlap-1]
                    r_new.weight = r.weight[:-overlap-1]
                    r.weight = r.weight[-overlap-1:]
                    r.startG1 = g1[g1pos[r.endG1]-overlap]
                    g1_runs[g1pos[r.startG1]].remove(r)
                else:
                    r_new.startG1 = g1[g1pos[r.startG1]+overlap+1]
                    r_new.weight = r.weight[overlap+1:]
                    r.weight = r.weight[:overlap+1]
                    r.endG1 = g1[g1pos[r.startG1]+overlap]
                    g1_runs[g1pos[r.endG1]].remove(r)
                r.startG2 = g2[j]
                log.info('Divided overlapping run in %s and %s' % (r_new, r))

                replaceByNew(g1_runs, g2_runs, g1pos[r_new.startG1],
                             g2pos[r_new.startG2], r, r_new)
                newRuns.add(r_new)

            elif g2pos[r.startG2] == j:
                if r.direction == '+':
                    g1_runs[g1pos[r.startG1]].remove(r)
                else:
                    g1_runs[g1pos[r.endG1]].remove(r)

            if len(g2) > j+1 and j < g2pos[r.endG2]:
                # run start cannot be larger than j
                log.info(('Run %s interfers with current run %s at ' +
                          'position G2:%s. Shifting.') % (r, m, g2[j]))
                r.startG2 = g2[j+1]
                if r.direction == '+':
                    r.startG1 = g1[g1pos[r.startG1]+1]
                    del r.weight[0]
                else:
                    r.endG1 = g1[g1pos[r.endG1]-1]
                    del r.weight[-1]
                log.info('Shifted run is now located at %s' % r)
                newRuns.add(r)
            elif r in newRuns:
                newRuns.remove(r)
    insertIntoRunList(newRuns, runList, alpha)


def mergeRuns(mod_g1, g1, g2, g1_runs, g2_runs, runList, alreadyMatched, alpha):
    g1pos = dict(zip(g1, range(len(g1))))
    g2pos = dict(zip(g2, range(len(g2))))

    newRuns = set()
    mod_g1 = list(mod_g1)
    for x in range(len(mod_g1)):
        g1i = mod_g1[x]
        i = g1pos[g1i]
        if len(g1) < i+2:
            continue

        # To understand this piece of code, one observation is important:
        # If r1 or r2 is already matched, then there exists only one combination
        # of possible merges. If r1 and r2 are both unmatched, several merges
        # are possible and all should be done.
        # After each merge between a matched and unmatched run, the newly
        # merged run must be completely matched, before further modification
        # points (mod_g1) can be processed.

        for r1, r2 in product(sorted(g1_runs[i].difference(g1_runs[i+1]),
                                     key=lambda x: x.getWeight(alpha), reverse=True),
                              sorted(g1_runs[i+1].difference(g1_runs[i]),
                                     key=lambda x: x.getWeight(alpha), reverse=True)):
            if r1.endG1 == g1[i] and r2.startG1 == g1[i+1] and \
                    r1.direction == r2.direction and \
                    r1.endG1[0] == r2.startG1[0] and \
                    r1.endG2[0] == r2.startG2[0] and \
                    ((r1.direction == '+' and
                      g2pos[r1.endG2] == g2pos[r2.startG2] - 1) or
                     (r1.direction == '-' and
                      g2pos[r2.endG2] == g2pos[r1.startG2] - 1)):

                log.info('Merge runs %s and %s.' % (r1, r2))
                if r1 in runList:
                    runList.remove(r1)
                if r2 in runList:
                    runList.remove(r2)
                if r1 in newRuns:
                    # ye-ah, this can happen too :/
                    newRuns.remove(r1)

                r2.startG1 = r1.startG1
                r2.weight = r1.weight + r2.weight
                if r1.direction == '+':
                    r2.startG2 = r1.startG2
                else:
                    r2.endG2 = r1.endG2
                log.info('Merged run is %s' % r2)
                replaceByNew(g1_runs, g2_runs, g1pos[r1.startG1],
                             g2pos[r1.startG2], r1, r2)
                if (r2 in alreadyMatched) ^ (r1 in alreadyMatched):
                    if r1 in alreadyMatched:
                        alreadyMatched.remove(r1)
                    # redo matching in case r1 xor r2 were not in matching before
                    insertIntoRunList(newRuns, runList, alpha)
                    return r2, set(mod_g1[x+1:])
                if r2 in alreadyMatched:
                    # actually, both are already matched
                    alreadyMatched.add(r2)
                    alreadyMatched.remove(r1)
                else:
                    # none is matched
                    newRuns.add(r2)

    insertIntoRunList(newRuns, runList, alpha)
    return None, []


def removeSingleGenes(genome, genome_runs):
    del_res = set()
    mod_res = set()
    i = 0
    while i < len(genome):
        if not genome_runs[i]:
            del_res.add(genome[i])
            mod_res.add(genome[i-1])
            del genome[i]
            del genome_runs[i]
        else:
            i += 1
    return del_res, mod_res


def findRandomRunSequence(g1, g2, dists, topXperCent, alpha):
    g2dists = dict()
    for g1i, x in list(dists.items()):
        for g2j, d in list(x.items()):
            if g2j not in g2dists:
                g2dists[g2j] = dict()
            g2dists[g2j][g1i] = d

    # copy g1, g2 and dists map, because we'll modify it. Also remove all genes
    # that do not contain edges.
    g1 = [x for x in g1 if x in dists and len(dists[x])]
    g2 = [x for x in g2 if x in g2dists and len(g2dists[x])]

    g1pos = dict(zip(g1, range(len(g1))))

    g1_runs, g2_runs, runs = getAllRuns(g1, g2, dists)
    log.info('Found %s runs.' % len(runs))
    # sort
    runList = sorted(runs, key=lambda x: x.getWeight(alpha))

    res = set()
    while runList:
        noOfAdjacencies = len([x for x in runList if x.getWeight(alpha) and
            x.getWeight(alpha) or 0])
        if noOfAdjacencies:
            randPos = randint(1, ceil(noOfAdjacencies * topXperCent))
        else:
            randPos = randint(1, ceil(len(runList) * topXperCent))
        log.info('From %s, select randomly among top %s run %s' %
                 (len(runList), int(ceil((noOfAdjacencies or len(runList)) * topXperCent)), runList[-randPos]))
        mx = runList.pop(-randPos)
        mod_g1 = set()
        while mx:
            res.add(mx)
            # update run list
            doMatching(g1, g2, g1_runs, g2_runs, mx, runList, alpha)
            del_g1, new_mod_g1 = removeSingleGenes(g1, g1_runs)
            if del_g1:
                log.info('Zombie genes removed from G1: %s' % ', '.join(map(str, del_g1)))
                # it can happen that a gene in mod_g1 has already been deleted
                # before being processed. This happens if there is a merge between
                # a matched and unmatched run. Then some genes remain unprocessed
                # while the merged run is re-matched. In this process, new genes
                # can be deleted. If one of the genes happens to be in mod_g1, it
                # should be deleted.
                for g in del_g1.intersection(mod_g1):
                    mod_g1.remove(g)

                g1pos = dict(zip(g1, range(len(g1))))
            # add new modification points
            mod_g1.update(new_mod_g1)

            del_g2, mod_g2 = removeSingleGenes(g2, g2_runs)
            if del_g2:
                log.info('Zombie genes removed from G2: %s' % ', '.join(map(str, del_g2)))
                for g2j in mod_g2:
                    for g1i, (d, _) in list(g2dists[g2j].items()):
                        if g1i in g1:
                            if d == '+':
                                mod_g1.add(g1i)
                            # d == DIRECTION_BOTH_STRANDS? both neighbors have to be added...
                            if d == '-':
                                mod_g1.add(g1[g1pos[g1i]-1])
            # merge runs
            mx, mod_g1 = mergeRuns(mod_g1, g1, g2, g1_runs, g2_runs,
                                   runList, res, alpha)

    if res:
        log.info('Matching finished. Longest run size is %s.' % (max(list(map(len, res)))))
    else:
        log.info('Matching finished, but no runs found. Empty input?')

    return (g1, g2, g1_runs, g2_runs, res)


def repeatMatching(g1, g2, g1_mod, g2_mod, g1_runs, g2_runs, dists, repMatching,
                   minCsSize, topXperCent, alpha):

    g1_mod_res = g1_mod
    g2_mod_res = g2_mod
    g1_runs_res = g1_runs
    g2_runs_res = g2_runs
    selectedRuns_res = list()

    g1pos = dict(zip(g1_mod, range(len(g1_mod))))
    g2pos = dict(zip(g2_mod, range(len(g2_mod))))

    noReps = repMatching

    while repMatching:
        for i in range(len(g1_runs)):
            run_set = g1_runs[i]
            if len(run_set) != 1:
                log.error('Expected run, set length of 1, but was told different: %s.' %
                          ', '.join(map(str, run_set)))
            run = next(run_set.__iter__())

            g1i = g1_mod[i]

            j = i-g1pos[run.startG1]
            if run.direction == '+':
                g2j = g2_mod[g2pos[run.startG2] + j]
            else:
                g2j = g2_mod[g2pos[run.endG2] - j]
            del dists[g1i][g2j]

            if not dists[g1i]:
                del dists[g1i]

        if not dists:
            log.info('Removed all edges in the input graph. Stopping iteration %s.' %
                     (noReps-repMatching+2))
            break

        g1_mod, g2_mod, g1_runs, g2_runs, selectedRuns = findRandomRunSequence(g1, g2, dists, topXperCent, alpha)
        checkMatching(g1_mod, g2_mod, g1_runs, g2_runs, selectedRuns, dists)

        log.info('Obtained %s adjacencies in matching of size %s from iteration %s.' %
                 (len(g1_mod) - len(selectedRuns), len(g1_mod), noReps-repMatching+2))

        # remove runs that fall below min length of minCsSize
        ff = lambda x: len(next(x.__iter__())) >= minCsSize
        g1_mod = [g1_mod[i] for i in range(len(g1_mod)) if ff(g1_runs[i])]
        g2_mod = [g2_mod[i] for i in range(len(g2_mod)) if ff(g2_runs[i])]
        g1_runs = list(filter(ff, g1_runs))
        g2_runs = list(filter(ff, g2_runs))
        selectedRuns = set([s for s in selectedRuns if len(s) >= minCsSize])

        # stop if no runs were found matching the criteria
        if not len(selectedRuns):
            log.info('No feasible runs found in matching round %s. Stopping iteration.' %
                     (noReps-repMatching+2))
            break

        log.info('%s feasible runs retained.' % len(selectedRuns))

        # reconciliate with result data
        g2pos = dict(zip(g2_mod, range(len(g2_mod))))
        g1pos = dict(zip(g1_mod, range(len(g1_mod))))
        g2pos_res = dict(zip(g2_mod_res, range(len(g2_mod_res))))
        g1pos_res = dict(zip(g1_mod_res, range(len(g1_mod_res))))

        g1_mod_new = sorted(set(g1_mod_res + g1_mod), key=lambda x: (x[0], x[1]))
        g2_mod_new = sorted(set(g2_mod_res + g2_mod), key=lambda x: (x[0], x[1]))
        g1_runs_new = list()
        g2_runs_new = list()

        for g1i in g1_mod_new:
            x = set()
            if g1i in g1pos_res:
                x.update(g1_runs_res[g1pos_res[g1i]])
            if g1i in g1pos:
                x.update(g1_runs[g1pos[g1i]])
            g1_runs_new.append(x)

        for g2j in g2_mod_new:
            x = set()
            if g2j in g2pos_res:
                x.update(g2_runs_res[g2pos_res[g2j]])
            if g2j in g2pos:
                x.update(g2_runs[g2pos[g2j]])
            g2_runs_new.append(x)

        g1_mod_res = g1_mod_new
        g2_mod_res = g2_mod_new
        g1_runs_res = g1_runs_new
        g2_runs_res = g2_runs_new

        selectedRuns_res.extend(selectedRuns)
        repMatching -= 1

    return (g1_mod_res, g2_mod_res, g1_runs_res, g2_runs_res, selectedRuns_res)


def printMatching(g1, g2, g1_runs, hasMultipleChromosomes, out):

    with open(out, 'w') as of:

        if hasMultipleChromosomes:
            of.write('Chr(G1)\tG1\tChr(G2)\tG2\tdirection\tedge weight\n')
        else:
            of.write('G1\tG2\tdirection\tedge weight\n')

        g2pos = dict(zip(g2, range(len(g2))))
        # g1pos = dict(zip(g1, range(len(g1))))

        cur_index = dict()
        for i in range(len(g1_runs)):
            run_set = g1_runs[i]
            for run in run_set:
                g1i = g1[i]
                j = 0
                if run in cur_index:
                    j = cur_index[run]
                if run.direction == '+':
                    g2j = g2[g2pos[run.startG2] + j]
                else:
                    g2j = g2[g2pos[run.endG2] - j]

                direction = run.direction == '+' and '1' or '-1'

                g1i1 = g1i[1] == -1 and 'TELOMERE_START' or g1i[1]
                g1i1 = g1i[1] == maxsize and 'TELOMERE_END' or g1i1
                g2j1 = g2j[1] == -1 and 'TELOMERE_START' or g2j[1]
                g2j1 = g2j[1] == maxsize and 'TELOMERE_END' or g2j1

                if hasMultipleChromosomes:
                    of.write('%s\t%s\t%s\t%s\t%s\t%s\n' %
                             (g1i[0], g1i1, g2j[0], g2j1, direction, run.weight[j]))
                else:
                    of.write('%s\t%s\t%s\t%s\n' % (g1i1, g2j1, direction, run.weight[j]))

                cur_index[run] = j+1


if __name__ == '__main__':
    from argparse import ArgumentParser as AP
    cli = AP(description='')
    cli.add_argument('-R', '--repeat-matching', type=int, metavar='N', default=0,
                     help='match N repetitions (default: 0)')
    cli.add_argument('-M', '--min-cs-size', type=int, metavar='N', default=1,
                     help='minimal cs size (default: 1)')  # should be conditional on -R
    cli.add_argument('-g', '--greedy', type=float, metavar='F', default=10.**-7)
    cli.add_argument('-e', '--edge_weight_threshold', type=float, default=0.0)
    cli.add_argument('-a', '--alpha', type=float, metavar='F', default=0.5)
    cli.add_argument('dist_file')
    args = cli.parse_args()
    repMatching = args.repeat_matching
    if repMatching > 0:
        repMatching -= 1

    log.basicConfig(filename=args.dist_file.rsplit('.')[0]+'.log', level=log.INFO,
                    format="%(levelname)s\t%(asctime)s\t++ %(message)s")

    multiChrom, g1, g2, dists = readDistsAndOrder(args.dist_file, args.edge_weight_threshold)
    g1_mod, g2_mod, g1_runs, g2_runs, selectedRuns = findRandomRunSequence(g1,
            g2, dists, args.greedy, args.alpha)
    checkMatching(g1_mod, g2_mod, g1_runs, g2_runs, selectedRuns, dists)

    # calculate number of breakpoints only from result of the first matching
    bkp = len(selectedRuns) - 1

    g1_mod, g2_mod, g1_runs, g2_runs, selectedRuns_new = repeatMatching(g1, g2,
            g1_mod, g2_mod, g1_runs, g2_runs, dists, repMatching,
            args.min_cs_size, args.greedy, args.alpha)

    selectedRuns.update(selectedRuns_new)

    # sum of weights of adjacencies
    wAdj = sum([r.getWeight(1) for r in selectedRuns])
    # sum of weights of all edges of the matching
    wEdg = sum([sum([x**2 for x in r.weight]) for r in selectedRuns])

    edg = sum(map(len, selectedRuns))

    matchfile = args.dist_file+'.matching'
    printMatching(g1_mod, g2_mod, g1_runs, multiChrom, matchfile)

    #
    # print matching scores
    #

    log.info(('FFAdj-MCS finished. Breakpoint distance between G1 and G2' +
              ' is %s with #edg = %s, adj(M) = %.3f and edg(M) = %.3f') %
             (bkp, edg, wAdj, wEdg))

    print('#bkp\t#edg\tadj\tedg')
    print('%s\t%s\t%.6f\t%.6f' % (bkp, edg, wAdj, wEdg))

