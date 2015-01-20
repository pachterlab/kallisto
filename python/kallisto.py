#! /usr/bin/python

import argparse
from collections import defaultdict, Counter
from Bio import SeqIO
import itertools
import time
import re
import os
import gzip

try:
    import cPickle as pickle
except:
    import pickle


def timed_print(string):
    print '[{0}] {1}'.format(time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), string)


# equiv_class = id, set of kmers, dict of counts by transcript
class EquivalenceClass:
    def __init__(self, ec_id, kmer_count = 0, counts = None ):
        self.id = ec_id
        self.kmer_count = kmer_count
        self.counts = counts or defaultdict(int)
    def __hash__(self):
        return hash(tuple(self.counts.items()))
    def min(self, arg):
        to_ret = EquivalenceClass(-1)
        is_self = True
        for key in self.counts:
            if key in arg.counts:
                if arg.counts[key] < self.counts[key]:
                    to_ret.counts[key] = arg.counts[key]
                    is_self = False
                else:
                    to_ret.counts[key] = self.counts[key]
            else:
                is_self = False
        if is_self:
            return self
        else:
            return to_ret


class memoized_min:
    def __init__(self, equiv_classes):
        self.equiv_classes = equiv_classes
        # set up ec_hash for recognizing equiv_classes after doing a min operation
        self.ec_hash = dict()
        for ec, ec_id in itertools.izip(self.equiv_classes.itervalues(), self.equiv_classes.iterkeys()):
#            self.ec_hash[tuple(ec.counts)] = ec_id
            self.ec_hash[tuple(ec.counts.iteritems())] = ec_id
        self.max_id = max([ec.id for ec in self.equiv_classes.itervalues()])
        self.memo = {}
    def __call__(self, *args):
        # with args[0] <= args[1], we can make use of the fact that our ecs are sorted to say that min(args[0], args[1]) <= args[0]
        if args[0] > args[1]:
            args = (args[1], args[0])
        # have we memoized these arguments? if so, great.
        if args in self.memo:
            return self.memo[args]
        else:
            # otherwise, actually compute min (again, makes use of the sorted order of the ecs. min is args[0] if anything)
            m = self.equiv_classes[args[0]].min(self.equiv_classes[args[1]])
            # if it's just the first guy then we're done
            if m.id == args[0]:
                m_id = m.id
            else:
                # see whether this is already an existing ec
                tup = tuple(m.counts.iteritems())
#                tup = tuple(m.counts)
                if tup in self.ec_hash:
                    m_id = self.ec_hash[tup]
                else:
                    # if not, add it (including to ec_hash)
                    self.max_id += 1
                    m_id = self.max_id
                    self.ec_hash[tup] = m_id
                    self.equiv_classes[m_id] = m
            m.id = m_id
            self.memo[args] = m_id
            return m_id


def kmers(seq, k, reverse_comp=False):
    for i in xrange(len(seq) - k + 1):
        yield seq[i:i+k]


def make_equiv_classes(file_name, k):
    kmer_ec = defaultdict(int) # holds the equiv class (id) for every k-mer. if we've never seen a k-mer before, then it'll get a new equiv class which is basically what equiv_classes[0] is
    equiv_classes = dict()
    equiv_classes[0] = EquivalenceClass(0)
    max_id = 0
    for trans_num, trans in enumerate(SeqIO.parse(file_name, "fasta")):
        if trans_num % 10000 == 0:
            timed_print(' - processing transcript {0}'.format(trans_num))
        for reverse in [0,1]:
            trans_dict = defaultdict(int)
            if reverse:
                for kmer in kmers(str(trans.reverse_complement().seq), k):
                    trans_dict[kmer] += 1
            else:
                for kmer in kmers(str(trans.seq), k):
                    trans_dict[kmer] += 1
            kmers_by_ec_and_count = defaultdict(list)
            for kmer, count in trans_dict.iteritems():
                kmers_by_ec_and_count[(kmer_ec[kmer],count)].append(kmer)
            for (ec_id,count) in sorted(kmers_by_ec_and_count.keys()):
                curr_kmers = kmers_by_ec_and_count[(ec_id,count)]
                max_id += 1
                equiv_classes[max_id] = EquivalenceClass(max_id, len(curr_kmers), dict(equiv_classes[ec_id].counts))
                equiv_classes[max_id].counts[2*trans_num + reverse] = count
                # remove these kmers from old equiv class and delete class if necessary
                if ec_id != 0:
                    equiv_classes[ec_id].kmer_count -= len(curr_kmers)
                    if equiv_classes[ec_id].kmer_count < 0:
                        raise Exception("EC with negative members?")
                    if equiv_classes[ec_id].kmer_count == 0:
                        del(equiv_classes[ec_id])
                # update kmer_ec
                for kmer in curr_kmers:
                    kmer_ec[kmer] = max_id

    return kmer_ec, equiv_classes


def process_reads(equiv_classes, kmer_ec, read_file, k, output_dir = None, thresh_frac = 0, dont_min = False, memo_min = None):
    if not memo_min:
        memo_min = memoized_min(equiv_classes)
    if output_dir:
        f = open('{0}/reads.kal'.format(output_dir), 'w')
    read_counts = defaultdict(int)

    if re.search('gz', read_file):
        handle1 = gzip.open(read_file.format('1'), 'r')
        handle2 = gzip.open(read_file.format('2'), 'r')
    else:
        handle1 = open(read_file.format('1'), 'r')
        handle2 = open(read_file.format('2'), 'r')

    if re.search('fastq', read_file) or re.search('fq', read_file):
        file1 = SeqIO.parse(handle1, "fastq")
        file2 = SeqIO.parse(handle2, "fastq")
    else:
        file1 = SeqIO.parse(handle1, "fasta")
        file2 = SeqIO.parse(handle2, "fasta")

    for read_num, (read1, read2) in enumerate(itertools.izip(file1, file2)):
        read2 = read2.reverse_complement()
        if dont_min:
            if output_dir:
                f.write('{0}\t'.format(read1.id))
            ecs = Counter([kmer_ec[kmer] for kmer in list(kmers(str(read1.seq),k)) + list(kmers(str(read2.seq),k))]).most_common()
            for ec, count in ecs:
                read_counts[ec] += count
            if output_dir:
                f.write(','.join(['{0} {1}'.format(ec, count) for ec, count in ecs]) + '\n')
        else:
            ecs = Counter([kmer_ec[kmer] for kmer in list(kmers(str(read1.seq),k)) + list(kmers(str(read2.seq),k)) if kmer in kmer_ec]).most_common()

            curr_ec = 0
            total_count = 0
            min_count = 0
            for ec_id, count in ecs:
                total_count += count
                if curr_ec == 0:
                    if ec_id != 0:
                        curr_ec = ec_id
                        min_count = count
                else:
                    new_ec = memo_min(curr_ec, ec_id)
                    if new_ec != 0:
                        curr_ec = new_ec
                        min_count += count
            if min_count < thresh_frac*total_count:
                curr_ec = 0

            if output_dir:
                f.write('{0}\t{1}\n'.format(read1.id, curr_ec))
            read_counts[curr_ec] += 1
        if read_num % 1000000 == 0:
            timed_print(' - processing read #{0}'.format(read_num))
    if output_dir:
        f.close()

    return read_counts, memo_min


def calc_eff_lens(transcriptome, read_len, frag_len_dist):
    lens = set([len(s) for s in transcriptome])
    eff_len = dict([(l, sum([frag_len_dist[fl]*(l - fl + 1) for fl in frag_len_dist if fl < l + 1])) for l in lens])

    return [eff_len[len(s)] for s in transcriptome]


def calc_trans_ec_weights(equiv_classes, kmer_ec, eff_lens, transcriptome, k = 0, read_len = 0, frag_len_dist = None, do_stat = False, memo_min = None):
    trans_ec_weight = [defaultdict(float) for i in range(len(transcriptome))]

    if not do_stat:
        for ec_id in equiv_classes.iterkeys():
            for trans, trans_count in equiv_classes[ec_id].counts.iteritems():
                if eff_lens[trans/2] > 0:
                    trans_ec_weight[trans/2][ec_id] = float(trans_count)/eff_lens[trans/2]
    else:
        if not memo_min:
            memo_min = memoized_min(equiv_classes)
        for trans_num, trans in enumerate(transcriptome):
            if trans_num % 10000 == 0:
                timed_print(' - processing transcript {0}'.format(trans_num))
            for rev in [0,1]:
                # get the ec for each k-mer in the transcript
                if rev:
                    ecs = [kmer_ec[kmer] for kmer in kmers(str(trans.reverse_complement().seq), k)]
                else:
                    ecs = [kmer_ec[kmer] for kmer in kmers(str(trans.seq), k)]
                # compute the mec for a read of length k starting at position i for all possible i
                mecs = []
                for i in xrange(len(trans) - read_len + 1):
                    curr_ec = ecs[i]
                    for ec in ecs[(i + 1):(i + read_len - k + 1)]:
                        curr_ec = memo_min(curr_ec, ec)
                    mecs.append(curr_ec)
                # record the total probability of seeing every possible mec
                counts = Counter([(l, mecs[p], mecs[p + l]) for l in frag_len_dist for p in xrange(len(mecs) - l)])
                for (l, mec1, mec2), count in counts.most_common():
                    trans_ec_weight[trans_num][memo_min(mec1,mec2)] += 1/2.0*count*frag_len_dist[l]
            for ec in trans_ec_weight[trans_num]:
                trans_ec_weight[trans_num][ec] /= eff_lens[trans_num]

    return trans_ec_weight, memo_min
    

def calculate_alpha(equiv_classes, read_counts, trans_ec_weights, num_iter, output_dir = None, alpha = None):
    # convert read counts to floats
    read_counts = dict(itertools.izip(read_counts.iterkeys(),map(float,read_counts.itervalues())))
    # some basic info
    num_trans = len(trans_ec_weights)
    total_mass = sum(read_counts.itervalues()) - read_counts[0]

    # initialize alphas if necessary
    if not alpha:
        alpha = [1.0/num_trans for x in range(num_trans)]

    if output_dir:
        bad_ecs = set()
    for t in range(num_iter):
        timed_print(' - iteration {0}'.format(t + 1))
        # for current alphas, calculate the denominator for every ec
        denom = defaultdict(float)
        new_alpha = num_trans*[0.0]
        for ec_id in read_counts.iterkeys():
            for trans, trans_count in equiv_classes[ec_id].counts.iteritems():
                denom[ec_id] += alpha[trans/2]*trans_ec_weights[trans/2][ec_id]
        # for every transcript, calculate total mass
        for trans in range(num_trans):
            for ec_id, ec_weight in trans_ec_weights[trans].iteritems():
                if ec_id in read_counts:
                    if denom[ec_id] == 0:
                        bad_ecs.add(ec_id)
                    else:
                        new_alpha[trans] += read_counts[ec_id]*(ec_weight*alpha[trans]/denom[ec_id])
        alpha = new_alpha

    if output_dir:
        f = open(output_dir + "/bad_ecs", "w")
        for ec_id in bad_ecs:
            f.write('ec_id: {0}\n'.format(ec_id))
        f.close()

        f = open('{0}/ec_splits.kal'.format(output_dir), 'w')
        of = open('{0}/ec_vecs.kal'.format(output_dir), 'w')
        denom = defaultdict(float)
        for ec_id in read_counts.iterkeys():
            for trans, trans_count in equiv_classes[ec_id].counts.iteritems():
                denom[ec_id] += alpha[trans/2]*trans_ec_weights[trans/2][ec_id]
        for ec_id in read_counts.iterkeys():
            first = True
            f.write('{0}\t'.format(ec_id))
            for trans, trans_count in equiv_classes[ec_id].counts.iteritems():
                if denom[ec_id] == 0:
                    continue
                if not first:
                    f.write(' || ')
                f.write('{0}, {1}'.format(trans, alpha[trans/2]*trans_ec_weights[trans/2][ec_id]/denom[ec_id]))
                first = False
            f.write('\n')

            first = True
            of.write('{0}\t'.format(ec_id))
            for trans, trans_count in equiv_classes[ec_id].counts.iteritems():
                if not first:
                    of.write(' || ')
                of.write('{0}, {1}'.format(trans, trans_count))
                first = False
            of.write('\n')
        f.close()
        of.close()

    return alpha


def read_transcriptome(file_name):
    return [trans for trans in SeqIO.parse(file_name, "fasta")]


def read_fld(file_name, thresh=0):
    with open(file_name) as f:
        fld = defaultdict(float, [(int(line.split()[0]), float(line.split()[1])) for line in f if float(line.split()[1]) > thresh])
    total = sum(list(fld.itervalues()))
    for l in fld:
        fld[l] /= total
    return fld


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Do Kallisto',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('trans_fasta', type=str, help='fasta file with transcript sequences')
    parser.add_argument('read_fasta', type=str, help='fasta file with read sequences')
    parser.add_argument('output_dir', type=str, help='directory to output things to')
    parser.add_argument('-i', '--index', type=str, help='file with fragment length distribution')
    parser.add_argument('-k', type=int, default=20, help='k-mer length')
    parser.add_argument('-n', '--num_iter', type=int, default=500, help='number of iterations for which to run the EM algorithm.')
    parser.add_argument('-thr', '--thresh_frac', type=float, default=0.0, help='minimum fraction of k-mers contributing to pseudoalignment')
    parser.add_argument('-fld', '--fld_file', type=str, help='file with fragment length distribution')
    parser.add_argument('-dm', '--dont_min', action="store_true", help='don\'t perform minimum operation on k-mer ECs (i.e. basically just do Sailfish)')
    parser.add_argument('-bs', '--bootstrap_samples', type=int, default=0, help='number of bootstrap samples to take')
    parser.add_argument('-bi', '--bootstrap_iter', type=int, default=100, help='number of EM iterations to do for each bootstrap sample')
    parser.add_argument('-st', '--do_stat', action="store_true", help='do things statistically')
    parser.add_argument('-so', '--same_orient', action="store_true", help='paired-end reads have the same orientation')
    parser.add_argument('-tew', '--tew_file', type=str, help='file giving EC weights for each element in the transcriptome')
    args = parser.parse_args()

    timed_print('Reading transcriptome')
    transcriptome = read_transcriptome(args.trans_fasta)

    if args.fld_file:
        frag_len_dist = read_fld(args.fld_file, thresh=10**-12)

        timed_print('Calculating effective transcript lengths')
        eff_lens = calc_eff_lens(transcriptome, 75, frag_len_dist)
        for i in range(len(eff_lens)):
            if eff_lens[i] < 0.1:
                eff_lens[i] = len(transcriptome[i].seq)
    else:
        frag_len_dist = None
        eff_lens = [max(1, len(s) - args.k + 1) for s in transcriptome]

    if not os.path.exists(args.output_dir):
        timed_print('Output directory created')
        os.mkdir(args.output_dir)

    timed_print('Creating equivalence classes')
    kmer_ec, equiv_classes = make_equiv_classes(args.trans_fasta, args.k)

    memo_min = None

    if args.tew_file:
        trans_ec_weights = []
        timed_print('Reading transcript-EC weight file')
        with open(args.tew_file) as f:
            for line in f:
                trans_ec_weights.append(defaultdict(float, [(int(s.split(',')[0]), float(s.split(',')[1])) for s in line.split("\t")[1].rstrip().split(" || ") if s]))
    else:
        if args.do_stat:
            timed_print('Calculating statistical transcript ec weights')
            trans_ec_weights, memo_min = calc_trans_ec_weights(equiv_classes, kmer_ec, eff_lens, transcriptome, args.k, 75, frag_len_dist, args.do_stat, memo_min)

    timed_print('Processing read file')
    read_counts, memo_min = process_reads(equiv_classes, kmer_ec, args.read_fasta, args.k, args.output_dir, args.thresh_frac, dont_min = args.dont_min)

    if not args.do_stat:
        timed_print('Calculating heuristic transcript ec weights')
        trans_ec_weights, memo_min = calc_trans_ec_weights(equiv_classes, kmer_ec, eff_lens, transcriptome, args.k, 75, frag_len_dist, args.do_stat, memo_min)

    with open('{0}/trans_ec_weights.kal'.format(args.output_dir), 'w') as f:
        for trans_num, ec_weights in enumerate(trans_ec_weights):
            first = True
            f.write('{0}\t'.format(trans_num))
            for ec_id, ec_weight in ec_weights.iteritems():
                if not first:
                    f.write(' || ')
                f.write('{0}, {1}'.format(ec_id, ec_weight))
                first = False
            f.write('\n')

    read_ecs = read_counts.keys()
    total_reads = sum(read_counts.values())
    read_probs = [float(n)/total_reads for n in read_counts.values()]
    read_counts = [read_counts]
    for i in range(args.bootstrap_samples):
        bs_read_counts = defaultdict(int, itertools.izip(read_ecs, np.random.multinomial(total_reads, read_probs, size=1)[0]))
        read_counts = read_counts + [bs_read_counts]

    with open('{0}/read_counts.kal'.format(args.output_dir), 'w') as f:
        for i in read_counts[0].keys():
            f.write("\t".join([str(i)] + [str(read_counts[j][i]) for j in range(len(read_counts))]) + "\n")


    timed_print('Running EM algorithm')
    alphas = [calculate_alpha(equiv_classes, read_counts[0], trans_ec_weights, args.num_iter, args.output_dir)]

    for i in range(args.bootstrap_samples):
        timed_print('Performing bootstrap iteration {0}'.format(i))
        alphas.append(calculate_alpha(equiv_classes, read_counts[i + 1], trans_ec_weights, args.bootstrap_iter, alpha = alphas[0]))

    # convert alphas to rhos
    rhos = []
    for i in range(len(alphas)):
        rhos.append([])
        for t, alpha_t in enumerate(alphas[i]):
            if eff_lens[t] == 0:
                if alpha_t > 0:
                    print('Error: transript {0} has zero effective length but non-zero abundance'.format(t)),
                    if i > 0:
                        print(' in bootstrap sample {0}'.format(i)),
                    print
                rhos[i].append(0)
            else:
                rhos[i].append(alpha_t/eff_lens[t])
        total = sum(rhos[i])
        rhos[i] = [rho/total for rho in rhos[i]]

    with open('{0}/results.kal'.format(args.output_dir), 'w') as f:
        f.write('name\teff_len\test_counts\ttpm' + '\t'.join(['tpm_{0}'.format(i) for i in range(1,len(alphas))]) + '\n')
        for t in range(len(alphas[0])):
            f.write('{0}\t{1}\t{2}\t'.format(transcriptome[t].id, eff_lens[t], alphas[0][t]))
            f.write("\t".join([str(10**6*rhos[j][t]) for j in range(len(alphas))]) + "\n")
