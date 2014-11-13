#! /usr/bin/python

from collections import defaultdict
from Bio import SeqIO
import itertools
import time

def timed_print( string ):
    print '[{0}] {1}'.format( time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()), string )

# equiv_class = id, set of kmers, dict of counts by transcript
class EquivalenceClass:
    """Equivalence class of a k-mer"""
    def __init__(self, ec_id, kmer_count = 0, counts = None ):
        self.id = ec_id
        self.kmer_count = kmer_count
        self.counts = counts or defaultdict(int)
    def __hash__( self ):
        return hash(tuple(self.counts.items()))
    def min( self, arg ):
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
    def __init__( self, equiv_classes ):
        self.equiv_classes = equiv_classes
        # set up ec_hash for recognizing equiv_classes after doing a min operation
        self.ec_hash = dict()
        for ec, ec_id in itertools.izip(self.equiv_classes.itervalues( ), self.equiv_classes.iterkeys( )):
            self.ec_hash[tuple(ec.counts)] = ec_id
        self.max_id = max([ec.id for ec in self.equiv_classes.itervalues()])
        self.memo = {}
    def __call__( self, *args ):
        ### maybe switch this to use try and except on key error?
        if args in self.memo:
            return self.memo[args]
        else:
            # take min
            m = self.equiv_classes[args[0]].min(self.equiv_classes[args[1]])
            # if it's just the first guy (we assume that args[0] < args[1], if anything) then we're done
            if m.id == args[0]:
                m_id = m.id
            else:
                if tuple(m.counts) in self.ec_hash:
                    m_id = self.ec_hash[tuple(m.counts)]
                else:
                    # if not, add it (including to ec_hash)
                    self.max_id += 1
                    m_id = self.max_id
                    self.ec_hash[tuple(m.counts)] = m_id
                    self.equiv_classes[m_id] = m
            m.id = m_id
            self.memo[args] = m_id
            return m_id


def kmers(seq,k):
    for i in xrange(len(seq)-k+1):
        yield seq[i:i+k]


#@profile
def make_equiv_classes( file_name, k ):
    kmer_ec = defaultdict(int) # holds the equiv class (id) for every k-mer. if we've never seen a k-mer before, then it'll get a new equiv class which is basically what equiv_classes[0] is
    equiv_classes = dict()
    equiv_classes[0] = EquivalenceClass(0)
    max_id = 0
    for trans_num, trans in enumerate(SeqIO.parse( file_name, "fasta" )):
        if trans_num % 1000 == 0:
            timed_print( 'processing transcript {0}'.format( trans_num ) )
        for reverse in [0,1]:
            trans_dict = defaultdict(int)
            if reverse:
                for kmer in kmers( str(trans.reverse_complement().seq), k ):
                    trans_dict[kmer] += 1
            else:
                for kmer in kmers( str(trans.seq), k ):
                    trans_dict[kmer] += 1
            kmers_by_ec_and_count = defaultdict(list)
            for kmer, count in trans_dict.iteritems():
                kmers_by_ec_and_count[(kmer_ec[kmer],count)].append(kmer)
            for (ec_id,count) in sorted(kmers_by_ec_and_count.keys()):
                curr_kmers = kmers_by_ec_and_count[(ec_id,count)]
                max_id += 1
                equiv_classes[max_id] = EquivalenceClass( max_id, len(curr_kmers), dict(equiv_classes[ec_id].counts) )
                equiv_classes[max_id].counts[2*trans_num + reverse] = count
                # remove these kmers from old equiv class and delete class if necessary
                if ec_id != 0:
                    equiv_classes[ec_id].kmer_count -= len(curr_kmers)
                    if( equiv_classes[ec_id].kmer_count < 0 ):
                        raise Exception("EC with negative members?")
                    if( equiv_classes[ec_id].kmer_count == 0 ):
                        del( equiv_classes[ec_id] )
                # update kmer_ec
                for kmer in curr_kmers:
                    kmer_ec[kmer] = max_id
#        if trans_num == 1000:
#            break
    return kmer_ec, equiv_classes


#@profile
def process_reads( equiv_classes, kmer_ec, read_file, k, memo_min = None ):
    if not memo_min:
        timed_print( 'creating ec_hash etc.' )
        memo_min = memoized_min( equiv_classes )
        timed_print( 'done' )
    read_counts = defaultdict(int)
#    for read_num, (read1, read2) in enumerate(itertools.izip( SeqIO.parse(read_file.format('1'),"fastq"), SeqIO.parse(read_file.format('2'),"fastq"))):
    f = open("uc001ugr.out",'w')
    for read_num, (read1, read2) in enumerate(itertools.izip( SeqIO.parse(read_file.format('1'),"fasta"), SeqIO.parse(read_file.format('2'),"fasta"))):
#        read2 = read2.reverse_complement()
        ecs = sorted(list(set([kmer_ec[kmer] for kmer in kmers(str(read1.seq),k) if kmer in kmer_ec] + [kmer_ec[kmer] for kmer in kmers(str(read2.seq),k) if kmer in kmer_ec])))
        if not ecs:
            curr_ec = 0
        else:
            curr_ec = ecs[0]
            for ec in ecs[1:]:
                curr_ec = memo_min( curr_ec, ec )
                if curr_ec == 0:
                    break
        read_counts[curr_ec] += 1
#        if 98064 not in equiv_classes[curr_ec].counts:
#            f.write( 'read #{0}:\n{1}\n{2}\n{3}\nMin = {4}\n\n'.format(read_num, read1.seq, read2.seq, '\n'.join(map(lambda i: 'EC #{0}: {1}'.format(i, repr(equiv_classes[i].counts)), ecs)), repr(equiv_classes[curr_ec].counts)) )
        if read_num % 100000 == 0:
            timed_print( 'processing read #{0}'.format(read_num) )
#        if read_num == 10000:
#            break
    return read_counts


def calculate_alphas( equiv_classes, read_counts, num_iter, lens, alpha = None ):
   # convert read counts to floats
   read_counts = dict(itertools.izip(read_counts.iterkeys(),map(float,read_counts.itervalues())))
   # invert ec=>transcript map
   transcript_ecs = defaultdict(dict)
   max_trans = 0
   total_mass = 0
   timed_print( "Starting ec=>transcript map inversion" )
   for ec_id in read_counts.iterkeys():
       total_mass += read_counts[ec_id]
       for trans, trans_count in equiv_classes[ec_id].counts.iteritems():
           max_trans = max(trans/2,max_trans)
           transcript_ecs[trans/2][ec_id] = trans_count
   num_trans = max_trans + 1
   total_mass = total_mass - read_counts[0]
   timed_print( "Done" )
   # initialize alphas if necessary
   if not alpha:
       alpha = [1.0/num_trans for x in range(num_trans)]
   # for current alphas, calculate the denominator for every ec
   for t in range(num_iter):
       timed_print( 'Iteration {0}'.format(t) )
       denom = defaultdict(int)
       masses = [0 for x in range(num_trans)]
       for ec_id in read_counts.iterkeys():
           for trans, trans_count in equiv_classes[ec_id].counts.iteritems():
               denom[ec_id] += alpha[trans/2]/lens[trans/2]*trans_count
       # for every transcript, calculate total mass
       for trans in range(num_trans):
           for ec_id, ec_weight in transcript_ecs[trans].iteritems():
               if denom[ec_id] == 0:
                   print( 'Arr, ec with zero stuff: {0}'.format( ec_id ))
               masses[trans] += read_counts[ec_id]*(ec_weight*alpha[trans]/(denom[ec_id]*lens[trans]))
       for trans in range(num_trans):
           alpha[trans] = masses[trans]/total_mass
   return alpha


def get_trans_lens( file_name ):
    lens = []
    for trans_num, trans in enumerate(SeqIO.parse( file_name, "fasta" )):
        lens.append( len(trans) )
    return lens


if __name__ == "__main__":
#    trans_file_name = "/home/bray/data/idx/Homo_sapiens.GRCh37.70.cdna.pc.fa"
#    read_file_name = "/home/bray/test_{0}.fastq"
    trans_file_name = "/scratch/bray/express_simdata/hg19-ucsc_transcripts.fa"
    read_file_name = "/scratch/bray/express_simdata/error_only/sim_30m_{0}.fa"
    k = 20
    num_iter = 100
    kmer_ec, equiv_classes = make_equiv_classes( trans_file_name, k )
    read_counts = process_reads( equiv_classes, kmer_ec, read_file_name, k )
    lens = get_trans_lens( trans_file_name )
    eff_lens = [l - 75 for l in eff_lens]
    alpha = calculate_alphas( equiv_classes, read_counts, num_iter, lens )

#    with open("alphas",'w') as f:
#        for i in range(len(alpha)):
#            f.write('{0}\n'.format(alpha[i]))
