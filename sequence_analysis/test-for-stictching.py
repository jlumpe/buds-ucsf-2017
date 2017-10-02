def readFastq(filename):
    sequences = []
    qualities = []
    with open (filename) as fh:
        while True:
            fh.readline() #skip name line
            seq = fh.readline().rstrip() #read base sequences
            fh.readline() #skip place holder line
            qual = fh.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

seqs1, quals1 = readFastq('R1_sample200x.fastq')
seqs3, quals3 = readFastq('R3_sample200x.fastq')

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

print "Reverse Complement of R3[17] is: %r" %reverseComplement(seqs3[16])
print "R1[17] is %r" %seqs1[16]

from Bio import pairwise2

alignments = pairwise2.align.globalxx(seqs1[16],reverseComplement(seqs3[16]))
#will need to change this part

from Bio.pairwise2 import format_alignment
print(format_alignment(*alignments[0]))
