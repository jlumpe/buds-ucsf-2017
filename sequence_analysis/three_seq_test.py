# removing n's from fastq
# 27 sept 2017

from Bio import SeqIO
import math


def call_filters(R1, R2, R3):
    """
    :param R1:
    :param R2:
    :param R3:
    :return:
    """
    filter_failures = []
    filter_passes = []

    for record1, record2, record3 in zip(R1, R2, R3):

        if  quality_filter(record1, 1) and quality_filter(record2, 1) and quality_filter(record3, 1) and n_filter(record2) and n_filter(record1) and n_filter(record3):
                filter_passes.append((record1, record2, record3))
        else:
            filter_failures.append((record1, record2, record3))
    return filter_passes, filter_failures


def n_filter(record):
    """
    :param record:
    :return:
    """

    if "N" not in record.seq:
        return True
    else:
        return False


def quality_filter(record, e_max):
    """
    quality score according to http://drive5.com/usearch/manual/exp_errs.html
    :param record:
    :param e_max:
    :return:
    """
    p = 0

    for phred in record.letter_annotations["phred_quality"]:
        p += (10**(-phred/10))

    if math.floor(p) < e_max:
        return True
    else:
        return False


if __name__ == "__main__":

    R1 = SeqIO.parse('R1_sample200x.fastq', "fastq")
    R2 = SeqIO.parse('R2_sample200x.fastq', "fastq")
    R3 = SeqIO.parse('R3_sample200x.fastq', "fastq")

    passes, failures = call_filters(R1, R2, R3)

    total = len(passes) + len(failures)

    print(len(passes)/total)