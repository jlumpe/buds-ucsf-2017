"""Analyze Illumina sequencing data for PUBS project."""

import re
from collections import namedtuple
from itertools import islice

import numpy as np
from Bio import SeqIO


# Regular expression for Illumina FASTQ header
ILLUMINA_HEADER_RE = re.compile(r'([^:]+):(\d+):([^:]+):(\d+):(\d+):(\d+):(\d+) (\d+):([YN]):(\d+):(\d+)')


# Named tuple for illumina header data
IlluminaHeader = namedtuple(
    'IlluminaHeader',
    'instrument run flowcell_id lane tile xpos ypos read is_filtered controlnum, samplenum',
    module=__name__,
)
IlluminaHeader.__doc__ = '''http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm'''


def zip_reads(files, format_, *, subsample=None, random=True):
    """
    Parse sequence reads from multiple files in parallel, yielding corresponding
    reads in tuples.

    :param files: Sequence of three sequence files, as file path strings or
        readable file-like objects
    :param str format_: Sequence file format as understood by
        :func:`Bio.SeqIO.parse`.
    :param subsample: If not None (default) subsample reads by this factor,
        ignoring the rest. Exact behavior depends on value of ``random``.
    :type subsample: int or float
    :param bool random: Only has an effect if ``subsample`` is not None. If
        True sample randomly, yielding on average only 1 out of this many
        read tuples. If False ``subsample`` should be an integer and will
        subsample in exact invervals of this many tuples.

    :returns: Iterator yielding N-tuples of :class:`Bio.SeqIO.SeqRecord`, where
        N is the length of ``files``.
    """

    # Use BioPython's parser func on each file, these are lazy and so return
    # iterators that only parse the next sequence on demand
    parsers = [SeqIO.parse(file, format_) for file in files]

    # Swap from N iterators yielding single records to one iterator yielding
    # N-tuples of records
    tuples = zip(*parsers)

    if subsample is None:
        # No subsampling - yield everything
        yield from tuples

    elif random:
        # Random subsampling
        from random import random
        p = 1 / subsample

        for tup in tuples:
            if random() < p:
                yield tup

    else:
        # Even slicing
        yield from islice(tuples, None, None, subsample)

    # Verify that all iterators are exhausted
    # (iterator returned by zip() stops when the first sub-iterator runs out)
    for parser in parsers:
        try:
            next(parser)

        except StopIteration:
            # Iterators raise StopIteration if next() called when they are empty
            pass

        else:
            # Wasn't empty
            from warnings import warn
            warn('Not all sequence files contain the same number of records')
            break


def default_context():
    """Create the default context dictionary.

    :rtype: dict
    """
    # TODO
    return dict()


def process_read_files(files, ctx=None, progress=False, **kwargs):
    """Parse and process reads from a triplet of files.

    :param files: Sequence of 3 sequence files, as file path strings or readable
        file-like objects.
    :param dict ctx: Context dictionary. If None will use default.
    :param bool progress: If True display a progress bar with :mod:`tqdm`.
    :param \\**kwargs: Additional keyword arguments to pass to
        :func:`.zip_reads`.

    :returns: ``(valid, filtered, failed)`` tuple. Each is a list of 2-tuples
        where the first element is the index of the read triplet within the
        files. Each triplet appears once in only one if the lists. ``valid``
        corresponds to reads which passed filtering and were processed
        successfully. The 2nd item for each pair is a dict of statistics
        calculated for the triplet. ``filtered`` corresponds to the filtered
        triplets and the 2nd item is a string containing the reason for being
        filtered. ``failed`` corresponds to triplets for which an exception was
        raised during processing, the 2nd item is the raised exception.
    """

    if len(files) != 3:
        raise ValueError('Must give three files')

    if ctx is None:
        ctx = default_context()

    triplets = zip_reads(files, 'fastq-illumina', **kwargs)

    if progress:
        # Wrap iterator in tqdm() to display progress
        from tqdm import tqdm
        triplets = tqdm(triplets)

    valid = []
    filtered = []
    failed = []

    # Process triplets
    for i, triplet in enumerate(triplets):
        try:
            stats = process_triplet(triplet, ctx)

        except Exception as exc:
            # Failed
            failed.append((i, exc))
            continue

        if stats is None:
            # Filtered
            # TODO
            filtered.append((i, '?'))

        else:
            # OK
            valid.append((i, stats))

    return valid, filtered, failed


def process_triplet(triplet, ctx):
    """Process a read triplet, apply filters, get mutation, barcode, and stats.

    :param tuple triplet: 3-tuple of :class:`Bio.SeqIO.SeqRecord` ojects for
        reads 1, 2, and 3.
    :param dict ctx: Context dictionary
    :returns: Dictionary of calculated values, or None if the triplet did not
        pass the filters.
    :rtype: dict
    """

    rec1, rec2, rec3 = triplet

    stats = dict()

    # TODO
    filters = []

    # Run filters
    for filter_ in filters:
        if not filter_(triplet):
            return None

    # Find mutation
    stats['mutation'] = find_mutation(r1, r3)

    # Find barcode sequence
    stats['barcode'] = get_barcode(r2)

    # TODO
    # Collect more statistics

    return stats


def find_mutation(read1, read3, ctx):
    """Find the mutation in the WT sequence from reads 1 and 3.

    :param read1: Read 1 of triplet (protein forward).
    :param read3: Read 3 of triplet (protein reverse).
    :param dict ctx: Context dictionary
    :returns: 2-tuple of ``(residue_index, new_codon)``. If no mutation was
        found, both tuple elements will be None.
    """
    # TODO


def get_barcode(read2, ctx):
    """Get the barcode sequence from read 2.

    :param read2: 2nd read of triplet.
    :type read2: Bio.SeqIO.SeqRecord
    :param dict ctx: Context dictionary
    :returns: Barcode sequence
    :rtype: str
    """
    # TODO


def parse_illumina_header(header):
    """Parse the header line from a read in an Illumia FASTQ file.

    Source: http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm

    :param str header: Header line text, first "@" character optional.

    :rtype: .IlluminaHeader
    """
    # Strip whitespace and remove leading @
    header = header.strip()
    if header[0] == '@':
        header = header[1:]

    # Parse w/ regex
    match = ILLUMINA_HEADER_RE.match(header)

    if match is None:
        raise ValueError('header not in expected format')

    # Unpack
    inst, run, fcell, lane, tile, x, y, read, filtered, control, sample = match.groups()

    try:
        return IlluminaHeader(
            instrument=inst,
            run=int(run),
            flowcell_id=fcell,
            lane=int(lane),
            tile=int(tile),
            xpos=int(x),
            ypos=int(y),
            read=int(read),
            is_filtered=filtered == 'Y',
            controlnum=int(control),
            samplenum=int(sample),
        )

    except ValueError as exc:
        # Can't parse int
        raise ValueError('Invalid header line: {}'.format(exc)) from None

