#!/usr/bin/env python3

"""Randomly subsample a set of FASTQ files.

This script uses a naive implementaiton that assumes each sequence in a FASTQ
file spans exactly four lines, with another sequence immediately thereafter.
I'm not actually sure if that is the case (e.g., can there be blank lines?).
"""

import sys
import os
import random
import argparse
from itertools import islice


parser = argparse.ArgumentParser(description='Randomly subsample a set of FASTQ files.')

parser.add_argument('factor', type=float, help='Factor by which to subsample the files')
parser.add_argument('files', type=str, nargs='+', help='FASTQ files to subsample')
parser.add_argument('-o', '--out', type=str, help='Directory to write output to')
parser.add_argument('-f', '--force', action='store_true', help='Overwite files in output directory')


def main(argv=None):

    args = parser.parse_args(argv)

    # Get output file paths
    outpaths = []
    for file in args.files:
        name, ext = os.path.splitext(os.path.basename(file))
        outname = '{}_sample{:.0f}x{}'.format(name, args.factor, ext)
        outpaths.append(os.path.join(args.out or '', outname))

    # Check output files do not exist
    if not args.force:
        for path in outpaths:
            if os.path.exists(path):
                print('Output file {} exists, refusing to overwrite'.format(path))

    # Open files
    instreams = [open(path) for path in args.files]
    outstreams = [open(path, 'wt') for path in outpaths]

    # Try creating a progress bar
    try:
        from tqdm import tqdm

    except ImportError:
        pbar = None

    else:
        pbar = tqdm()

    # Probability to sample each set of reads
    p = 1 / args.factor

    sampled = 0
    total = 0

    while True:

        # Read four lines from each file
        line_groups = [list(islice(fobj, 4)) for fobj in instreams]

        # Done if no lines read
        if not any(line_groups):
            break

        # Check for partial or missing reads
        for file, lines in zip(args.files, line_groups):
            if not lines:
                print('Reached end of file {} before others'.format(file))
                sys.exit(1)

            if len(lines) < 4:
                print('File {} ends in middle of read'.format(file))
                sys.exit(1)

        # Sample?
        if random.random() < p:
            sampled += 1

            # Write output
            for fobj, lines in zip(outstreams, line_groups):
                for line in lines:
                    fobj.write(line)

        # Update progress
        total += 1
        if pbar is not None:
            pbar.update()

    # Close stuff
    if pbar is not None:
        pbar.close()

    for fobj in instreams:
        fobj.close()

    for fobj in outstreams:
        fobj.close()

    print('Read {} sequences, sampled {} ({:.1f}x)'.format(total, sampled, total/sampled))


if __name__ == '__main__':
    main()

