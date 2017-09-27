"""Analyze Illumina sequencing data for PUBS project."""

import re
from collections import namedtuple


# Regular expression for Illumina FASTQ header
ILLUMINA_HEADER_RE = re.compile(r'@([^:]+):(\d+):(\d+):(\d+):(\d+)#(\d+)/(\d+)')


IlluminaHeader = namedtuple(
    'IlluminaHeader',
    'instrument_name lane tile xcoord ycoord index read',
    module=__name__,
)
IlluminaHeader.__doc__ = '''\
NamedTuple of header information for Illuma sequencing read.

.. attribute:: instrument_name

    Unique instrument name (``str``).

.. attribute:: lane

    Flowcell lane (``int``).

.. attribute:: tile

    Tile number within flowcell lane (``int``).

.. attribute:: xcoord

    X coordinate of cluster within the tile (``int``).

.. attribute:: ycoord

    Y coordinate of cluster within the tile (``int``).

.. attribute:: index

    Index number for a multiplexed sample (``str``). Not sure if integer, leave
    as string for now.

.. attribute:: read

    Index of read within pair/triplet/etc. (``int``, one-indexed).
'''


def parse_illumina_header(header):
    """Parse the header line from a read in an Illumia FASTQ file.

    :param str header: Header line text, including first "@" character.

    :rtype: .IlluminaHeader
    """
    try:
        match = ILLUMINA_HEADER_RE.match(header)

        if match is None:
            raise ValueError('header not in expected format')

        instrument, lane, tile, xcoord, ycoord, index, read = match.groups()

        return IlluminaHeader(
            instrument,
            int(lane),
            int(tile),
            int(xcoord),
            int(ycoord),
            index,
            int(read),
        )

    except ValueError as exc:
        # Doesn't match regex or can't parse int
        raise ValueError('Invalid header line: {}'.format(exc)) from None

