import pytest

import seqanalysis as sa


def test_parse_header():
    """Test parse_illumina_header function."""

    # Example from Illumina help page
    header = '@SIM:1:FCX:1:15:6329:1045 1:N:0:2'
    expected = ('SIM', 1, 'FCX', 1, 15, 6329, 1045, 1, False, 0, 2)

    # Basic parsing
    parsed = sa.parse_illumina_header(header)
    assert isinstance(parsed, sa.IlluminaHeader)
    assert parsed == expected

    # Without @
    assert sa.parse_illumina_header(header[1:]) == expected

    # Invalid format - doesn't match regex
    with pytest.raises(ValueError):
        sa.parse_illumina_header('@this:is:not:a:proper:header')

    # Invalid format - not integers
    with pytest.raises(ValueError):
        sa.parse_illumina_header('@SIM:1:FCX:1:15:notanint:1045 1:N:0:2')

    # Invalid format - not Y/N
    with pytest.raises(ValueError):
        sa.parse_illumina_header('@SIM:1:FCX:1:15:6329:1045 1:Q:0:2')

