import pytest

import seqanalysis as sa


def test_parse_header():
    """Test parse_illumina_header function."""

    # Wikipedia example
    header = sa.parse_illumina_header('@HWUSI-EAS100R:6:73:941:1973#0/1')
    assert isinstance(header, sa.IlluminaHeader)
    assert header == ('HWUSI-EAS100R', 6, 73, 941, 1973, '0', 1)

    # Invalid format - doesn't match regex
    with pytest.raises(ValueError):
        sa.parse_illumina_header('@this:is:not:a#proper/header')

    # Invalid format - not integers
    with pytest.raises(ValueError):
        sa.parse_illumina_header('@HWUSI-EAS100R:6:NOTANINT:941:1973#0/1')

