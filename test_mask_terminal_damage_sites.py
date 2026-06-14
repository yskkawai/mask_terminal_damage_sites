"""
Tests for mask_terminal_damage_sites.py

These tests build small synthetic BAM files (and a small SNP file) on the
fly, run the masking script's core function (`process_bam`) on them, and
check that the resulting BAM was modified exactly as expected.
"""

import pysam
import pytest

import mask_terminal_damage_sites as mts


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

CHROM_NAME = "chr1"
CHROM_LEN = 1000

# 20-base read sequence used for all test reads.
SEQ_20 = "ACGTACGTACGTACGTACGT"
QUAL_20 = [30] * 20


def make_read(qname, ref_start, cigartuples, is_reverse, seq=SEQ_20,
               qual=None, with_md_nm=True):
    """Build a pysam.AlignedSegment for the test BAM."""
    if qual is None:
        qual = list(QUAL_20)

    a = pysam.AlignedSegment()
    a.query_name = qname
    a.query_sequence = seq
    a.flag = 16 if is_reverse else 0
    a.reference_id = 0
    a.reference_start = ref_start
    a.mapping_quality = 60
    a.cigartuples = cigartuples
    a.query_qualities = pysam.qualitystring_to_array(
        "".join(chr(q + 33) for q in qual)
    )
    if with_md_nm:
        # MD/NM values themselves are not biologically meaningful here;
        # we only care whether they get removed by the masking code.
        a.set_tag("NM", 0)
        a.set_tag("MD", str(len(seq)))
    return a


@pytest.fixture
def header():
    return {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM_NAME, "LN": CHROM_LEN}],
    }


@pytest.fixture
def snp_file(tmp_path):
    """SNP file (Eigenstrat-style: ID CHR GENPOS POS A1 A2).

    - chr1:101  C/T  -> ct_sites  (used by readA, readB, readF)
    - chr1:310  C/T  -> ct_sites  (used by readE, middle of read)
    - chr1:150  G/A  -> ga_sites  (used by readC, readD)
    """
    path = tmp_path / "test.snp"
    path.write_text(
        "rs1 1 0.0 101 C T\n"
        "rs2 1 0.0 310 C T\n"
        "rs3 1 0.0 150 G A\n"
    )
    return str(path)


@pytest.fixture
def input_bam(tmp_path, header):
    """Build a small coordinate-sorted BAM with the test reads.

    All reads are 20bp, mask_length will be the default of 5.

    readA_ct_fwd : forward, 20M @ ref_start=100 (1-based 101-120)
                   -> qpos0 maps to chr1:101 (C/T SNP), terminal,
                      forward -> should be masked (N at qpos0)

    readB_ct_rev : reverse, 20M @ ref_start=100 (1-based 101-120)
                   -> qpos0 maps to chr1:101 (C/T SNP), terminal,
                      reverse -> uses ga_sites -> NOT masked

    readE_ct_mid : forward, 20M @ ref_start=300 (1-based 301-320)
                   -> qpos9 maps to chr1:310 (C/T SNP), but qpos9 is
                      NOT terminal (mask_length=5) -> NOT masked

    readF_softclip : forward, 5S15M @ ref_start=100
                   -> qpos0-4 are soft-clipped (no ref pos) -> skipped
                   -> qpos5 maps to chr1:101 (C/T SNP) but qpos5 is
                      NOT terminal -> NOT masked; read unchanged

    readC_ga_rev : reverse, 20M @ ref_start=130 (1-based 131-150)
                   -> qpos19 maps to chr1:150 (G/A SNP), terminal,
                      reverse -> should be masked (N at qpos19)

    readD_ga_fwd : forward, 20M @ ref_start=130 (1-based 131-150)
                   -> qpos19 maps to chr1:150 (G/A SNP), terminal,
                      forward -> uses ct_sites -> NOT masked
    """
    path = tmp_path / "input.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        bam.write(make_read("readA_ct_fwd", 100, [(0, 20)], is_reverse=False))
        bam.write(make_read("readB_ct_rev", 100, [(0, 20)], is_reverse=True))
        bam.write(make_read("readE_ct_mid", 300, [(0, 20)], is_reverse=False))
        bam.write(make_read("readF_softclip", 100, [(4, 5), (0, 15)], is_reverse=False))
        bam.write(make_read("readC_ga_rev", 130, [(0, 20)], is_reverse=True))
        bam.write(make_read("readD_ga_fwd", 130, [(0, 20)], is_reverse=False))
    return str(path)


def _read_dict(bam_path):
    """Return {query_name: AlignedSegment} for all reads in a BAM."""
    out = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            out[read.query_name] = read
    return out


@pytest.fixture
def output_reads(tmp_path, input_bam, snp_file):
    """Run process_bam and return {query_name: AlignedSegment}."""
    output_bam = tmp_path / "output.bam"
    stats = mts.process_bam(
        input_bam=input_bam,
        snp_file=snp_file,
        output_bam=str(output_bam),
        mask_length=5,
        keep_tags=False,
        build_index=False,
    )
    reads = _read_dict(str(output_bam))
    return reads, stats


# ---------------------------------------------------------------------------
# normalize_chrom
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("raw,expected", [
    ("1", "1"),
    ("chr1", "1"),
    ("X", "X"),
    ("chrX", "X"),
    ("23", "X"),
    ("Y", "Y"),
    ("chrY", "Y"),
    ("24", "Y"),
    ("MT", "MT"),
    ("chrMT", "MT"),
    ("M", "MT"),
    ("chrM", "MT"),
    ("90", "MT"),
])
def test_normalize_chrom(raw, expected):
    assert mts.normalize_chrom(raw) == expected


# ---------------------------------------------------------------------------
# load_snp_sites
# ---------------------------------------------------------------------------

def test_load_snp_sites_eigenstrat_style(tmp_path):
    p = tmp_path / "x.snp"
    p.write_text(
        "rs123 1 0.0 100000 C T\n"
        "rs456 1 0.0 100050 G A\n"
        "rs789 2 0.0 100100 A C\n"  # ignored: not CT/GA pair
    )
    ct_sites, ga_sites = mts.load_snp_sites(str(p))
    assert ("1", 100000) in ct_sites
    assert ("1", 100050) in ga_sites
    assert ("2", 100100) not in ct_sites
    assert ("2", 100100) not in ga_sites


def test_load_snp_sites_bim_style(tmp_path):
    p = tmp_path / "x.bim"
    p.write_text(
        "1 rs123 0 100000 C T\n"
        "1 rs456 0 100050 G A\n"
    )
    ct_sites, ga_sites = mts.load_snp_sites(str(p))
    assert ("1", 100000) in ct_sites
    assert ("1", 100050) in ga_sites


def test_load_snp_sites_chrom_normalization(tmp_path):
    p = tmp_path / "x.snp"
    p.write_text(
        "rsA chr1 0.0 100 C T\n"
        "rsB 23 0.0 200 G A\n"   # Eigenstrat code for X
        "rsC chrY 0.0 300 C T\n"
    )
    ct_sites, ga_sites = mts.load_snp_sites(str(p))
    assert ("1", 100) in ct_sites
    assert ("X", 200) in ga_sites
    assert ("Y", 300) in ct_sites


# ---------------------------------------------------------------------------
# End-to-end masking behaviour
# ---------------------------------------------------------------------------

def test_ct_snp_forward_read_terminal_base_masked(output_reads):
    reads, _ = output_reads
    r = reads["readA_ct_fwd"]
    assert r.query_sequence[0] == "N"
    assert r.query_qualities[0] == 0
    # rest of the read untouched
    assert r.query_sequence[1:] == SEQ_20[1:]
    assert list(r.query_qualities[1:]) == QUAL_20[1:]


def test_ct_snp_reverse_read_not_masked(output_reads):
    reads, _ = output_reads
    r = reads["readB_ct_rev"]
    assert r.query_sequence == SEQ_20
    assert list(r.query_qualities) == QUAL_20


def test_ga_snp_reverse_read_terminal_base_masked(output_reads):
    reads, _ = output_reads
    r = reads["readC_ga_rev"]
    assert r.query_sequence[19] == "N"
    assert r.query_qualities[19] == 0
    # rest of the read untouched
    assert r.query_sequence[:19] == SEQ_20[:19]
    assert list(r.query_qualities[:19]) == QUAL_20[:19]


def test_ga_snp_forward_read_not_masked(output_reads):
    reads, _ = output_reads
    r = reads["readD_ga_fwd"]
    assert r.query_sequence == SEQ_20
    assert list(r.query_qualities) == QUAL_20


def test_middle_of_read_snp_not_masked(output_reads):
    reads, _ = output_reads
    r = reads["readE_ct_mid"]
    assert r.query_sequence == SEQ_20
    assert list(r.query_qualities) == QUAL_20


def test_softclip_region_not_masked(output_reads):
    reads, _ = output_reads
    r = reads["readF_softclip"]
    assert r.query_sequence == SEQ_20
    assert list(r.query_qualities) == QUAL_20


def test_md_nm_tags_removed_for_masked_reads(output_reads):
    reads, _ = output_reads

    masked = reads["readA_ct_fwd"]
    assert not masked.has_tag("MD")
    assert not masked.has_tag("NM")

    masked2 = reads["readC_ga_rev"]
    assert not masked2.has_tag("MD")
    assert not masked2.has_tag("NM")


def test_md_nm_tags_kept_for_unmasked_reads(output_reads):
    reads, _ = output_reads
    for qname in ("readB_ct_rev", "readD_ga_fwd", "readE_ct_mid", "readF_softclip"):
        r = reads[qname]
        assert r.has_tag("MD")
        assert r.has_tag("NM")


def test_summary_stats(output_reads):
    _, stats = output_reads
    assert stats["n_ct_sites"] == 2
    assert stats["n_ga_sites"] == 1
    assert stats["n_reads_processed"] == 6
    assert stats["n_bases_masked"] == 2
    assert stats["n_bases_masked_fwd"] == 1
    assert stats["n_bases_masked_rev"] == 1
    assert stats["n_reads_tags_removed"] == 2


# ---------------------------------------------------------------------------
# keep-tags option
# ---------------------------------------------------------------------------

def test_keep_tags_option(tmp_path, input_bam, snp_file):
    output_bam = tmp_path / "output_keep.bam"
    mts.process_bam(
        input_bam=input_bam,
        snp_file=snp_file,
        output_bam=str(output_bam),
        mask_length=5,
        keep_tags=True,
        build_index=False,
    )
    reads = _read_dict(str(output_bam))

    # Still masked...
    r = reads["readA_ct_fwd"]
    assert r.query_sequence[0] == "N"
    assert r.query_qualities[0] == 0
    # ...but tags kept.
    assert r.has_tag("MD")
    assert r.has_tag("NM")


# ---------------------------------------------------------------------------
# pass-through of unmapped/secondary/supplementary reads
# ---------------------------------------------------------------------------

def test_unmapped_read_passed_through(tmp_path, header, snp_file):
    input_path = tmp_path / "input_unmapped.bam"
    with pysam.AlignmentFile(str(input_path), "wb", header=header) as bam:
        unmapped = make_read("unmapped_read", 0, None, is_reverse=False,
                              with_md_nm=False)
        unmapped.is_unmapped = True
        unmapped.cigartuples = None
        unmapped.reference_id = -1
        unmapped.reference_start = -1
        bam.write(unmapped)

    output_path = tmp_path / "output_unmapped.bam"
    stats = mts.process_bam(
        input_bam=str(input_path),
        snp_file=snp_file,
        output_bam=str(output_path),
        mask_length=5,
    )

    reads = _read_dict(str(output_path))
    assert reads["unmapped_read"].query_sequence == SEQ_20
    assert stats["n_reads_processed"] == 0
