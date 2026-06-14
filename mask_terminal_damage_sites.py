#!/usr/bin/env python3
"""
mask_terminal_damage_sites.py

Ancient DNA damage masking utility.

Hard-masks SNP-dependent, strand-dependent terminal damage sites in BAM
reads. The idea is analogous to ``pileupCaller --singleStrandMode``, but
restricted to the terminal ``mask_length`` bases (in read/query coordinates)
of each read.

Masking rule
------------
For each read base within the terminal ``mask_length`` bp (in *read* /
query coordinates) that has a corresponding reference position:

* If the read is mapped to the forward strand (``not read.is_reverse``) and
  the reference position is a C/T (or T/C) SNP (as listed in the SNP file),
  the base is set to ``N`` and its quality is set to 0.

* If the read is mapped to the reverse strand (``read.is_reverse``) and the
  reference position is a G/A (or A/G) SNP, the base is set to ``N`` and its
  quality is set to 0.

Only the affected query bases are modified -- the rest of the read (and the
read itself) is kept as-is. Supplementary, secondary and unmapped reads are
passed through unchanged.

Because SEQ/QUAL are modified, the existing ``MD`` and ``NM`` tags may
become inconsistent with the new sequence. By default these tags are
removed from modified reads; use ``--keep-tags`` to keep them as-is.
"""

import argparse
import logging
import sys

import pysam


# ---------------------------------------------------------------------------
# Chromosome name normalization
# ---------------------------------------------------------------------------

# Eigenstrat-style numeric codes used for human sex chromosomes / mtDNA.
_EIGENSTRAT_CODE_MAP = {
    "23": "X",
    "24": "Y",
    "90": "MT",
}


def normalize_chrom(chrom):
    """Normalize a chromosome name from a BAM file or a SNP file.

    Handles, among others:

    * ``"1"``, ``"chr1"``                       -> ``"1"``
    * ``"X"``, ``"chrX"``, ``"23"``             -> ``"X"``
    * ``"Y"``, ``"chrY"``, ``"24"``             -> ``"Y"``
    * ``"MT"``, ``"chrMT"``, ``"M"``, ``"chrM"``, ``"90"`` -> ``"MT"``

    Parameters
    ----------
    chrom : str
        Raw chromosome name as found in the BAM header or SNP file.

    Returns
    -------
    str
        Normalized chromosome name.
    """
    c = chrom.strip()
    if c[:3].lower() == "chr":
        c = c[3:]
    c = c.upper()
    if c == "M":
        c = "MT"
    c = _EIGENSTRAT_CODE_MAP.get(c, c)
    return c


# ---------------------------------------------------------------------------
# SNP file parsing
# ---------------------------------------------------------------------------

def _looks_like_chrom_token(token):
    """Return True if ``token`` normalizes to a value that looks like a
    plausible (human) chromosome name (1-22, X, Y, MT)."""
    norm = normalize_chrom(token)
    if norm in ("X", "Y", "MT"):
        return True
    if norm.isdigit():
        n = int(norm)
        if 1 <= n <= 22:
            return True
    return False


def load_snp_sites(snp_file):
    """Parse an Eigenstrat/PLINK-BIM-like SNP file.

    Two column orderings are accepted (both whitespace separated, >= 6
    columns)::

        SNP_ID  CHR  GENPOS  PHYSPOS  A1  A2     (Eigenstrat .snp style)
        CHR     SNP_ID  GENPOS  PHYSPOS  A1  A2  (PLINK .bim style)

    Only the chromosome, physical position (1-based) and the two alleles
    are used.

    Parameters
    ----------
    snp_file : str
        Path to the SNP file.

    Returns
    -------
    (ct_sites, ga_sites) : tuple of dict
        ``ct_sites`` maps ``(normalized_chrom, pos_1based) -> True`` for
        C/T (or T/C) SNPs. ``ga_sites`` maps the same key space for G/A
        (or A/G) SNPs. All other allele combinations are ignored.
    """
    ct_sites = {}
    ga_sites = {}

    with open(snp_file) as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 6:
                logging.warning(
                    "SNP file line %d: fewer than 6 columns, skipping: %r",
                    lineno, line,
                )
                continue

            f0, f1, _genpos, pos_str, a1, a2 = fields[:6]

            # Determine which of the first two columns is the chromosome.
            if _looks_like_chrom_token(f1):
                # Eigenstrat-like: SNP_ID CHR GENPOS PHYSPOS A1 A2
                chrom = f1
            elif _looks_like_chrom_token(f0):
                # PLINK BIM-like: CHR SNP_ID GENPOS PHYSPOS A1 A2
                chrom = f0
            else:
                logging.warning(
                    "SNP file line %d: could not identify chromosome "
                    "column, skipping: %r", lineno, line,
                )
                continue

            try:
                pos = int(pos_str)
            except ValueError:
                logging.warning(
                    "SNP file line %d: non-integer position %r, skipping",
                    lineno, pos_str,
                )
                continue

            alleles = {a1.upper(), a2.upper()}
            key = (normalize_chrom(chrom), pos)

            if alleles == {"C", "T"}:
                ct_sites[key] = True
            elif alleles == {"G", "A"}:
                ga_sites[key] = True
            # other allele combinations (e.g. indels, A/C, ...) are ignored

    return ct_sites, ga_sites


# ---------------------------------------------------------------------------
# Per-read masking
# ---------------------------------------------------------------------------

def mask_read_bases(read, mask_length, ct_sites, ga_sites, norm_chrom):
    """Hard-mask terminal SNP-dependent damage sites in a single read.

    The read's ``query_sequence`` / ``query_qualities`` are modified
    in-place when at least one base needs masking.

    Parameters
    ----------
    read : pysam.AlignedSegment
        The read to (potentially) modify. Must be mapped.
    mask_length : int
        Number of bases at each end of the read (in query/read
        coordinates) to consider for masking.
    ct_sites, ga_sites : dict
        Dictionaries as returned by :func:`load_snp_sites`.
    norm_chrom : str
        Normalized chromosome name of ``read``.

    Returns
    -------
    int
        Number of bases masked in this read (0 if none).
    """
    seq = read.query_sequence
    if seq is None:
        return 0

    read_length = len(seq)
    if read_length == 0:
        return 0

    quals = read.query_qualities
    if quals is None:
        # No quality information available (e.g. SEQ present but QUAL "*").
        # We cannot meaningfully set quality to 0, so skip masking for
        # this read rather than fabricate quality values.
        logging.debug(
            "Read %s has no query qualities; skipping masking",
            read.query_name,
        )
        return 0
    quals = list(quals)
    seq_list = list(seq)

    # Site dictionary depends on the strand the read is mapped to.
    site_dict = ga_sites if read.is_reverse else ct_sites

    n_masked = 0

    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
        if query_pos is None or ref_pos is None:
            # Deletion/N (no query base) or insertion/soft-clip (no
            # reference position) -- not maskable.
            continue

        is_terminal = (
            query_pos < mask_length
            or query_pos >= read_length - mask_length
        )
        if not is_terminal:
            continue

        ref_pos_1based = ref_pos + 1  # get_aligned_pairs is 0-based
        if (norm_chrom, ref_pos_1based) in site_dict:
            seq_list[query_pos] = "N"
            quals[query_pos] = 0
            n_masked += 1

    if n_masked == 0:
        return 0

    # NOTE: assigning to query_sequence resets query_qualities, so the
    # sequence must be set first and the (updated) qualities afterwards.
    read.query_sequence = "".join(seq_list)
    read.query_qualities = quals

    return n_masked


def remove_md_nm_tags(read):
    """Remove ``MD`` and ``NM`` tags from ``read`` if present.

    Returns
    -------
    bool
        True if at least one of the tags was present and removed.
    """
    removed = False
    if read.has_tag("MD"):
        read.set_tag("MD", None)
        removed = True
    if read.has_tag("NM"):
        read.set_tag("NM", None)
        removed = True
    return removed


# ---------------------------------------------------------------------------
# Main processing loop
# ---------------------------------------------------------------------------

def process_bam(input_bam, snp_file, output_bam, mask_length=5,
                 keep_tags=False, build_index=False):
    """Process ``input_bam``, masking terminal damage sites, and write
    the result to ``output_bam``.

    Parameters
    ----------
    input_bam : str
        Path to the input BAM file.
    snp_file : str
        Path to the SNP file (Eigenstrat/PLINK-BIM-like).
    output_bam : str
        Path to the output BAM file.
    mask_length : int
        Number of terminal bases (in read coordinates) to consider for
        masking at each end of the read.
    keep_tags : bool
        If True, do not remove MD/NM tags from modified reads.
    build_index : bool
        If True, build a ``.bai`` index for ``output_bam`` after writing.

    Returns
    -------
    dict
        Summary statistics (also emitted via the logging module).
    """
    ct_sites, ga_sites = load_snp_sites(snp_file)
    logging.info("Loaded %d C/T SNP site(s)", len(ct_sites))
    logging.info("Loaded %d G/A SNP site(s)", len(ga_sites))

    n_reads_processed = 0
    n_bases_masked = 0
    n_bases_masked_fwd = 0
    n_bases_masked_rev = 0
    n_reads_tags_removed = 0

    with pysam.AlignmentFile(input_bam, "rb") as infile:
        with pysam.AlignmentFile(output_bam, "wb", template=infile) as outfile:
            for read in infile.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    outfile.write(read)
                    continue

                n_reads_processed += 1

                norm_chrom = normalize_chrom(read.reference_name)

                n_masked = mask_read_bases(
                    read, mask_length, ct_sites, ga_sites, norm_chrom
                )

                if n_masked > 0:
                    n_bases_masked += n_masked
                    if read.is_reverse:
                        n_bases_masked_rev += n_masked
                    else:
                        n_bases_masked_fwd += n_masked

                    if not keep_tags:
                        if remove_md_nm_tags(read):
                            n_reads_tags_removed += 1

                outfile.write(read)

    if build_index:
        pysam.index(output_bam)

    stats = {
        "n_ct_sites": len(ct_sites),
        "n_ga_sites": len(ga_sites),
        "n_reads_processed": n_reads_processed,
        "n_bases_masked": n_bases_masked,
        "n_bases_masked_fwd": n_bases_masked_fwd,
        "n_bases_masked_rev": n_bases_masked_rev,
        "n_reads_tags_removed": n_reads_tags_removed,
    }

    logging.info("Processed reads (primary, mapped): %d", n_reads_processed)
    logging.info("Total masked bases: %d", n_bases_masked)
    logging.info("  - masked bases on forward-mapped reads: %d", n_bases_masked_fwd)
    logging.info("  - masked bases on reverse-mapped reads: %d", n_bases_masked_rev)
    logging.info("Reads with MD/NM tags removed: %d", n_reads_tags_removed)

    return stats


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_arg_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Hard-mask SNP-dependent, strand-dependent terminal damage "
            "sites (the terminal N bp of each read) in a BAM file."
        )
    )
    parser.add_argument(
        "--input", required=True, metavar="BAM",
        help="Input BAM file (coordinate-sorted).",
    )
    parser.add_argument(
        "--snp-file", required=True, metavar="SNP",
        help="SNP file in Eigenstrat (.snp) or PLINK (.bim) style "
             "(whitespace-separated, >= 6 columns).",
    )
    parser.add_argument(
        "--output", required=True, metavar="BAM",
        help="Output BAM file.",
    )
    parser.add_argument(
        "--mask-length", type=int, default=5, metavar="N",
        help="Number of bases at each end of the read (in read/query "
             "coordinates) to consider for masking. Default: 5.",
    )
    parser.add_argument(
        "--keep-tags", action="store_true",
        help="Keep MD/NM tags on modified reads (they may become "
             "inconsistent with the new SEQ/QUAL). Default: remove them "
             "from modified reads.",
    )
    parser.add_argument(
        "--index", action="store_true",
        help="Build a BAM index (.bai) for the output BAM after writing.",
    )
    parser.add_argument(
        "--log-level", default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity. Default: INFO.",
    )
    return parser


def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
        stream=sys.stderr,
    )

    if args.mask_length < 0:
        parser.error("--mask-length must be >= 0")

    process_bam(
        input_bam=args.input,
        snp_file=args.snp_file,
        output_bam=args.output,
        mask_length=args.mask_length,
        keep_tags=args.keep_tags,
        build_index=args.index,
    )


if __name__ == "__main__":
    main()
