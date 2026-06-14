# mask_terminal_damage_sites

A Python script that hard-masks (`N` + quality 0) SNP-dependent,
strand-dependent terminal damage sites in the first/last N bp of each
read in a BAM file, as a damage-mitigation step for ancient DNA (aDNA)
analysis.

It applies a similar idea to `pileupCaller --singleStrandMode`, but
restricted to the terminal N bp of each read.

## Overview

- Forward-mapped reads: bases within the terminal N bp that overlap a
  **C/T (or T/C) SNP** position are replaced with `N`
- Reverse-mapped reads: bases within the terminal N bp that overlap a
  **G/A (or A/G) SNP** position are replaced with `N`
- Terminal positions are determined by query position in the read
  sequence (soft-clip / insertion / deletion regions are excluded)
- Only the affected bases are rewritten; the read itself is never dropped
- `MD` / `NM` tags are removed from modified reads (can be kept with
  `--keep-tags`)
- Supplementary / secondary / unmapped reads are passed through unchanged

## Requirements

- Python 3
- [pysam](https://github.com/pysam-developers/pysam)
- pytest (for running the tests)

```bash
pip install pysam pytest
```

## Usage

```bash
python mask_terminal_damage_sites.py \
  --input sample.bam \
  --snp-file 1240k.snp \
  --output sample.masked.bam \
  --mask-length 5 \
  --index
```

### Options

| Option | Description | Default |
| --- | --- | --- |
| `--input` | Input BAM (coordinate-sorted) | required |
| `--snp-file` | SNP file (Eigenstrat .snp or PLINK .bim format) | required |
| `--output` | Output BAM | required |
| `--mask-length` | Number of terminal bases (at each end of the read) to consider for masking | 5 |
| `--keep-tags` | Keep `MD`/`NM` tags instead of removing them | tags removed |
| `--index` | Build a `.bai` index for the output BAM | not built |
| `--log-level` | Logging level (`DEBUG`/`INFO`/`WARNING`/`ERROR`) | `INFO` |

## SNP file format

Both Eigenstrat (`SNP_ID CHR GENPOS POS A1 A2`) and PLINK BIM
(`CHR SNP_ID GENPOS POS A1 A2`) column orders are supported.

```text
rs3094315 1 0.02012958 817186 G A
rs12124819 1 0.02024177 841166 A G
```

SNPs whose allele pair is `{C,T}` are loaded into `ct_sites`, and SNPs
whose allele pair is `{G,A}` are loaded into `ga_sites`. All other allele
combinations are ignored.

### Chromosome name variants

Differences in chromosome naming between the BAM and SNP file
(`1` / `chr1`, `X` / `chrX` / `23`, `MT` / `chrMT` / `M` / `chrM` / `90`,
etc.) are normalized internally.

## Log output

The following summary statistics are written to stderr:

- Number of C/T SNP sites and G/A SNP sites loaded
- Number of reads processed
- Total number of masked bases (forward / reverse breakdown)
- Number of reads with `MD`/`NM` tags removed

## Tests

```bash
pytest test_mask_terminal_damage_sites.py -v
```

Using small synthetic BAM/SNP files, the tests verify:

- Only the terminal base of a forward read is masked at a C/T SNP
- Reverse reads are not masked at C/T SNPs
- Only the terminal base of a reverse read is masked at a G/A SNP
- Forward reads are not masked at G/A SNPs
- SNPs in the middle of a read are not masked
- Soft-clipped regions are excluded from masking
- `MD`/`NM` tags are removed/kept correctly

## Notes

- Because SEQ/QUAL are rewritten, the existing `MD`/`NM` tags may become
  inconsistent. By default they are removed from modified reads.
- Reads without quality information (`QUAL="*"`) are skipped for masking.
