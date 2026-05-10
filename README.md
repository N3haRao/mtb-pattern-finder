# M. tuberculosis Genome Pattern Finder

This is my submission for the Rock Lab technical interview task. The project searches for 24 specific DNA sequence patterns in the *Mycobacterium tuberculosis* H37Rv reference genome and maps every hit to annotated genomic features.

---

## Project Description

The *M. tuberculosis* H37Rv genome is about 4.4 million nucleotides long and is encoded as a single circular chromosome. Given a set of 24 seven-character patterns (where `N` means "any nucleotide"), I search the entire genome using compiled regular expressions, then cross-reference every hit against a GFF3 annotation file to report which annotated feature (gene, CDS, tRNA, etc.) each hit falls in.

The result is a detailed, tab-separated output listing every match position, the matched sequence, and the overlapping annotated region IDs.

---

## Repository Layout

```
task/
  pattern_finder.py          # Main analysis script (start here)
  genome.txt                 # M. tuberculosis H37Rv genome (FASTA)
  annotation.txt             # Genomic feature annotations (GFF3)
  pattern_results.txt        # Full output: every match with region annotations
  pattern_results_summary.txt  # Condensed summary statistics
  analysis_writeup.md        # My writeup: tricky parts, decisions, and next steps
  task.md                    # The original task description from Dr. DeJesus 
  requirements.txt           # Python dependencies (standard library only)
```

---

## How to Run

I only use the Python standard library, so there is nothing to install beyond Python 3.6+.

```bash
git clone https://github.com/N3haRao/mtb-pattern-finder.git
cd task
python3 pattern_finder.py
```

The script reads `genome.txt` and `annotation.txt` from the working directory, prints progress to the console, and writes:
- `pattern_results.txt` - complete results (~7 MB, one line per match)
- A summary to the console

---

## Input Files

### genome.txt - FASTA format

This is the complete *M. tuberculosis* H37Rv genome sequence (NCBI accession NC\_000962.3), 4,411,709 nucleotides long. The file has a single header line starting with `>` and then the sequence in 70-character-wide wrapped lines.

### annotation.txt - GFF3 format

This is the standard GFF3 annotation file for the same genome. The relevant columns are:

| Column | Content |
|--------|---------|
| 3 | Feature type (gene, CDS, mRNA, exon, tRNA, ...) |
| 4 | Start coordinate (1-based, inclusive) |
| 5 | End coordinate (1-based, inclusive) |
| 7 | Strand (+ or -) |
| 9 | Semicolon-separated attributes, including `ID=...` and `Name=...` |

I only keep features that have an `ID` attribute, as instructed.

---

## The Patterns

All 24 patterns are 7 nucleotides long. The first two characters are `N` (wildcard), meaning any of A, T, G, or C. My script translates each `N` to `[ATGC]` and uses a compiled regex for the search.

```
NNAGAAG   NNAGAAT   NNAGAAA   NNGGAAG   NNAGAAC   NNGGAAA
NNAGCAT   NNAGGAG   NNAGGAT   NNAGCAA   NNGGAAC   NNGGAAT
NNAGCAG   NNAGGAA   NNAGGAC   NNGGGAG   NNGGGAT   NNGGGAA
NNAGCAC   NNGGGAC   NNGGCAT   NNGGCAG   NNGGCAA   NNGGCAC
```

---

## Results at a Glance

| Metric | Value |
|--------|-------|
| Patterns searched | 24 |
| Total matches found | 88,372 |
| Patterns with at least one match | 24 (100%) |
| Annotated regions with at least one hit | 12,220 / 12,265 (99.6%) |
| Most frequent pattern | NNGGCAC (8,333 matches) |
| Least frequent pattern | NNAGAAT (1,053 matches) |

### Pattern match counts

| Pattern | Matches |
|---------|---------|
| NNGGCAC | 8,333 |
| NNGGCAG | 7,352 |
| NNGGCAA | 6,621 |
| NNAGCAG | 6,102 |
| NNAGCAC | 6,060 |
| NNGGCAT | 5,954 |
| NNGGGAT | 4,043 |
| NNGGAAC | 3,798 |
| NNGGGAC | 3,485 |
| NNGGAAG | 3,211 |
| NNAGGAC | 3,018 |
| NNGGGAA | 2,897 |
| NNAGCAA | 2,849 |
| NNGGAAT | 2,830 |
| NNAGAAC | 2,738 |
| NNGGGAG | 2,727 |
| NNAGCAT | 2,576 |
| NNAGGAT | 2,481 |
| NNGGAAA | 2,359 |
| NNAGGAG | 2,263 |
| NNAGGAA | 2,222 |
| NNAGAAG | 2,167 |
| NNAGAAA | 1,233 |
| NNAGAAT | 1,053 |

---

## Key Design Decisions

**Regex over brute force.** For each pattern I compile a single regular expression (replacing `N` with `[ATGC]`) and use Python's `re.finditer`. This delegates matching to a C-level DFA implementation and runs the full genome scan in a few seconds per pattern.

**Load the whole genome into memory.** At 4.4 MB the sequence is tiny by current standards, so I read it once into a string rather than streaming. This avoids repeated disk I/O across 24 pattern searches.

**Both strands.** DNA is double-stranded. I search both the forward sequence and its reverse complement so that patterns occurring on the minus strand are not missed.

**1-based coordinates throughout.** GFF3 is 1-based. I convert Python's 0-based `re.match.start()` results immediately on output so that all reported positions are directly comparable to the annotation file.

**Overlap detection for annotated features.** For each match I do a linear scan over all annotated regions and report every region whose interval overlaps the 7-bp match window. Multiple overlapping annotations (e.g., gene, mRNA, CDS, exon all at the same locus) are all reported.

---

## Edge Cases I Handled

- **Reverse complement strand**: patterns that occur on the minus strand would be missed if I only search the forward sequence. I generate the reverse complement of the whole genome and search it as well.
- **Overlapping pattern matches**: `re.finditer` does not return overlapping matches by default. I use `re.finditer` with `pos` stepping to catch all overlapping occurrences.
- **Patterns spanning region boundaries**: a 7-bp hit can start in one annotated region and end in another. My overlap test (`not (end < region.start or start > region.end)`) catches this correctly.
- **Multiple annotations at the same position**: the GFF3 file has gene, mRNA, exon, and CDS records for the same locus. I report all of them.
- **Uppercase normalization**: I convert the genome sequence to uppercase on load to avoid any case-sensitivity issues.
- **Annotation lines without ID tags**: the task says to skip these, so I do.
- **Malformed GFF3 lines**: comment lines (`#`) and lines with fewer than 9 tab-separated fields are skipped gracefully.

---

## Output Format

`pattern_results.txt` contains one section per pattern:

```
Pattern: NNAGAAG
Total matches: 2167

Position    Sequence    Strand    Regions
580         GCAGAAG     +         RVBD0001; transcript:RVBD0001; CDS:RVBD0001
1911        TGAGAAG     +         none
3915        GAAGAAG     +         RVBD0003; transcript:RVBD0003; CDS:RVBD0003
```

---

## Writeup

See `analysis_writeup.md` for my full discussion of what was tricky, what decisions I made, and what I would improve with more time.

---

## Dependencies

The core analysis script (`pattern_finder.py`) uses only the Python standard library:

- `re` - regular expression matching
- `sys` - error handling and exit codes
- `collections.defaultdict` - grouping results by region
- `typing` - type annotations for readability
