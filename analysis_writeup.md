# Analysis Writeup

## Overview

I searched the complete *M. tuberculosis* H37Rv reference genome (4,411,709 nucleotides) for 24 specific 7-nucleotide patterns in which the first two characters are wildcards (`N` = A, T, G, or C). I found a total of 88,372 matches across all patterns and mapped each one to any overlapping annotated genomic feature from the GFF3 annotation file. Of the 12,265 annotated regions that have an ID tag, 12,220 (99.6%) contain at least one pattern match.

---

## What Was Tricky

**1. The double-stranded nature of DNA**

The most important thing I had to think about was that DNA is double-stranded. The genome file gives me the forward (plus) strand sequence, but genes and other features can be transcribed from either strand. If I only searched the forward strand I would miss every pattern occurrence that sits on the minus strand. I handled this by generating the reverse complement of the full genome sequence and running the same regex search on it. For minus-strand matches, I then convert the match position back to forward-strand 1-based coordinates so all results are in the same coordinate system as the GFF3 annotation file.

**2. Overlapping matches**

Python's `re.finditer` does not return overlapping matches by default. After finding a match ending at position `i`, it advances to position `i + 1` (past the end of the match). For 7-bp patterns this means consecutive matches shifted by just 1 or 2 bp could be silently dropped. I handle this by using a `while` loop with `re.search(sequence, pos)` and advancing `pos` by only 1 after each match (not by the full pattern length), which ensures every overlapping occurrence is captured.

**3. Coordinate system mismatch**

GFF3 annotation files use 1-based, inclusive coordinates. Python's `re` module returns 0-based match positions. I had to be careful to add 1 to every match start position before reporting it and before doing the interval overlap check against the annotation regions, otherwise I would report positions that are off by one relative to the reference coordinates.

**4. Multiple overlapping annotations at the same locus**

The annotation file has separate records for gene, mRNA, exon, and CDS at the same genomic locus. A single 7-bp pattern match can therefore overlap four or more distinct annotation entries simultaneously. I report all of them, because each layer of annotation contains different biological information (gene identity, transcript structure, coding sequence boundaries). The output file makes this clear with semicolon-separated lists of overlapping region IDs.

**5. Patterns spanning region boundaries**

A 7-bp match can start inside one annotated region and end inside a neighboring one, or start in an intergenic gap and end inside a gene. My overlap test (`not (query_end < region.start or query_start > region.end)`) correctly identifies partial overlaps, so these boundary-crossing matches are reported against all regions they touch.

---

## Decisions I Made

**Using compiled regular expressions.** I translate each `N` to `[ATGC]` and compile the pattern with `re.compile`. Python's regex engine is implemented in C and uses automata-based matching, which is much faster than anything I could write in pure Python. Searching 4.4 million nucleotides per pattern takes well under a second.

**Loading the full genome into memory.** At 4.4 MB the sequence is tiny. Loading it once into a Python string lets me do 24 regex searches (plus 24 more on the reverse complement) without any disk I/O overhead after the initial load. On any modern system this is a straightforward choice.

**Keeping all annotation feature types.** The task says to ignore entries without an ID tag, but not to filter by feature type. I keep gene, mRNA, exon, CDS, tRNA, rRNA, and ncRNA records because they each carry distinct biological meaning. A researcher looking at the results might care whether a match is inside a coding sequence versus just inside the broader gene boundary.

**Linear scan for region overlap detection.** For each of the 88,000+ matches I scan all ~12,000 annotation regions to find overlaps. This is O(matches x regions), roughly 1 billion simple comparisons total. Each comparison is just two integer subtractions, so in practice this runs in a few seconds. For this dataset size the simplicity is worth it.

---

## What I Would Do Differently With More Time

**Interval tree for overlap queries.** Replacing the linear scan with an interval tree (e.g. the `intervaltree` Python package, or a hand-rolled centered interval tree) would reduce each overlap query from O(regions) to O(log regions + hits). For this 12,000-region dataset the speedup would be modest, but it would matter a lot for larger eukaryotic genomes with millions of annotated features.

**Statistical enrichment testing.** The raw counts I report show where the patterns occur, but not whether they are significantly enriched or depleted in specific feature types (genes vs. intergenic regions, coding vs. regulatory). I would add a simple hypergeometric or chi-squared test comparing observed versus expected hit rates in each feature category.

**BED and bedGraph output.** Producing a BED file (one line per match, with chromosome, start, end, score, and strand) would let anyone load the results directly into a genome browser like IGV or UCSC. A bedGraph file of per-nucleotide match density would make it easy to visualize hotspots.

**Parallel processing.** The 24 pattern searches are completely independent of each other, so I could run them in parallel using `multiprocessing.Pool` or `concurrent.futures`. This would reduce the wall-clock time roughly proportionally to the number of CPU cores.

**Strand-specific region lookup.** Currently I report all overlapping regions regardless of strand. A more precise analysis would only report a region as a hit if the pattern strand matches the region strand (or if the user explicitly wants both). This matters for understanding whether the pattern is in the sense or antisense orientation relative to a gene.
