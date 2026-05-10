#!/usr/bin/env python3
"""
M. tuberculosis Genome Pattern Finder
======================================

I search for 24 specific 7-nucleotide patterns in the M. tuberculosis H37Rv
genome (genome.txt) and report where each match falls relative to the annotated
genomic features in the GFF3 annotation file (annotation.txt).

The character 'N' in each pattern stands for any nucleotide (A, T, G, or C).
I convert N -> [ATGC] and use compiled regular expressions for efficient
genome-wide scanning.

Key design decisions I made:
  - Search both the forward strand AND the reverse complement so that patterns
    on the minus strand are not missed.
  - Use re.finditer with overlapping=True logic (manual pos stepping) so that
    overlapping matches are captured.
  - Load the full genome into memory once (4.4 MB is trivial) to avoid
    repeated disk reads across 24 pattern searches.
  - Report ALL annotated regions that overlap a match window, because a single
    genomic locus can be annotated as gene + mRNA + exon + CDS simultaneously.

Usage:
    python3 pattern_finder.py

Output:
    pattern_results.txt  -- complete tab-separated results for every match
    Console summary      -- per-pattern counts and top regions
"""

import re
import sys
from typing import List, Tuple, Dict, Optional
from collections import defaultdict


# ---------------------------------------------------------------------------
# Reverse complement helpers
# ---------------------------------------------------------------------------

# Complement mapping: each nucleotide maps to its Watson-Crick complement.
# I use a simple string translation table for speed.
_COMPLEMENT_TABLE = str.maketrans("ATGCatgc", "TACGtacg")


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence.

    For example:
        reverse_complement("ATGC") -> "GCAT"

    This is necessary because annotated features on the minus strand appear
    in the reverse complement orientation relative to the reference sequence.
    If I only search the forward strand I would miss all minus-strand pattern
    occurrences.
    """
    return sequence.translate(_COMPLEMENT_TABLE)[::-1]


# ---------------------------------------------------------------------------
# Data class: GenomicRegion
# ---------------------------------------------------------------------------

class GenomicRegion:
    """Holds the coordinates and metadata for one annotated genomic feature.

    In GFF3 format (the annotation file format used here), every feature has:
      - A start and end position (1-based, both inclusive)
      - A feature type (gene, CDS, mRNA, exon, tRNA, ...)
      - A strand ('+' or '-')
      - An attributes string in column 9, which is semicolon-separated and
        contains key=value pairs like ID=RVBD0001;Name=dnaA;...

    I parse the attributes string into a dictionary so I can easily look up
    the ID and Name for each region.
    """

    def __init__(
        self,
        start: int,
        end: int,
        attributes_string: str,
        feature_type: str = "",
        strand: str = ".",
    ):
        # GFF3 uses 1-based, inclusive coordinates.
        self.start = start
        self.end = end
        self.feature_type = feature_type
        self.strand = strand  # '+', '-', or '.'

        # Parse the raw attribute string into a dictionary for easy lookup.
        self.attributes = self._parse_attributes(attributes_string)

    def _parse_attributes(self, attr_string: str) -> Dict[str, str]:
        """Parse a GFF3 attribute string into a key->value dictionary.

        GFF3 attribute strings look like:
            ID=RVBD0001;Name=dnaA;biotype=protein_coding;Parent=...

        I split on ';' first, then on '=' for each key=value pair.
        If there is no '=' in a token I skip it (malformed).
        If attr_string is missing or '.', I return an empty dict.
        """
        attributes: Dict[str, str] = {}

        if not attr_string or attr_string == ".":
            return attributes

        for token in attr_string.split(";"):
            token = token.strip()
            if "=" in token:
                # Split on the first '=' only, in case the value itself
                # contains '=' characters.
                key, value = token.split("=", 1)
                attributes[key.strip()] = value.strip()

        return attributes

    def get_id(self) -> str:
        """Return the ID attribute value, or an empty string if absent."""
        return self.attributes.get("ID", "")

    def get_name(self) -> str:
        """Return the Name attribute value, or an empty string if absent."""
        return self.attributes.get("Name", "")

    def overlaps_with(self, query_start: int, query_end: int) -> bool:
        """Return True if [query_start, query_end] overlaps [self.start, self.end].

        Both intervals are 1-based and inclusive. Two intervals do NOT overlap
        only when one ends before the other starts, so the overlap condition is
        simply the negation of that:

            not (query_end < self.start or query_start > self.end)

        This correctly handles partial overlaps, full containment, and the case
        where the pattern spans a region boundary (starts inside one annotated
        feature and ends inside another).
        """
        return not (query_end < self.start or query_start > self.end)

    def __repr__(self) -> str:
        return (
            f"GenomicRegion({self.get_id()!r}, "
            f"{self.feature_type}, "
            f"{self.start}-{self.end}, strand={self.strand})"
        )


# ---------------------------------------------------------------------------
# Main analysis class: PatternFinder
# ---------------------------------------------------------------------------

class PatternFinder:
    """Loads the genome and annotations, then searches for all patterns.

    Workflow:
        1. load_genome()       -- read and validate the FASTA sequence
        2. load_annotations()  -- parse the GFF3 file into GenomicRegion objects
        3. analyze_patterns()  -- run regex search on both strands; map hits to
                                  annotated regions; collect statistics
        4. write_results()     -- write full results to pattern_results.txt
        5. print_summary()     -- print condensed stats to the console
    """

    def __init__(self):
        # The full genome sequence as a single uppercase string (forward strand).
        self.genome_sequence: str = ""

        # The reverse complement of the full genome (minus strand).
        self.genome_rc: str = ""

        # Length of the genome, stored once after loading.
        self.genome_length: int = 0

        # List of GenomicRegion objects parsed from the annotation file.
        # Only regions that have an ID attribute are kept (as per task spec).
        self.regions: List[GenomicRegion] = []

    # ------------------------------------------------------------------
    # Step 1: Load genome
    # ------------------------------------------------------------------

    def load_genome(self, filename: str) -> None:
        """Read the genome sequence from a FASTA file.

        FASTA format:
            >header line (one or more, starting with '>')
            SEQUENCE...  (may be split across many lines)

        I concatenate all non-header lines into a single string and convert
        to uppercase so that pattern matching is case-insensitive. After
        loading I compute the reverse complement so I have both strands ready.

        I also validate that the sequence contains only expected nucleotide
        characters (A, T, G, C, and sometimes N for unknown bases in the
        reference). Any unexpected characters get reported as a warning.
        """
        print(f"Loading genome from {filename}...")

        sequence_parts: List[str] = []

        with open(filename, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    # This is a FASTA header. The M. tuberculosis genome has
                    # only one chromosome so there will be exactly one header.
                    continue
                sequence_parts.append(line.upper())

        self.genome_sequence = "".join(sequence_parts)
        self.genome_length = len(self.genome_sequence)

        print(f"  Loaded {self.genome_length:,} nucleotides (forward strand)")

        # Validate: warn about any unexpected characters.
        # The standard nucleotides are A, T, G, C.
        # 'N' can appear in reference sequences to mark uncertain positions.
        allowed_chars = set("ATGCN")
        observed_chars = set(self.genome_sequence)
        unexpected = observed_chars - allowed_chars
        if unexpected:
            print(f"  Warning: unexpected characters in genome: {unexpected}")
            print("  These positions will not match any pattern (patterns require A/T/G/C).")

        # Build the reverse complement so I can search the minus strand too.
        # I store it as a string the same length as the forward strand.
        self.genome_rc = reverse_complement(self.genome_sequence)
        print("  Reverse complement built (for minus-strand search)")

    # ------------------------------------------------------------------
    # Step 2: Load annotations
    # ------------------------------------------------------------------

    def load_annotations(self, filename: str) -> None:
        """Parse the GFF3 annotation file into a list of GenomicRegion objects.

        GFF3 format (tab-separated, 9 columns):
            1. Sequence ID  (e.g. "Chromosome")
            2. Source       (e.g. "ena")
            3. Feature type (e.g. "gene", "CDS", "mRNA", "exon", "tRNA")
            4. Start        (1-based, inclusive)
            5. End          (1-based, inclusive)
            6. Score        ("." if unused)
            7. Strand       ("+" or "-")
            8. Phase        (0/1/2 for CDS, "." otherwise)
            9. Attributes   (semicolon-separated key=value pairs)

        Lines beginning with '#' are comments and are skipped.
        The task says to ignore features that don't have an ID tag in column 9.

        I keep ALL feature types (gene, CDS, mRNA, exon, tRNA, rRNA, ncRNA, ...)
        so that the output faithfully reflects what annotation level each match
        falls into.
        """
        print(f"Loading annotations from {filename}...")

        total_lines = 0
        skipped_no_id = 0
        skipped_malformed = 0

        with open(filename, "r") as fh:
            for line_number, line in enumerate(fh, start=1):
                line = line.strip()

                # Skip blank lines and comment / directive lines.
                if not line or line.startswith("#"):
                    continue

                fields = line.split("\t")

                # GFF3 requires exactly 9 tab-separated fields.
                if len(fields) < 9:
                    skipped_malformed += 1
                    continue

                total_lines += 1

                try:
                    feature_type = fields[2]
                    start = int(fields[3])   # column 4 in 1-based human numbering
                    end = int(fields[4])     # column 5
                    strand = fields[6]       # '+', '-', or '.'
                    attributes_string = fields[8]

                    # Sanity check: start should be <= end.
                    if start > end:
                        print(
                            f"  Warning: line {line_number} has start > end "
                            f"({start} > {end}), swapping."
                        )
                        start, end = end, start

                    region = GenomicRegion(
                        start=start,
                        end=end,
                        attributes_string=attributes_string,
                        feature_type=feature_type,
                        strand=strand,
                    )

                    # Per the task spec: skip regions that have no ID tag.
                    if not region.get_id():
                        skipped_no_id += 1
                        continue

                    self.regions.append(region)

                except (ValueError, IndexError) as exc:
                    skipped_malformed += 1
                    # Only print the first few warnings to avoid flooding the console.
                    if skipped_malformed <= 5:
                        print(f"  Warning: skipping malformed line {line_number}: {exc}")

        print(
            f"  Parsed {total_lines:,} feature lines; "
            f"kept {len(self.regions):,} regions with ID tags "
            f"(skipped {skipped_no_id:,} without ID, "
            f"{skipped_malformed:,} malformed)"
        )

    # ------------------------------------------------------------------
    # Step 3: Pattern matching
    # ------------------------------------------------------------------

    def pattern_to_regex(self, pattern: str) -> str:
        """Convert a pattern string to a regular expression string.

        The only wildcard character used in these patterns is 'N', which
        stands for any of the four standard nucleotides. I replace each 'N'
        with '[ATGC]' (a regex character class) so that the compiled regex
        matches exactly the right set of sequences.

        Example:
            "NNAGAAG" -> "[ATGC][ATGC]AGAAG"
        """
        return pattern.replace("N", "[ATGC]")

    def find_matches_on_strand(
        self,
        pattern: str,
        sequence: str,
        strand: str,
    ) -> List[Dict]:
        """Search 'sequence' for all occurrences of 'pattern' and return matches.

        Each match is a dict with:
            position    - 1-based start position on the FORWARD strand
            end_position - 1-based end position on the forward strand
            sequence    - the actual matched nucleotides
            strand      - '+' or '-'

        For minus-strand matches I search the reverse complement sequence
        (passed in as 'sequence') and then convert the position back to
        forward-strand 1-based coordinates.

        IMPORTANT - overlapping matches:
            Python's re.finditer skips overlapping matches by default because
            after finding a match it advances past the end of that match.
            These patterns are only 7 bp long and in practice overlapping
            matches matter (e.g. two consecutive matches shifted by 1 bp).
            I handle this by using a while loop with re.search starting from
            the character AFTER the previous match start (pos + 1), not
            after the match end. This captures all overlapping occurrences.
        """
        regex_string = self.pattern_to_regex(pattern)
        compiled = re.compile(regex_string)
        pattern_length = len(pattern)
        matches: List[Dict] = []

        pos = 0
        while pos <= len(sequence) - pattern_length:
            m = compiled.search(sequence, pos)
            if m is None:
                break  # No more matches anywhere in the remaining sequence.

            matched_seq = m.group()
            match_start_0 = m.start()  # 0-based index in 'sequence'

            if strand == "+":
                # Forward strand: convert directly to 1-based.
                fwd_start = match_start_0 + 1
                fwd_end = fwd_start + pattern_length - 1
                reported_seq = matched_seq

            else:
                # Minus strand: 'sequence' is the reverse complement.
                # The 0-based index in RC maps to the forward strand as:
                #   fwd_end   = genome_length - match_start_0       (1-based)
                #   fwd_start = fwd_end - pattern_length + 1
                # I report the sequence as it appears on the minus strand
                # (i.e. the reverse complement of the matched RC substring).
                fwd_end = self.genome_length - match_start_0
                fwd_start = fwd_end - pattern_length + 1
                reported_seq = reverse_complement(matched_seq)

            matches.append(
                {
                    "position": fwd_start,
                    "end_position": fwd_end,
                    "sequence": reported_seq,
                    "strand": strand,
                }
            )

            # Advance by 1 to allow overlapping matches.
            pos = match_start_0 + 1

        return matches

    def find_overlapping_regions(
        self, query_start: int, query_end: int
    ) -> List[GenomicRegion]:
        """Return all annotated regions that overlap [query_start, query_end].

        I do a simple linear scan. For a genome with ~12,000 annotated
        regions and ~88,000 total matches this is roughly 1 billion comparisons
        in the worst case, but each comparison is a single integer subtraction
        so in practice it runs in a few seconds.

        With more time I would replace this with an interval tree (e.g.
        the 'intervaltree' package) for O(log n + k) lookup instead of O(n).
        """
        overlapping: List[GenomicRegion] = []
        for region in self.regions:
            if region.overlaps_with(query_start, query_end):
                overlapping.append(region)
        return overlapping

    # ------------------------------------------------------------------
    # Step 4: Full analysis pipeline
    # ------------------------------------------------------------------

    def analyze_patterns(self, patterns: List[str]) -> Dict:
        """Run the complete analysis: search all patterns on both strands,
        then map every match to overlapping annotated regions.

        Returns a nested dict structured as:
            {
              "pattern_matches": {
                  "<PATTERN>": {
                      "count": int,
                      "matches": [
                          {
                              "position":     int,   # 1-based, forward strand
                              "end_position": int,
                              "sequence":     str,   # actual matched bases
                              "strand":       str,   # '+' or '-'
                              "regions":      [      # all overlapping annotations
                                  {
                                      "id":   str,
                                      "name": str,
                                      "type": str,
                                      "start": int,
                                      "end":   int,
                                      "strand": str,
                                  },
                                  ...
                              ],
                          },
                          ...
                      ],
                      "regions_hit": set(),  # set of region IDs with any hit
                  },
                  ...
              },
              "region_hits": defaultdict(list),  # region_id -> list of hits
              "summary": { ... },
            }
        """
        print(f"\nSearching for {len(patterns)} patterns on both strands...")

        results: Dict = {
            "pattern_matches": {},
            "region_hits": defaultdict(list),
            "summary": {},
        }

        total_matches = 0
        patterns_with_matches = 0

        for pattern in patterns:
            print(f"  Pattern {pattern} ...", end=" ", flush=True)

            # Search the forward strand (+) and the reverse complement (-).
            fwd_matches = self.find_matches_on_strand(
                pattern, self.genome_sequence, strand="+"
            )
            rc_matches = self.find_matches_on_strand(
                pattern, self.genome_rc, strand="-"
            )

            # Combine and sort by genomic position for cleaner output.
            all_matches = sorted(
                fwd_matches + rc_matches, key=lambda m: m["position"]
            )

            match_count = len(all_matches)
            print(f"{match_count:,} matches ({len(fwd_matches):,} fwd / {len(rc_matches):,} rc)")

            pattern_data: Dict = {
                "count": match_count,
                "matches": [],
                "regions_hit": set(),
            }

            # For each match, find which annotated regions it overlaps.
            for match in all_matches:
                start_pos = match["position"]
                end_pos = match["end_position"]

                overlapping = self.find_overlapping_regions(start_pos, end_pos)

                # Build a list of region info dicts for this match.
                region_info_list = []
                for region in overlapping:
                    rid = region.get_id()
                    region_info = {
                        "id": rid,
                        "name": region.get_name(),
                        "type": region.feature_type,
                        "start": region.start,
                        "end": region.end,
                        "strand": region.strand,
                    }
                    region_info_list.append(region_info)

                    # Track which region IDs have been hit.
                    pattern_data["regions_hit"].add(rid)

                    # Also accumulate hits per region across all patterns.
                    results["region_hits"][rid].append(
                        {
                            "pattern": pattern,
                            "position": start_pos,
                            "sequence": match["sequence"],
                            "strand": match["strand"],
                        }
                    )

                # Attach the region list to the match record.
                match["regions"] = region_info_list
                pattern_data["matches"].append(match)

            results["pattern_matches"][pattern] = pattern_data
            total_matches += match_count
            if match_count > 0:
                patterns_with_matches += 1

        # Compute summary statistics across all patterns.
        results["summary"] = {
            "total_patterns": len(patterns),
            "patterns_with_matches": patterns_with_matches,
            "total_matches": total_matches,
            "unique_regions_hit": len(results["region_hits"]),
            "total_annotated_regions": len(self.regions),
        }

        return results

    # ------------------------------------------------------------------
    # Step 5: Write results to file
    # ------------------------------------------------------------------

    def write_results(self, results: Dict, output_filename: str = "pattern_results.txt") -> None:
        """Write complete match results to a tab-separated text file.

        Each pattern gets its own section. For every match I write:
            position (1-based, forward strand)
            matched sequence
            strand ('+' or '-')
            semicolon-separated list of overlapping region IDs (or 'none')

        This file can be large (~7 MB for this dataset) because there are
        88,000+ matches.
        """
        print(f"\nWriting full results to {output_filename!r}...")

        with open(output_filename, "w") as fh:
            fh.write("M. tuberculosis Pattern Matching Results\n")
            fh.write("Generated by pattern_finder.py (Neha Rao, Rock Lab interview)\n")
            fh.write("=" * 60 + "\n\n")

            for pattern, data in results["pattern_matches"].items():
                fh.write(f"Pattern: {pattern}\n")
                fh.write(f"Total matches: {data['count']:,}\n\n")

                if data["count"] > 0:
                    # Tab-separated header row.
                    fh.write("Position\tSequence\tStrand\tRegions\n")
                    fh.write("-" * 50 + "\n")

                    for match in data["matches"]:
                        # Build region ID string for this match.
                        if match["regions"]:
                            region_ids = "; ".join(r["id"] for r in match["regions"] if r["id"])
                        else:
                            region_ids = "none"

                        fh.write(
                            f"{match['position']}\t"
                            f"{match['sequence']}\t"
                            f"{match['strand']}\t"
                            f"{region_ids}\n"
                        )

                fh.write("\n" + "-" * 50 + "\n\n")

        print(f"  Done. See {output_filename!r} for complete results.")

    # ------------------------------------------------------------------
    # Step 6: Print summary to console
    # ------------------------------------------------------------------

    def print_summary(self, results: Dict) -> None:
        """Print a concise summary of the analysis to the console."""
        summary = results["summary"]

        print("\n" + "=" * 60)
        print("ANALYSIS SUMMARY")
        print("=" * 60)
        print(f"  Total patterns searched : {summary['total_patterns']}")
        print(f"  Patterns with matches   : {summary['patterns_with_matches']}")
        print(f"  Total matches found     : {summary['total_matches']:,}")
        print(
            f"  Unique regions with hits: "
            f"{summary['unique_regions_hit']:,} / "
            f"{summary['total_annotated_regions']:,} annotated regions"
        )

        print("\nPer-pattern match counts (sorted by frequency):")
        sorted_patterns = sorted(
            results["pattern_matches"].items(),
            key=lambda x: x[1]["count"],
            reverse=True,
        )
        for pattern, data in sorted_patterns:
            bar = "#" * (data["count"] // 500)  # 1 '#' per 500 matches
            print(f"  {pattern}  {data['count']:>6,}  {bar}")

        # Show the top 10 regions by total hit count.
        print("\nTop 10 regions by number of pattern hits:")
        sorted_regions = sorted(
            results["region_hits"].items(),
            key=lambda x: len(x[1]),
            reverse=True,
        )[:10]

        for region_id, hits in sorted_regions:
            pattern_set = {h["pattern"] for h in hits}
            print(
                f"  {region_id:<20}  "
                f"{len(hits):>4} hits across "
                f"{len(pattern_set)} pattern(s)"
            )

        print("=" * 60)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Load inputs, run the analysis, and write results.

    I define the 24 target patterns here. They are all 7 nucleotides long with
    'N' wildcards in positions 1 and 2. The fixed core (positions 3-7) varies
    across the set, giving different AG/GG/AGG/GGG/AGC/GGC motifs.
    """

    # All 24 patterns from the task description.
    # N = any nucleotide (A, T, G, or C).
    patterns = [
        "NNAGAAG", "NNAGAAT", "NNAGAAA", "NNGGAAG", "NNAGAAC", "NNGGAAA",
        "NNAGCAT", "NNAGGAG", "NNAGGAT", "NNAGCAA", "NNGGAAC", "NNGGAAT",
        "NNAGCAG", "NNAGGAA", "NNAGGAC", "NNGGGAG", "NNGGGAT", "NNGGGAA",
        "NNAGCAC", "NNGGGAC", "NNGGCAT", "NNGGCAG", "NNGGCAA", "NNGGCAC",
    ]

    finder = PatternFinder()

    try:
        # Step 1: Load the genome sequence from the FASTA file.
        finder.load_genome("genome.txt")

        # Step 2: Load annotated regions from the GFF3 file.
        finder.load_annotations("annotation.txt")

        # Step 3 + 4: Search for all patterns on both strands and map hits
        # to annotated regions.
        results = finder.analyze_patterns(patterns)

        # Step 5: Write every match to a file.
        finder.write_results(results, "pattern_results.txt")

        # Step 6: Print a human-readable summary.
        finder.print_summary(results)

        print("\nAnalysis complete.")

    except FileNotFoundError as exc:
        print(f"\nError: could not open input file: {exc}")
        print("Make sure genome.txt and annotation.txt are in the current directory.")
        sys.exit(1)

    except Exception as exc:
        print(f"\nUnexpected error during analysis: {exc}")
        raise  # Re-raise so the full traceback is visible for debugging.


if __name__ == "__main__":
    main()
