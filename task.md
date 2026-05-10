# Rock Lab Coding Task

This is the original task description I received from Dr. DeJesus at the Rockefeller University Rock Lab.

## Intro

The task will be briefly described below. I will not describe the task in all possible detail or mention all possible issues, on purpose. I am also only describing the input files minimally, as I want to see if you're familiar with these kinds of files and know what information they contain and how you could use it. Basically, I want you to think critically about the task and be able to anticipate and address things that would be relevant for someone doing this type of task, or "edge-cases" that could arise.

When looking at your code I'll be mainly assessing these things (in order of importance):

1. **The different cases, potential issues, and details you considered, and how (and if) you handled them.** This is by far the most important one.
2. **How efficient the code is.** Primarily in execution time and time complexity. It doesn't have to be super-fast, but if it takes days to run then that's not great.
3. **How "nice" the code is.** This is vague but it includes things like readability, reusability, etc.

**Also, please do not use Claude, ChatGPT, GitHub co-pilot, or any other AI assistants to write the code for you or to think about the problem for you. You can use AI for looking up syntax or debugging an error if you really need to (like you would use Stack Overflow), but please do not discuss the task with it.  The whole point is to see your design decisions, and how you thought about the problem.**

## Submission

When you are done, please upload your source code, output files, and any other necessary files here: [link]

You have a little over 48 hours to finish the task.

## Task

I attached two files here. The file called `genome.txt` contains the DNA sequence of *M. tuberculosis*.

The file called `annotation.txt` is an annotation of *M. tuberculosis*, describing where regions of interest appear within the genome. This file is tab-separated. Columns 4 and 5 contain the start and end coordinates of these regions. Column 9 contains a semicolon-separated list of features/tags associated with these regions. The tag called `ID` (e.g. `ID=XXX`) would have the unique ID or the name for that region. You can ignore cases that don't have this ID tag. You can potentially use other tags you see, if you think they would be useful.

The task is to take these patterns below and find their location in the DNA sequence (provided in the genome file). The character `N` in the patterns stands for any character (i.e. N = A, T, G, or C).

### Patterns

```
NNAGAAG
NNAGAAT
NNAGAAA
NNGGAAG
NNAGAAC
NNGGAAA
NNAGCAT
NNAGGAG
NNAGGAT
NNAGCAA
NNGGAAC
NNGGAAT
NNAGCAG
NNAGGAA
NNAGGAC
NNGGGAG
NNGGGAT
NNGGGAA
NNAGCAC
NNGGGAC
NNGGCAT
NNGGCAG
NNGGCAA
NNGGCAC
```

In addition, you should also provide in what, if any, of the regions defined in the annotation file a pattern occurs in.

Finally, include a short writeup (a paragraph or two) describing what was tricky, what decisions you made, and what, if anything, you would do differently with more time.

Ideally you would write this in Python, but feel free to use any language you prefer.

---

## Checklist: How I addressed each requirement

| Requirement | Where addressed |
|-------------|-----------------|
| Find all 24 patterns in the genome | `pattern_finder.py` - `find_matches_on_strand()` |
| N = any nucleotide (A, T, G, C) | N replaced with `[ATGC]` regex character class |
| Report location of each pattern | `pattern_results.txt` - 1-based position per match |
| Report which annotated regions each match falls in | `find_overlapping_regions()` in `pattern_finder.py` |
| Ignore annotation entries without ID tag | Enforced in `load_annotations()` |
| Writeup: tricky parts, decisions, next steps | `analysis_writeup.md` |
| Python implementation | Yes, standard library only |
| Handle edge cases (strand, overlaps, boundaries) | See `analysis_writeup.md` and code comments |
