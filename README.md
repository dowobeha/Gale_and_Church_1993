# A Program for Aligning Sentences in Bilingual Corpora

The accompanying program `align_regions.c` was originally published in the Appendix of [Gale & Church (1993)](http://aclweb.org/anthology/J/J93/J93-1004.pdf) in the journal Computational Linguistics, Volume 19, Number 1, March 1993.

A small number of very minor changes have been made to the program in order to enable compilation using modern C compilers:
* Two include statements were changed, replacing deprecated headers with their modern equivalents.* The original program also included code where a struct was passed as a parameter to a function where the function required a pointer to a struct; such code has been changed to explicitly take the address of the relevant struct object.
* Whitespace has been altered in some places to improve readability.

This repository also includes two small sample text files and a Makefile, demonstrating how to compile and run the program. The code was transcribed from the Appendix by Lane Schwartz.


The code was written by William A. Gale and Kenneth W. Church, with Michael D. Riley, and was presented with the following note:

> The following code is the core of align. It is a C language program that inputs two 
> text files, with one token (word) per line. The text files contain a number of delimiter
> tokens. There are two types of delimiter tokens: "hard" and "soft." The hard regions
> (e.g., paragraphs) may not be changed, and there must be equal numbers of them in
> the two input files. The soft regions (e.g., sentences) may be deleted (1-0), inserted (0-
> 1), substituted (1-1), contracted (2-1), expanded (1-2), or merged (2-2) as necessary so
> that the output ends up with the same number of soft regions. The program generates
> two output files. The two output files contain an equal number of soft regions, each
> on a line. If the -v command line option is included, each soft region is preceded by
> its probability score.
