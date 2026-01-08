# A Program for Aligning Sentences in Bilingual Corpora

The accompanying program `align_regions.c` was originally published in the Appendix of [Gale & Church (1993)](http://aclweb.org/anthology/J/J93/J93-1004.pdf) in the journal Computational Linguistics, Volume 19, Number 1, March 1993.

A small number of very minor changes have been made to the program in order to enable compilation using modern C compilers:
* Two include statements were changed, replacing deprecated headers with their modern equivalents.
* One include statement ("unistd.h") was added; the original code assumed that the getopt function would be provided by the compiler, which is no longer necessarily the case.
* The original program also included code where a struct was passed as a parameter to a function where the function required a pointer to a struct; such code has been changed to explicitly take the address of the relevant struct object.
* Whitespace has been altered in some places to improve readability.

In order to enable the code to compile under C23 and later, the following changes have also been made:
* K&R [old-style function definitions](https://www.gnu.org/software/c-intro-and-ref/manual/html_node/Old_002dStyle-Function-Definitions.html) have been refactored to use the [function parameter variable](https://www.gnu.org/software/c-intro-and-ref/manual/html_node/Function-Parameter-Variables.html) declaration style that has been required since C89, and which is now enforced as of C23.

This repository also includes two small sample text files and a Makefile, demonstrating how to compile and run the program. The code was transcribed from the Appendix by Lane Schwartz.


The code was written by [William A. Gale](https://linguistlist.org/issues/13/13-2047.html) and [Kenneth W. Church](http://researcher.watson.ibm.com/researcher/view.php?person=us-kwchurch), with [Michael D. Riley](http://research.google.com/pubs/author125.html), and was presented with the following note:

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
