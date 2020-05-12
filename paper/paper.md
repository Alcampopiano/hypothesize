---
title: 'Hypothesize: Robust Statistics for Python'
tags:
  - Python
  - R
  - statistics
  - statistical analysis
  - bootstrapping
  - trimmed mean
  - data analysis
  - data science
  - social science
  - hypothesis testing
authors:
  - name: Allan Campopiano
    orcid: 0000-0002-3280-4447
    affiliation: 1
  - name: Rand R. Wilcox
    orcid: 0000-0002-2524-2976
    affiliation: 2
    
affiliations:
  - name: Halton Catholic District School Board
    index: 1
  - name: University of Southern California
    index: 2
date: 08 May 2020
bibliography: paper.bib
---

# Summary

Hypothesize is a robust statistics library for Python that is used for comparing groups and
measuring associations. Robust methods, in particular those based on the trimmed mean [@20000755025] 
and/or bootstrapping [@bradley1993introduction], routinely outperform traditional statistical approaches 
in terms of power and accuracy. This is especially true when dealing with
distributions that produce outliers [@wilcox1998many; @wilcox2013introduction].

Hypothesize is based on Rand R. Wilcox's collection of [R functions](https://dornsife.usc.edu/labs/rwilcox/software/)
which contains hundreds of robust methods developed since the 1960's. 
Hypothesize brings many of these functions into the Python library ecosystem with the goal
of making robust computations as easy as possible for researchers. 
Hypothesize keeps the barrier-to-entry low, for example,

 - To easily incorporate Hypothesize with standard data processing tools
 [see @reback2020pandas; also @mckinney-proc-scipy-2010], all top-level 
 functions take a Pandas DataFrame/Series as input and return a Python Dictionary.
 
 - The API maps cleanly onto features of the user's statistical design. 
 This makes it easier to find and discover the set of appropriate functions for a
 given use case.
 
 - All top-level functions can be run directly in the browser alongside the documentation via 
[Google Colab Notebooks](https://colab.research.google.com/notebooks/intro.ipynb) 
(no local installation required).

# Acknowledgements

The authors would like to thank Lisa Collimore, James Desjardins, 
and the Halton Catholic District School Board for their support
of this project.

# References