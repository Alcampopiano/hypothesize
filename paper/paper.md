---
title: 'Hypothesize: Robust Statistics for Python'
tags:
  - Python
  - statistics
  - R
  - bootstrapping
  - trimmed mean
  - data science
authors:
  - name: Allan Campopiano
    orcid: 0000-0002-9623-3401
    affiliation: 1
  - name: Rand R. Wilcox
    orcid: 0000-0002-5223-6168
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
measuring associations. Robust methods, in particular those based on the trimmed mean [@tuk:1962] 
and bootstrapping [@efron:1993], routinely outperform traditional statistical approaches 
in terms of power and accuracy. This is especially true when dealing with
distributions that produce outliers [@wil:1998; @wil:2012].

Hypothesize is based on Rand R. Wilcox's collection of [R functions](https://dornsife.usc.edu/labs/rwilcox/software/)
which contains hundreds of robust methods developed since 1960's. 
Hypothesize brings many of these functions into the Python library ecosystem with the goal
of making robust computations as easy as possible for researchers. 
Examplesof how Hypothesize keeps the barrier-to-entry low are as follows:

 - To easily incorporate Hypothesize with standard data processing tools,
 all top-level functions take a Pandas DataFrame/Series as input and return a Python Dictionary.
 
 - The API maps cleanly onto to important features of the user's statistical design. 
 This promotes discoverability since for any given design there are often several 
 functions from which to choose.
 
 - All top-level functions can be run directly in the browser alongside the documentation via 
[Google Colab Notebooks](https://colab.research.google.com/notebooks/intro.ipynb) 
(no installation required).

# Acknowledgements

# References