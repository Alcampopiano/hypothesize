# Hypothesize <a href="https://Alcampopiano.github.io/hypothesize/"><img align="right" src="https://github.com/Alcampopiano/hypothesize/blob/master/docs/docs/img/vp_inv.png" height="50"></img></a>

[![status](https://joss.theoj.org/papers/caf4095b3cdcc3adbb0252c995d59926/status.svg)](https://joss.theoj.org/papers/caf4095b3cdcc3adbb0252c995d59926)
![tests](https://github.com/Alcampopiano/hypothesize/workflows/tests/badge.svg)
[![PyPI version](https://img.shields.io/pypi/v/hypothesize?style=flat-square)](https://pypi.org/project/hypothesize/)
[![PyPI - Downloads](https://img.shields.io/pypi/dw/hypothesize?style=flat-square)](https://pypistats.org/packages/hypothesize)
[![license](https://img.shields.io/pypi/l/hypothesize?style=flat-square)](https://github.com/Alcampopiano/hypothesize/blob/master/LICENSE)

A Python package for hypothesis testing using robust statistics

## Basic Example

### A robust measure of association with winsorized correlation
[<img src="https://deepnote.com/buttons/launch-in-deepnote-white-small.svg">](https://deepnote.com/launch?name=wincor&url=https://github.com/Alcampopiano/hypothesize/blob/master/examples/wincor.ipynb
)

```python
from hypothesize.measuring_associations import wincor
from hypothesize.utilities import create_example_data

# creating an example DataFrame with columns "cell_1" and "cell_2"
df=create_example_data(2)

results=wincor(df.cell_1, df.cell_2)

# returning the correlation, number of observations, p-value, and winsorized covariance
print(results)
{'cor': 0.11, 'nval': 50, 'sig': 0.44, 'wcov': 0.01}
```

## Documentation
:book: Please visit the [Hypothesize documentation site](https://Alcampopiano.github.io/hypothesize/).
Note that each statistical test in the can be launched 
directly in [Deepnote's](deepnote.com) hosted notebook environment—complete with sample data 
(as shown in the example above 👆). 

## Citing Hypothesize

[![status](https://joss.theoj.org/papers/caf4095b3cdcc3adbb0252c995d59926/status.svg)](https://joss.theoj.org/papers/caf4095b3cdcc3adbb0252c995d59926)

If you use Hypothesize in academic work, please use the following citation:

Campopiano, A., & Wilcox, R. R. (2020). Hypothesize: Robust Statistics for Python. 
Journal of Open Source Software, 5(50), 2241, https://doi.org/10.21105/joss.02241

BibTex:

```bib
@article{Campopiano2020,
  doi = {10.21105/joss.02241},
  url = {https://doi.org/10.21105/joss.02241},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {50},
  pages = {2241},
  author = {Allan Campopiano and Rand R. Wilcox},
  title = {Hypothesize: Robust Statistics for Python},
  journal = {Journal of Open Source Software}
}
```
