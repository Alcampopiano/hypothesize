# Hypothesize
![tests](https://github.com/Alcampopiano/hypothesize/workflows/tests/badge.svg)
[![PyPI version](https://img.shields.io/pypi/v/hypothesize?style=flat-square)](https://pypi.org/project/hypothesize/)
[![PyPI - Downloads](https://img.shields.io/pypi/dw/hypothesize?style=flat-square)](https://pypistats.org/packages/hypothesize)
[![license](https://img.shields.io/pypi/l/hypothesize?style=flat-square)](https://github.com/Alcampopiano/hypothesize/blob/master/LICENSE)
[![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/hypothesize_notebook_for_colab.ipynb)

<img src="https://github.com/Alcampopiano/hypothesize/blob/master/docs/docs/img/vp_inv.png?raw=true" alt="drawing" width="150"/>

A Python package for hypothesis testing using robust statistics.

:book: Please visit the [Hypothesize documentation site](https://Alcampopiano.github.io/hypothesize/) for more details.

## Basic Example
### A robust measure of association using winsorized correlation

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