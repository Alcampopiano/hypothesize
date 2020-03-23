# Hypothesize
[![PyPI version](https://img.shields.io/pypi/v/hypothesize?style=flat-square)](https://pypi.org/project/hypothesize/)
[![PyPI - Downloads](https://img.shields.io/pypi/dw/hypothesize?style=flat-square)](https://pypi.org/project/hypothesize/)
[![license](https://img.shields.io/pypi/l/hypothesize?style=flat-square)](https://github.com/Alcampopiano/hypothesize/blob/master/LICENSE)
[![Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/hypothesize_notebook_for_colab.ipynb)

A Python package for comparing groups and measuring associations using robust statistics.

This package is a port of Rand R. Wilcox's R library [WRS](https://dornsife.usc.edu/labs/rwilcox/software/). The functions in this repository, as well as the issues of robustness, are described in his book "[Introduction to Robust Estimation and Hypothesis Testing](https://play.google.com/store/books/details?id=8f8nBb4__EYC&gl=ca&hl=en-CA&source=productsearch&utm_source=HA_Desktop_US&utm_medium=SEM&utm_campaign=PLA&pcampaignid=MKTAD0930BO1&gclid=CjwKCAiA44LzBRB-EiwA-jJipJzyqx9kwNMq5MMU7fG2RrwBK9F7sirX4pfhS8wO7k9Uz_Sqf2P28BoCYzcQAvD_BwE&gclsrc=aw.ds)".

Please visit the [Hypothesize documentation site](https://Alcampopiano.github.io/hypothesize/).

### :warning: This repository is still in the early stages of development

## Installation
The Hypothesize package is available on [PyPI](https://pypi.org/project/hypothesize/). To install it, simply type:

```python
pip install hypothesize

```

## Expand the following topics to see examples

<details>
<summary><strong> <font size="3">How to compare two groups </font></strong></summary>
    
#### Load data from a CSV
    
```python
import pandas as pd

df=pd.read_csv("/home/allan/two_groups_data.csv")

df.head()
```
    
    
|    |   Group_1 |   Group_2 |
|---:|----------:|----------:|
|  0 | 0.0446518 |  0.90675  |
|  1 | 0.763458  |  0.291555 |
|  2 | 0.71039   |  0.59828  |
|  3 | 0.175208  |  0.268073 |
|  4 | 0.957819  |  0.222688 |
    
    
#### Import the desired function and pass in the data for each group
- This example uses the bootstrapped-t method with 20% trimmed means
- The output is a dictionary containing the results (95% confidence interval, p_value, test statistics, etc...)
    

```python
from hypothesize.compare_groups_with_single_factor import yuenbt

results=yuenbt(df.Group_1, df.Group_2)

print(results['ci'])
```
 
<p>

[-0.3115715617702292, 0.10636703554225341]

<p>
    
</details>

<details>
 <summary> <strong> <font size="3">How to compare groups in a factorial design</font></strong></summary>
    
#### Load data from a CSV
    
```python
import pandas as pd

df=pd.read_csv("/home/allan/two_way_data.csv")

df.head() 
```

|    |   cell_1_1 |   cell_1_2 |   cell_1_3 |   cell_2_1 |   cell_2_2 |   cell_2_3 |
|---:|-----------:|-----------:|-----------:|-----------:|-----------:|-----------:|
|  0 |  0.0446518 |   0.90675  |   0.795696 |  0.519486  |   0.333636 |  0.232153  |
|  1 |  0.763458  |   0.291555 |   0.84158  |  0.0339891 |   0.511235 |  0.732503  |
|  2 |  0.71039   |   0.59828  |   0.110407 |  0.898072  |   0.769496 |  0.0484005 |
|  3 |  0.175208  |   0.268073 |   0.888728 |  0.287442  |   0.100153 |  0.210394  |
|  4 |  0.957819  |   0.222688 |   0.834161 |  0.599158  |   0.655308 |  0.203486  |
    
#### Import the desired function and pass in the data
- This example uses a 2-by-3 design
- One approach is to use a set of linear contrasts that will test all main effects and interactions
- Then, the bootstrap-t method and the 20% trimmed mean can be used
- CIs are adjusted to control for FWE for each family of tests (factor A, factor B, and the interactions)
- All pairwise contrasts are created internally using the `con2way` function
- The results are a dictionary of DataFrames that contain various statistics for each factor and the interactions
    
```python
from hypothesize.compare_groups_with_two_factors import bwmcp

results=bwmcp(J=2, K=3, x=df)
```
<p>

```python
results['factor_A']  
    
```
<p>    
    
|    |   con_num |    psihat |       se |     test |   crit_value |   p_value |
|---:|----------:|----------:|---------:|---------:|-------------:|----------:|
|  0 |         0 | 0.0393584 | 0.169849 | 0.231726 |      3.35959 |  0.941569 |   
    
    
<p>

```python
results['factor_B']  
    
```
<p>    
    
|    |   con_num |     psihat |       se |       test |   crit_value |   p_value |
|---:|----------:|-----------:|---------:|-----------:|-------------:|----------:|
|  0 |         0 | -0.104506  | 0.126135 | -0.828529  |       2.4329 |  0.452421 |
|  1 |         1 | -0.0931364 | 0.151841 | -0.613382  |       2.4329 |  0.552588 |
|  2 |         2 |  0.01137   | 0.135392 |  0.0839783 |       2.4329 |  0.923205 |
    
<p>

```python
results['factor_AB']  
    
```
<p>    
    
|    |   con_num |     psihat |       se |      test |   crit_value |   p_value |
|---:|----------:|-----------:|---------:|----------:|-------------:|----------:|
|  0 |         0 | -0.100698  | 0.126135 | -0.798336 |       2.3771 |  0.410684 |
|  1 |         1 | -0.037972  | 0.151841 | -0.250078 |       2.3771 |  0.804674 |
|  2 |         2 |  0.0627261 | 0.135392 |  0.463291 |       2.3771 |  0.659432 |
    
    
    
</details>

<details>
 <summary> <strong> <font size="3">How to compute a robust correlation</font></strong></summary>
    
#### Load data from a CSV
    
```python
import pandas as pd

df=pd.read_csv("/home/allan/two_groups_data.csv")

df.head() 
```
|    |   Group_1 |   Group_2 |
|---:|----------:|----------:|
|  0 | 0.0446518 |  0.90675  |
|  1 | 0.763458  |  0.291555 |
|  2 | 0.71039   |  0.59828  |
|  3 | 0.175208  |  0.268073 |
|  4 | 0.957819  |  0.222688 |
    
    
#### Import the desired function and pass in the data for each group
- One approach is to winsorize the x and y data
- A heteroscedastic method for testing zero correlation is also provided in this package but not shown here 
 - Please see the function `corb` which uses the percentile bootstrap to compute a 1-alpha CI and p_value for any correlation   
- The output is a dictionary containing various statistics (the winsorized correlation, winsorized covariance, etc...)

```python
from hypothesize.measuring_associations import wincor

results=wincor(df.Group_1, df.Group_2)

print(results['wcor'])
```

<p>
    
-0.05690314435050796
    
</details>
 
<p>
    
---
    
Note the following:
- The amount of trimming, alpha, and other common parameters can be explicitly specified, otherwise defaults are used implicitly (20% trimming, alpha=.05, etc...)
- The API is subject to change as the package is still in the early stages of development