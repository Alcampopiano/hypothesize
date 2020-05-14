# Function Reference

Hypothesize exposes the following top-level functions for comparing 
groups and measuring associations. The function names, code, and descriptions
are kept generally consistent with Wilcox's WRS package. 
If you want to learn more about the theory and research
behind any given function here, see Wilcox's books, especially
[Introduction to Robust Estimation and Hypothesis Testing](https://play.google.com/store/books/details?id=8f8nBb4__EYC&gl=ca&hl=en-CA&source=productsearch&utm_source=HA_Desktop_US&utm_medium=SEM&utm_campaign=PLA&pcampaignid=MKTAD0930BO1&gclid=CjwKCAiA44LzBRB-EiwA-jJipJzyqx9kwNMq5MMU7fG2RrwBK9F7sirX4pfhS8wO7k9Uz_Sqf2P28BoCYzcQAvD_BwE&gclsrc=aw.ds).

---
Jump to:

**[Comparing groups with a single factor](function_guide.md#comparing-groups-with-a-single-factor)**

**[Comparing groups with two factors](function_guide.md#comparing-groups-with-two-factors)**

**[Measuring associations](function_guide.md#measuring-associations)**

**[Other important functions](function_guide.md#other-important-functions)**

---

## Comparing groups with a single factor
Statistical tests analogous to a 1-way ANOVA or T-tests. 
That is, group tests that have a single factor.
---

## Independent groups

### l2drmci

`#!py l2drmci(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

Compute a bootstrap confidence interval for a
measure of location associated with the distribution of x-y. 
That is, compare x and y by looking at all possible difference scores
in random samples of `x` and `y`. `x` and `y` are possibly dependent.

_Parameters:_

**x: Pandas Series**

Data for group one

**y: Pandas Series**

Data for group two

**est: function**

Measure of location (currently only `trim_mean` is supported)

***args: list/value**

Parameter(s) for measure of location (e.g., .2)

**pairwise_drop_na: bool**

If True, treat data as dependent and remove any row with missing data. If False,
remove missing data for each group seperately (cannot deal with unequal sample sizes)

**alpha: float**

Alpha level (default is .05)

**nboot: int**

Number of bootstrap samples (default is 2000)

**seed: bool**

Random seed for reprodicible results (default is `False`)

_Returns:_

Dictionary of results

**ci: list** 

Confidence interval

**p_value: float** 

p-value


<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/l2drmci.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### linconb

`#!py linconb(x, con, tr=.2, alpha=.05, nboot=599, seed=False)`

Compute a 1-alpha confidence interval for a set of d linear contrasts
involving trimmed means using the bootstrap-t bootstrap method.
Independent groups are assumed. CIs are adjusted to control FWE 
(p values are not adjusted).

_Parameters:_

**x: DataFrame**

Each column represents a group of data

**con: array**

`con` is a J (number of columns) by d (number of contrasts)
matrix containing the contrast coefficents of interest.
All linear constrasts can be created automatically by using the function [con1way](J)
(the result of which can be used for `con`).

**tr: float**

Proportion to trim (default is .2)

**alpha: float**

Alpha level (default is .05)

**nboot: int**

Number of bootstrap samples (default is 2000)

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictionary of results

**con: array**

Contrast matrix

**crit: float**

Critical value

**n: list**

Number of observations for each group

**psihat: DataFrame**

Difference score and CI for each contrast

**test: DataFrame**

Test statistic, standard error, and p-value for each contrast

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/linconb.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### pb2gen

`#!py pb2gen(x, y, est, *args, alpha=.05, nboot=2000, seed=False)`

Compute a bootstrap confidence interval for the
the difference between any two parameters corresponding to two
independent groups.

_Parameters:_

**x: Pandas Series**

Data for group one

**y: Pandas Series**

Data for group two

**est: function**

Measure of location (currently only `trim_mean` is supported)

***args: list/value**

Parameter(s) for measure of location (e.g., .2)

**alpha: float**

Alpha level (default is .05)

**nboot: int**

Number of bootstrap samples (default is 2000)

**seed: bool**

Random seed for reprodicible results (default is `False`)

_Return:_

Dictionary of results

**ci: list** 

Confidence interval

**est_1: float**

Estimated value (based on `est`) for group one

**est_2: float**

Estimated value (based on `est`) for group two

**est_dif: float**

Estimated difference between group one and two

**n1: int**

Number of observations in group one

**n2: int**

Number of observations in group two

**p_value: float** 

p-value

**variance: float**

Variance of group one and two

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/pb2gen.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### tmcppb

`#!py tmcppb(x, est, *args, con=None, bhop=False, alpha=.05, nboot=None, seed=False)`

Multiple comparisons for J independent groups using trimmed means and
the percentile bootstrap method. Rom’s method is used to control the 
probability of one or more type I errors. For C > 10 hypotheses, 
or when the goal is to test at some level other than .05 and .01, 
Hochberg’s method is used. Setting the argument `bhop` to `True` uses the
Benjamini–Hochberg method instead.

_Parameters:_

**x: Pandas DataFrame**

Each column represents a group of data

**est: function**

Measure of location (currently only `trim_mean` is supported)

***args: list/value**

Parameter(s) for measure of location (e.g., .2)

**con: array**

`con` is a J (number of columns) by d (number of contrasts)
matrix containing the contrast coefficents of interest.
All linear constrasts can be created automatically by using the function [con1way](J)
(the result of which can be used for `con`). The default is `None` and in this
case all linear contrasts are created automatically.
 
**bhop: bool**

If `True`, the Benjamini–Hochberg method is used to control FWE

**alpha: float**

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples (default is 2000)

**seed: bool**

Random seed for reproducible results. Default is `False`.

_Return:_

Dictionary of results

**con: array**

Contrast matrix

**num_sig: int**

Number of statistically significant results

**output: DataFrame**

Difference score, p-value, critical value, and CI for each contrast

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/tmcppb.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### yuenbt

`#!py yuenbt(x, y, tr=.2, alpha=.05, nboot=599, seed=False)`

Compute a 1-alpha confidence interval for the difference between
the trimmed means corresponding to two independent groups.
The bootstrap-t method is used. During the bootstrapping, 
the absolute value of the test statistic is used (the "two-sided method").

_Parameters:_

**x: Pandas Series**

Data for group one

**y: Pandas Series**

Data for group two

**tr: float**

Proportion to trim (default is .2)

**alpha: float**

Alpha level (default is .05)

**nboot: int**

Number of bootstrap samples (default is 599)

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictionary of results

**ci: list** 

Confidence interval

**est_dif: float**

Estimated difference between group one and two

**est_1: float**

Estimated value (based on `est`) for group one

**est_2: float**

Estimated value (based on `est`) for group two

**p_value: float** 

p-value

**test_stat: float**

Test statistic

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/yuenbt.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

## Dependent groups

### bootdpci

`#!py bootdpci(x, est, *args, nboot=None, alpha=.05, dif=True, BA=False, SR=False)`

Use percentile bootstrap method, compute a .95 confidence interval 
for the difference between a measure of location or scale 
when comparing two dependent groups.

The argument `dif` defaults to `True` indicating 
that difference scores will be used, in which case Hochberg’s 
method is used to control FWE. If `dif` is `False`, measures of 
location associated with the marginal distributions are used 
instead. 

If `dif` is `False` and `BA` is `True`, the bias adjusted 
estimate of the generalized p-value is recommended.
Using `BA`=`True` (when `dif`=`False`) 
is recommended when comparing groups 
with M-estimators and MOM, but it is not necessary when 
comparing 20% trimmed means (Wilcox & Keselman, 2002). 

The so-called the SR method, which is a slight 
modification of Hochberg's (1988) "sequentially rejective" 
method can be applied to control FWE, especially when 
comparing one-step M-estimators or M-estimators.

_Parameters:_

**x: Pandas DataFrame**

Each column represents a group of data

**est: function**

Measure of location (currently only `trim_mean` is supported)

***args: list/value**

Parameter(s) for measure of location (e.g., .2)

**alpha: float**

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is `None`
in which case `nboot` will be chosen for you 
based on the number of contrasts.

**dif: bool**

When `True`, use difference scores, otherwise use marginal distributions

**BA: bool**

When `True`, use the bias adjusted estimate of the 
generalized p-value is applied (e.g., when `dif` is `False`)

**SR: bool**

When `True`, use the modified "sequentially rejective", especially when 
comparing one-step M-estimators or M-estimators

_Return:_

Dictionary of results

**con: array**

Contrast matrix

**num_sig: int**

Number of statistically significant results

**output: DataFrame**

Difference score, p-value, critical value, and CI for each contrast

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/bootdpci.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### rmmcppb

`#!py rmmcppb(x, est, *args,  alpha=.05, con=None, dif=True, nboot=None, BA=False,
    hoch=False, SR=False, seed=False)`

Use a percentile bootstrap method to compare dependent groups.
By default, compute a .95 confidence interval for all linear contrasts
specified by con, a J-by-C matrix, where C is the number of
contrasts to be tested, and the columns of `con` are the
contrast coefficients. If con is not specified, 
all pairwise comparisons are done.

If `est` is the function `onestep` or `mom` (these are not implemeted yet),
method SR can be used to control the probability of at least one Type I error.
Otherwise, Hochberg's method is used.

If `dif` is `False` and `BA` is `True`, the bias adjusted 
estimate of the generalized p-value is recommended.
Using `BA`=`True` (when `dif`=`False`) 
is recommended when comparing groups 
with M-estimators and MOM, but it is not necessary when 
comparing 20% trimmed means (Wilcox & Keselman, 2002). 

Hochberg's sequentially rejective method can be used and is used 
if n>=80.

_Parameters:_

**x: Pandas DataFrame**

Each column represents a group of data

**est: function**

Measure of location (currently only `trim_mean` is supported)

***args: list/value**

Parameter(s) for measure of location (e.g., .2)

**alpha: float**

Alpha level (default is .05)

**con: array**

`con` is a J (number of columns) by d (number of contrasts)
matrix containing the contrast coefficents of interest.
All linear constrasts can be created automatically by using the function [con1way](J)
(the result of which can be used for `con`). The default is `None` and in this
case all linear contrasts are created automatically.
 
**dif: bool**

When `True`, use difference scores, otherwise use marginal distributions

**nboot: int**

Number of bootstrap samples. Default is `None`
in which case `nboot` will be chosen for you 
based on the number of contrasts.

**BA: bool**

When `True`, use the bias adjusted estimate of the 
generalized p-value is applied (e.g., when `dif` is `False`)

**hoch: bool**

When `True`, Hochberg's sequentially rejective method can be used and is used 
if n>=80.

**SR: bool**

When `True`, use the modified "sequentially rejective", especially when 
comparing one-step M-estimators or M-estimators.

**seed: bool**

Random seed for reprodicible results (default is `False`)

_Return:_

Dictionary of results

**con: array**

Contrast matrix

**num_sig: int**

Number of statistically significant results

**output: DataFrame**

Difference score, p-value, critical value, and CI for each contrast

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/rmmcppb.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### lindepbt

`#!py lindepbt(x, tr=.2, con=None, alpha=.05, nboot=599, dif=True, seed=False)`

Multiple comparisons on trimmed means with FWE controlled with Rom's method
Using a bootstrap-t method.


_Parameters:_

**x: Pandas DataFrame**

Each column in the data represents a different group 

**tr: float**

Proportion to trim (default is .2)

**con: array**

`con` is a J (number of groups) by d (number of contrasts) 
matrix containing the contrast coefficents of interest.
All linear constrasts can be created automatically by using the function [con1way](J)
(the result of which can be used for `con`). The default is `None` and in this 
case all linear contrasts are created automatically.

**alpha: float**

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples (default is 2000)

**dif: bool**

When `True`, use difference scores, otherwise use marginal distributions

**seed: bool**

Random seed for reprodicible results (default is `False`)

_Return:_

Dictionary of results

**con: array**

Contrast matrix

**num_sig: int**

Number of observations for each group

**psihat: DataFrame**

Difference score and CI for each contrast

**test: DataFrame**

Test statistic, p-value, critical value, and standard error
for each contrast

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/lindepbt.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### ydbt

`#!py ydbt(x, y, tr=.2, alpha=.05, nboot=599, side=True, seed=False)`

Using the bootstrap-t method,
compute a .95 confidence interval for the difference between
the marginal trimmed means of paired data.
By default, 20% trimming is used with 599 bootstrap samples.

_Parameters:_

**x: Pandas Series**

Data for group one

**y: Pandas Series**

Data for group two

**tr: float**

Proportion to trim (default is .2)

**alpha: float**

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples (default is 2000)

**side: bool**
When `True` the function returns a symmetric CI and a p value, 
otherwise the function returns equal-tailed CI (no p value)

**seed: bool**

Random seed for reprodicible results (default is `False`)

_Return:_

Dictionary of results

**ci: list**

Confidence interval

**dif: float**

Difference between group one and two

**p_value: float**

p-value

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/ydbt.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

## Comparing groups with two factors
Statistical tests analogous to a 2-way ANOVA. 
That is, group tests that have a two factors.
---

## Dependent groups

### wwmcppb

`#!py wwmcppb(J, K, x,  est, *args,  alpha=.05, dif=True,
            nboot=None, BA=True, hoch=True, seed=False)`

Do all multiple comparisons for a within-by-within design 
using a percentile bootstrap method.A sequentially rejective 
method is used to control alpha.
Hochberg's method can be used and is if n>=80.

_Parameters:_

**J: int**

Number of J levels associated with Factor A

**K: int**

Number of K levels associated with Factor B

**x: Pandas DataFrame**

Each column represents a cell in the factorial design. For example,
a 2x3 design would correspond to a DataFrame with 6 columns 
(levels of Factor A x levels of Factor B).

Order your columns according to the following pattern
(traversing each row in a matrix): 

 - the first column contains data for level 1 of Factor A
 and level 1 of Factor B
 
 - the second column contains data for level 1 of Factor A
 and level 2 of Factor B
 
 - column `K` contains the data for level 1 of Factor A 
 and level `K` of Factor B
 
 - column `K` + 1 contains the data for level 2 of Factor A
 and level 1 of Factor B
 
 - and so on ... 
 
**est: function**

Measure of location (currently only `trim_mean` is supported)
 
***args: list/value**

Parameter(s) for measure of location (e.g., .2)

**alpha: float**

Alpha level (default is .05)

**dif: bool**

When `True`, use difference scores, otherwise use marginal distributions

**nboot: int**

Number of bootstrap samples (default is 599)

**BA: bool**

When `True`, use the bias adjusted estimate of the 
generalized p-value is applied (e.g., when `dif` is `False`)

**hoch: bool**

When `True`, Hochberg's sequentially 
rejective method can be used to control FWE

**seed: bool**

Random seed for reprodicible results (default is `False`)

_Return:_

Dictionary of results

The following results are returned for factor A, factor B, 
and the interaction. See the keys `'factor_A'`, `'factor_A'`, and `'factor_AB'`, 
respectively.

**con: array**

Contrast matrix

**num_sig: int**

Number of statistically significant results

**output: DataFrame**

Difference score, p-value, critical value, and CI for each contrast

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/wwmcppb.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### wwmcpbt

`#!py wwmcpbt(J, K, x, tr=.2, alpha=.05, nboot=599, seed=False)`

Do multiple comparisons for a within-by-within design.
using a bootstrap-t method and trimmed means.
All linear contrasts relevant to main effects and interactions
are tested. With trimmed means FWE is
controlled with Rom's method.
    
_Parameters:_

**J: int**

Number of J levels associated with Factor A

**K: int**

Number of K levels associated with Factor B

**x: Pandas DataFrame**

Each column represents a cell in the factorial design. For example,
a 2x3 design would correspond to a DataFrame with 6 columns 
(levels of Factor A x levels of Factor B).

Order your columns according to the following pattern
(traversing each row in a matrix): 

 - the first column contains data for level 1 of Factor A
 and level 1 of Factor B
 
 - the second column contains data for level 1 of Factor A
 and level 2 of Factor B
 
 - column `K` contains the data for level 1 of Factor A 
 and level `K` of Factor B
 
 - column `K` + 1 contains the data for level 2 of Factor A
 and level 1 of Factor B
 
 - and so on ... 

**tr: float**

Proportion to trim (default is .2)

**alpha: float**

Alpha level (default is .05)

**nboot: int**

Number of bootstrap samples (default is 599)

**seed: bool**

Random seed for reprodicible results (default is `False`)

_Return:_

Dictionary of results

The following results are returned for factor A, factor B, 
and the interaction. See the keys `'factor_A'`, `'factor_A'`, and `'factor_AB'`, 
respectively.

**con: array**

Contrast matrix

**num_sig: int**

Number of statistically significant results

**psihat: DataFrame**

Difference score and CI for each contrast

**test: DataFrame**

Test statistic, p-value, critical value, and standard error for each contrast

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/wwmcpbt.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

## Mixed designs

These designs are also known as "split-plot" or "between-within" desgins.
Hypothesize follws the common convention that assumes that the between-subjects
factor is factor A and the within-subjects conditions are Factor B. For example,
in a 2x3 mixed design, factor A has two levels. For each of these levels,
there are 3 within-subjects conditions.

!!! danger "Make sure your DataFrame corresponds to a between-within design, not the other way around"
    For example, in a 2x3 mixed design, the first 3 columns correspond to the first level of Factor A. The last 3
    columns correspond to the second level of factor A.

### bwamcp

`#!py bwamcp(J, K, x, tr=.2, alpha=.05, pool=False)`

All pairwise comparisons among levels of Factor A
in a mixed design using trimmed means. The `pool`
option allows you to pool dependent groups across 
Factor A for each level of Factor B.

_Parameters:_

**J: int**

Number of J levels associated with Factor A

**K: int**

Number of K levels associated with Factor B

**x: Pandas DataFrame**

Each column represents a cell in the factorial design. For example,
a 2x3 design would correspond to a DataFrame with 6 columns 
(levels of Factor A x levels of Factor B).

Order your columns according to the following pattern
(traversing each row in a matrix): 

 - the first column contains data for level 1 of Factor A
 and level 1 of Factor B
 
 - the second column contains data for level 1 of Factor A
 and level 2 of Factor B
 
 - column `K` contains the data for level 1 of Factor A 
 and level `K` of Factor B
 
 - column `K` + 1 contains the data for level 2 of Factor A
 and level 1 of Factor B
 
 - and so on ... 

**tr: float**

Proportion to trim (default is .2)

**alpha: float**

Alpha level (default is .05)

**pool: bool**

If `True`, pool dependent groups together (default is `False`).
Otherwise generate pairwise contrasts 
across factor A for each level of factor B.

_Return:_

Dictionary of results

**n: list**

Number of observations for each group

**psihat: DataFrame**

Difference score and CI, amd p-value for each contrast

**test: DataFrame**

Test statistic, critical value, standard error, and degrees of freedom for each contrast

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/bwamcp.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### bwbmcp

`#!py bwbmcp(J, K, x, tr=.2, con=None, alpha=.05,
           dif=True, pool=False, hoch=False)`

All pairwise comparisons among levels of Factor B
in a mixed design using trimmed means. The `pool`
option allows you to pool dependent groups across 
Factor A for each level of Factor B.

Rom's method is used to control for FWE (when alpha is 0.5, .01, 
or when number of comparisons are > 10).
Hochberg's method can also be used. Note that CIs are adjusted based on the 
corresponding critical p-value after controling for FWE.

    
_Parameters:_

**J: int**

Number of J levels associated with Factor A

**K: int**

Number of K levels associated with Factor B

**x: Pandas DataFrame**

Each column represents a cell in the factorial design. For example,
a 2x3 design would correspond to a DataFrame with 6 columns 
(levels of Factor A x levels of Factor B).

Order your columns according to the following pattern
(traversing each row in a matrix): 

 - the first column contains data for level 1 of Factor A
 and level 1 of Factor B
 
 - the second column contains data for level 1 of Factor A
 and level 2 of Factor B
 
 - column `K` contains the data for level 1 of Factor A 
 and level `K` of Factor B
 
 - column `K` + 1 contains the data for level 2 of Factor A
 and level 1 of Factor B
 
 - and so on ... 

**tr: float**

Proportion to trim (default is .2)

**con: array**

`con` is a K by d (number of contrasts)
matrix containing the contrast coefficents of interest.
All linear constrasts can be created automatically by using the function [con1way](K)
(the result of which can be used for `con`).

**alpha: float**

Alpha level (default is .05)

**dif: bool**

When `True`, use difference scores, otherwise use marginal distributions

**pool: bool**

If `True`, pool dependent groups together (default is `False`).
Otherwise generate pairwise contrasts 
across factor A for each level of factor B.

**hoch: bool**

When `True`, Hochberg's sequentially 
rejective method can be used to control FWE

_Return:_

Dictionary or List of Dictionaries depending on `pool` parameter. If `pool`
is set to False, all pairwise comparisons for Factor B 
are computed and returned as elements in a list corresponding to 
each level of Factor A.

**con: array**

Contrast matrix

**n: int**

Number of observations for Factor B

**num_sig: int**

Number of statistically significant results

**psihat: DataFrame**

Difference score between group X and group Y, and CI
for each contrast

**test: DataFrame**

Test statistic, p-value, critical value, and standard 
error for each contrast

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/bwbmcp.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### bwmcp

`#!py bwmcp(J, K, x, alpha=.05, tr=.2, nboot=599, seed=False)`

A bootstrap-t for multiple comparisons among
for all main effects and interactions in a between-by-within design.
The analysis is done by generating bootstrap samples and
using an appropriate linear contrast.

_Parameters:_

**J: int**

Number of J levels associated with Factor A

**K: int**

Number of K levels associated with Factor B

**x: Pandas DataFrame**

Each column represents a cell in the factorial design. For example,
a 2x3 design would correspond to a DataFrame with 6 columns 
(levels of Factor A x levels of Factor B).

Order your columns according to the following pattern
(traversing each row in a matrix): 

 - the first column contains data for level 1 of Factor A
 and level 1 of Factor B
 
 - the second column contains data for level 1 of Factor A
 and level 2 of Factor B
 
 - column `K` contains the data for level 1 of Factor A 
 and level `K` of Factor B
 
 - column `K` + 1 contains the data for level 2 of Factor A
 and level 1 of Factor B
 
 - and so on ...
 
**alpha: float**

Alpha level (default is .05)

**tr: float**

Proportion to trim (default is .2)

**nboot: int**

Number of bootstrap samples (default is 500)

**seed: bool**

Random seed for reprodicible results (default is `False`)

_Return:_

Dictionary of results

**contrast_coef: dict**

Dictionary of arrays where each value corresponds to the contrast matrix
for factor A, factor B, and the interaction

**factor_A: DataFrame**

Difference score, standard error, test statistic, 
critical value, and p-value for each contrast relating to Factor A


**factor_B: DataFrame**

Difference score, standard error, test statistic, 
critical value, and p-value for each contrast relating to Factor B

**factor_AB: DataFrame**

Difference score, standard error, test statistic, 
critical value, and p-value for each contrast relating to the interaction

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/bwmcp.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### bwimcp

`#!py bwimcp(J, K, x, tr=.2, alpha=.05)`

Multiple comparisons for interactions
in a split-plot design.
The analysis is done by taking difference scores
among all pairs of dependent groups and
determining which of
these differences differ across levels of Factor A
using trimmed means. FWE is controlled via Hochberg's 
method. For MOM or M-estimators 
(possibly not implemented yet), use spmcpi which 
uses a bootstrap method

_Parameters:_

**J: int**

Number of J levels associated with Factor A

**K: int**

Number of K levels associated with Factor B

**x: Pandas DataFrame**

Each column represents a cell in the factorial design. For example,
a 2x3 design would correspond to a DataFrame with 6 columns 
(levels of Factor A x levels of Factor B).

Order your columns according to the following pattern
(traversing each row in a matrix): 

 - the first column contains data for level 1 of Factor A
 and level 1 of Factor B
 
 - the second column contains data for level 1 of Factor A
 and level 2 of Factor B
 
 - column `K` contains the data for level 1 of Factor A 
 and level `K` of Factor B
 
 - column `K` + 1 contains the data for level 2 of Factor A
 and level 1 of Factor B
 
 - and so on ...
 
**tr: float**

Proportion to trim (default is .2)

**alpha: float**

Alpha level (default is .05)

_Return:_

Dictionary of results

**con: array**

Contrast matrix

**output: DataFrame**

Difference score, p-value, and critical value for each contrast relating to the interaction

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/bwimcp.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### bwmcppb

`#!py bwmcppb(J, K, x, est, *args, alpha=.05,
            nboot=500, bhop=True, seed=True)`

(note: this is for trimmed means only depite the `est` arg. 
This will be fixed eventually. Use `trim_mean` from SciPy)

A percentile bootstrap for multiple comparisons
for all main effects and interactions
The analysis is done by generating bootstrap samples and
using an appropriate linear contrast.
    
Uses Rom's method to control FWE. Setting the 
argument `bhop` to `True` uses the Benjamini–Hochberg 
method instead.

_Parameters:_

**J: int**

Number of J levels associated with Factor A

**K: int**

Number of K levels associated with Factor B

**x: Pandas DataFrame**

Each column represents a cell in the factorial design. For example,
a 2x3 design would correspond to a DataFrame with 6 columns 
(levels of Factor A x levels of Factor B).

Order your columns according to the following pattern
(traversing each row in a matrix): 

 - the first column contains data for level 1 of Factor A
 and level 1 of Factor B
 
 - the second column contains data for level 1 of Factor A
 and level 2 of Factor B
 
 - column `K` contains the data for level 1 of Factor A 
 and level `K` of Factor B
 
 - column `K` + 1 contains the data for level 2 of Factor A
 and level 1 of Factor B
 
 - and so on ...
 
**est: function**

Measure of location (currently only `trim_mean` is supported)

***args: list/value**

Parameter(s) for measure of location (e.g., .2)

**alpha: float**

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples (default is 500)

**bhop: bool**

When `True`, use the Benjamini–Hochberg 
method to control FWE

**seed: bool**

Random seed for reprodicible results (default is `False`)

_Return:_

Dictionary of DataFrames for each Factor and the interaction.
See the keys `'factor_A'`, `'factor_B'`, and `'factor_AB'`

Each DataFrame contains the difference score, p-value, 
critical value, and CI for each contrast.

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/bwmcppb.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### spmcpa

`#!py spmcpa(J, K, x, est, *args,
           avg=False, alpha=.05, nboot=None, seed=False)`

All pairwise comparisons among levels of Factor A 
in a mixed design. A sequentially rejective 
method is used to control FWE. The `avg` option
controls whether or not to average data across levels
of Factor B prior to performing the statistical test. 
If `False`, contrasts are created to test across Factor A
for each level of Factor B.
    
_Parameters:_

**J: int**

Number of J levels associated with Factor A

**K: int**

Number of K levels associated with Factor B

**x: Pandas DataFrame**

Data for group one

**est: function**

Measure of location (currently only `trim_mean` is supported)

***args: list/value**

Parameter(s) for measure of location (e.g., .2)

**avg: bool**

If `False`, contrasts are created to test across Factor A
for each level of Factor B (default is `False`)

**alpha: float**

Alpha level (default is .05)

**nboot: int**

Number of bootstrap samples 
(default is `None` in which case the 
number is chosen based on the number of contrasts). 

**seed: bool**

Random seed for reprodicible results (default is `False`)

_Return:_

Dictionary of results

**con: array**

Contrast matrix

**num_sig: int**

Number of statistically significant results

**output: DataFrame**

Difference score, p-value, critical value, and CI for each contrast

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/spmcpa.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### spmcpb

`#!py spmcpb(J, K, x, est, *args, dif=True, 
alpha=.05, nboot=599, seed=False)`

All pairwise comparisons among levels of Factor B
in a split-plot design. A sequentially rejective 
method is used to control FWE.

If `est` is `onestep` or `mom` (not be implemeted yet),
method SR is used to control the probability of at 
least one Type I error. Otherwise, Hochberg is used.

_Parameters:_

**J: int**

Number of J levels associated with Factor A

**K: int**

Number of K levels associated with Factor B

**x: Pandas DataFrame**

Data for group one

**est: function**

Measure of location (currently only `trim_mean` is supported)

***args: list/value**

Parameter(s) for measure of location (e.g., .2)

**dif: bool**

When `True`, use difference scores, otherwise use marginal distributions

**alpha: float**

Alpha level (default is .05)

**nboot: int**

Number of bootstrap samples (default is 599)

**seed: bool**

Random seed for reprodicible results (default is `False`)

_Return:_

Dictionary of results

**con: array**

Contrast matrix

**num_sig: int**

Number of statistically significant results

**output: DataFrame**

Difference score, p-value, critical value, and CI for each contrast

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/spmcpb.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### spmcpi

`#!py spmcpi(J, K, x, est, *args, alpha=.05, 
nboot=None, SR=False, seed=False)`

Multiple comparisons for interactions
in a split-plot design.
The analysis is done by taking difference scores
among all pairs of dependent groups and
determining which of
these differences differ across levels of Factor A.

The so-called the SR method, which is a slight 
modification of Hochberg's (1988) "sequentially rejective" 
method can be applied to control FWE, especially when 
comparing one-step M-estimators or M-estimators.

_Parameters:_

**J: int**

Number of J levels associated with Factor A

**K: int**

Number of K levels associated with Factor B

**x: Pandas DataFrame**

Data for group one

**est: function**

Measure of location (currently only `trim_mean` is supported)

***args: list/value**

Parameter(s) for measure of location (e.g., .2)

**alpha: float**

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples (default is `None`
in which case the number is 
chosen based on the number of contrasts)

**SR: bool**

When `True`, use the slight 
modification of Hochberg's (1988) "sequentially rejective" 
method to control FWE

**seed: bool**

Random seed for reprodicible results (default is `False`)

_Return:_

Dictionary of results

**con: array**

Contrast matrix

**num_sig: int**

Number of statistically significant results

**output: DataFrame**

Difference score, p-value, critical value, and CI for each contrast

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/spmcpi.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

## Measuring associations

For statistical tests and measurements that 
include robust correlations and tests of independence.
Note that regression functions will be added here eventually.

---

### corb

`#!py corb(corfun, x, y, alpha, nboot, *args, seed=False)`

Compute a 1-alpha confidence interval for a 
correlation using percentile bootstrap method
The function `corfun` is any function that returns a
correlation coefficient. The functions pbcor and
wincor follow this convention. When using 
Pearson's correlation, and when n<250, use
lsfitci instead (not yet implemented).

_Parameters:_

**corfun: function**

corfun is any function that returns a correlation coefficient

**x: Pandas Series**

Data for group one

**y: Pandas Series**

Data for group two

**alpha: float**

Alpha level (default is .05)

**nboot: int**

Number of bootstrap samples

***args: list/value**

List of arguments to corfun (e.g., .2)

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictionary of results

**ci: list**

Confidence interval

**cor: float**

Correlation estimate

**p_value: float**

p-value

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/corb.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### pball

`#!py pball(x, beta=.2)`

Compute the percentage bend correlation matrix 
for all pairs of columns in `x`. This function also 
returns the two-sided significance level for all pairs 
of variables, plus a test of zero correlation
among all pairs.

_Parameters:_

**x: Pandas DataFrame**

Each column represents a variable to use in the correlations

**beta: float**

`0 < beta < .5`. Beta is analogous to trimming in 
other functions and related to the measure of 
dispersion used in the percentage bend
calculation.

_Return:_

Dictionary of results

**H: float**

The test statistic $H$.Reject null if $H > \chi^2_{1−\alpha}$ , 
the 1−α quantile.

**H_p_value: float**

p-value corresponding to the test that all correlations are equal to zero

**p_value: array**

p-value matrix corresponding to each pairwise correlation

**pbcorm: array**

Correlation matrix

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/pball.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### pbcor

`#!py pbcor(x, y, beta=.2)`

Compute the percentage bend 
correlation between `x` and `y`

_Parameters:_

**x: Pandas Series**

Data for group one

**y: Pandas Series**

Data for group two

**beta: float**

`0 < beta < .5`. Beta is analogous to trimming in 
other functions and related to the measure of 
dispersion used in the percentage bend
calculation.

_Return:_

Dictionary of results

**cor: float**

Correlation

**nval: int**

Number of observations

**p_value**

p-value

**test: float**

Test statistic

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/pbcor.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### winall

`#!py winall(x, tr=.2)`

Compute the Winsorized correlation and covariance matrix 
for all pairs of columns in `x`. This function also 
returns the two-sided significance level for all pairs 
of variables, plus a test of zero correlation
among all pairs.

_Parameters:_

**x: Pandas DataFrame**

Each column represents a variable to use in the correlations

**tr: float**

Proportion to winsorize (default is .2)

_Return:_

Dictionary of results

**center: array**

Trimmed mean for each group

**p_value: array**

p-value array corresponding to the pairwise correlations

**wcor: array**

Winsorized correlation matrix

**wcov: array**

Winsorized covariance matrix

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/winall.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

### wincor

`#!py wincor(x, y, tr=.2)`

Compute the winsorized correlation between `x` and `y`.
This function also returns the winsorized covariance.

_Parameters:_

**x: Pandas Series**

Data for group one

**y: Pandas Series**

Data for group two

**tr: float**

Proportion to winsorize (default is .2)

_Return:_

Dictionary of results

**cor: float**

Winsorized correlation

**nval: int**

Number of observations

**sig: float**

p-value

**wcov: float**

Winsorized covariance

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/wincor.ipynb" 
target="_blank" class="button">Try this example yourself in Colab!</a>

## Other important functions

The following functions are sometimes required by Hypothesize
as input arguments. They are also potentially useful on their own.

---

### trim_mean

Calculate the sample mean after removing a proportion of values from each tail.
This is Scipy's implementation of the trimmed mean.

`#!py trim_mean(x, tr)`

_Parameters:_

**x: array or DataFrame**

Array or DataFrame of observations

**tr: float**

Proportion to trim

_Return:_

**float or list**

The trimmed mean(s)

---

### con1way

`#!py con1way(J)`

Return all linear contrasts for J groups

_Parameters:_

**J: int**

Number of groups

_Return:_

**array**

Array of contrasts where the rows correspond to groups and the columns are the contrasts to be used

---

### con2way

`#!py con2way(J, K)`

Return all linear contrasts for Factor A, Factor B, and the interaction

_Parameters:_

**J: int**

Number of levels for Factor A

**K: int**

Number of levels for Factor B

_Return:_

**list of arrays**

Each item in the list contains the contrasts for Factor A, Factor B, and the interaction, in that order.
For each array, the rows correspond to groups and the columns are the contrasts to be used

### create_example_data

`#!py create_example_data(design_values, missing_data_proportion=0, save_array_path=None, seed=False)`

Create a Pandas DataFrame of random data with a certain number of columns which
correspond to a design of a given shape (e.g., 1-way, two groups, 2-way design).
There is also an option to randomly add a proportion of null values.
The purpose of this function is to make it easy to demonstrate and test the package.

_Parameters:_

**design_values: int or list**

An integer or list indicating the design shape. For example, `[2,3]` indicates
a 2-by-3 design and will produce a six column DataFrame 
with appropriately named columns.

**missing_data_proportion: float**

Proportion of randomly missing data

**save_array_path: str (default is `None`)**

Save each group as an array for loading into R by specifying a path (e.g. , `'/home/allan/'`).
If left unset (i.e., `None`), no arrays will be saved.

**seed: bool**

Set random seed for reproducible results
