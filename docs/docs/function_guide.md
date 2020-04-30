# Function Reference

Hypothesize exposes the following top-level functions for comparing 
groups and measuring associations. The function names,  code, and descriptions
are kept generally consistent with Wilcox's WRS package. 
If you want to learn more about the theory and research
behind any given function here, see Wilcox's books, especially
[Introduction to Robust Estimation and Hypothesis Testing](https://play.google.com/store/books/details?id=8f8nBb4__EYC&gl=ca&hl=en-CA&source=productsearch&utm_source=HA_Desktop_US&utm_medium=SEM&utm_campaign=PLA&pcampaignid=MKTAD0930BO1&gclid=CjwKCAiA44LzBRB-EiwA-jJipJzyqx9kwNMq5MMU7fG2RrwBK9F7sirX4pfhS8wO7k9Uz_Sqf2P28BoCYzcQAvD_BwE&gclsrc=aw.ds).

---

## Comparing groups with a single factor

<a name="comp_single_ind"></a>
## Independent groups

### l2drmci

---

`l2drmci(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/l2drmci.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### linconb

---

`linconb(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/linconb.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### pb2gen

---

`pb2gen(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/pb2gen.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### tmcppb

---

`tmcppb(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/tmcppb.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### yuenbt

---

`yuenbt(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/yuenbt.ipynb" 
target="_blank">Try this function right now in Colab!</a>

<a name="comp_single_dep"></a>
## Dependent groups

### bootdpci

---

`bootdpci(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/bootdpci.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### rmmcppb

---

`rmmcppb(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/rmmcppb.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### lindepbt

---

`lindepbt(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/lindepbt.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### ydbt

---

`ydbt(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/ydbt.ipynb" 
target="_blank">Try this function right now in Colab!</a>

## Comparing groups with two factors

<a name="comp_double_dep"></a>
## Dependent groups

### wwmcppb

---

`wwmcppb(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/wwmcppb.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### wwmcpbt

---

`wwmcpbt(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/wwmcpbt.ipynb" 
target="_blank">Try this function right now in Colab!</a>

<a name="comp_double_mix"></a>
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

---

`bwamcp(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/bwamcp.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### bwbmcp

---

`bwbmcp(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/bwbmcp.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### bwcmp

---

`bwcmp(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/bwcmp.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### bwimcp

---

`bwimcp(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/bwimcp.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### bwmcppb

---

`bwmcppb(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/bwmcppb.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### spmcpa

---

`spmcpa(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/spmcpa.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### spmcpb

---

`spmcpb(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/spmcpb.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### spmcpi

---

`spmcpi(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/spmcpi.ipynb" 
target="_blank">Try this function right now in Colab!</a>

## Measuring associations

These functions deal with correlational statistics (and eventually regression).
For some of these functions, the input data are given as a Pandas Series for `x` and for `y'.

### corb

---

`corb(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/corb.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### pball

---

`pball(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/pball.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### pbcor

---

`pbcor(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/pbcor.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### winall

---

`winall(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/winall.ipynb" 
target="_blank">Try this function right now in Colab!</a>

### wincor

---

`wincor(x, y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False)`

---
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

Alpha level. Default is .05.

**nboot: int**

Number of bootstrap samples. Default is 2000.

**seed: bool**

Random seed for reprodicible results. Default is `False`.

_Return:_

Dictonary of results

<a class="btn btn-info btn-lg btn-block" 
href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/wincor.ipynb" 
target="_blank">Try this function right now in Colab!</a>

