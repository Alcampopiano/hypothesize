# Overview

The benefits of using robust methods for hypothesis testing 
have been known for the last half century. 
They have been shown to subtantially increase power and accuracy when compared to 
traditional approaches.
The issues of robustness and the functions in this library are described in detail in Rand R. Wilcox's book 
[Introduction to Robust Estimation and Hypothesis Testing](https://play.google.com/store/books/details?id=8f8nBb4__EYC&gl=ca&hl=en-CA&source=productsearch&utm_source=HA_Desktop_US&utm_medium=SEM&utm_campaign=PLA&pcampaignid=MKTAD0930BO1&gclid=CjwKCAiA44LzBRB-EiwA-jJipJzyqx9kwNMq5MMU7fG2RrwBK9F7sirX4pfhS8wO7k9Uz_Sqf2P28BoCYzcQAvD_BwE&gclsrc=aw.ds).

The code and function names in Hypothesize are based on Wilcox's R functions in the [WRS](somwhere) package. 
Hypothesize simply brings many of these helpful and well-studied robust methods to the Python ecosystem. 
In addition, Hypothesize provides a user-friendly API and package structure 
as well as one-click, [ready-to-run examples](function_guide.md) for every top-level 
function.

## Hypothesize is easy to use

Hypothesize's API is friendly and 
consistent, making it easy for you to discover 
and use robust functions that are appropriate for 
your statistical design.

### Package Structure

Hypothesize organizes functions
based on the statistical design. The following visualizations show
how the package is structured and how
this is reflected in practice when importing from the library:

```mermaid
graph TB
linkStyle default interpolate basis
A[Hypothesize]
A --> B(compare groups with single factor)
A --> C(compare groups with two factors)
A --> D(measure associations)

B --> F(f<sub>1</sub>)
B --> G(f<sub>2</sub>)
B --- H(f<sub>n</sub>)

C --> F1(f<sub>1</sub>)
C --> G2(f<sub>2</sub>)
C --> H3(f<sub>n</sub>)

D --> F5(f<sub>1</sub>)
D --> G6(f<sub>2</sub>)
D --> H7(f<sub>n</sub>)
```
	
---
![Screenshot](img/package_import_viz.gif)

---
## Hypothesize is flexible and powerful

Depending on the statistical test, 
Hypothesize allows you to specify parameters to control the
following (not an exhaustive list):

- the estimator
    - trimmed mean
    - winsorized correlation
    - percentage bend correlation
    - mean
    - median
    
    !!! note "Hypothesize focuses on methods that relate to the trimmed mean"
        However, in places where Hypothesize accepts a function as an input,
        there is a good chance that other estimators will also work
        (e.g., one-step M-estimator, MOM)


- family-wise error (FWE) using sequentially rejective methods
- how group differences are calculated (e.g., marginal, pairwise, all combinations of differences)
- linear contrasts
- the alpha level
- the random seed (for reprodicible results when using bootstrap-based tests)

<br>

Visit the [tutorial section](basic_tutorial.md) and the 
[function documentation](function_guide.md) for complete examples
using Hypothesize.