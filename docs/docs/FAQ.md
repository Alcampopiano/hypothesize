# Frequently asked questions

No attempt is made to fully explain the following
concepts, but hopefully this gets
you started. The Internet has plenty of resources on these topics
if you would like to learn more.

## What is a trimmed mean?

The trimmed mean involves calculating the sample mean after
removing a proportion of values from each
tail of the distribution. In symbols the trimmed mean is expressed as
follows:

$$
\bar{X}_t = \frac{X_{(g+1)}\,+,...,+\,X_{(n-g)}}{n-2g}
$$

where $X_1, \,X_2,\,...\,,X_n$ is a random sample and
$X_{(1)}, \le X_{(2)}\,,...,\,\le X_{(n)}$ are the observations in
ascending order. The proportion to trim is $\gamma\,(0\lt \gamma \lt.5)$
and $g = [ \gamma n ]$ rounded down to the nearest integer.

## What is bootstrapping?

In the context of hypothesis testing and generally speaking,
bootstrapping involves taking many random samples (with replacement)
from the data at hand in order to estimate a sampling
distribution of interest. This is in contrast to traditional methods
which assume the shape of the particular sampling distribution under study.
Once we have an emprically derived sampling distribution,
obtaining CIs and p values is relatively straightforward.

## What is a contrast matrix?

A contrast matrix is an array of 1's and 0's
that indicates how conditions (or groups) are
to be compared. Hypothesize uses these matrices
internally to make certain calculations more
convienient and requires them to be specified 
as inputs to some functions.

The rows in a contrast matrix correspond to
the groups in your design (in order). The columns
correspond to each linear contrasts that you want to
make. That is, they indicate which conditions are being compared.

!!! success "Not a fan of contrast matrices?"
    Don't worry, Hypothesize can generate all pairwise
    contrasts automatically, even for factorial designs.
    See the functions [con1Way]() and [con2way](); however,
    it is useful to understand how to read a contrast matrix.

### Here are some examples of contrasts matrices for different designs

=== "design with 2 groups"
    | contrast 1 |
    |------------|
    |  1         |
    |  -1        |
    
=== "design with 3 groups"
    | contrast 1 | contrast 2 | contrast 3 | 
    |------------|------------|------------|
    |  1         |   1        |    0       | 
    | -1         |   0        |    1       | 
    |  0         |  -1        |    -1      | 

=== "2x2 design"
    **Factor A**
    
    | contrast 1 | 
    |------------|
    |  1         |  
    |  1         |  
    | -1         |  
    | -1         |  
    
    **Factor B**
    
    | contrast 1 | 
    |------------|
    |  1         |  
    | -1         |  
    |  1         |  
    | -1         | 
    
    **Interaction**
    
    | contrast 1 | 
    |------------|
    |  1         |  
    | -1         |  
    | -1         |  
    |  1         | 
    
=== "2x3 design"
    **Factor A**
    
    | contrast 1 |   
    |------------|
    |  1         |  
    |  1         |  
    |  1         |  
    | -1         |  
    | -1         |  
    | -1         |  
        
    **Factor B**
    
    | contrast 1 | contrast 2 | contrast 3 | 
    |------------|------------|------------|
    |  1         |   1        |    0       | 
    | -1         |   0        |    1       | 
    |  0         |  -1        |    -1      | 
    |  1         |   1        |    0       | 
    | -1         |   0        |    1       | 
    |  0         |  -1        |    -1      | 
    
    **Interactions**
    
    | contrast 1 | contrast 2 | contrast 3 | 
    |------------|------------|------------|
    |  1         |   1        |    0       | 
    | -1         |   0        |    1       | 
    |  0         |  -1        |    -1      | 
    |  -1        |  -1        |     0      | 
    |  1         |   0        |   -1       | 
    |  0         |   1        |    1       | 