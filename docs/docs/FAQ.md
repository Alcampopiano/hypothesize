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

First, it is helpful to imagine your
design arranged into a JxK matrix. 

$$
A=\begin{bmatrix} 
a_{1,1} & a_{1,2} & ... & a_{1,K} \\ 
a_{2,1} & a_{2,2} & ... & a_{2,K} \\
a_{J,1} & a_{J,2} & ... & a_{J,K}
\end{bmatrix}
$$

A contrast matrix specifies which cells (or elements) in the above
design are to be compared. The rows in a contrast matrix
correspond to the cells in your design. The columns correspond
to the contrasts that you wish to make.
    
### Examples of contrast matrices for different designs

Matrix notation is used to explain which cells are
being compared, followed by the corresponding 
contrast matrix.

=== "design with 2 groups"
    
    ${A_{1,1} - A_{1,2}}$
    
    | contrast 1 |
    |------------|
    |  1         |
    |  -1        |
    
=== "design with 3 groups"

    1. $\Large{A_{1,1} - A_{1,2}}$  
    2. $\Large{A_{1,1} - A_{1,3}}$  
    3. $\Large{A_{1,2} - A_{1,3}}$  

    | contrast 1 | contrast 2 | contrast 3 | 
    |------------|------------|------------|
    |  1         |   1        |    0       | 
    | -1         |   0        |    1       | 
    |  0         |  -1        |    -1      | 

=== "2x2 design"
    **Factor A**
    
    $\Large{(A_{1,1} + A_{1,2})-(A_{2,1} + A_{2,2})}$  
    
    | contrast 1 | 
    |------------|
    |  1         |  
    |  1         |  
    | -1         |  
    | -1         |  
    
    **Factor B**
    
    $\Large{(A_{1,1} + A_{2,1})-(A_{1,2} + A_{2,2})}$  
    
    | contrast 1 | 
    |------------|
    |  1         |  
    | -1         |  
    |  1         |  
    | -1         | 
    
    **Interaction**
    
    $\Large{(A_{1,1} + A_{2,2})-(A_{1,2} + A_{2,1})}$  
    
    That is, the difference of the differences

    | contrast 1 | 
    |------------|
    |  1         |  
    | -1         |  
    | -1         |  
    |  1         | 
    
=== "2x3 design"
    **Factor A**
    
    $\Large{(A_{1,1} + A_{1,2} + A_{1,3})-(A_{2,1} + A_{2,2} + A_{2,3})}$  
    
    | contrast 1 |   
    |------------|
    |  1         |  
    |  1         |  
    |  1         |  
    | -1         |  
    | -1         |  
    | -1         |  
        
    **Factor B**
    
    1. $\Large{(A_{1,1} + A_{2,1})-(A_{1,2} + A_{2,2})}$  
    -  $\Large{(A_{1,1} + A_{2,1})-(A_{1,3} + A_{2,3})}$   
    -  $\Large{(A_{1,2} + A_{2,2})-(A_{1,3} + A_{2,3})}$    
    
    | contrast 1 | contrast 2 | contrast 3 | 
    |------------|------------|------------|
    |  1         |   1        |    0       | 
    | -1         |   0        |    1       | 
    |  0         |  -1        |    -1      | 
    |  1         |   1        |    0       | 
    | -1         |   0        |    1       | 
    |  0         |  -1        |    -1      | 
    
    **Interactions**
    
    1. $\Large{(A_{1,1} + A_{2,2})-(A_{1,2} + A_{2,1})}$  
    -  $\Large{(A_{1,1} + A_{2,3})-(A_{1,3} + A_{2,1})}$   
    -  $\Large{(A_{1,2} + A_{2,3})-(A_{1,3} + A_{2,2})}$  
    
    | contrast 1 | contrast 2 | contrast 3 | 
    |------------|------------|------------|
    |  1         |   1        |    0       | 
    | -1         |   0        |    1       | 
    |  0         |  -1        |    -1      | 
    |  -1        |  -1        |     0      | 
    |  1         |   0        |   -1       | 
    |  0         |   1        |    1       | 
    
    
!!! success "Not a fan of contrast matrices?"
    Don't worry, Hypothesize can generate all linear
    contrasts automatically (see functions [con1Way]()
    and [con2way]()). However, it is useful to 
    understand this concept so that you know
    which comparisons are being made and 
    how to specify your own if necessary.
    
<br>