# Basic Tutorial

The following tutorial demonstrates how to perform a 
robust hypothesis test using 20% trimmed means and 
the bootstrap-t test. The data correspond to a 
2 (between-subjects) x 3 (within-subjects) factorial design. 

### Getting your data into Hypothesize

In Hypothesize, input data are always specified as a Pandas DataFrame or Series. 
In this example, we have a 2x3 factorial design so the data would take the form of 
a six-column DataFrame (i.e., J levels x K levels). Using Pandas you can read your data into Python and 
use one of the appropriate functions from Hypothesize. In this case we will use the function `bwmcp`
but there are [many others](function_guide.md) to choose from.

!!! note ""What about my column names?""
    Don't worry, Hypothesize doesn't make use of your column names. 
    Feel free to name them however you like!


```python
import pandas as pd

df=pd.read_csv('my_data.csv')

df.head() 
```

| cell_1_1   |   cell_1_2 |   cell_1_3 |   cell_2_1 |   cell_2_2 |   cell_2_3 |
|------------|------------|------------|------------|------------|------------|
|  0.04      |   0.90     |   0.79     |  0.51      |   0.33     |  0.23      |
|  0.76      |   0.29     |   0.84     |  0.03      |   0.5      |  0.73      |
|  0.71      |   0.59     |   0.11     |  0.89      |   0.76     |  0.04      |
|  0.17      |   0.26     |   0.88     |  0.28      |   0.1      |  0.21      |
|  0.95      |   0.22     |   0.83     |  0.59      |   0.65     |  0.20      |
    
```python
from hypothesize.compare_groups_with_two_factors import bwmcp

results=bwmcp(J=2, K=3, x=df)
```

### Examining your results

The results are returned as a Python Dictionary containing simple Python objects
 or DataFrames (when the results are best given as a matrix). For example, here are the 
 previously computed results for the interaction returned as a DataFrame.

```python
results['factor_AB']
```
    
|   con_num |     psihat |       se |      test |   crit_value |   p_value |
|---------- |----------- |--------- |---------- |------------- |---------- |
|         0 | -0.100698  | 0.126135 | -0.798336 |       2.3771 |  0.410684 |
|         1 | -0.037972  | 0.151841 | -0.250078 |       2.3771 |  0.804674 |
|         2 |  0.0627261 | 0.135392 |  0.463291 |       2.3771 |  0.659432 |

<br>

<a href="https://colab.research.google.com/github/Alcampopiano/hypothesize/blob/master/examples/hypothesize_notebook_for_colab.ipynb" 
target="_blank" class="button">Try more examples in Colab!</a>
