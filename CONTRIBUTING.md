# Feedback and contribution

Feedback, bug reports, and contributions are welcome via the 
[Hypothesize GitHub Repository](http://github.com/Alcampopiano/hypothesize/).

## How to contribute new functions to Hypothesize

A great way to contribute would be to choose a function from the 
[WRS](https://dornsife.usc.edu/labs/rwilcox/software/) that does not yet exist in
Hypothesize and convert it to Python. There is a current wish list 
[here](https://github.com/Alcampopiano/hypothesize/issues/2)
but certainly any WRS function would be a welcome addition to the library. A list of the currently available
functions in Hypothesize can be found in the documentation's
[function reference](https://alcampopiano.github.io/hypothesize/function_guide/).

#### Create example data to be used in R and Python

It is helpful to be able to create some example data that can be used in both R and Python. 
One way to do this is to use Hypothesize's 
[create_example_data](https://alcampopiano.github.io/hypothesize/function_guide/#create_example_data) function. 
It will generate a DataFrame of random data (to be used in Python) as 
well save Numpy arrays that can be read into R with the
[RcppCNPy](https://cran.r-project.org/web/packages/RcppCNPy/index.html) 
and [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) libraries.

#### IDE for R and Python

It is convenient to use the same IDE when converting functions from R to Python.
One suggestion is to use PyCharm's 
[r-language-for-intellij](https://plugins.jetbrains.com/plugin/6632-r-language-for-intellij)
Plugin. This makes it possible to have an interpreter and editor for 
both languages in the same IDE. Like so:

<img src="https://github.com/Alcampopiano/hypothesize/blob/master/docs/docs/img/ide_pycharm.png?raw=true" alt="drawing"/>

Of course there are many ways that one might go about converting WRS functions to Python. 
These are merely suggestions. 

### Setting up your Git environment

1. Install the latest version of Hypothesize locally using 
    
    ```
    $ pip install git+https://github.com/Alcampopiano/hypothesize/
    ```

2. Fork the repository on GitHub and clone the fork to you local
machine. For more details on forking see the [GitHub
Documentation](https://help.github.com/en/articles/fork-a-repo).
    
    ```
    $ git clone https://github.com/YOUR-USERNAME/hypothesize.git
    ```

3. Create a sync to the original upstream repository by creating a so-called 
[remote](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/configuring-a-remote-for-a-fork).

    ```
    $ git remote add upstream https://github.com/Alcampopiano/hypothesize.git
    $ git checkout master
    $ git pull upstream master
    ```

Now you will have all of the updates in the master branch of your local fork.
Note that git will complain if you've committed changes to your local master
branch that are not on the upstream repository. This is one reason why it's good practice to avoid
working directly on your master branch.

### Commiting new code to Hypothesize

1. Create a new local branch and commit changes to your remote branch:

    ```
    $ git checkout -b <branch-name>
    ```
    
    With this branch checked-out, make the desired changes to the package.
    When you are happy with your changes, you can commit them to a remote branch by running
    
    ```
    $ git add <modified-file>
    $ git commit -m "Some descriptive message about your change"
    $ git push origin <branch-name>
    ```

2. Write a unit test for your code (optional)

    Hypothesize uses `pytest` for unit testing. The strategy currently used for testing
    is to pickle results that are assumed to be correct and compare those
    against fresh results from the modified code (see the
    [tests](https://github.com/Alcampopiano/hypothesize/tree/master/hypothesize/tests) folder for examples).
    If you would like to write a test for your new code, you may follow the strategy 
    described above or come up with another way to test your code. To run the test suite,
    first navigate to the "tests" directory then use the `pytest` command from your terminal.

3. Submit a pull request (PR) to merge your new branch to Hypothesize's master branch

    For details on creating a PR see GitHub documentation [Creating a pull
    request](https://help.github.com/en/articles/creating-a-pull-request). 

