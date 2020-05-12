# Feedback and Contribution

Feedback, bug reports, and contributions are welcome via the 
[Hypothesize GitHub Repository](http://github.com/Alcampopiano/hypothesize/).

## How to Contribute Code to Hypothesize

### Working with R and Python, side-by-side

Once you have choosen a function that you like from the [WRS](https://dornsife.usc.edu/labs/rwilcox/software/)
collection, create some example data that can be used in both languages
with Hypothesize's `create_example_data`. It will save numpy arrays 
that can be read into R using the [RcppCNPy](https://cran.r-project.org/web/packages/RcppCNPy/index.html) 
and [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) libraries. 
This makes it easy to ensure the results in R match those in Python. 

In terms of IDE setup, I use PyCharm's 
[r-language-for-intellij](https://plugins.jetbrains.com/plugin/6632-r-language-for-intellij)
Plugin. This allows me to have an interpreter and editor for 
both languages in the same IDE. Like so:

<img src="https://github.com/Alcampopiano/hypothesize/blob/master/docs/docs/img/ide_pycharm.png?raw=true" alt="drawing"/>

### Setting up your Git Environment

Install the latest version of Hypothesize locally using 
```
$ pip install git+https://github.com/Alcampopiano/hypothesize/
```
Next step is to fork the repository on GitHub and clone the fork to you local
machine. For more details on forking see the [GitHub
Documentation](https://help.github.com/en/articles/fork-a-repo).

```
$ git clone https://github.com/YOUR-USERNAME/hypothesize.git
```

You can have a single clone of the repository that points to both your fork and
the main package repository. These pointers to GitHub are called "remotes".
On your local clone you should run

```
$ git remote add upstream https://github.com/Alcampopiano/hypothesize.git
$ git checkout master
$ git pull upstream master
```

And then you'll have all the updates in the master branch of your local fork.
Note that git will complain if you've committed changes to your local master
branch that are not on upstream (this is one reason it's good practice to **never**
work directly on your master branch).

### Creating a Branch

Once your local environment is up-to-date, you can create a new git branch which will
contain your contribution:
```
$ git checkout -b <branch-name>
```

With this branch checked-out, make the desired changes to the package.
When you are happy with your changes, you can commit them to your branch by running

```
$ git add <modified-file>
$ git commit -m "Some descriptive message about your change"
$ git push origin <branch-name>
```

Finally you will need to submit a pull request (PR) on GitHub asking to merge
your example branch into Hypothesize's master branch. For details on creating a PR see GitHub
documentation [Creating a pull
request](https://help.github.com/en/articles/creating-a-pull-request). You can
add more details about your example in the PR such as motivation for the
example or why you thought it would be a good addition.  You will get feedback
in the PR discussion if anything needs to be changed. To make changes continue
to push commits made in your local example branch to origin and they will be
automatically shown in the PR. 

Hopefully your PR will be answered in a timely manner and your contribution will
help others in the future.

### Testing your Changes

Hypothesize uses `pytest` for unit testing. The strategy currently used for testing
is to pickle results that are assumed to be correct and compare those
against fresh results from the modified code (see the "test_data" folder).

If you write a new function, consider writing a test for it. 
You may follow the strategy described above if you like, or come up with another
way to test your code.

When you submit a pull request, the continuous integration test suite will
run tests to validate the correctness of certain functions. It is a good 
idea run these existing tests locally to see if they pass. 
Feel free to pickle new results and save them in the "test_data" folder
if this is needed for testing purposes.
