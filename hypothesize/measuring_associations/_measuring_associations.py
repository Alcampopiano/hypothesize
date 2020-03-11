__all__ = ["wincor", "corb"]

import numpy as np
from scipy.stats.mstats import winsorize
from scipy.stats import t

def wincor(x, y, tr=.2):

    """
    Compute the winsorized correlation between x and y.
    This function also returns the winsorized covariance
    Pairwise deletion of missing values is performed.

    x is a vector, or it can be a matrix with two columns when y=NULL

    :param x: x is an array (n by 1)
    :param y: y is an array (n by 1)
    :param tr: amount of winsorization
    :return: winsorized correlation and winsorized covariance, p_value, and number of rows
    """

    m1 = np.c_[x, y] # cbind
    m1 = m1[~np.isnan(m1).any(axis=1)]
    nval = m1.shape[0]
    x = m1[:, 0]
    y = m1[:, 1]
    g = np.floor(tr * len(x))
    xvec = winsorize(x, limits=(tr,tr))
    yvec = winsorize(y, limits=(tr,tr))
    wcor = np.corrcoef(xvec, yvec)[0,1]
    wcov = np.cov(xvec, yvec)[0,1]
    test = wcor * np.sqrt((len(x) - 2) / (1. - wcor ** 2))
    sig = 2 * (1 - t.cdf(abs(test), len(x) - 2 * g - 2))

    return wcor, wcov, sig, nval

def corb(corfun, x, y, alpha, nboot, *args, seed=False):

    """
    Compute a 1-alpha confidence interval for a correlation using percentile bootstrap method

    The function corfun is any function that returns a
    correlation coefficient. The functions pbcor and
    wincor follow this convention.

    When using Pearson's correlation, and when n<250, use
    lsfitci instead (not yet implemented).

    :param x: x is an array (n by 1)
    :param y: y is an array (n by 1)
    :param corfun: corfun is any function that returns a correlation coefficient
    :param alpha: alpha level
    :param nboot: number of bootstrap samples
    :param seed: random seed for reproducibility
    :param args: list of arguments to corfun (e.g., args=[.2])
    :return: CI, p_value, correlation estimate
    """

    m1 = np.c_[x, y] # cbind
    m1 = m1[~np.isnan(m1).any(axis=1)]
    nval = m1.shape[0]
    x = m1[:, 0]
    y = m1[:, 1]
    est = corfun(x, y, *args)[0]

    if seed:
        np.random.seed(seed)

    data_inds = np.random.choice(len(x), size=(nboot, len(x)))
    bvec = np.array([corbsub(row_inds, x, y, corfun, *args) for row_inds in data_inds])

    ihi = int(np.floor((1 - alpha / 2) * nboot + .5))
    ilow = int(np.floor((alpha / 2) * nboot + .5))
    bsort = sorted(bvec)
    corci = [bsort[ilow], bsort[ihi]]
    phat = sum(bvec < 0) / nboot
    sig =  2 * min(phat, 1 - phat)

    return corci, sig, est

def corbsub(isub, x, y, corfun, *args):

    """
    Compute correlation for x[isub] and y[isub]
    isub is a vector of length n,
    a bootstrap sample from the sequence of integers
    0, 1, 2, 3, ..., n

    This function is used by other functions when computing
    bootstrap estimates.

    corfun is some correlation function
    """

    corbsub_results = corfun(x[isub], y[isub], *args)[0]

    return corbsub_results