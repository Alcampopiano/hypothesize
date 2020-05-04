__all__ = ["wincor", "pbcor", "corb", "pball", "winall"]

import numpy as np
from scipy.stats.mstats import winsorize
from scipy.stats import t, chi2, trim_mean
from hypothesize.utilities import pandas_to_arrays

def wincor(x, y, tr=.2):

    """
    Compute the winsorized correlation between x and y.
    This function also returns the winsorized covariance
    Pairwise deletion of missing values is performed.

    :param x: Pandas Series
    :param y: Pandas Series
    :param tr: amount of winsorization
    :return: winsorized correlation and winsorized covariance, p_value, and number of rows
    """

    if type(x) is not np.ndarray:
        x, y=pandas_to_arrays([x, y])

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

    res={'cor': wcor, 'wcov': wcov, 'sig': sig, 'nval': nval}

    return res

def pbcor(x, y, beta=.2):

    """
    Compute the percentage bend correlation between x and y.

    :param x: data array
    :param y: data array
    :param beta: bending constant for omega sub N
    :return: percentage bend correlation
    """

    if type(x) is not np.ndarray:
        x, y = pandas_to_arrays([x, y])

    if len(x) != len(y):
        raise Exception("The arrays do not have equal lengths")

    m1 = np.c_[x, y] # cbind
    m1 = m1[~np.isnan(m1).any(axis=1)]
    nval = m1.shape[0]
    x = m1[:, 0]
    y = m1[:, 1]
    temp = np.sort(abs(x - np.median(x)))
    omhatx = temp[int(np.floor((1 - beta) * len(x)))-1]
    temp = np.sort(abs(y - np.median(y)))
    omhaty = temp[int(np.floor((1 - beta) * len(y)))-1]

    a = (x - pbos(x, beta)) / omhatx
    b = (y - pbos(y, beta)) / omhaty

    a = np.where(a <= -1, -1, a)
    a = np.where(a >= 1, 1, a)
    b = np.where(b <= -1, -1, b)
    b = np.where(b >= 1, 1, b)

    pbcor_result = sum(a * b) / np.sqrt(sum(a ** 2) * sum(b ** 2))
    test = pbcor_result * np.sqrt((len(x) - 2) / (1 - pbcor_result ** 2))
    sig = 2 * (1 - t.cdf(abs(test), len(x) - 2))

    res = {'cor': pbcor_result, 'test': test, 'p_value': sig, 'nval': nval}
    return res

def pbos(x, beta=.2):

    """
    Compute the one-step percentage bend measure of location

    :param x:
    :param beta:
    :return:
    """

    temp = np.sort(abs(x - np.median(x)))
    omhatx = temp[int(np.floor((1 - beta) * len(x)))-1]
    psi = (x - np.median(x)) / omhatx
    i1 = len(psi[psi < -1])
    i2 = len(psi[psi > 1])

    sx = np.where(psi < -1, 0, x)
    sx = np.where(psi > 1, 0, sx)

    pbos_result = (sum(sx) + omhatx * (i2 - i1)) / (len(x) - i1 - i2)

    return  pbos_result

def corb(corfun, x, y, alpha, nboot, *args, seed=False):

    """
    Compute a 1-alpha confidence interval for a correlation using percentile bootstrap method

    The function corfun is any function that returns a
    correlation coefficient. The functions pbcor and
    wincor follow this convention.

    When using Pearson's correlation, and when n<250, use
    lsfitci instead (not yet implemented).

    :param x: Series
    :param y: Series
    :param corfun: corfun is any function that returns a correlation coefficient
    :param alpha: alpha level
    :param nboot: number of bootstrap samples
    :param seed: random seed for reproducibility
    :param args: list of arguments to corfun (e.g., args=[.2])
    :return: CI, p_value, correlation estimate
    """

    x, y=pandas_to_arrays([x, y])


    m1 = np.c_[x, y] # cbind
    m1 = m1[~np.isnan(m1).any(axis=1)]
    nval = m1.shape[0]
    x = m1[:, 0]
    y = m1[:, 1]
    est = corfun(x, y, *args)['cor']#[0]

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

    #return corci, sig, est
    return {'ci': corci, 'p_value': sig, 'cor': est}

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

    corbsub_results = corfun(x[isub], y[isub], *args)['cor']#[0]

    return corbsub_results

def pball(x, beta=.2):

    """
    Compute the percentage bend correlation matrix for the
    data in the n by p matrix m.

    This function also returns the two-sided significance level
    for all pairs of variables, plus a test of zero correlations
    among all pairs. (See chapter 6 for details)

    :param x: Pandas DataFrame
    :param beta:
    :return:
    """

    m=x.values
    ncol=m.shape[1]

    pbcorm=np.zeros([ncol, ncol])
    temp=np.ones([ncol, ncol])
    siglevel=np.full([ncol, ncol], np.nan)
    #cmat = np.zeros([ncol, ncol])

    for i in range(ncol):
        for j in range(i,ncol):
            if i < j:
                pbc = pbcor(m[:, i], m[:, j], beta)
                pbcorm[i, j] = pbc['cor']
                temp[i, j] = pbcorm[i, j]
                temp[j, i] = pbcorm[i, j]
                siglevel[i, j] = pbc['p_value']
                siglevel[j, i] = siglevel[i, j]


    tstat = pbcorm * np.sqrt((m.shape[0] - 2) / (1 - pbcorm ** 2))
    cmat = np.sqrt((m.shape[0] - 2.5) * np.log(1 + tstat ** 2 / (m.shape[0] - 2)))
    bv = 48 * (m.shape[0] - 2.5) ** 2
    cmat = \
    cmat + (cmat ** 3 + 3 * cmat) / bv - (4 * cmat ** 7 + 33 * cmat ** 5 + 240 ** cmat ** 3 + 855 * cmat) / \
        (10 * bv ** 2 + 8 * bv * cmat ** 4 + 1000 * bv)

    H = np.sum(cmat ** 2)
    df = ncol * (ncol - 1) / 2
    h_siglevel = 1 - chi2.cdf(H, df)

    results={"pbcorm": temp, "p_value": siglevel,
             "H":H, "H_p_value": h_siglevel}

    return results

def winall(x, tr=.2):

    """
    Compute the Winsorized correlation and covariance matrix for the
    data in the n by p matrix m.

    This function also returns the two-sided significance level
    for all pairs of variables, plus a test of zero correlations
    among all pairs. (See chapter 6 for details.)

    :param x: Pandas DataFrame
    :param tr:
    :return:
    """

    m = x.values
    ncol = m.shape[1]

    wcor = np.ones([ncol, ncol])
    wcov = np.zeros([ncol, ncol])
    siglevel = np.full([ncol, ncol], np.nan)

    for i in range(ncol):
        #ip = i
        for j in range(i,ncol):
            val = wincor(m[:, i], m[:, j], tr)
            wcor[i, j] = val['cor']
            wcor[j, i] = wcor[i, j]

            if i == j:
                wcor[i, j] = 1

            wcov[i, j] = val['cor']
            wcov[j, i] = wcov[i, j]

            if i != j:
                siglevel[i, j] = val['sig']
                siglevel[j, i] = siglevel[i, j]

    m=m[~np.isnan(m).any(axis=1)]
    cent=trim_mean(m, tr)

    return {"wcor": wcor, "wcov": wcov, "center": cent, "p_value": siglevel}
