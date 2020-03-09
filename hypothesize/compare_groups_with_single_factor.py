import numpy as np
from scipy.stats import trim_mean
from hypothesize.utilities import winvar, trimse

def yuenbt(x, y, tr=.2, alpha=.05, nboot=599, seed=False):

    """
    Compute a 1-alpha confidence interval for the difference between
    the trimmed means corresponding to two independent groups.
    The bootstrap-t method is used. During the bootstrapping, the absolute value of the test
    statistic is used (the "two-sided method").

    :param x: group one data (array)
    :param y: group two data (array)
    :param tr: proportion to trim
    :param alpha: alpha level
    :param nboot: number of bootstrap samples
    :param seed: seed value to set for reproducible results
    :return: CI, test_stat, p_value, est_x, est_y, est_dif
    """
    if seed:
        np.random.seed(seed)

    ci=[]
    x=x[~np.isnan(x)]
    y=y[~np.isnan(y)]

    xcen = x - trim_mean(x, tr)
    ycen = y - trim_mean(y, tr)

    test_stat = (trim_mean(x, tr) - trim_mean(y, tr)) / \
           np.sqrt(trimse(x, tr = tr) ** 2 + trimse(y, tr = tr) ** 2)

    datax = np.random.choice(xcen, size=(nboot, len(x)))
    datay = np.random.choice(ycen, size=(nboot, len(y)))

    top = trim_mean(datax, .2, axis=1) - trim_mean(datay, .2, axis=1)

    #botx = list(map(lambda row: trimse(row,.2), datax))
    botx = np.array([trimse(x) for x in datax])
    boty = np.array([trimse(x) for x in datay])
    tval = top / np.sqrt(botx ** 2 + boty ** 2)
    tval = abs(tval)
    tval = sorted(tval)
    icrit = int(np.floor((1 - alpha) * nboot + .5))
    #ibot = int(np.floor(alpha * nboot / 2 + .5))
    #itop = int(np.floor((1 - alpha / 2) * nboot + .5))
    se = np.sqrt((trimse(x, tr)) ** 2 + (trimse(y, tr)) ** 2)
    ci.append(trim_mean(x, tr) - trim_mean(y, tr) - tval[icrit] * se)
    ci.append(trim_mean(x, tr) - trim_mean(y, tr) + tval[icrit] * se)
    p_value = sum(np.abs(test_stat) <= np.abs(tval)) / nboot
    est_x = trim_mean(x,tr)
    est_y = trim_mean(y, tr)
    est_dif = est_x - est_y


    return ci, test_stat, p_value, est_x, est_y, est_dif



