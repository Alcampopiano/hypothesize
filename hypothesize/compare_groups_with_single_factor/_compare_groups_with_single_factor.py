__all__ = ["yuenbt", "linconb"]

import numpy as np
import pandas as pd
from scipy.stats import trim_mean
from hypothesize.utilities import trimse, lincon, \
    trimparts, trimpartt, pandas_to_arrays

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
    :return: dict of CI, test_stat, p_value, est_x, est_y, est_dif
    """

    x, y=pandas_to_arrays([x, y])

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

    results = {'ci': ci, 'test_stat': test_stat, 'p_value': p_value,
               'est_x': est_x, 'est_y': est_y, 'est_dif': est_dif}

    return results

def linconb(x, con, tr, alpha, nboot, seed=False):

    """
    Compute a 1-alpha confidence interval for a set of d linear contrasts
    involving trimmed means using the bootstrap-t bootstrap method.
    Independent groups are assumed.

    CIs are adjusted to control FWE (p values are not adjusted)

    x[1] contains the data for the first group, x[2] the data
    for the second group, etc. len(x)= the number of groups = J

    Missing values are automatically removed.

    con is a J by d matrix containing the contrast coefficents of interest.
    If unspecified, all pairwise comparisons are performed.
    For example, con[:,0]=[1,1,-1,-1,0,0] and con[:,1]=[1,-1,0,0,1,-1]
    will test two contrasts: (1) the sum of the first two trimmed means is
    equal to the sum of the second two, and (2) the difference between
    the first two is equal to the difference between the trimmed means of
    groups 5 and 6.

    The default number of bootstrap samples is nboot=599

    :param x: 1 x number of groups array. x[0] contains the data for the first group, and so on
    :param tr: amount of trimming
    :param con: contrast matrix (see con1way)
    :param alpha: alpha level
    :param nboot: number of bootstrap samples
    :param seed: random seed for reproducibility
    :return: n for each group, psihat, test statistic, critical value, contrast matrix
    """

    x=pandas_to_arrays(x)

    J = len(x)
    x = np.asarray([j[~np.isnan(j)] for j in x])
    #Jm = J - 1
    #d = (J ** 2 - J) / 2

    if con.shape[0] != len(x):
      raise Exception("The number of groups does not match the number of contrast coefficients.")

    bvec = np.zeros([nboot, J, 2])

    if seed:
        np.random.seed(seed)

    nsam = [len(xi) for xi in x]
    for j in range(J):

        xcen = x[j] - trim_mean(x[j], tr)
        data = np.random.choice(xcen, size=(nboot, len(x[j])))

        for i, row in enumerate(data):
            bvec[i,j,:]=trimparts(row, tr)

    m1 = bvec[:,:,0].T
    m2 = bvec[:,:, 1].T
    boot = np.zeros([con.shape[1], nboot])
    for d in range(con.shape[1]):
        top = np.asarray([trimpartt(row, con[:,d]) for row in m1.T])
        consq = con[:, d] ** 2
        bot = np.asarray([trimpartt(row,consq) for row in m2.T])
        boot[d,:] = np.abs(top) / np.sqrt(bot)

    testb=np.asarray([max(row) for row in boot.T])
    ic = int(np.floor((1 - alpha) * nboot) -1) # one less than R
    testb = np.sort(testb)
    psihat = np.zeros([con.shape[1], 4])
    test = np.zeros([con.shape[1], 4])

    for d in range(con.shape[1]):
        test[d, 0] = d
        psihat[d, 0] = d
        testit = lincon(x, np.array([con[:,d]]).T, tr, alpha) # column slice of contrast matrix
        #test[d, 1]=testit['test'][0, 1]
        test[d, 1]=testit['test']['test'][0]
        #pval = np.mean((abs(testit['test'][0, 1]) < boot[d,:]))
        pval = np.mean((abs(testit['test']['test'][0]) < boot[d,:]))
        test[d, 3] = pval
        print(testit['test'])
        print(testit['psihat'])
        # psihat[d, 2] = testit['psihat'][0, 1] - testb[ic] * testit['test'][0, 3]
        # psihat[d, 3] = testit['psihat'][0, 1] + testb[ic] * testit['test'][0, 3]
        # psihat[d, 1] = testit['psihat'][0, 1]
        psihat[d, 2] = testit['psihat']['psihat'][0] - testb[ic] * testit['test']['se'][0]
        psihat[d, 3] = testit['psihat']['psihat'][0] + testb[ic] * testit['test']['se'][0]
        psihat[d, 1] = testit['psihat']['psihat'][0]
        #test[d, 2] = testit['test'][0, 3]
        test[d, 2] = testit['test']['se'][0]



    psihat_col_names=['contrast_index', 'psihat', 'ci_low', 'ci_up']
    test_col_names = ['contrast_index', 'test', 'se', 'p_value']

    psihat = pd.DataFrame(psihat, columns=psihat_col_names)
    test=pd.DataFrame(test, columns=test_col_names)

    return {'n': nsam, 'psihat': psihat,  'test': test, 'crit': testb[ic], 'con': con}



