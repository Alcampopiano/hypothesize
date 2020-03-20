import numpy as np
from scipy.stats.mstats import winsorize
from scipy.stats import trim_mean
from scipy.stats import t
from hypothesize.measuring_associations import wincor
import pandas as pd

def con1way(J):

    """
    Create contrast coefficients for all pairwise comparisons given a 1-way design
    :param J: The number of levels (i.e., groups)
    :return con: array of contrast coefficients
    """

    Ja = (J ** 2 - J) // 2
    con=np.zeros([J, Ja])

    ind=0
    for j in range(J):
        for k in range(J):
            if j<k:
                con[j, ind] = 1
                con[k, ind] = 0-1
                ind += 1

    return con

def con2way(J,K):

    """
    For a J by K design, create the contrast coefficients for all pairwaise
    comparisons for each main effect and all interactions

    :param J: The number of levels for the first factor (factor A)
    :param K: The number of levels for the second factor (factor B)
    :return conA, conB, conAB: arrays containing contrast coefficients for each factor and the interaction
    """

    Ja =(J ** 2 - J) // 2
    Ka =(K ** 2 - K) // 2
    JK = J * K
    conA=np.zeros([JK, Ja])

    ic = 0
    for j in range(J):
        for jj in range(J):
            if j < jj:
                mat=np.zeros([J,K])
                mat[j,:] = 1
                mat[jj,:] = 0-1
                conA[:, ic] = mat.reshape(conA.shape[0]).T
                ic+=1

    conB=np.zeros([JK, Ka])
    ic = 0
    for k in range(K):
        for kk in range(K):
            if k < kk:
                mat = np.zeros([J, K])
                mat[:, k] = 1
                mat[:, kk] = 0-1
                conB[:, ic] = mat.reshape(conB.shape[0]).T
                ic+=1


    conAB = np.zeros([JK, Ka * Ja])
    ic = 0
    for j in range(J):
        for jj in range(J):
            if j < jj:
                for k in range(K):
                    for kk in range(K):
                        if k < kk:
                            mat=np.zeros([J,K])
                            mat[j, k] = 1
                            mat[j, kk] = 0-1
                            mat[jj, k] = 0-1
                            mat[jj, kk] = 1
                            ic = ic + 1

                        conAB[:, ic-1] = mat.reshape(conAB.shape[0]).T

    return conA, conB, conAB

def winvar(x, tr=.2):
    """
    Compute the gamma Winsorized variance for the data in the vector x.
    tr is the amount of Winsorization which defaults to .2.
    Nan values are removed.

    :param x:
    :param tr:
    :return:
    """

    y=winsorize(x, limits=(tr,tr))
    wv = np.var(y, ddof=1)

    # x=x[~np.isnan(x)]
    # y=np.sort(x)
    # n=len(x)
    # ibot = int(np.floor(tr * n))
    # itop = len(x) - ibot -1
    # xbot = y[ibot]
    # xtop = y[itop]
    # y = np.where(y <= xbot, xbot, y)
    # y = np.where(y >= xtop, xtop, y)
    # wv = np.var(y, ddof=1) # DF to be consistent with Wilcox/R

    return wv

def trimse(x, tr=.2):

    """
    Estimate the standard error of the gamma trimmed mean
    The default amount of trimming is tr=.2.

    :param x:
    :param tr:
    :return:
    """
    #x=x[~np.isnan(x)]
    trimse_result = np.sqrt(winvar(x, tr)) / ((1 - 2 * tr) * np.sqrt(len(x)))

    return trimse_result

def lincon(x, con=None, tr=.2, alpha=.05, seed=False):

    """
    A heteroscedastic test of d linear contrasts using trimmed means.

    It is assumed all groups are independent.

    con is a J by d matrix containing the contrast coefficients that are used. If con==None,
    groups are pooled together.

    Missing values are automatically removed.

    :param x: array of data
    :param alpha: alpha value
    :param con: contrast matrix or None (for pooled groups)
    :param tr: trim proportion
    :param seed: random seed
    :return:
    """

    if tr >= .5:
      raise Exception("Use the function medpb to compare medians (it may not be implemented in this package")

    flag=True if alpha == .05 or alpha == .01 else False

    J = len(x)
    sam = []
    h = np.zeros(J)
    w = np.zeros(J)
    xbar = np.zeros(J)

    for j in range(J):
      x[j] = x[j][~np.isnan(x[j])]
      sam.append(len(x[j]))
      h[j] = len(x[j]) - 2 * np.floor(tr * len(x[j]))
      # h is the number of observations in the jth group after trimming.
      w[j] = ((len(x[j]) - 1) * winvar(x[j], tr)) / (h[j] * (h[j] - 1))
      xbar[j] = trim_mean(x[j], tr)

    if type(con) is not np.ndarray:
        CC = (J ** 2 - J) // 2
        psihat=np.zeros([CC,8])
        test=np.full([CC,6], np.nan)

        jcom=0
        for j in range(J):
            for k in range(J):
                if j<k:
                    test[jcom, 2] = abs(xbar[j] - xbar[k]) / np.sqrt(w[j] + w[k])
                    sejk = np.sqrt(w[j] + w[k])
                    test[jcom, 4] = sejk
                    psihat[jcom, 0] = j
                    psihat[jcom, 1] = k
                    test[jcom, 0] = j
                    test[jcom, 1] = k
                    psihat[jcom, 2] = (xbar[j] - xbar[k])
                    df = (w[j] + w[k]) ** 2 / (w[j] ** 2 / (h[j] - 1) + w[k] ** 2 / (h[k] - 1))
                    test[jcom, 5] = df
                    psihat[jcom, 5] = 2 * (1 - t.cdf(test[jcom, 2], df))
                    psihat[jcom, 6] = xbar[j]
                    psihat[jcom, 7] = xbar[k]

                    if CC>20:
                        flag=False

                    if flag:

                        if alpha==.05:
                            crit=smmcrit(df,CC)

                        if alpha==.01:
                            crit=smmcrit01(df,CC)

                    elif not flag or CC>28:
                        dfvec = np.repeat(df, CC)
                        crit=smmvalv2(dfvec, alpha, seed=seed)

                    test[jcom, 3] = crit
                    psihat[jcom, 3] = (xbar[j] - xbar[k]) - crit * sejk
                    psihat[jcom, 4] = (xbar[j] - xbar[k]) + crit * sejk
                    jcom+=1

        results_test = pd.DataFrame({'group_x': test[:, 0],
                                     'group_y': test[:,1],
                                     'test': test[:,2],
                                     'crit': test[:, 3],
                                     'se': test[:, 4],
                                     'df': test[:, 5]
                                     })

        results_psihat = pd.DataFrame({'group_x': psihat[:, 0],
                                       'group_y': psihat[:, 1],
                                       'psihat': psihat[:, 2],
                                       'ci_low': psihat[:,3],
                                       'ci_up': psihat[:, 4],
                                       'p_value': psihat[:, 5],
                                       'est_1': psihat[:, 6],
                                       'est_2': psihat[:, 7]
                                       })

        return {'test': results_test, 'psihat': results_psihat}

    elif type(con) is np.ndarray:

        if con.shape[0] != len(x):
            raise Exception("The number of groups does not match the number of contrast coefficients.")

        psihat = np.zeros([con.shape[1], 5])
        test = np.zeros([con.shape[1], 5])

        for d in range(con.shape[1]):
            psihat[d, 0] = d
            psihat[d, 1] = sum(con[:, d] * xbar)
            sejk = np.sqrt(sum(con[:, d] ** 2 * w))
            test[d, 0] = d
            test[d, 1] = sum(con[:, d] * xbar) / sejk
            df = (sum(con[:, d] ** 2 * w)) ** 2 / sum(con[:, d] ** 4 * w ** 2 / (h - 1))

            if flag:

                if alpha == .05:
                    crit = smmcrit(df, con.shape[1])

                if alpha == .01:
                    crit = smmcrit01(df, con.shape[1])

            if not flag:
                dfvec = np.repeat(df, con.shape[1])
                crit = smmvalv2(dfvec, alpha, seed)

            test[d, 2] = crit
            test[d, 3] = sejk
            test[d, 4] = df
            psihat[d, 2] = psihat[d, 1] - crit * sejk
            psihat[d, 3] = psihat[d, 1] + crit * sejk
            psihat[d, 4] = 2 * (1 - t.cdf(abs(test[d, 1]), df))


        results_test = pd.DataFrame({'con_num': test[:, 0],
                                     'test': test[:,1],
                                     'crit': test[:,2],
                                     'se': test[:, 3],
                                     'df': test[:, 4]})

        results_psihat = pd.DataFrame({'con_num': psihat[:, 0],
                                       'psihat': psihat[:, 1],
                                       'ci_low': psihat[:,2],
                                       'ci_up': psihat[:, 3],
                                       'p_value': psihat[:, 4]
                                       })

        results = {'n': sam, 'test': results_test, 'psihat': results_psihat}

        return results

def smmcrit(nuhat, C):

    """
    Determine the .95 quantile of the C-variate Studentized maximum
    modulus distribution using linear interpolation on inverse
    degrees of freedom
    If C=1, this function returns the .975 quantile of Student's t
    distribution.

    :param nuhat: df
    :param C: number of contrasts
    :return: critical value
    """

    if C - round(C) != 0:
        raise Exception("The number of contrasts, C, must be an  integer")

    if C >= 29:
        raise Exception("C must be less than or equal to 28")

    if C <= 0:
        raise Exception("C must be greater than or equal to 1")

    if nuhat < 2:
        raise Exception("The degrees of freedom must be greater than or equal to 2")

    if C == 1:
        smmcrit_result = t.ppf(.975, nuhat)

    if C >= 2:
        C = C - 2 # 1 less than R to get proper column indexing
        m1 = np.zeros([20, 27])

        m1[0, :] = [
        5.57,
        6.34,
        6.89,
        7.31,
        7.65,
        7.93,
        8.17,
        8.83,
        8.57,
        8.74,
        8.89,
        9.03,
        9.16,
        9.28,
        9.39,
        9.49,
        9.59,
        9.68,
        9.77,
        9.85,
        9.92,
        10.00,
        10.07,
        10.13,
        10.20,
        10.26,
        10.32
        ]
        m1[1, :] = [
        3.96,
        4.43,
        4.76,
        5.02,
        5.23,
        5.41,
        5.56,
        5.69,
        5.81,
        5.92,
        6.01,
        6.10,
        6.18,
        6.26,
        6.33,
        6.39,
        6.45,
        6.51,
        6.57,
        6.62,
        6.67,
        6.71,
        6.76,
        6.80,
        6.84,
        6.88,
        6.92
        ]
        m1[2, :] = [
        3.38,
        3.74,
        4.01,
        4.20,
        4.37,
        4.50,
        4.62,
        4.72,
        4.82,
        4.89,
        4.97,
        5.04,
        5.11,
        5.17,
        5.22,
        5.27,
        5.32,
        5.37,
        5.41,
        5.45,
        5.49,
        5.52,
        5.56,
        5.59,
        5.63,
        5.66,
        5.69
        ]
        m1[3, :] = [
        3.09,
        3.39,
        3.62,
        3.79,
        3.93,
        4.04,
        4.14,
        4.23,
        4.31,
        4.38,
        4.45,
        4.51,
        4.56,
        4.61,
        4.66,
        4.70,
        4.74,
        4.78,
        4.82,
        4.85,
        4.89,
        4.92,
        4.95,
        4.98,
        5.00,
        5.03,
        5.06
        ]
        m1[4, :] = [
        2.92,
        3.19,
        3.39,
        3.54,
        3.66,
        3.77,
        3.86,
        3.94,
        4.01,
        4.07,
        4.13,
        4.18,
        4.23,
        4.28,
        4.32,
        4.36,
        4.39,
        4.43,
        4.46,
        4.49,
        4.52,
        4.55,
        4.58,
        4.60,
        4.63,
        4.65,
        4.68
        ]
        m1[5, :] = [
        2.80,
        3.06,
        3.24,
        3.38,
        3.49,
        3.59,
        3.67,
        3.74,
        3.80,
        3.86,
        3.92,
        3.96,
        4.01,
        4.05,
        4.09,
        4.13,
        4.16,
        4.19,
        4.22,
        4.25,
        4.28,
        4.31,
        4.33,
        4.35,
        4.38,
        4.39,
        4.42
        ]
        m1[6, :] = [
        2.72,
        2.96,
        3.13,
        3.26,
        3.36,
        3.45,
        3.53,
        3.60,
        3.66,
        3.71,
        3.76,
        3.81,
        3.85,
        3.89,
        3.93,
        3.96,
        3.99,
        4.02,
        4.05,
        4.08,
        4.10,
        4.13,
        4.15,
        4.18,
        4.19,
        4.22,
        4.24
        ]
        m1[7, :] = [
        2.66,
        2.89,
        3.05,
        3.17,
        3.27,
        3.36,
        3.43,
        3.49,
        3.55,
        3.60,
        3.65,
        3.69,
        3.73,
        3.77,
        3.80,
        3.84,
        3.87,
        3.89,
        3.92,
        3.95,
        3.97,
        3.99,
        4.02,
        4.04,
        4.06,
        4.08,
        4.09
        ]
        m1[8, :] = [
        2.61,
        2.83,
        2.98,
        3.10,
        3.19,
        3.28,
        3.35,
        3.41,
        3.47,
        3.52,
        3.56,
        3.60,
        3.64,
        3.68,
        3.71,
        3.74,
        3.77,
        3.79,
        3.82,
        3.85,
        3.87,
        3.89,
        3.91,
        3.94,
        3.95,
        3.97,
        3.99
        ]
        m1[9, :] = [
        2.57,
        2.78,
        2.93,
        3.05,
        3.14,
        3.22,
        3.29,
        3.35,
        3.40,
        3.45,
        3.49,
        3.53,
        3.57,
        3.60,
        3.63,
        3.66,
        3.69,
        3.72,
        3.74,
        3.77,
        3.79,
        3.81,
        3.83,
        3.85,
        3.87,
        3.89,
        3.91
        ]
        m1[10, :] = [
        2.54,
        2.75,
        2.89,
        3.01,
        3.09,
        3.17,
        3.24,
        3.29,
        3.35,
        3.39,
        3.43,
        3.47,
        3.51,
        3.54,
        3.57,
        3.60,
        3.63,
        3.65,
        3.68,
        3.70,
        3.72,
        3.74,
        3.76,
        3.78,
        3.80,
        3.82,
        3.83
        ]
        m1[11, :] = [
        2.49,
        2.69,
        2.83,
        2.94,
        3.02,
        3.09,
        3.16,
        3.21,
        3.26,
        3.30,
        3.34,
        3.38,
        3.41,
        3.45,
        3.48,
        3.50,
        3.53,
        3.55,
        3.58,
        3.59,
        3.62,
        3.64,
        3.66,
        3.68,
        3.69,
        3.71,
        3.73
        ]
        m1[12, :] = [
        2.46,
        2.65,
        2.78,
        2.89,
        2.97,
        3.04,
        3.09,
        3.15,
        3.19,
        3.24,
        3.28,
        3.31,
        3.35,
        3.38,
        3.40,
        3.43,
        3.46,
        3.48,
        3.50,
        3.52,
        3.54,
        3.56,
        3.58,
        3.59,
        3.61,
        3.63,
        3.64
        ]
        m1[13, :] = [
        2.43,
        2.62,
        2.75,
        2.85,
        2.93,
        2.99,
        3.05,
        3.11,
        3.15,
        3.19,
        3.23,
        3.26,
        3.29,
        3.32,
        3.35,
        3.38,
        3.40,
        3.42,
        3.44,
        3.46,
        3.48,
        3.50,
        3.52,
        3.54,
        3.55,
        3.57,
        3.58
        ]
        m1[14, :] = [
        2.41,
        2.59,
        2.72,
        2.82,
        2.89,
        2.96,
        3.02,
        3.07,
        3.11,
        3.15,
        3.19,
        3.22,
        3.25,
        3.28,
        3.31,
        3.33,
        3.36,
        3.38,
        3.39,
        3.42,
        3.44,
        3.46,
        3.47,
        3.49,
        3.50,
        3.52,
        3.53
        ]
        m1[15, :] = [
        2.38,
        2.56,
        2.68,
        2.77,
        2.85,
        2.91,
        2.97,
        3.02,
        3.06,
        3.09,
        3.13,
        3.16,
        3.19,
        3.22,
        3.25,
        3.27,
        3.29,
        3.31,
        3.33,
        3.35,
        3.37,
        3.39,
        3.40,
        3.42,
        3.43,
        3.45,
        3.46
        ]
        m1[16, :] = [
        2.35,
        2.52,
        2.64,
        2.73,
        2.80,
        2.87,
        2.92,
        2.96,
        3.01,
        3.04,
        3.07,
        3.11,
        3.13,
        3.16,
        3.18,
        3.21,
        3.23,
        3.25,
        3.27,
        3.29,
        3.30,
        3.32,
        3.33,
        3.35,
        3.36,
        3.37,
        3.39
        ]
        m1[17, :] = [
        2.32,
        2.49,
        2.60,
        2.69,
        2.76,
        2.82,
        2.87,
        2.91,
        2.95,
        2.99,
        3.02,
        3.05,
        3.08,
        3.09,
        3.12,
        3.14,
        3.17,
        3.18,
        3.20,
        3.22,
        3.24,
        3.25,
        3.27,
        3.28,
        3.29,
        3.31,
        3.32
]
        m1[18, :] = [
        2.29,
        2.45,
        2.56,
        2.65,
        2.72,
        2.77,
        2.82,
        2.86,
        2.90,
        2.93,
        2.96,
        2.99,
        3.02,
        3.04,
        3.06,
        3.08,
        3.10,
        3.12,
        3.14,
        3.16,
        3.17,
        3.19,
        3.20,
        3.21,
        3.23,
        3.24,
        3.25
        ]
        m1[19, :] = [
        2.24,
        2.39,
        2.49,
        2.57,
        2.63,
        2.68,
        2.73,
        2.77,
        2.79,
        2.83,
        2.86,
        2.88,
        2.91,
        2.93,
        2.95,
        2.97,
        2.98,
        3.01,
        3.02,
        3.03,
        3.04,
        3.06,
        3.07,
        3.08,
        3.09,
        3.11,
        3.12
        ]

        if nuhat >= 200:
            smmcrit_result = m1[19, C]

        if nuhat < 200:
            nu = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 24, 30, 40, 60, 200])
            temp = abs(nu - nuhat)
            find = np.argsort(temp)

        if temp[find[0]] == 0:
            smmcrit_result = m1[find[0], C]

        if temp[find[0]] != 0:
            if nuhat > nu[find[0]]:
                smmcrit_result = m1[find[0], C] - (1 / nu[find[0]] - 1 / nuhat) * (m1[find[0], C] - m1[find[0] + 1, C]) /\
                          (1 / nu[find[0]] - 1 / nu[find[0] + 1])

        if nuhat < nu[find[0]]:
            smmcrit_result = m1[find[0] - 1, C] - (1 / nu[find[0] - 1] - 1 / nuhat) * (m1[find[0] - 1, C] - m1[find[0], C]) /\
                      (1 / nu[find[0] - 1] - 1 / nu[find[0]])

    return smmcrit_result

def smmcrit01(nuhat, C):

    """
    Determine the .99 quantile of the C-variate Studentized maximum
    modulus distribution using linear interpolation on inverse
    degrees of freedom
    If C=1, this function returns the .975 quantile of Student's t
    distribution.

    :param nuhat: df
    :param C: number of contrasts
    :return: critical value
    """

    if C - round(C) != 0:
        raise Exception("The number of contrasts, C, must be an  integer")

    if C >= 29:
        raise Exception("C must be less than or equal to 28")

    if C <= 0:
        raise Exception("C must be greater than or equal to 1")

    if nuhat < 2:
        raise Exception("The degrees of freedom must be greater than or equal to 2")

    if C == 1:
        smmcrit01_result = t.ppf(.995, nuhat)

    if C >= 2:

        C = C - 2 # 1 less than R to get proper column indexing
        m1 = np.zeros([20, 27])

        m1[0, :] = [
            12.73,
            14.44,
            15.65,
            16.59,
            17.35,
            17.99,
            18.53,
            19.01,
            19.43,
            19.81,
            20.15,
            20.46,
            20.75,
            20.99,
            20.99,
            20.99,
            20.99,
            20.99,
            22.11,
            22.29,
            22.46,
            22.63,
            22.78,
            22.93,
            23.08,
            23.21,
            23.35
        ]
        m1[1, :] = [
            7.13,
            7.91,
            8.48,
            8.92,
            9.28,
            9.58,
            9.84,
            10.06,
            10.27,
            10.45,
            10.61,
            10.76,
            10.90,
            11.03,
            11.15,
            11.26,
            11.37,
            11.47,
            11.56,
            11.65,
            11.74,
            11.82,
            11.89,
            11.97,
            12.07,
            12.11,
            12.17
        ]
        m1[2, :] = [
            5.46,
            5.99,
            6.36,
            6.66,
            6.89,
            7.09,
            7.27,
            7.43,
            7.57,
            7.69,
            7.80,
            7.91,
            8.01,
            8.09,
            8.17,
            8.25,
            8.32,
            8.39,
            8.45,
            8.51,
            8.57,
            8.63,
            8.68,
            8.73,
            8.78,
            8.83,
            8.87
        ]
        m1[3,:]=[
            4.70,
            5.11,
            5.39,
            5.63,
            5.81,
            5.97,
            6.11,
            6.23,
            6.33,
            6.43,
            6.52,
            6.59,
            6.67,
            6.74,
            6.81,
            6.87,
            6.93,
            6.98,
            7.03,
            7.08,
            7.13,
            7.17,
            7.21,
            7.25,
            7.29,
            7.33,
            7.36
        ]
        m1[4,:]=[
            4.27,
            4.61,
            4.85,
            5.05,
            5.20,
            5.33,
            5.45,
            5.55,
            5.64,
            5.72,
            5.79,
            5.86,
            5.93,
            5.99,
            6.04,
            6.09,
            6.14,
            6.18,
            6.23,
            6.27,
            6.31,
            6.34,
            6.38,
            6.41,
            6.45,
            6.48,
            6.51
        ]
        m1[5,:]=[
            3.99,
            4.29,
            4.51,
            4.68,
            4.81,
            4.93,
            5.03,
            5.12,
            5.19,
            5.27,
            5.33,
            5.39,
            5.45,
            5.50,
            5.55,
            5.59,
            5.64,
            5.68,
            5.72,
            5.75,
            5.79,
            5.82,
            5.85,
            5.88,
            5.91,
            5.94,
            5.96
        ]
        m1[6,:]=[
            3.81,
            4.08,
            4.27,
            4.42,
            4.55,
            4.65,
            4.74,
            4.82,
            4.89,
            4.96,
            5.02,
            5.07,
            5.12,
            5.17,
            5.21,
            5.25,
            5.29,
            5.33,
            5.36,
            5.39,
            5.43,
            5.45,
            5.48,
            5.51,
            5.54,
            5.56,
            5.59
        ]
        m1[7,:]=[
            3.67,
            3.92,
            4.10,
            4.24,
            4.35,
            4.45,
            4.53,
            4.61,
            4.67,
            4.73,
            4.79,
            4.84,
            4.88,
            4.92,
            4.96,
            5.01,
            5.04,
            5.07,
            5.10,
            5.13,
            5.16,
            5.19,
            5.21,
            5.24,
            5.26,
            5.29,
            5.31
        ]
        m1[8,:]=[
            3.57,
            3.80,
            3.97,
            4.09,
            4.20,
            4.29,
            4.37,
            4.44,
            4.50,
            4.56,
            4.61,
            4.66,
            4.69,
            4.74,
            4.78,
            4.81,
            4.84,
            4.88,
            4.91,
            4.93,
            4.96,
            4.99,
            5.01,
            5.03,
            5.06,
            5.08,
            5.09
        ]
        m1[9,:]=[
            3.48,
            3.71,
            3.87,
            3.99,
            4.09,
            4.17,
            4.25,
            4.31,
            4.37,
            4.42,
            4.47,
            4.51,
            4.55,
            4.59,
            4.63,
            4.66,
            4.69,
            4.72,
            4.75,
            4.78,
            4.80,
            4.83,
            4.85,
            4.87,
            4.89,
            4.91,
            4.93
        ]
        m1[10,:]=[
            3.42,
            3.63,
            3.78,
            3.89,
            3.99,
            4.08,
            4.15,
            4.21,
            4.26,
            4.31,
            4.36,
            4.40,
            4.44,
            4.48,
            4.51,
            4.54,
            4.57,
            4.59,
            4.62,
            4.65,
            4.67,
            4.69,
            4.72,
            4.74,
            4.76,
            4.78,
            4.79
        ]
        m1[11,:]=[
            3.32,
            3.52,
            3.66,
            3.77,
            3.85,
            3.93,
            3.99,
            4.05, # bug in WRS
            4.10,
            4.15,
            4.19,
            4.23,
            4.26,
            4.29,
            4.33,
            4.36,
            4.39,
            4.41,
            4.44,
            4.46,
            4.48,
            4.50,
            4.52,
            4.54,
            4.56,
            4.58,
            4.59
        ]
        m1[12,:]=[
            3.25,
            3.43,
            3.57,
            3.67,
            3.75,
            3.82,
            3.88,
            3.94,
            3.99,
            4.03,
            4.07,
            4.11,
            4.14,
            4.17,
            4.19,
            4.23,
            4.25,
            4.28,
            4.29,
            4.32,
            4.34,
            4.36,
            4.38,
            4.39,
            4.42,
            4.43,
            4.45
        ]
        m1[13,:]=[
            3.19,
            3.37,
            3.49,
            3.59,
            3.68,
            3.74,
            3.80,
            3.85,
            3.89,
            3.94,
            3.98,
            4.01,
            4.04,
            4.07,
            4.10,
            4.13,
            4.15,
            4.18,
            4.19,
            4.22,
            4.24,
            4.26,
            4.28,
            4.29,
            4.31,
            4.33,
            4.34
        ]
        m1[14,:]=[
            3.15,
            3.32,
            3.45,
            3.54,
            3.62,
            3.68,
            3.74,
            3.79,
            3.83,
            3.87,
            3.91,
            3.94,
            3.97,
            3.99,
            4.03,
            4.05,
            4.07,
            4.09,
            4.12,
            4.14,
            4.16,
            4.17,
            4.19,
            4.21,
            4.22,
            4.24,
            4.25
        ]
        m1[15,:]=[
            3.09,
            3.25,
            3.37,
            3.46,
            3.53,
            3.59,
            3.64,
            3.69,
            3.73,
            3.77,
            3.80,
            3.83,
            3.86,
            3.89,
            3.91,
            3.94,
            3.96,
            3.98,
            4.00,
            4.02,
            4.04,
            4.05,
            4.07,
            4.09,
            4.10,
            4.12,
            4.13
        ]
        m1[16,:]=[
            3.03,
            3.18,
            3.29,
            3.38,
            3.45,
            3.50,
            3.55,
            3.59,
            3.64,
            3.67,
            3.70,
            3.73,
            3.76,
            3.78,
            3.81,
            3.83,
            3.85,
            3.87,
            3.89,
            3.91,
            3.92,
            3.94,
            3.95,
            3.97,
            3.98,
            4.00,
            4.01
        ]
        m1[17,:]=[
            2.97,
            3.12,
            3.22,
            3.30,
            3.37,
            3.42,
            3.47,
            3.51,
            3.55,
            3.58,
            3.61,
            3.64,
            3.66,
            3.68,
            3.71,
            3.73,
            3.75,
            3.76,
            3.78,
            3.80,
            3.81,
            3.83,
            3.84,
            3.85,
            3.87,
            3.88,
            3.89
        ]
        m1[18,:]=[
            2.91,
            3.06,
            3.15,
            3.23,
            3.29,
            3.34,
            3.38,
            3.42,
            3.46,
            3.49,
            3.51,
            3.54,
            3.56,
            3.59,
            3.61,
            3.63,
            3.64,
            3.66,
            3.68,
            3.69,
            3.71,
            3.72,
            3.73,
            3.75,
            3.76,
            3.77,
            3.78
        ]
        m1[19,:]=[
            2.81,
            2.93,
            3.02,
            3.09,
            3.14,
            3.19,
            3.23,
            3.26,
            3.29,
            3.32,
            3.34,
            3.36,
            3.38,
            3.40,
            3.42, # bug in WRS
            3.44, # bug in WRS
            3.45,
            3.47,
            3.48,
            3.49,
            3.50,
            3.52,
            3.53,
            3.54,
            3.55,
            3.56,
            3.57
        ]

        if nuhat >= 200:
            smmcrit01_result = m1[19, C]

        if nuhat < 200:
            nu = np.array([2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 20, 24, 30, 40, 60, 200])
            temp = abs(nu - nuhat)
            find = np.argsort(temp)

        if temp[find[0]] == 0:
            smmcrit01_result = m1[find[0], C]

        if temp[find[0]] != 0:
            if nuhat > nu[find[0]]:
                smmcrit01_result = m1[find[0], C] - (1 / nu[find[0]] - 1 / nuhat) * (m1[find[0], C] - m1[find[0] + 1, C]) /\
                          (1 / nu[find[0]] - 1 / nu[find[0] + 1])

        if nuhat < nu[find[0]]:
            smmcrit01_result = m1[find[0] - 1, C] - (1 / nu[find[0] - 1] - 1 / nuhat) * (m1[find[0] - 1, C] - m1[find[0], C]) /\
                      (1 / nu[find[0] - 1] - 1 / nu[find[0]])


    return smmcrit01_result

def smmvalv2(dfvec, alpha, seed, iters=20000):

    if seed:
        np.random.seed(seed)

    J = len(dfvec)
    z=np.empty((iters, J))
    z[:] = np.nan

    for j in range(J):
        z[:, j] = abs(np.random.standard_t(dfvec[j], size=iters))

    vals = np.asarray([np.max(x) for x in z])
    vals = np.sort(vals)
    ival = round((1 - alpha) * iters)
    qval = vals[ival]

    return qval

def trimparts(x, tr):

    """
    Compute the trimmed mean, effective sample size, and squared standard error.

    :param x: the data array
    :param tr: amount to be trimmed
    :return:
    """

    tm = trim_mean(x, tr)
    h1 = len(x) - 2 * np.floor(tr * len(x))
    sqse = (len(x) - 1) * winvar(x, tr) / (h1 * (h1 - 1))
    trimparts_results=[tm, sqse]

    return trimparts_results

def trimpartt(x, con):
    return sum(con * x)

def covmtrim(x, tr):

    p=len(x)
    n=len(x[0])
    h=len(x[0]) - 2 * np.floor(tr * len(x[0]))
    covest = np.zeros([p, p])
    covest[0, 0] = (n - 1) * winvar(x[0], tr) / (h * (h - 1))

    for j in range(1,p):

        covest[j, j] = (n - 1) * winvar(x[j], tr) / (h * (h - 1))

        for k in range(j):

            covest[j, k] = (n - 1) * \
                wincor(x[j], x[k], tr)['wcov'] / (h * (h - 1))

            covest[k, j] = covest[j, k]


    return covest

def lindep(x, con, cmat, tr):

    J = len(x)
    xbar = np.zeros([1, J])

    for j in range(J):
        xbar[0,j]=trim_mean(x[j], tr)

    ncon=con.shape[1]
    psihat = np.zeros([ncon,4])
    w = cmat

    for d in range(ncon):
        psihat[d,0]=d
        psihat[d,1]=np.sum(con[:,d]*xbar)
        sejk = np.sqrt(con[:,d].T @ w @ con[:,d])
        psihat[d,2]=sejk
        psihat[d,3]=psihat[d,1]/sejk

    res=pd.DataFrame(psihat, columns=['con_num', 'psihat', 'se', 'test'])

    return res

def array_to_list_of_arrays(x):
    pass

def list_of_arrays_to_array(x):
    pass

def create_random_2_factor_data(jth_ns, k_levels):

    # [20,15,11]
    # 4
    # creates 3x4 design

    x=[]
    for n in jth_ns:

        data = [np.random.rand(n) for _ in range(k_levels)]
        #data = [i for i in data]
        x+=data



    for i, xi in enumerate(x):
        np.save(f'/home/allan/test_{i + 1}.npy', xi)
        #x[np.random.randint(0, len(x))]=np.nan

    return x

def remove_nans_across_dependent_groups(x, J, K):


    ind_low=0
    ind_up=K
    all_data=[]
    for j in range(J):

        x_slice=np.c_[x[ind_low:ind_up]]
        x_slice = x_slice.T[~np.isnan(x_slice.T).any(axis=1)]
        xx = [i for i in x_slice.T]
        all_data += xx
        ind_low+=K
        ind_up+=K


    return all_data
