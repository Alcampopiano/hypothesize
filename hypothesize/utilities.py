import numpy as np
from scipy.stats.mstats import winsorize
from scipy.stats import trim_mean
from scipy.stats import t
import pandas as pd
# pd.set_option('display.max_rows', 500)
# pd.set_option('display.max_columns', 500)
# pd.set_option('display.width', 1000)

def con1way(J):

    """
    Return all linear contrasts for J groups


    :param J: int
    Number of groups

    :return:

    array:

    Array of contrasts where the rows correspond to groups and the columns are the contrasts to be used
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
    :param J: int
    Number of levels for Factor A

    :param K: int
    Number of levels for Factor B

    :return:

    list of arrays:

    Each item in the list contains the contrasts for Factor A, Factor B, and the interaction, in that order.
    For each array, the rows correspond to groups and the columns are the contrasts to be used

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

    #results={"conA": conA, "conB": conB, "conAB": conAB}

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

def trimci(x, tr=.2, alpha=.05, null_value=0):

    """
    Compute a 1-alpha confidence interval for the trimmed mean
    The default amount of trimming is tr=.2

    :param x: 1-D array
    :param tr:
    :param alpha:
    :param null_value: The p-value returned by this function is based on the value
        specified by the argument null_value, which defaults to 0
    :return:
    """

    x=x[~np.isnan(x)]
    se = np.sqrt(winvar(x, tr)) / ((1 - 2 * tr) * np.sqrt(len(x)))
    trimci_res = np.zeros(2)
    df = len(x) - 2 * np.floor(tr * len(x)) - 1
    trimci_res[0] = trim_mean(x, tr) - t.ppf(1 - alpha / 2, df) * se
    trimci_res[1] = trim_mean(x, tr) + t.ppf(1 - alpha / 2, df) * se
    test = (trim_mean(x, tr) - null_value) / se
    sig = 2 * (1 - t.cdf(abs(test), df))

    results={"ci": trimci_res, "estimate": trim_mean(x,tr),
             "test_stat": test, "se": se, "p_value": sig,
             "n": len(x)}

    return results

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

    #type(con) is not np.ndarray
    if con is None:

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

        #return {'contrast_matrix': con, 'n': sam, 'test': results_test, 'psihat': results_psihat}

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

    results = {'con': con, 'n': sam, 'test': results_test, 'psihat': results_psihat}

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

    from hypothesize.measuring_associations import wincor

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

def yuen(x, y, tr=.2, alpha=.05):

    """
   Perform Yuen's test for trimmed means on the data in x and y.
   The default amount of trimming is 20%
   Missing values are automatically removed.

   A confidence interval for the trimmed mean of x minus the
   the trimmed mean of y is computed and returned in yuen['ci'].
   The p-value is returned in yuen['p_value']

   x, y: The data for the two groups are stored in x and y
   tr=.2: indicates that the default amount of trimming is .2
          tr=0 results in using the sample mean

   For an omnibus test with more than two independent groups,
   use t1way (may not be implemented yet).

    :param x:
    :param y:
    :param tr:
    :param alpha:
    :return:
    """

    if tr ==.5:
        raise Exception("Using tr=.5 is not allowed; use a method designed for medians "
                        "(they may not be implemented yet")
    if tr>.25:
        raise Warning("with tr>.25 type I error control might be poor")

    x=x[~np.isnan(x)]
    y=y[~np.isnan(y)]

    h1 = len(x) - 2 * np.floor(tr * len(x))
    h2 = len(y) - 2 * np.floor(tr * len(y))
    q1 = (len(x) - 1) * winvar(x, tr) / (h1 * (h1 - 1))
    q2 = (len(y) - 1) * winvar(y, tr) / (h2 * (h2 - 1))
    df = (q1 + q2) ** 2 / ((q1 ** 2 / (h1 - 1)) + (q2 ** 2 / (h2 - 1)))
    crit = t.ppf(1 - alpha / 2, df)
    dif = trim_mean(x, tr) - trim_mean(y, tr)
    low = dif - crit * np.sqrt(q1 + q2)
    up = dif + crit * np.sqrt(q1 + q2)
    test = abs(dif / np.sqrt(q1 + q2))
    yuen_results = 2 * (1 - t.cdf(test, df))


    results={'n1': len(x), 'n2': len(y),
             'est_1': trim_mean(x, tr), 'est_2': trim_mean(y, tr),
             'ci': [low, up], 'p_value': yuen_results,
             'dif': dif, 'se': np.sqrt(q1 + q2),
             'test_stat': test, 'crit': crit,
             'df': df}

    return results

def yuend(x, y, tr=.2, alpha=.05):

    """
     Compare the trimmed means of two dependent random variables
     using the data in x and y.
     The default amount of trimming is 20%

     Any pair with a missing value is eliminated

     A confidence interval for the trimmed mean of x minus the
     the trimmed mean of y is computed and returned in yuend['ci'].
     The significance level is returned in yuend['p_value']

     For inferences based on difference scores, use trimci

    :param x:
    :param y:
    :param tr:
    :param alpha:
    :return:
    """

    from hypothesize.measuring_associations import wincor

    if type(x) is not np.ndarray:
        x, y=pandas_to_arrays([x, y])

    m = np.c_[x, y] # cbind
    m = m[~np.isnan(m).any(axis=1)]
    x = m[:,0]
    y = m[:, 1]

    h1 = len(x) - 2 * np.floor(tr * len(x))
    q1 = (len(x) - 1) * winvar(x, tr)
    q2 = (len(y) - 1) * winvar(y, tr)
    q3 = (len(x) - 1) * wincor(x, y, tr)['wcov']

    df = h1 - 1
    se = np.sqrt((q1 + q2 - 2 * q3) / (h1 * (h1 - 1)))
    crit = t.ppf(1 - alpha / 2, df)
    dif = trim_mean(x, tr) - trim_mean(y, tr)
    low = dif - crit * se
    up = dif + crit * se
    test = dif / se
    yuend_res = 2 * (1 - t.cdf(abs(test), df))

    keys=['ci', 'p_value', 'est1', 'est2',
          'dif', 'se', 'teststat', 'n', 'df']

    vals=[[low, up], yuend_res, trim_mean(x,tr), trim_mean(x,tr),
          dif, se, test, len(x), df]

    return dict(zip(keys,vals))

def pandas_to_arrays(obj):

    if type(obj) is pd.core.frame.DataFrame:
        x=[a for a in obj.values.T]

        return x

    elif type(obj) is list:

        x=[]
        for s in obj:
            x.append(s.values)

        return x

def create_example_data(design_values, missing_data_proportion=0, save_array_path=None, seed=False):

    """
    Create a Pandas DataFrame of random data with a certain number of columns which
    correspond to a design of a given shape (e.g., 1-way, two groups, 2-way design).
    There is also an option to randomly add a proportion of null values.
    The purpose of this function is to make it easy to demonstrate and test the package.

    :param design_values: list or int

    An integer or list indicating the design shape. For example, `[2,3]` indicates
    a 2-by-3 design and will produce a six column DataFrame
    with appropriately named columns

    :param missing_data_proportion: float

    Proportion of randomly missing data

    :param save_array_path: str (default is `None`)

    Save each group as an array for loading into R by specifying a path (e.g. , `'/home/allan/'`).
    If left unset (i.e., `None`), no arrays will be saved.

    :param seed: bool

    Set random seed for reproducible results
    """

    if seed:
        np.random.seed(seed)

    if type(design_values) is not list:
        design_values=[design_values]

    design_values=[i for i in design_values if i != 1]

    cell_str=[]

    if len(design_values) == 1:
        J=design_values[0]

        for j in range(1, J + 1):
            cell_str.append(f'cell_{j}')

        x = np.random.rand(50, J)

    elif len(design_values) == 2:

        J, K=design_values
        for j in range(1,J+1):
            for k in range(1,K+1):
                cell_str.append(f'cell_{j}_{k}')

        x = np.random.rand(50, J*K)

    elif len(design_values) == 3:

        J,K,L=design_values
        for l in range(1,L+1):
            for j in range(1,J+1):
                for k in range(1, K+1):
                    cell_str.append(f'cell_{j}_{k}_{l}')


        x = np.random.rand(50, J*K*L)


    if missing_data_proportion:
        x.ravel()[np.random.choice(x.size, int(round(missing_data_proportion*x.size)), replace=False)] = np.nan

    x=pd.DataFrame(x, columns=cell_str)

    if save_array_path:

        for i, xi in enumerate(x.values.T):
            #xi = xi[~np.isnan(xi)]
            np.save(f'{save_array_path}test_{i + 1}.npy', xi)
            # x[np.random.randint(0, len(x))]=np.nan

    return x

def remove_nans_based_on_design(x, design_values, design_type):

    """
    Remove nans in a way that considers whether or not a factor is within or between subjects.
    That is, for all within subjects cells at a give jth level, remove entire rows where
    any column contains a nan.

    :param x: list of arrays
    :param design_values: Only needed for mixed factorial designs. A list or integer indicating the design shape.
        For example, [2,3] would indicate a 2-by-3 design (2 J levels, 3 K levels).
    :param design_type: One of the following:

        'dependent_groups', 'independent_groups', 'between_within', 'between_between', 'within_within'

    :return:
    """

    if type(design_values) is not list:
        design_values=[design_values]

    design_values=[i for i in design_values if i != 1]

    design_types=['dependent_groups', 'independent_groups', 'between_within', 'between_between', 'within_within']

    if design_type not in design_types:
        print("Please specify an available design type from the following list:")
        print(design_types)
        raise Exception("invalid design type")

    if design_type in ('dependent_groups', 'within_within'):

        #x=np.array([np.c_[x[0:]]])
        x=np.c_[x[0:]]
        x = x.T[~np.isnan(x.T).any(axis=1)]
        x_clean = [i for i in x.T]

    elif design_type in ('independent_groups', 'between_between'):
        x_clean = [i[~np.isnan(i)] for i in x]

    elif design_type == 'between_within':

        J, K = design_values

        ind_low = 0
        ind_up = K
        x_clean = []
        for j in range(J):
            x_slice = np.c_[x[ind_low:ind_up]]
            x_slice = x_slice.T[~np.isnan(x_slice.T).any(axis=1)]
            xx = [i for i in x_slice.T]
            x_clean += xx
            ind_low += K
            ind_up += K


    return x_clean

def bptdpsi(x,con):

    return np.sum(con*x)

def rmanogsub(isub, x, est, *args):
    tsub=est(x[isub], *args)
    return tsub

def rmmcp(x, con=None, tr=.2, alpha=.05, dif=True, hoch=True):

    """
    MCP on trimmed means with FWE controlled with Hochberg's method
    will use Rom's method if alpha=.05 or .01 and number of tests is <=10

    Note: confidence intervals are adjusted based on the corresponding critical p-value.

    :param x: list of arrays
    :param con:
    :param tr:
    :param alpha:
    :param dif:
    :param hoch:
    :return:
    """

    from hypothesize.measuring_associations import wincor

    flagcon = False
    #x = np.vstack(x).T
    J = x.shape[1]
    xbar=np.zeros(J)
    x=x[~np.isnan(x).any(axis=1)]
    nval=x.shape[0]
    nrow=x.shape[0]
    h1 = nrow - 2 * np.floor(tr * nrow)
    df=h1-1

    for j in range(J):
        xbar[j]=trim_mean(x[:,j], .2)

    if con is not None: #np.sum(con**2 != 0): # is not none?
        CC=con.shape[1]

    if con is None: #sum(con**2)==0: # is none?
        CC = (J ** 2 - J) / 2

    ncon=CC

    if alpha==.05:
        dvec=np.array([.05,
          .025,
          .0169,
          .0127,
          .0102,
          .00851,
          .0073,
          .00639,
          .00568,
          .00511])

        if ncon > 10:
            avec=.05/np.arange(ncon,11+1)[::-1]
            dvec = np.concatenate([dvec, avec])

    if alpha==.01:
        dvec=np.array([.01,
        .005,
        .00334,
        .00251,
        .00201,
        .00167,
        .00143,
        .00126,
        .00112,
        .00101])

        if ncon > 10:
            avec=.01/np.arange(ncon,11+1)[::-1]
            dvec = np.concatenate([dvec, avec])

    if hoch:
        dvec = alpha / np.arange(1,ncon+1)

    if alpha != .05 and alpha != .01:
        dvec = alpha / np.arange(1, ncon + 1)

    if con is None: #sum(con**2)==0:
        flagcon=True
        psihat=np.zeros([int(CC),5])
        test=np.full([int(CC),6], np.nan)#np.zeros([int(CC),6])

        temp1=np.empty(0) #np.zeros(J*J) #np.array([0])
        jcom=0
        for j in range(J):
            for k in range(J):
                if j<k:
                    q1 = (nrow - 1) * winvar(x[:, j], tr)
                    q2 = (nrow - 1) * winvar(x[:, k], tr)
                    q3 = (nrow - 1) * wincor(x[:, j], x[:, k], tr)['wcov']
                    sejk = np.sqrt((q1 + q2 - 2 * q3) / (h1 * (h1 - 1)))

                    if not dif:
                        test[jcom, 5] = sejk
                        test[jcom, 2] = (xbar[j] - xbar[k]) / sejk
                        temp1=np.append(temp1, 2 * (1 - t.cdf(abs(test[jcom, 2]), df)))
                        test[jcom, 3] = temp1[jcom]
                        psihat[jcom, 0] = j
                        psihat[jcom, 1] = k
                        test[jcom, 0] = j
                        test[jcom, 1] = k
                        psihat[jcom, 2] = (xbar[j] - xbar[k])

                    elif dif:
                        dv = x[:, j] - x[:, k]
                        test[jcom, 5] = trimse(dv, tr)
                        temp = trimci(dv, alpha=alpha / CC, tr=tr)
                        test[jcom, 2] = temp['test_stat']
                        temp1=np.append(temp1, temp['p_value'])
                        test[jcom, 3] = temp1[jcom]
                        psihat[jcom, 0] = j
                        psihat[jcom, 1] = k
                        test[jcom, 0] = j
                        test[jcom, 1] = k
                        psihat[jcom, 2] = trim_mean(dv, tr)
                        psihat[jcom, 3] = temp['ci'][0]
                        psihat[jcom, 4] = temp['ci'][1]

                    jcom+=1

        if hoch:
            dvec = alpha / np.arange(1, ncon+1)

        temp2=(-temp1).argsort()
        zvec = dvec[0:int(ncon)]
        sigvec = (test[temp2, 3] >= zvec)

        if np.sum(sigvec) < ncon:
            dd = ncon - np.sum(sigvec)
            ddd = np.sum(sigvec) + 1

        test[temp2, 4] = zvec

        if not dif:
            psihat[:, 3] = psihat[:, 2] - t.ppf(1 - test[:, 4] / 2, df) * test[:, 5]
            psihat[:, 4] = psihat[:, 2] + t.ppf(1 - test[:, 4] / 2, df) * test[:, 5]

    elif con is not None: #sum(con ** 2) > 0:

        if con.shape[0] != x.shape[1]:
            raise Exception("The number of groups does not match the number "
                            "of contrast coefficients.")

        ncon = con.shape[1]
        psihat = np.zeros([ncon, 4])
        test = np.zeros([ncon, 5])
        temp1=np.empty(0)

        for d in range(ncon):
            psihat[d,0]=d

            if not dif:
                psihat[d, 1] = np.sum(con[:, d] * xbar)
                sejk=0

                for j in range(J):
                    for k in range(J):
                        djk = (nval - 1) * wincor(x[:, j], x[:, k], tr)['wcov'] / (h1 * (h1 - 1))
                        sejk = sejk + con[j, d] * con[k, d] * djk

                sejk = np.sqrt(sejk)
                test[d, 0] = d
                test[d, 1] = np.sum(con[:, d] * xbar) / sejk
                test[d, 4] = sejk
                temp1 = np.append(temp1, 2 * (1 - t.cdf(abs(test[d, 1]), df)))

            elif dif:

                for j in range(J):
                    if j==0:
                        dval = con[j, d] * x[:,j]

                    elif j>0:
                        dval=dval + con[j,d] * x[:,j]

                temp1=np.append(temp1, trimci(dval, tr=tr)['p_value'])
                test[d, 0] = d
                test[d, 1] = trimci(dval, tr=tr)['test_stat']
                test[d, 4] = trimse(dval, tr=tr)
                psihat[d, 1] = trim_mean(dval, tr)

        test[:, 2] = temp1
        temp2 = (-temp1).argsort()
        zvec = dvec[0:ncon]
        sigvec = (test[temp2, 2] >= zvec)

        if sum(sigvec) < ncon:
            dd = ncon - sum(sigvec)
            ddd = sum(sigvec) + 1

        test[temp2, 3] = zvec
        psihat[:, 2] = psihat[:, 1] - t.ppf(1 - test[:, 3] / 2, df) * test[:, 4]
        psihat[:, 3] = psihat[:, 1] + t.ppf(1 - test[:, 3] / 2, df) * test[:, 4]

    num_sig=test.shape[0]

    if flagcon:
        ior = (-test[:,4]).argsort()

        for j in range(test.shape[0]):
            if test[ior[j], 3] <= test[ior[j], 4]:
                break
            else:
                num_sig=num_sig - 1

    elif not flagcon:
        ior = (-test[:, 3]).argsort()

        for j in range(test.shape[0]):
            if test[ior[j], 2] <= test[ior[j], 3]:
                break
            else:
                num_sig=num_sig - 1

    return {"n": nval, "test": test, "psihat": psihat,
            "con": con, "num_sig": num_sig}

def trimcibt(x, tr=.2, alpha=.05, nboot=None, seed=False):

    """
    Compute a 1-alpha confidence interval for the trimmed mean
    using a bootstrap percentile t method.

    The default amount of trimming is tr=.2.

    During the bootstrapping, the absolute value of the test
    statistic is used (the "two-sided method").

    :param x:
    :param tr:
    :param alpha:
    :param nboot:
    :param seed:
    :return:
    """

    if seed:
        np.random.seed(seed)

    x=x[~np.isnan(x)]
    test=trim_mean(x,tr)/trimse(x,tr)
    data=np.random.choice(x, size=(nboot, len(x)))
    data=data-trim_mean(x,tr)
    top=np.array([trim_mean(row, tr) for row in data])
    bot = np.array([trimse(row, tr) for row in data])
    tval=np.sort(abs(top/bot))
    #tval=abs(tval)
    #tval=np.sort(tval)
    icrit = round((1 - alpha) * nboot) - 1 # one less for python
    #ibot = round(alpha * nboot / 2) - 1 # one less for python
    #itop = nboot - ibot - 2 # -2 for python

    ci_low_ci_up = [trim_mean(x, tr) - tval[icrit] * trimse(x, tr),
                     trim_mean(x, tr) + tval[icrit] * trimse(x, tr)]

    p_value = (np.sum(np.abs(test) <= np.abs(tval))) / nboot

    return {"estimate": trim_mean(x, tr), "ci": ci_low_ci_up,
            "test_stat": test, "p_value": p_value, "n": len(x)}

def llocv2(x, est, *args):

    """

    :param x:
    :param est:
    :param args:
    :return:
    """

    if type(x) is not list:
        val=est(x,*args)

    elif type(x) is list:
        val=np.full(len(x), np.nan)

        for i in range(len(x)):
            val[i]=est(x[i], *args)

    return {'center': val}

def linhat(x, con, est, *args):

    """

    :param x:
    :param con:
    :param est:
    :param args:
    :return:
    """

    psihat = np.full([con.shape[1]], np.nan)
    xbar = llocv2(x, est, *args)['center']
    for i in range(con.shape[1]):
        psihat[i] = np.sum(con[:, i] * xbar)

    return psihat

def mkdocstrings_to_pycharm_docstrings(mkdocstr):

    """
    Convert docstrings from mkdocs site to what pycharm uses

    :param mkdocstr: docstring from a function on mkdocs site
    :return: pycharm-style docstring
    """

    import re

    mkdocstr=re.sub('_(P|p)arameter(s|):_\\n', '', mkdocstr)
    mkdocstr = re.sub('_(R|r)eturn(s|):_\\n', ':return:', mkdocstr)
    mkdocstr = re.sub('(\*\*|\*\*\*)\\b', ':param ', mkdocstr)
    mkdocstr = re.sub('\*\*($|\\s)', '', mkdocstr)

    ind=mkdocstr.find(':return:')

    start_str=mkdocstr[:ind]
    end_str=mkdocstr[ind:]
    end_str=end_str.replace(':param ', '')
    final_str=start_str+end_str

    return print(final_str)

