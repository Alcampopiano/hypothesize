__all__ = ["bwmcp", "bwamcp"]

import numpy as np
import pandas as pd
from scipy.stats import trim_mean
from hypothesize.utilities import con2way, lindep, covmtrim, lincon
from hypothesize.utilities import remove_nans_across_dependent_groups, pandas_to_arrays
import more_itertools as mit
# np.set_printoptions(linewidth=300)

def bwmcp(J, K, x, alpha=.05, tr=.2, nboot=599, seed=False):

    """
    A bootstrap-t for multiple comparisons among
    for all main effects and interactions in a between-by-within design.
    The analysis is done by generating bootstrap samples and
    using an appropriate linear contrast.

    The variable x is a Pandas DataFrame where the first column
    contains the data for the first level of both factors: level 1,1.
    The second column contains the data for level 1 of the
    first factor and level 2 of the second: level 1,2.
    x.iloc[:,K] is the data for level 1,K. x.iloc[:,K+1] is the data for level 2,1.
    x.iloc[:, 2K] is level 2,K, etc.

    :param J: Number of groups in Factor A
    :param K: Number of groups in Factor B
    :param x: Pandas DataFrame
    :param alpha: significance level
    :param nboot: number of bootstrap samples
    :param tr: amount to trim
    :param seed: random seed for reproducibility
    :return: results dictionary

    """

    x=pandas_to_arrays(x)

    conA, conB, conAB=con2way(J,K)
    p=J*K
    v=np.zeros([p,p])
    data=[xi-trim_mean(xi, tr) for xi in x]
    x=remove_nans_across_dependent_groups(x,J,K)

    ilow=0
    iup=K
    for j in range(J):
        v[ilow:iup,ilow:iup]=covmtrim(x[ilow:iup], tr)
        ilow+=K
        iup+=K

    A = lindep(x, conA, cmat = v, tr = tr)
    B = lindep(x, conB, cmat=v, tr=tr)
    AB = lindep(x, conAB, cmat=v, tr=tr)

    testA = []
    testB = []
    testAB = []
    aboot = np.empty([nboot, conA.shape[1]])
    bboot = np.empty([nboot, conB.shape[1]])
    abboot = np.empty([nboot, conAB.shape[1]])

    if seed:
        np.random.seed(seed)

    for ib in range(nboot):

        bsam=[]
        ilow=0
        iup=K
        for j in range(J):

            nv=len(x[ilow])
            bdat=np.random.randint(nv, size=nv)

            for k in range(ilow,iup):
                bsam.append(data[k][bdat])

            ilow+=K
            iup+=K

        ilow=0
        iup=K
        for j in range(J):
            v[ilow:iup, ilow:iup] = covmtrim(bsam[ilow:iup], tr)
            ilow += K
            iup += K

        temp = abs(lindep(bsam, conA, cmat = v, tr = tr)['test'])
        aboot[ib,:] = temp
        testA.append(max(temp))
        temp = abs(lindep(bsam, conB, cmat = v, tr = tr)['test'])
        bboot[ib,:] = temp
        testB.append(max(temp))
        temp = abs(lindep(bsam, conAB, cmat = v, tr = tr)['test'])
        abboot[ib,:] = temp
        testAB.append(max(temp))

    pbA = []
    pbB = []
    pbAB = []
    for j in range(aboot.shape[1]):
        pbA.append(np.mean((abs(A['test'][j]) < aboot[:, j])))

    for j in range(bboot.shape[1]):
        pbB.append(np.mean((abs(B['test'][j]) < bboot[:, j])))

    for j in range(abboot.shape[1]):
        pbAB.append(np.mean((abs(AB['test'][j]) < abboot[:, j])))

    critA = sorted(testA)
    critB = sorted(testB)
    critAB = sorted(testAB)
    ic = int(np.floor((1-alpha) * nboot)) - 1
    critA = critA[ic]
    critB = critB[ic]
    critAB = critAB[ic]

    A['crit_value'] = critA
    A=pd.concat([A, pd.Series(pbA, name='p_value')], axis=1)

    B['crit_value'] = critB
    B=pd.concat([B, pd.Series(pbB, name='p_value')], axis=1)

    AB['crit_value'] = critAB
    AB=pd.concat([AB, pd.Series(pbAB, name='p_value')], axis=1)

    res={'factor_A': A, 'factor_B': B, 'factor_AB': AB,
         'contrast_coef': {'conA': conA, 'conB': conB, 'conAB': conAB}}

    return res

def bwamcp(J, K, x, tr=.2, alpha=.05, pool=False):


    """
    All pairwise comparisons among levels of Factor A
    in a mixed design using trimmed means.

    The variable x is a Pandas DataFrame where the first column
    contains the data for the first level of both factors: level 1,1.
    The second column contains the data for level 1 of the
    first factor and level 2 of the second: level 1,2.
    x.iloc[:,K] is the data for level 1,K. x.iloc[:,K+1] is the data for level 2,1.
    x.iloc[:, 2K] is level 2,K, etc.

    :param J: number of levels for factor A
    :param K: number of levels for factor B
    :param x: Pandas DataFrame
    :param tr: amount of trimming
    :param alpha: alpha level
    :param pool: if "True", pool dependent groups together.
        Otherwise generate pairwise contrasts across factor A for each level of factor B

    :return:
    """

    x = pandas_to_arrays(x)
    x = remove_nans_across_dependent_groups(x, J, K)

    if pool:
        data = [np.concatenate(x[i:i + K]) for i in range(0, J * K, K)]
        results = lincon(data, con=None, tr=tr, alpha=alpha)

    elif not pool:

        MJK = K * (J ** 2 - J) // 2
        c=np.zeros([J*K,MJK])
        n_idioms=J-1
        idioms=[]
        K_mult=K
        for i in range(n_idioms):
            tmp=np.concatenate([[1], np.repeat(0, K_mult-1), [-1]])
            idioms.append(tmp)
            K_mult*=2

        col_ind=0
        for idiom in idioms:
            num_rep_idiom=len(list(mit.windowed(c[:,0], n=len(idiom))))

            row_start=0
            for _ in range(num_rep_idiom):
                c[row_start:row_start+len(idiom),col_ind]=idiom
                row_start+=1
                col_ind+=1

        results = lincon(x, con=c, tr=tr, alpha=alpha)

    return results












