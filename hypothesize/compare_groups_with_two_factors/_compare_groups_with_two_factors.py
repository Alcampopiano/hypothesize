__all__ = ["bwmcp"]

import numpy as np
import pandas as pd
from scipy.stats import trim_mean
from hypothesize.utilities import con2way, lindep, covmtrim
from hypothesize.utilities import remove_nans_across_dependent_groups

def bwmcp(J, K, x, alpha, nboot, tr=.2, seed=False):

    """
    A bootstrap-t for multiple comparisons among
    for all main effects and interactions in a between-by-within design.
    The analysis is done by generating bootstrap samples and
    using an appropriate linear contrast.

    The variable x is assumed to contain the raw
    data stored in an array.
    x[[1]] contains the data
    for the first level of both factors: level 1,1.
    x[[2]] is assumed to contain the data for level 1 of the
    first factor and level 2 of the second: level 1,2
    x[[K]] is the data for level 1,K
    x[[K+1]] is the data for level 2,1, x[[2K]] is level 2,K, etc.

    :param J: Number of groups in Factor A
    :param K: Number of groups in Factor B
    :param x: data array
    :param con: contrast matrix
    :param alpha: significance level
    :param nboot: number of bootstrap samples
    :param tr: amount to trim
    :param seed: random seed for reproducibility
    :return:


    for dev:
    create fake data
    x=np.asarray([np.random.rand(20) for i in range(6)])
    np.save('/home/allan/temp.npy', x)

    In R:
    install.packages("RcppCNPy")
    library("RcppCNPy")
    x=npyLoad("/home/allan/temp.npy")
    x=listm(t(x))

    """

    conA, conB, conAB=con2way(J,K)
    p=J*K
    v=np.zeros([p,p])
    data=[xi-trim_mean(xi, tr) for xi in x]
    x=remove_nans_across_dependent_groups(x,J,K)

    ilow=0
    iup=K
    for j in range(J):
        v[ilow:iup,ilow:iup]=covmtrim(x[ilow:iup], tr)
        ilow+=iup
        iup+=iup

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
        ilow=0
        iup=K
        bsam=[]

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
            ilow += iup
            iup += iup

        temp = abs(lindep(bsam, conA, cmat = v, tr = tr)['test'])
        aboot[ib,:] = temp
        testA.append(max(temp))
        temp = abs(lindep(bsam, conB, cmat = v, tr = tr)['test'])
        bboot[ib,:] = temp
        testB.append(max(temp))
        temp = abs(lindep(bsam, conAB, cmat = v, tr = tr)['test'])
        abboot[ib,:] = temp
        testAB.append(max(temp))

    # YOU ARE HERE









