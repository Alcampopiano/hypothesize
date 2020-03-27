__all__ = ["bwmcp", "bwamcp", "bwbmcp", "bwimcp"]

import numpy as np
import pandas as pd
from scipy.stats import trim_mean
from hypothesize.utilities import con2way, lindep, covmtrim, \
    lincon, winvar,  trimse
from hypothesize.measuring_associations import wincor
from hypothesize.utilities import remove_nans_based_on_design, pandas_to_arrays
import more_itertools as mit
from scipy.stats import t

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
    x=remove_nans_based_on_design(x, design_values=[J,K], design_type='between_within')

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
    x = remove_nans_based_on_design(x, design_values=[J,K], design_type='between_within')

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

def bwbmcp(J, K, x, tr=.2, con=None, alpha=.05,
           dif=True, pool=False, hoch=False):

    """
    All pairwise comparisons among levels of Factor B
    in a split-plot design using trimmed means.

    Data can be pooled for each level
    of Factor B. This function calls rmmcp.

    The variable x is a Pandas DataFrame where the first column
    contains the data for the first level of both factors: level 1,1.
    The second column contains the data for level 1 of the
    first factor and level 2 of the second: level 1,2.
    x.iloc[:,K] is the data for level 1,K. x.iloc[:,K+1] is the data for level 2,1.
    x.iloc[:, 2K] is level 2,K, etc.

    :param J:
    :param K:
    :param x:
    :param tr:
    :param con:
    :param alpha:
    :param dif:
    :param pool:
    :param hoch:
    :return:
    """

    x = pandas_to_arrays(x)
    x = remove_nans_based_on_design(x, design_values=[J,K], design_type='between_within')

    if pool:
        data = [np.concatenate(x[i:i+J*K+1:K]) for i in range(K)]
        results=rmmcp(data, con=con, tr=tr,alpha=alpha,dif=dif,hoch=hoch)

    else:

        results=[]
        j_ind=0
        for j in range(J):
            data=x[j_ind:j_ind+K]
            tests = rmmcp(data, con=con, tr=tr, alpha=alpha, dif=dif, hoch=hoch)
            results.append(tests)
            j_ind+=K

    return results

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

    from hypothesize.utilities import trimci

    flagcon = False
    x = np.vstack(x).T
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

def bwimcp():
    pass






































