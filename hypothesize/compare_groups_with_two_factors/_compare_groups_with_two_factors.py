__all__ = ["bwmcp", "bwamcp", "bwbmcp", "bwimcp", "spmcpa", "spmcpb", "spmcpi", "wwmcppb"]

import numpy as np
import pandas as pd
from scipy.stats import trim_mean
from hypothesize.utilities import con2way, lindep, covmtrim, \
    lincon, winvar, trimse, yuen, con1way, bptdpsi, rmanogsub
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
        tests=rmmcp(data, con=con, tr=tr,alpha=alpha,dif=dif,hoch=hoch)

        col_names_test = ["group_x", "group_y", "test", "p_value", "p_crit", "se"]
        col_names_psihat = ["group_x", "group_y", "psihat", "ci_lower", "ci_upper"]
        test_df = pd.DataFrame(tests['test'], columns=col_names_test)
        psihat_df = pd.DataFrame(tests['psihat'], columns=col_names_psihat)

        results = {"test": test_df, 'psihat': psihat_df, 'n': tests['n'], 'con': tests['con'],
                 'num_sig': tests['num_sig']}


    else:

        results=[]
        j_ind=0
        for j in range(J):
            data=x[j_ind:j_ind+K]
            tests = rmmcp(data, con=con, tr=tr, alpha=alpha, dif=dif, hoch=hoch)

            col_names_test=["group_x", "group_y", "test", "p_value", "p_crit", "se"]
            col_names_psihat = ["group_x", "group_y", "psihat", "ci_lower", "ci_upper"]
            test_df=pd.DataFrame(tests['test'], columns=col_names_test)
            psihat_df=pd.DataFrame(tests['psihat'], columns=col_names_psihat)
            tests={"test": test_df, 'psihat': psihat_df, 'n': tests['n'], 'con': tests['con'], 'num_sig': tests['num_sig']}

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

def bwimcp(J, K, x, tr=.2, alpha=.05):

    """
    Multiple comparisons for interactions
    in a split-plot design.
    The analysis is done by taking difference scores
    among all pairs of dependent groups and
    determining which of
    these differences differ across levels of Factor A
    using trimmed means.

    FWE is controlled via Hochberg's method
    To adjusted p-values, use the function p.adjust

    For MOM or M-estimators (possibly not implemented yet),
    use spmcpi which uses a bootstrap method

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
    :param alpha:
    :return:
    """

    x=pandas_to_arrays(x)
    x=remove_nans_based_on_design(x, [J, K], 'between_within')

    MJ = (J ** 2 - J) // 2
    MK = (K ** 2 - K) // 2
    JMK = J * MK
    MJMK = MJ * MK
    Jm = J - 1

    output = np.zeros([MJMK, 7])

    m = np.array(np.arange(J*K)).reshape(J,K)

    ic=0
    test=np.array([])
    for j in range(J):
        for jj in range(J):
            if j < jj:
                for k in range(K):
                    for kk in range(K):
                        if k<kk:
                            output[ic, 0]=j
                            output[ic, 1]=jj
                            output[ic, 2]=k
                            output[ic, 3]=kk
                            x1 = x[m[j, k]] - x[m[j, kk]]
                            x2 = x[m[jj, k]] - x[m[jj, kk]]
                            temp = yuen(x1, x2)
                            output[ic, 4] = trim_mean(x1, tr) - trim_mean(x2, tr)
                            test=np.append(test, temp['p_value'])
                            output[ic, 5] = test[ic]


                            ic+=1

    ncon = len(test)
    dvec = alpha / np.arange(1, ncon+1)
    temp2 = (-test).argsort()
    zvec = dvec[0:ncon]
    #sigvec = (test[temp2] >= zvec)
    output[temp2, 6] = zvec
    output[:, 6] = output[:, 6]


    col_names=["A_x", "A_y", "B_x", "B_y", "psihat", "p_value", "p_crit"]
    results=pd.DataFrame(output, columns=col_names)

    return results

def rmmcppbd(x,  est, *args, alpha=.05, con=None,
             nboot=None, hoch=True, seed=False):

    """
      Use a percentile bootstrap method to compare dependent groups
      based on difference scores.
      By default,
      compute a .95 confidence interval for all linear contrasts
      specified by con, a J by C matrix, where  C is the number of
      contrasts to be tested, and the columns of con are the
      contrast coefficients.
      If con is not specified, all pairwise comparisons are done.

      nboot is the bootstrap sample size. If not specified, a value will
      be chosen depending on the number of contrasts there are.

      A sequentially rejective method is used to control alpha.
      If n>=80, hochberg's method is used.

    :param x:
    :param y:
    :param alpha:
    :param con:
    :param est:
    :param nboot:
    :param hoch:
    :param seed:
    :return:
    """

    x = x[~np.isnan(x).any(axis=1)]
    J=x.shape[1]
    n=x.shape[0]
    if n>=80:
        hoch=True

    #Jm=J-1
    if con is None:
        con=con1way(J)

    d = con.shape[1]
    if not nboot:

      if d <= 10:
          nboot = 3000

      elif d <= 6:
          nboot = 2000

      elif d <= 4:
          nboot = 1000

      else:
          nboot=5000

    connum=d
    xx=x@con

    if seed:
        np.random.seed(seed)

    psihat=np.zeros([connum, nboot])
    data=np.random.randint(n, size=(nboot,n))

    # wilcox's implementation in R is a bit more complicated,
    # I have simplified. Hopefully correctly.
    for ib in range(nboot):
        psihat[:,ib]=est(xx[data[ib,:], :], *args)

    test = np.full(connum, np.nan)
    icl = round(alpha * nboot // 2) #+ 1
    icu = nboot - icl -  2 #- 1
    cimat=np.full([connum, 2], np.nan)

    for ic in range(connum):

      test[ic] =(sum(psihat[ic, :] > 0) + .5 * sum(psihat[ic, :] == 0)) / nboot
      test[ic] = min(test[ic], 1 - test[ic])
      temp = np.sort(psihat[ic, :])
      cimat[ic, 0] = temp[icl]
      cimat[ic, 1] = temp[icu]

    test = 2 * test
    ncon = con.shape[1]

    if alpha == .05:
      dvec =[.025,
          .025,
          .0169,
          .0127,
          .0102,
          .00851,
          .0073,
          .00639,
          .00568,
          .00511]

      if ncon > 10:
        avec = .05 / np.arange(11, ncon+1)
        dvec = np.append(dvec, avec)

    elif alpha == .01:
      dvec =[.005,
          .005,
          .00334,
          .00251,
          .00201,
          .00167,
          .00143,
          .00126,
          .00112,
          .00101]

      if ncon > 10:
        avec = .01 / np.arange(11,ncon+1)
        dvec = np.append(dvec, avec)

    else:
      dvec = alpha / np.arange(1,ncon+1)
      dvec[1] = alpha / 2

    if hoch:
      dvec = alpha / (2 * np.arange(1,ncon+1))

    dvec = 2 * dvec
    temp2 = (-test).argsort()
    ncon = con.shape[1]
    zvec = dvec[:ncon]
    output=np.zeros([connum, 6])

    tmeans=est(xx,*args)
    output[temp2, 3] = zvec

    for ic in range(ncon):
      output[ic, 1] = tmeans[ic]
      output[ic, 0] = ic
      output[ic, 2] = test[ic]
      output[ic, 4:6] = cimat[ic,:]

    num_sig = np.sum(output[:, 2] <= output[:, 3])

    return {"output": output, "con": con, "num_sig": num_sig}

def rmmcppb(x,  est, *args,  alpha=.05, con=None,
            dif=True, nboot=None, BA=False,
            hoch=False, SR=False, seed=False):

    """
    Use a percentile bootstrap method to  compare dependent groups.
    By default,
    compute a .95 confidence interval for all linear contrasts
    specified by con, a J-by-C matrix, where  C is the number of
    contrasts to be tested, and the columns of con are the
    contrast coefficients.
    If con is not specified, all pairwise comparisons are done.

    If est=onestep or mom (may not be implemeted yet),
    method SR (see my book on robust methods)
    is used to control the probability of at least one Type I error.
    Otherwise, Hochberg is used.

    dif=True indicates that difference scores are to be used
    dif=False indicates that measure of location associated with
    marginal distributions are used instead.

    nboot is the bootstrap sample size. If not specified, a value will
    be chosen depending on the number of contrasts there are.

    A sequentially rejective method is used to control alpha using method SR.

    Argument BA: When using dif=False, BA=True uses a correction term
    when computing a p-value.

    :param x:
    :param y:
    :param alpha:
    :param con:
    :param est:
    :param dif:
    :param nboot:
    :param BA:
    :param hoch:
    :param SR:
    :param seed:
    :return:
    """

    if hoch:
        SR=False

    if SR:
        raise Exception("onestep and mom estimators are not yet implemented"
                        "and only these can be used with SR method. Please set SR to False for now.")

    if dif:
        print("analysis is being done on difference scores",
              "each confidence interval has probability coverage of 1-alpha.")

        temp=rmmcppbd(x,est, *args, alpha=alpha,con=con,
                      nboot=nboot,hoch=True)

        return {'output': temp['output'],
                'con':  temp['con'], "num_sig": temp['num_sig']}

    else:
        print("dif=False so using marginal distributions")

        if not BA:
            print("If and when MOM and/or onestep estimators are implemeted, "
                  "it is suggested to use BA=True and hoch=T")

        J=x.shape[1]
        xcen=np.full([x.shape[0], x.shape[1]], np.nan)
        for j in range(J):
            xcen[:, j] = x[:, j] - est(x[:, j], *args)

        if con is None:
            con=con1way(J)

        d=con.shape[1]

        if nboot is None:
            if d<4:
                nboot=1000
            elif d>4:
                nboot=5000

        n=x.shape[0]
        connum=con.shape[1]

        if seed:
            np.random.seed(seed)

        xbars=est(x,*args)

        psidat=np.zeros(connum)
        for ic in range(connum):
            psidat[ic]=np.sum(con[:,ic] * xbars)

        psihat=np.zeros([connum, nboot])
        psihatcen=np.zeros([connum, nboot])
        bvec=np.full([nboot,J], np.nan)
        bveccen = np.full([nboot, J], np.nan)
        data=np.random.randint(n,size=(nboot,n))
        for ib in range(nboot):
            bvec[ib,:] = est(x[data[ib,:],:], *args)
            bveccen[ib, :] = est(xcen[data[ib, :], :], *args)

        test=np.full(connum, np.nan)
        bias=np.full(connum, np.nan)

        for ic in range(connum):
            psihat[ic,:]=[bptdpsi(row, con[:, ic]) for row in bvec]
            psihatcen[ic,:] = [bptdpsi(row, con[:,ic]) for row in bveccen]
            bias[ic] = np.sum((psihatcen[ic,:] > 0)) / nboot - .5
            ptemp =(np.sum(psihat[ic,:] > 0) + .5 * np.sum(psihat[ic,:] == 0)) / nboot

            if BA:
                test[ic] = ptemp - .1 * bias[ic]

            if not BA:
                test[ic] = ptemp

            test[ic] = np.min([test[ic], 1 - test[ic]])
            test[ic] = np.max([test[ic], 0])  # bias corrected might be less than zero

        test=2*test
        ncon=con.shape[1]
        dvec=alpha/np.arange(1,ncon+1)

        if SR:

            if alpha == .05:

                dvec =[.025,
                .025,
                .0169,
                .0127,
                .0102,
                .00851,
                .0073,
                .00639,
                .00568,
                .00511]

                dvecba = [.05,
                .025,
                .0169,
                .0127,
                .0102,
                .00851,
                .0073,
                .00639,
                .00568,
                .00511]

                if ncon > 10:
                    avec = .05 / np.arange(11,ncon+1)
                    dvec = np.append(dvec, avec)

            elif alpha == .01:

                dvec =[.005,
                .005,
                .00334,
                .00251,
                .00201,
                .00167,
                .00143,
                .00126,
                .00112,
                .00101]

                dvecba =[.01,
                .005,
                .00334,
                .00251,
                .00201,
                .00167,
                .00143,
                .00126,
                .00112,
                .00101]

                if ncon > 10:
                    avec = .01 / np.arange(11,ncon+1)
                    dvec = np.append(dvec, avec)


            else:

                dvec = alpha / np.arange(1,ncon+1)
                dvecba = dvec
                dvec[1] = alpha

        if hoch:
            dvec=alpha/np.arange(1,ncon+1)

        dvecba=dvec
        temp2 = (-test).argsort()
        zvec = dvec[:ncon]

        if BA:
            zvec = dvecba[:ncon]

        output=np.zeros([connum, 6])
        tmeans=est(x, *args)

        output[temp2, 3] = zvec
        for ic in range(ncon):
            output[ic, 1] = np.sum(con[:, ic] * tmeans)
            output[ic, 0] = ic
            output[ic, 2] = test[ic]
            temp = np.sort(psihat[ic, :])
            icl = round(alpha * nboot / 2) #+ 1
            icu = nboot - icl - 1 #nboot - (icl - 1)
            output[ic, 4] = temp[icl]
            output[ic, 5] = temp[icu]

    num_sig = output.shape[0]
    ior = (-output[:, 2]).argsort()
    for j in range(output.shape[0]):
        if output[ior[j], 2] <= output[ior[j], 3]:
            break
        else:
            num_sig = num_sig - 1

    results={"output": output, "con": con, "num_sig": num_sig}

    return results

def spmcpa(J, K, x, est, *args,
           avg=False, alpha=.05, nboot=None, seed=False):

    """
    All pairwise comparisons among levels of Factor A
    in a mixed design using trimmed means.

    The variable x is a Pandas DataFrame where the first column
    contains the data for the first level of both factors: level 1,1.
    The second column contains the data for level 1 of the
    first factor and level 2 of the second: level 1,2.
    x.iloc[:,K] is the data for level 1,K. x.iloc[:,K+1] is the data for level 2,1.
    x.iloc[:, 2K] is level 2,K, etc.
    :return:
    """

    x=pandas_to_arrays(x)
    x=remove_nans_based_on_design(x, design_values=[J,K], design_type='between_within')

    if seed:
        np.random.seed(seed)

    nvec=[len(nj) for nj in x[:J*K:K]]

    if avg:
        con=con1way(J)

    elif not avg:

        MJK = K * (J ** 2 - J) // 2
        con=np.zeros([J*K,MJK])
        n_idioms=J-1
        idioms=[]
        K_mult=K

        for i in range(n_idioms):
            tmp=np.concatenate([[1], np.repeat(0, K_mult-1), [-1]])
            idioms.append(tmp)
            K_mult*=2

        col_ind=0
        for idiom in idioms:
            num_rep_idiom=len(list(mit.windowed(con[:,0], n=len(idiom))))

            row_start=0
            for _ in range(num_rep_idiom):
                con[row_start:row_start+len(idiom),col_ind]=idiom
                row_start+=1
                col_ind+=1

    d=con.shape[1]

    if not nboot:
        if d<=4:
            nboot=1000
        else:
            nboot=5000

    xx=x.copy()
    bloc=np.full([J, nboot], np.nan)
    mvec=np.array([])

    ik=0
    for j in range(J):
        x=np.full([nvec[j], K], np.nan)

        for k in range(K):
            x[:, k] = xx[ik]

            if not avg:
                mvec=np.append(mvec, est(xx[ik], *args))

            ik += 1

        tempv = est(x, *args)

        data=np.random.randint(nvec[j], size=(nboot, nvec[j]))
        bvec=np.full([nboot, K], np.nan)

        for k in range(K):
            temp = x[:, k]
            bvec[:, k] = [rmanogsub(data_row, temp, est, *args) for data_row in data]

        if avg:
            mvec=np.append(mvec, np.mean(tempv))
            bloc[j,:] = np.mean(bvec, axis=1)

        elif not avg:

            if j==0:
                bloc=bvec.copy()

            elif j>0:
                bloc=np.c_[bloc, bvec]

    if avg:
        bloc=bloc.T

    connum=d
    psihat=np.zeros([connum, nboot])
    test=np.full(connum, np.nan)

    for ic in range(connum):

      psihat[ic, :] = [bptdpsi(row, con[:, ic]) for row in bloc]
      test[ic] = (np.sum(psihat[ic, :] > 0) + .5 * np.sum(psihat[ic, :] == 0)) / nboot
      test[ic] = np.min([test[ic], 1 - test[ic]])

    ncon=con.shape[1]

    if alpha == .05:

        dvec = [.025,
                .025,
                .0169,
                .0127,
                .0102,
                .00851,
                .0073,
                .00639,
                .00568,
                .00511]

        if ncon > 10:
            avec = .05 / np.arange(11, ncon + 1)
            dvec = np.append(dvec, avec)

    elif alpha == .01:

        dvec = [.005,
                .005,
                .00334,
                .00251,
                .00201,
                .00167,
                .00143,
                .00126,
                .00112,
                .00101]

        if ncon > 10:
            avec = .01 / np.arange(11, ncon + 1)
            dvec = np.append(dvec, avec)


    else:

        dvec = alpha / np.arange(1, ncon + 1)
        dvec[0] = alpha/2

    temp2=(-test).argsort()
    zvec=dvec[:ncon]
    output=np.zeros([connum,6])

    tmeans=mvec
    output[temp2,3]=zvec
    for ic in range(ncon):
        output[ic, 1] = np.sum(con[:, ic] * tmeans)
        output[ic, 0] = ic
        output[ic, 2] = test[ic]
        temp = np.sort(psihat[ic, :])
        icl = int(round(dvec[ncon-1] * nboot)) #+ 1
        icu = nboot - icl - 1 #(icl - 1)
        output[ic, 4] = temp[icl]
        output[ic, 5] = temp[icu]

    output[:, 2] = 2 * output[:, 2]
    output[:, 3] = 2 * output[:, 3]
    num_sig = np.sum(output[:, 2] <= output[:, 3])

    col_names=["con_num", "psihat", "p_value",
               "p_crit", "ci_lower", "ci_upper"]
    output=pd.DataFrame(output, columns=col_names)

    results={"output": output, "con": con, "num_sig": num_sig}

    return results

def spmcpb(J, K, x, est, *args, dif=True, alpha=.05, nboot=599, seed=False):

    """
    All pairwise comparisons among levels of Factor B
    in a split-plot design.

    The variable x is a Pandas DataFrame where the first column
    contains the data for the first level of both factors: level 1,1.
    The second column contains the data for level 1 of the
    first factor and level 2 of the second: level 1,2.
    x.iloc[:,K] is the data for level 1,K. x.iloc[:,K+1] is the data for level 2,1.
    x.iloc[:, 2K] is level 2,K, etc.

    If dif=True, the analysis is done based on all pairs
    of difference scores.
    Otherwise, marginal measures of location are used.

    :param J:
    :param K:
    :param x:
    :param est:
    :param dif:
    :param alpha:
    :param nboot:
    :param seed:
    :return:
    """

    x=pandas_to_arrays(x)
    x=remove_nans_based_on_design(x, design_values=[J,K], design_type='between_within')

    #JK=J*K

    if seed:
        np.random.seed(seed)

    #nvec=np.array([])
    x_mat=np.empty(shape=(0,K))
    jj=0
    kk=K
    for j in range(J):
        #nvec=np.append(nvec, len(x[j]))
        x_mat=np.concatenate([x_mat, np.r_[x[jj:kk]].T], axis=0)
        jj+=K
        kk+=K

    temp=rmmcppb(x_mat, est, *args, nboot=nboot, dif=dif, alpha=alpha)

    col_names=['con_num', 'psihat', 'p_value', 'p_crit', 'ci_lower', 'ci_upper']
    test_res=pd.DataFrame(temp['output'], columns=col_names)
    results={"output": test_res, 'con': temp['con'],
             'num_sig': temp['num_sig']}


    return results

def spmcpi(J, K, x, est, *args, alpha=.05, nboot=None, SR=False, seed=False):

    """
    Multiple comparisons for interactions
    in a split-plot design.
    The analysis is done by taking difference scores
    among all pairs of dependent groups and
    determining which of
    these differences differ across levels of Factor A.

    :return:
    """

    x=pandas_to_arrays(x)
    x=remove_nans_based_on_design(x, design_values=[J,K], design_type='between_within')

    if seed:
        np.random.seed(seed)

    nvec=[len(nj) for nj in x[:J*K:K]]
    MK = (K ** 2 - K) // 2
    MJK = K * (J ** 2 - J) // 2
    con = np.zeros([J * K, MJK])
    n_idioms = J - 1
    idioms = []
    K_mult = K

    for i in range(n_idioms):
        tmp = np.concatenate([[1], np.repeat(0, K_mult - 1), [-1]])
        idioms.append(tmp)
        K_mult *= 2

    col_ind = 0
    for idiom in idioms:
        num_rep_idiom = len(list(mit.windowed(con[:, 0], n=len(idiom))))

        row_start = 0
        for _ in range(num_rep_idiom):
            con[row_start:row_start + len(idiom), col_ind] = idiom
            row_start += 1
            col_ind += 1

    d=con.shape[1]

    if not nboot:
        if d<=4:
            nboot=1000
        else:
            nboot=5000


    xx=x.copy()
    bloc=np.full([J, nboot], np.nan)
    mvec=np.array([])
    it=0
    for j in range(1,J+1):
        x=np.full([nvec[j-1], MK], np.nan)
        im=0
        for k in range(K):
            for kk in range(K):
                if k<kk:
                    kp = j * K + k - K
                    kpp = j * K + kk - K
                    x[:, im] = xx[kp] - xx[kpp]
                    mvec=np.append(mvec, est(x[:,im], *args))
                    it+=1
                    im+=1

        data=np.random.randint(nvec[j-1], size=(nboot, nvec[j-1]))
        bvec=np.full([nboot, MK], np.nan)

        for k in range(MK):
            temp = x[:, k]
            bvec[:, k] = [rmanogsub(data_row, temp, est, *args) for data_row in data]

        if j==1:
            bloc=bvec.copy()

        elif j>1:
            bloc=np.c_[bloc, bvec]

    connum=d
    psihat=np.zeros([connum, nboot])
    test=np.full(connum, np.nan)
    for ic in range(connum):

      psihat[ic, :] = [bptdpsi(row, con[:, ic]) for row in bloc]
      test[ic] = (np.sum(psihat[ic, :] > 0) + .5 * np.sum(psihat[ic, :] == 0)) / nboot
      test[ic] = np.min([test[ic], 1 - test[ic]])

    ncon=con.shape[1]
    dvec = alpha / np.arange(1,ncon+1)

    if SR:
        raise Exception("onestep and mom estimators are not yet implemented"
                        "and only these can be used with SR method. Please set SR to False for now.")

        # THIS CODE IS UNREACHABLE UNTIL SR CONDITION IS ALLOWED
        # if alpha == .05:
        #
        #     dvec = [.025,
        #             .025,
        #             .0169,
        #             .0127,
        #             .0102,
        #             .00851,
        #             .0073,
        #             .00639,
        #             .00568,
        #             .00511]
        #
        #     if ncon > 10:
        #         avec = .05 / np.arange(11, ncon + 1)
        #         dvec = np.append(dvec, avec)
        #
        # elif alpha == .01:
        #
        #     dvec = [.005,
        #             .005,
        #             .00334,
        #             .00251,
        #             .00201,
        #             .00167,
        #             .00143,
        #             .00126,
        #             .00112,
        #             .00101]
        #
        #     if ncon > 10:
        #         avec = .01 / np.arange(11, ncon + 1)
        #         dvec = np.append(dvec, avec)
        #
        # else:
        #
        #     dvec = alpha / np.arange(1, ncon + 1)
        #     dvec[0] = alpha/2

    temp2=(-test).argsort()
    zvec=dvec[:ncon]
    output=np.zeros([connum,6])

    tmeans=mvec
    output[temp2,3]=zvec
    for ic in range(ncon):
        output[ic, 1] = np.sum(con[:, ic] * tmeans)
        output[ic, 0] = ic
        output[ic, 2] = test[ic]
        temp = np.sort(psihat[ic, :])
        icl = int(round(dvec[ncon-1] * nboot)) #+ 1
        icu = nboot - icl - 1 #(icl - 1)
        output[ic, 4] = temp[icl]
        output[ic, 5] = temp[icu]

    output[:, 2] = 2 * output[:, 2]

    if SR:
        output[:, 3] = 2 * output[:, 3]

    num_sig = np.sum(output[:, 2] <= output[:, 3])

    col_names=["con_num", "psihat", "p_value",
               "p_crit", "ci_lower", "ci_upper"]
    output=pd.DataFrame(output, columns=col_names)

    results={"output": output, "con": con, "num_sig": num_sig}

    return results

def wwmcppb(J, K, x,  est, *args,  alpha=.05, dif=True,
            nboot=None, BA=True, hoch=True, seed=False):

    """
    Do all multiple comparisons for a within-by-within design using a percentile bootstrap method

    The variable x is a Pandas DataFrame where the first column
    contains the data for the first level of both factors: level 1,1.
    The second column contains the data for level 1 of the
    first factor and level 2 of the second: level 1,2.
    x.iloc[:,K] is the data for level 1,K. x.iloc[:,K+1] is the data for level 2,1.
    x.iloc[:, 2K] is level 2,K, etc.

    If dif=True, the analysis is done based on all pairs
    of difference scores.
    Otherwise, marginal measures of location are used.

    :param x:
    :param est:
    :param args:
    :param alpha:
    :param dif:
    :param nboot:
    :param BA:
    :param hoch:
    :param seed:
    :return:
    """

    x=pandas_to_arrays(x)
    x=remove_nans_based_on_design(x, design_values=[J,K], design_type='within_within')
    x_mat=np.r_[x].T

    conA, conB, conAB = con2way(J, K)
    A = rmmcppb(x_mat, est, *args, con = conA, alpha = alpha, dif = dif,
                nboot = nboot,BA = BA, hoch = hoch, seed = seed)

    B = rmmcppb(x_mat, est, *args, con = conB, alpha = alpha, dif = dif,
                nboot = nboot,BA = BA, hoch = hoch, seed = seed)

    AB = rmmcppb(x_mat, est, *args, con = conAB, alpha = alpha, dif = dif,
                nboot = nboot,BA = BA, hoch = hoch, seed = seed)

    col_names=['con_num', 'psihat', 'p_value', 'p_crit', 'ci_lower', 'ci_upper']

    [X.update({"output": pd.DataFrame(X['output'], columns=col_names)})
        for X in [A,B,AB]]

    results={'factor_A': A, 'factor_B': B, "factor_AB": AB}

    return results






































