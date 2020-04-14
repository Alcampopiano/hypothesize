__all__ = ["bwmcp", "bwamcp", "bwbmcp", "bwimcp",
           "spmcpa", "spmcpb", "spmcpi", "wwmcppb",
           "wwmcpbt", "bwmcppb"]

import numpy as np
import pandas as pd
from scipy.stats import trim_mean
from hypothesize.utilities import con2way, lindep, covmtrim, lincon, yuen, con1way, bptdpsi, rmanogsub, \
    rmmcp, linhat, remove_nans_based_on_design, pandas_to_arrays
from hypothesize.compare_groups_with_single_factor import rmmcppb, lindepbt
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

    if con is None:
        col_names_test = ["group_x", "group_y", "test", "p_value", "p_crit", "se"]
        col_names_psihat = ["group_x", "group_y", "psihat", "ci_lower", "ci_upper"]

    else:
        col_names_test = ["con_num", "test", "p_value", "p_crit", "se"]
        col_names_psihat = ["con_num", "psihat", "ci_lower", "ci_upper"]

    if pool:
        data = [np.concatenate(x[i:i+J*K+1:K]) for i in range(K)]
        data=np.vstack(data).T
        tests=rmmcp(data, con=con, tr=tr,alpha=alpha,dif=dif,hoch=hoch)

        test_df = pd.DataFrame(tests['test'], columns=col_names_test)
        psihat_df = pd.DataFrame(tests['psihat'], columns=col_names_psihat)

        results = {"test": test_df, 'psihat': psihat_df, 'n': tests['n'], 'con': tests['con'],
                 'num_sig': tests['num_sig']}


    else:

        results=[]
        j_ind=0
        for j in range(J):
            data=x[j_ind:j_ind+K]
            data = np.vstack(data).T
            tests = rmmcp(data, con=con, tr=tr, alpha=alpha, dif=dif, hoch=hoch)

            test_df=pd.DataFrame(tests['test'], columns=col_names_test)
            psihat_df=pd.DataFrame(tests['psihat'], columns=col_names_psihat)
            tests={"test": test_df, 'psihat': psihat_df, 'n': tests['n'], 'con': tests['con'], 'num_sig': tests['num_sig']}

            results.append(tests)
            j_ind+=K


    return results

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

    :param J:
    :param K:
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

def wwmcpbt(J, K, x, tr=.2, alpha=.05, nboot=599, seed=False):

    """
    Do multiple comparisons for a within-by-within design.
    using a bootstrap-t method and trimmed means.
    All linear contrasts relevant to main effects and interactions
    are tested.

    :param J:
    :param K:
    :param x:
    :param tr:
    :param alpha:
    :param nboot:
    :param seed:
    :return:
    """

    print("ask wilcox if dif is supposed to be a argument here")

    x = pandas_to_arrays(x)
    x = remove_nans_based_on_design(x, design_values=[J, K], design_type='within_within')
    x_mat = np.r_[x].T

    conA, conB, conAB = con2way(J, K)

    A = lindepbt(x_mat, tr=tr, con=conA, alpha=alpha, nboot=nboot, seed=seed)
    B = lindepbt(x_mat, tr=tr, con=conB, alpha=alpha, nboot=nboot, seed=seed)
    AB = lindepbt(x_mat, tr=tr, con=conAB, alpha=alpha, nboot=nboot, seed=seed)

    col_names_test=['con_num', 'test', 'p_value', 'p_crit', 'se']
    col_names_psihat=['con_num', 'psihat', 'ci_lower', 'ci_upper']

    [X.update({"test": pd.DataFrame(X['test'], columns=col_names_test)})
        for X in [A,B,AB]]

    [X.update({"psihat": pd.DataFrame(X['psihat'], columns=col_names_psihat)})
        for X in [A,B,AB]]

    results={'factor_A': A, 'factor_B': B, "factor_AB": AB}

    return results

def bwmcppb(J, K, x, est, *args, alpha=.05,
            nboot=500, bhop=True, seed=True):

    """
    A percentile bootstrap for multiple comparisons
    for all main effects and interactions
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
    :param est: estimator (e.g., trim_mean)
    :param args: positional arguments for estimator (e.g., .2 for trim_mean)
    :param alpha: significance level
    :param nboot: number of bootstrap samples
    :param bhop:
    :param seed: random seed for reproducibility
    :return: results dictionary
    """

    print("ask wilcox if this is onlly for trimmed means "
          "depite the optional est arg",
          "I did, he said this is taken care of in 5th edition")

    x = pandas_to_arrays(x)
    x = remove_nans_based_on_design(x, design_values=[J, K], design_type='between_within')

    conA, conB, conAB = con2way(J, K)

    A = bwmcppb_sub(J, K, x, est, *args, con=conA,
        alpha = alpha, nboot = nboot, bhop = bhop, seed = seed)

    B = bwmcppb_sub(J, K, x, est, *args, con=conB,
                    alpha=alpha, nboot=nboot, bhop=bhop, seed=seed)

    AB = bwmcppb_sub(J, K, x, est, *args, con=conAB,
                    alpha=alpha, nboot=nboot, bhop=bhop, seed=seed)

    # col_names_test=['con_num', 'test', 'p_value', 'p_crit', 'se']
    # col_names_psihat=['con_num', 'psihat', 'ci_lower', 'ci_upper']
    #
    # [X.update({"test": pd.DataFrame(X['test'], columns=col_names_test)})
    #     for X in [A,B,AB]]
    #
    # [X.update({"psihat": pd.DataFrame(X['psihat'], columns=col_names_psihat)})
    #     for X in [A,B,AB]]

    results={'factor_A': pd.DataFrame(A), 'factor_B': pd.DataFrame(B), "factor_AB": pd.DataFrame(AB)}

    return results

def bwmcppb_sub(J, K, x, est, *args, con=None, alpha=.05,
                nboot=500, bhop=True, seed=True):

    """

    :param J:
    :param K:
    :param x:
    :param est:
    :param args:
    :param con:
    :param alpha:
    :param nboot:
    :param bhop:
    :param seed:
    :return:
    """

    nvec=max([len(j) for j in x])
    ncon=con.shape[1]
    #xx=x.copy()

    if seed:
        np.random.seed(seed)

    bsam = []
    aboot=np.full([nboot, ncon], np.nan)
    tvec = linhat(x, con, est, *args)

    ilow = 0
    iup = K
    for ib in range(nboot):

        for j in range(J):

            nv = len(x[ilow])
            bdat=np.random.randint(nv, size=nv)

            for k in range(ilow,iup):
                bsam.append(x[k][bdat])

            ilow = ilow + K
            iup = iup + K

        ilow = 0
        iup = K

        aboot[ib, :] = linhat(bsam, con, est, *args)
        bsam=[]

    pbA=np.full(aboot.shape[1], np.nan)

    for j in range(aboot.shape[1]):
        pbA[j] = np.mean(aboot[:, j] > 0)
        pbA[j] = 2 * np.min([pbA[j], 1 - pbA[j]])

    if not bhop:
        dvec=alpha/np.arange(1, ncon+1)

    else:
        dvec=(ncon - np.arange(1, ncon+1) + 1) \
             * alpha / ncon

    outputA=np.zeros([ncon, 6])

    test = pbA
    temp2 = (-test).argsort()
    zvec = dvec[:ncon]
    outputA[temp2, 3] = zvec
    icl = int(np.round(dvec[-1] * nboot / 2)) # + 1
    icu = nboot - icl -3 # - 1
    outputA[:, 1] = tvec

    for ic in range(ncon):
        outputA[ic, 0] = ic
        outputA[ic, 2] = test[ic]
        temp = np.sort(aboot[:, ic])
        outputA[ic, 4] = temp[icl]
        outputA[ic, 5] = temp[icu]


    return {"con_num": outputA[:, 0],
            "psihat": outputA[:, 1],
            "p_value": outputA[:,2],
            "p_crit": outputA[:,3],
            "ci_lower": outputA[:,4],
            "ci_upper": outputA[:,5]
            }


































