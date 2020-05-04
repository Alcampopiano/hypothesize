__all__ = ["yuenbt", "pb2gen", "linconb", "rmmcppb",
           "lindepbt", "bootdpci", "ydbt", "tmcppb", "l2drmci"]

import numpy as np
import pandas as pd
from scipy.stats import trim_mean
from hypothesize.utilities import yuend, trimse, lincon, trimparts, trimpartt, pandas_to_arrays, \
    con1way, con2way, bptdpsi, rmmcp, trimcibt, remove_nans_based_on_design

def yuenbt(x, y, tr=.2, alpha=.05, nboot=599, seed=False):

    """
    Compute a 1-alpha confidence interval for the difference between
    the trimmed means corresponding to two independent groups.
    The bootstrap-t method is used. During the bootstrapping, the absolute value of the test
    statistic is used (the "two-sided method").

    :param x: group one data; Pandas Series
    :param y: group two data; Pandas Series
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

def linconb(x, con, tr=.2, alpha=.05, nboot=599, seed=False):

    """
    Compute a 1-alpha confidence interval for a set of d linear contrasts
    involving trimmed means using the bootstrap-t bootstrap method.
    Independent groups are assumed.

    CIs are adjusted to control FWE (p values are not adjusted)

    x is a Pandas DataFrame where each column represents a group of the data.

    Missing values are automatically removed.

    con is a J by d matrix containing the contrast coefficents of interest.
    If unspecified, all pairwise comparisons are performed.
    For example, con[:,0]=[1,1,-1,-1,0,0] and con[:,1]=[1,-1,0,0,1,-1]
    will test two contrasts: (1) the sum of the first two trimmed means is
    equal to the sum of the second two, and (2) the difference between
    the first two is equal to the difference between the trimmed means of
    groups 5 and 6.

    The default number of bootstrap samples is nboot=599

    :param x: Pandas DataFrame
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
        #print(testit['test'])
        #print(testit['psihat'])
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

def rmmcppb(x,  est, *args,  alpha=.05, con=None,
            dif=True, nboot=None, BA=False,
            hoch=False, SR=False, seed=False):

    """
    Use a percentile bootstrap method to compare dependent groups.
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

    called_directly=False
    if type(x) is pd.core.frame.DataFrame:
        called_directly=True
        x=x.dropna().values

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

        if called_directly:

            col_names = ['con_num', 'psihat', 'p_value', 'p_crit', 'ci_lower', 'ci_upper']

            return {'output': pd.DataFrame(temp['output'], columns=col_names),
                    'con':  temp['con'], "num_sig": temp['num_sig']}

        else:

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

    if called_directly:
        col_names=['con_num', 'psihat', 'p_value', 'p_crit', 'ci_lower', 'ci_upper']
        results={"output": pd.DataFrame(output, columns=col_names), "con": con, "num_sig": num_sig}
        print(results)

    else:
        results={"output": output, "con": con, "num_sig": num_sig}


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

def lindepbt(x, tr=.2, con=None, alpha=.05, nboot=599, dif=True, seed=False):

    """
    MCP on trimmed means with FWE controlled with Rom's method
    Using a bootstrap-t method.

    dif=T, difference scores are used. And for linear contrasts a simple
    extension is used.

    dif=F, hypotheses are tested based on the marginal trimmed means.


    :param x:
    :param tr:
    :param con:
    :param alpha:
    :param nboot:
    :param dif:
    :param seed:
    :return:
    """

    called_directly=False
    if type(x) is pd.DataFrame:
        x = pandas_to_arrays(x)
        x = remove_nans_based_on_design(x, design_values=len(x), design_type='dependent_groups')
        x = np.r_[x].T
        called_directly=True

    from hypothesize.measuring_associations import wincor

    if seed:
        np.random.seed(seed)

    if con is None:
        con=con2way(1,x.shape[1])[1] # all pairwise
        ncon = con.shape[1]

    else:
        ncon = con.shape[1]

    x = x[~np.isnan(x).any(axis=1)]
    n=x.shape[0]
    J=x.shape[1]
    nval=x.shape[0]
    h1 = nval - 2 * np.floor(tr * nval)
    #df=h1-1
    xbar=trim_mean(x, tr)

    if alpha == .05:

        dvec = [.05,
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

        dvec = [.01,
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


    psihat=np.zeros([ncon,4])
    test = np.zeros([ncon, 5])
    temp1=np.array([])

    for d in range(ncon):
        psihat[d, 0] = d

        if not dif:
            psihat[d, 1] = np.sum(con[:, d] * xbar)
            sejk = 0

            for j in range(J):
                for k in range(J):
                    djk = (nval - 1) * wincor(x[:, j], x[:, k], tr)['wcov'] / (h1 * (h1 - 1))
                    sejk = sejk + con[j, d] * con[k, d] * djk

            sejk = np.sqrt(sejk)
            test[d, 0] = d
            test[d, 1] = np.sum(con[:, d] * xbar) / sejk
            test[d, 4] = sejk

            data=np.random.randint(n, size=(nboot, n))
            xcen = np.full([x.shape[0], x.shape[1]], np.nan)
            for j in range(J):
                xcen[:, j] = x[:, j] - trim_mean(x[:, j], tr)

            bvec=[lindep_sub(data_row, xcen, con[:,d], tr=tr)
                  for data_row in data]

            bsort = np.sort(np.abs(bvec))
            ic = round((1 - alpha) * nboot) - 1 # correct for python with the "- 1"?
            psihat[d, 2] = psihat[d, 1] - bsort[ic] * test[d, 4]
            psihat[d, 3] = psihat[d, 1] + bsort[ic] * test[d, 4]
            p_value = np.mean(np.abs(test[d, 1]) <= np.abs(bvec))
            temp1 = np.append(temp1, p_value)

        elif dif:

            for j in range(J):
                if j==0:
                    dval=con[j,d] * x[:,j]

                elif j>0:
                    dval=dval+con[j,d] * x[:,j]

            temp = trimcibt(dval,tr=tr,alpha=alpha,nboot=nboot,seed=seed)
            temp1 = np.append(temp1, temp['p_value'])
            test[d, 0] = d
            test[d, 1]=temp['test_stat'] ## missing in R?
            test[d, 4] = trimse(dval, tr=tr)
            psihat[d, 1] = trim_mean(dval, tr)
            psihat[d, 2] = temp['ci'][0]
            psihat[d, 3] = temp['ci'][1]

    test[:, 2] = temp1
    temp2 =  (-temp1).argsort()
    zvec = dvec[:ncon]
    test[temp2, 3] = zvec

    # if flagcon
    num_sig = np.sum(test[:, 2] <= test[:, 3])

    if called_directly:

        test=pd.DataFrame(test, columns=["con_num", "test", "p_value", "p_crit", "se"])
        psihat=pd.DataFrame(psihat, columns=["con_num", "psihat", "ci_lower", "ci_upper"])


    return {'test': test, 'psihat': psihat, 'con': con, 'num_sig': num_sig}

def lindep_sub(data, x, con = None, tr = .2):

    con = con.reshape(len(con), 1) # make 2D col vector
    res = rmmcp(x[data,:], con=con, tr=tr, dif=False)['test'][:, 1]

    return res[0]

def pb2gen(x, y, est, *args, alpha=.05, nboot=2000, seed=False):

    """
    Compute a bootstrap confidence interval for the
    the difference between any two parameters corresponding to two
    independent groups.

    :param x: Series
    :param y: Series
    :param est:
    :param args:
    :param alpha:
    :param nboot:
    :param seed:
    :return:
    """

    x, y = pandas_to_arrays([x, y])

    x=x[~np.isnan(x)]
    y=y[~np.isnan(y)]

    if seed:
        np.random.seed(seed)


    datax = np.random.choice(x, size=(nboot, len(x)))
    datay = np.random.choice(y, size=(nboot, len(y)))

    bvecx=est(datax, *args, axis=1)
    bvecy = est(datay, *args, axis=1)

    bvec = np.sort(bvecx - bvecy)
    low = round((alpha / 2) * nboot) #+ 1
    up = nboot - low - 2
    temp = np.sum(bvec < 0) / nboot + np.sum(bvec == 0) / (2 * nboot)
    sig_level = 2 * (min(temp, 1 - temp))
    se = np.var(bvec)

    results={'est_1': est(x,*args),
             'est_2': est(y,*args),
             'est_dif': est(x, *args) - est(y, *args),
             'ci': [bvec[low], bvec[up]],
             'p_value': sig_level,
             'variance': se,
             'n1': len(x),
             'n2': len(y)}

    return results

def bootdpci(x, est, *args, nboot=None, alpha=.05,
             dif=True, BA=False, SR=False):

    """
    Use percentile bootstrap method,
    compute a .95 confidence interval for the difference between
    a measure of location or scale
    when comparing two dependent groups.

    :param x: Series
    :param est:
    :param args:
    :param nboot:
    :param alpha:
    :param dif:
    :param BA:
    :param SR:
    :return:
    """

    # replace with actual estimators when implemented
    if SR and est not in ('onestep', 'mom'):
        SR=False
        print("setting SR to False. SR=True should apparently "
              "only be used with onestep or mom")

    ## in R
    # okay=False
    # if est in (onestep, mom):
    #     okay=True
    #
    # if not okay:
    #     SR=False

    results=rmmcppb(x, est, *args, nboot=nboot,alpha=alpha,
                   SR=SR, dif=dif, BA=BA)

    col_names = ['con_num', 'psihat', 'p_value', 'p_crit', 'ci_lower', 'ci_upper']
    results.update({'output': pd.DataFrame(results['output'], columns=col_names)})

    return results

def ydbt(x, y, tr=.2, alpha=.05, nboot=599, side=True, seed=False):

    """
      Using the bootstrap-t method,
      compute a .95 confidence interval for the difference between
      the marginal trimmed means of paired data.
      By default, 20% trimming is used with nboot=599 bootstrap samples.

      side=False returns equal-tailed ci (no p value)
      side=True returns symmetric ci and a p value

    :param x: Series
    :param y: Series
    :param tr:
    :param alpha:
    :param nboot:
    :param side:
    :param seed:
    :return:
    """

    x = pandas_to_arrays([x, y])
    x=remove_nans_based_on_design(x, 2, 'dependent_groups')
    x,y=[x[0], x[1]]

    if seed:
        np.random.seed(seed)

    data = np.random.randint(len(x), size=(nboot, len(x)))

    xcen = x - trim_mean(x, tr)
    ycen = y - trim_mean(y, tr)

    bvec=[tsub(row, xcen, ycen, tr) for row in data]

    dotest = yuend(x, y, tr=tr)

    estse = dotest['se']
    p_value = np.nan
    dif = trim_mean(x, tr) - trim_mean(y, tr)
    ci=[]

    if not side:
        print('p_value is only returned when side=True')
        ilow = round((alpha / 2) * nboot) -1
        ihi = nboot - ilow - 2
        bsort = np.sort(bvec)
        ci.append(dif - bsort[ihi] * estse)
        ci.append(dif - bsort[ilow + 1] * estse)

    else:
        bsort = np.sort(np.abs(bvec))
        ic = round((1 - alpha) * nboot)-1
        ci.append(dif - bsort[ic] * estse)
        ci.append(dif + bsort[ic] * estse)
        p_value = (np.sum(np.abs(dotest['teststat']) <= np.abs(bvec))) / nboot


    return {'ci': ci, 'dif': dif, 'p_value': p_value}

def tsub(isub, x, y, tr):

    """
    Compute test statistic for trimmed means
    when comparing dependent groups.
    By default, 20% trimmed means are used.
    isub is an array of length n of random integers
    to control bootstrap sampling.

    This function is used by ydbt

    :param isub:
    :param x:
    :param y:
    :param tr:
    :return:
    """

    tsub_res = yuend(x[isub], y[isub], tr = tr)['teststat']

    return tsub_res

def tmcppb(x, est, *args, con=None, bhop=False, alpha=.05, nboot=None, seed=False):

    """
    Multiple comparisons for  J independent groups using trimmed means

    A percentile bootstrap method with Rom's method is used.

    est is the measure of location and defaults
    to the trimmed mean (currently only trimmed mean is implemented)

    :param x:
    :param est:
    :param args:
    :param con:
    :param bhop:
    :param alpha:
    :param nboot:
    :param seed:
    :return:
    """

    x=pandas_to_arrays(x)
    x=remove_nans_based_on_design(x, len(x), 'independent_groups')
    J=len(x)

    mvec = [est(i, *args) for i in x]

    if con is None:
        con=con1way(J)

    ncon=con.shape[1]

    if not nboot:
      nboot = 5000
      if J <= 8:
        nboot = 4000
      elif J <= 3:
        nboot = 2000

    if not bhop:

        if alpha == .05:
            dvec=[.05,
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
                dvec = [dvec, avec]

        elif alpha == .01:
            dvec =[.01,
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
                dvec = [dvec, avec]

        else: #not (alpha != .05 or alpha != .01):
            dvec = alpha / np.arange(1,ncon+1)

    else:
        dvec = (ncon - np.arange(1,ncon+1) + 1) * alpha / ncon

    if seed:
        np.random.seed(seed)

    bvec=np.full([J,nboot], np.nan)
    for i, j in enumerate(x):
        data = np.random.choice(j, size=(nboot, len(j)))
        bvec[i,:]=[est(row, *args) for row in data]

    bcon=con.T @ bvec
    tvec=con.T @ mvec
    test=np.full(ncon, np.nan)
    for d in range(ncon):
        tv = np.sum(bcon[d,:] == 0) / nboot
        test[d] = np.sum(bcon[d, :] > 0) / nboot + .5 * tv
        if test[d] > .5:
            test[d] = 1 - test[d]

    output=np.full([ncon,6], np.nan)
    test=2*test
    temp2=(-test).argsort()
    zvec = dvec[:ncon]
    output[temp2, 3] = zvec
    icl = int(np.round(dvec[-1] * nboot / 2) + 1) - 1
    icu = nboot - icl - 3

    for ic in range(ncon):
        output[ic, 1] = tvec[ic]
        output[ic, 0] = ic
        output[ic, 2] = test[ic]
        temp = np.sort(bcon[ic, :])
        output[ic, 4] = temp[icl]
        output[ic, 5] = temp[icu]


    num_sig = np.sum(output[:, 2] <= output[:, 3])
    cols=["con_num","psihat", "p_value", "p_crit", "ci_lower", "ci_upper"]
    output=pd.DataFrame(output, columns=cols)

    results={'output': output, 'con': con, 'num_sig': num_sig}

    return results

def l2drmci(x,y, est, *args, pairwise_drop_na=True, alpha=.05, nboot=2000, seed=False):

    """
      Compute a bootstrap confidence interval for a
      measure of location associated with
      the distribution of x-y. That is, compare x and y by looking at all possible difference scores
      in random samples of x and y.

      est indicates which measure of location
      will be used (currently only trimmed mean is implemented)

      x and y are possibly dependent


    :param x:
    :param y:
    :param est:
    :param args:
    :param pairwise_drop_na: if True,
        treat data as dependent and remove any row with missing data. If False,
        remove missing data for each group seperately (cannot deal with unequal sample sizes)

    :param alpha:
    :param nboot:
    :param seed:
    :return:
    """

    x, y = pandas_to_arrays([x, y])

    if pairwise_drop_na:
        m1 = np.c_[x, y]  # cbind
        x = m1[~np.isnan(m1).any(axis=1)]

    else:
        x = x[~np.isnan(x)]
        y = y[~np.isnan(y)]

        if len(x) != len(y):
             raise Exception("With unequal sample sizes, you might consider wmwpb "
                    "(currently not implemented)")

        else:
            x = np.c_[x, y]  # cbind

    if seed:
        np.random.seed(seed)

    data = np.random.choice(x.shape[0], size=(nboot, len(x)))

    bvec=np.full(nboot, np.nan)
    for i in range(nboot):
        bvec[i] = \
             loc2dif(x[data[i,:], 0], x[data[i,:], 1], est, *args,
                     drop_na=pairwise_drop_na)

    bvec=np.sort(bvec)
    low = int(np.round((alpha / 2) * nboot) + 1) -1
    up = nboot - low -2
    temp = np.sum(bvec < 0) / nboot + np.sum(bvec == 0) / (2 * nboot)
    sig_level = 2 * (np.min([temp, 1 - temp]))
    ci=[bvec[low], bvec[up]]

    results=dict(zip(['ci', 'p_value'], [ci, sig_level]))

    return results

def loc2dif(x,y, est, *args, drop_na=True):

    """
    Compute a measure of location associated with the
    distribution of x-y, the typical difference between two randomly sampled values.
    The measure of location is indicated by the argument
    est.

    x and y are paired data or independent variables having the same length.
    If x and y have different lengths, use the function wmwloc (not currently implemented)

    Advantage of this estimator: relatively high efficiency even under normality versus
    using sample means.

    :param x:
    :param y:
    :param est:
    :param args:
    :param drop_na:
    :return:
    """

    if drop_na:
        m1 = np.c_[x, y]  # cbind
        m1 = m1[~np.isnan(m1).any(axis=1)]
        x, y = [m1[:,0], m1[:,1]]

    else:
        x=x[~np.isnan(x)]
        y=y[~np.isnan(y)]

    temp=np.subtract.outer(x,y).reshape(len(x)*len(y))
    val=est(temp, *args)

    return val











