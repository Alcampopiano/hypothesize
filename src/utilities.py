import numpy as np

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
                        print(k, kk)
                        if k < kk:
                            mat=np.zeros([J,K])
                            mat[j, k] = 1
                            mat[j, kk] = 0-1
                            mat[jj, k] = 0-1
                            mat[jj, kk] = 1
                            ic = ic + 1

                        conAB[:, ic-1] = mat.reshape(conAB.shape[0]).T

    return conA, conB, conAB

