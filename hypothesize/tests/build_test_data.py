from hypothesize.measuring_associations import *
from hypothesize.compare_groups_with_single_factor import *
from hypothesize.compare_groups_with_two_factors import *
from hypothesize.utilities import create_example_data, trim_mean, con1way, con2way
import numpy as np
import pickle

alpha=.05
nboot=100
tr=.2
beta=.2

def pkl_l2drmci():

    np.random.seed(42)
    df = create_example_data(2)
    results = l2drmci(df.cell_1, df.cell_2, trim_mean, tr)
    pickle.dump(results, open("hypothesize/tests/test_data/l2drmci.pkl", "wb"))

def pkl_linconb():

    np.random.seed(42)
    df = create_example_data(3)
    results = linconb(df, con1way(3))
    pickle.dump(results, open("hypothesize/tests/test_data/linconb.pkl", "wb"))

def pkl_pb2gen():

    np.random.seed(42)
    df = create_example_data(2)
    results = pb2gen(df.cell_1, df.cell_2, trim_mean, tr)
    pickle.dump(results, open("hypothesize/tests/test_data/pb2gen.pkl", "wb"))

def pkl_tmcppb():

    np.random.seed(42)
    df = create_example_data(3)
    results = tmcppb(df, trim_mean, tr)
    pickle.dump(results, open("hypothesize/tests/test_data/tmcppb.pkl", "wb"))

def pkl_yuenbt():

    np.random.seed(42)
    df = create_example_data(2)
    results = yuenbt(df.cell_1, df.cell_2)
    pickle.dump(results, open("hypothesize/tests/test_data/yuenbt.pkl", "wb"))

def pkl_bootdpci():

    np.random.seed(42)
    df = create_example_data(3)
    results = bootdpci(df, trim_mean, tr)
    pickle.dump(results, open("hypothesize/tests/test_data/bootdpci.pkl", "wb"))

def pkl_rmmcppb():

    np.random.seed(42)
    df = create_example_data(3)
    results = rmmcppb(df, trim_mean, tr)
    pickle.dump(results, open("hypothesize/tests/test_data/rmmcppb.pkl", "wb"))

def pkl_lindepbt():

    np.random.seed(42)
    df = create_example_data(3)
    results = lindepbt(df)
    pickle.dump(results, open("hypothesize/tests/test_data/lindepbt.pkl", "wb"))

def pkl_ydbt():

    np.random.seed(42)
    df = create_example_data(2)
    results = ydbt(df.cell_1, df.cell_2)
    pickle.dump(results, open("hypothesize/tests/test_data/ydbt.pkl", "wb"))

def pkl_wwmcppb():

    np.random.seed(42)
    df = create_example_data(6)
    results = wwmcppb(2, 3, df, trim_mean, tr)
    pickle.dump(results, open("hypothesize/tests/test_data/wwmcppb.pkl", "wb"))

def pkl_wwmcpbt():

    np.random.seed(42)
    df = create_example_data(6)
    results = wwmcpbt(2, 3, df, tr)
    pickle.dump(results, open("hypothesize/tests/test_data/wwmcpbt.pkl", "wb"))

def pkl_bwamcp():

    np.random.seed(42)
    df = create_example_data(6)
    results = bwamcp(2, 3, df)
    pickle.dump(results, open("hypothesize/tests/test_data/bwamcp.pkl", "wb"))

def pkl_bwbmcp():

    np.random.seed(42)
    df = create_example_data(6)
    results = bwbmcp(2, 3, df)
    pickle.dump(results, open("hypothesize/tests/test_data/bwbmcp.pkl", "wb"))

def pkl_bwmcp():

    np.random.seed(42)
    df = create_example_data(6)
    results = bwmcp(2, 3, df)
    pickle.dump(results, open("hypothesize/tests/test_data/bwmcp.pkl", "wb"))

def pkl_bwimcp():

    np.random.seed(42)
    df = create_example_data(6)
    results = bwimcp(2, 3, df)
    pickle.dump(results, open("hypothesize/tests/test_data/bwimcp.pkl", "wb"))

def pkl_bwmcppb():

    np.random.seed(42)
    df = create_example_data(6)
    results = bwmcppb(2, 3, df, trim_mean, tr)
    pickle.dump(results, open("hypothesize/tests/test_data/bwmcppb.pkl", "wb"))

def pkl_spmcpa():

    np.random.seed(42)
    df = create_example_data(6)
    results = spmcpa(2, 3, df, trim_mean, tr)
    pickle.dump(results, open("hypothesize/tests/test_data/spmcpa.pkl", "wb"))

def pkl_spmcpb():

    np.random.seed(42)
    df = create_example_data(6)
    results = spmcpb(2, 3, df, trim_mean, tr)
    pickle.dump(results, open("hypothesize/tests/test_data/spmcpb.pkl", "wb"))

def pkl_spmcpi():

    np.random.seed(42)
    df = create_example_data(6)
    results = spmcpi(2, 3, df, trim_mean, tr)
    pickle.dump(results, open("hypothesize/tests/test_data/spmcpi.pkl", "wb"))

def pkl_corb():

    np.random.seed(42)
    df = create_example_data(2)
    results = corb(wincor, df.cell_1, df.cell_2, alpha, nboot, tr)
    pickle.dump(results, open("hypothesize/tests/test_data/corb.pkl", "wb"))

def pkl_pball():

    np.random.seed(42)
    df = create_example_data(3)
    results = pball(df)
    pickle.dump(results, open("hypothesize/tests/test_data/pball.pkl", "wb"))

def pkl_pbcor():

    np.random.seed(42)
    df = create_example_data(2)
    results = pbcor(df.cell_1, df.cell_2)
    pickle.dump(results, open("hypothesize/tests/test_data/pbcor.pkl", "wb"))

def pkl_winall():

    np.random.seed(42)
    df = create_example_data(3)
    results = winall(df)
    pickle.dump(results, open("hypothesize/tests/test_data/winall.pkl", "wb"))

def pkl_wincor():

    np.random.seed(42)
    df = create_example_data(2)
    results = wincor(df.cell_1, df.cell_2)
    pickle.dump(results, open("hypothesize/tests/test_data/wincor.pkl", "wb"))
