from hypothesize.measuring_associations import *
from hypothesize.compare_groups_with_single_factor import *
from hypothesize.compare_groups_with_two_factors import *
from hypothesize.utilities import create_example_data, trim_mean, con1way
import numpy as np
import pandas as pd
from pandas._testing import assert_frame_equal
import pickle
import os

alpha=.05
nboot=100
tr=.2
beta=.2

try:
    os.chdir('hypothesize/tests')
except:
    pass

def run_all_pkl_funcs():

    from hypothesize.tests import build_test_data

    for i in dir(build_test_data):
        item = getattr(build_test_data,i)
        if callable(item) and i.startswith('pkl'):
            item()

def build_truth_list(expected_results):

    truth_list=[]

    if type(expected_results) is list:

        for item in expected_results:
            nested_truth_list=build_truth_list(item)
            truth_list.append(nested_truth_list)

    elif type(expected_results) is dict:

        for k in expected_results:

            if type(expected_results[k]) is dict:
                nested_truth_list=[True] * len(expected_results[k])

                truth_list.append(nested_truth_list)
            else:
                truth_list.append(True)

    return truth_list

def check_dict_items_equality(expected_results, actual_results):

    actual_truth=[]

    if type(expected_results) is list:
        for exp_item, act_item in zip(expected_results, actual_results):
            nested_truth = check_dict_items_equality(exp_item, act_item)
            actual_truth.append(nested_truth)

    elif type(expected_results) is dict:

        for k in expected_results:

            if type(expected_results[k]) is np.ndarray:

                # truth=True if not np.testing.assert_array_equal(expected_results[k], actual_results[k]) \
                #     else False

                truth=True if not np.testing.assert_allclose(expected_results[k], actual_results[k]) \
                    else False

                actual_truth.append(truth)

            elif type(expected_results[k]) is pd.DataFrame:

                # truth=True if not assert_frame_equal(expected_results[k], actual_results[k]) \
                #     else False

                truth=True if not assert_frame_equal(expected_results[k], actual_results[k], check_less_precise=True) \
                    else False

                actual_truth.append(truth)

            elif type(expected_results[k]) is dict:
                nested_truth=check_dict_items_equality(expected_results[k], actual_results[k])
                actual_truth.append(nested_truth)

            else:

                if expected_results[k] is None and actual_results[k] is None: \
                    truth = True
                else:
                    truth=True if not np.testing.assert_almost_equal(expected_results[k], actual_results[k]) \
                        else False

                actual_truth.append(truth)

    return actual_truth

def test_l2drmci():

    np.random.seed(42)
    df = create_example_data(2)
    results = l2drmci(df.cell_1, df.cell_2, trim_mean, tr)
    expected = pickle.load(open("test_data/l2drmci.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    #assert results == expected
    assert actual_truth == expected_truth

def test_linconb():

    np.random.seed(42)
    df = create_example_data(3)
    results = linconb(df, con1way(3))
    expected = pickle.load(open("test_data/linconb.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_pb2gen():

    np.random.seed(42)
    df = create_example_data(2)
    results = pb2gen(df.cell_1, df.cell_2, trim_mean, tr)
    expected = pickle.load(open("test_data/pb2gen.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    #assert results == expected
    assert actual_truth == expected_truth

def test_tmcppb():

    np.random.seed(42)
    df = create_example_data(3)
    results = tmcppb(df, trim_mean, tr)
    expected = pickle.load(open("test_data/tmcppb.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_yuenbt():

    np.random.seed(42)
    df = create_example_data(2)
    results = yuenbt(df.cell_1, df.cell_2)
    expected = pickle.load(open("test_data/yuenbt.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    #assert results == expected
    assert actual_truth == expected_truth

def test_bootdpci():

    np.random.seed(42)
    df = create_example_data(3)
    results = bootdpci(df, trim_mean, tr)
    expected = pickle.load(open("test_data/bootdpci.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_rmmcppb():

    np.random.seed(42)
    df = create_example_data(3)
    results = rmmcppb(df, trim_mean, tr)
    expected = pickle.load(open("test_data/rmmcppb.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_lindepbt():

    np.random.seed(42)
    df = create_example_data(3)
    results = lindepbt(df)
    expected = pickle.load(open("test_data/lindepbt.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_ydbt():

    np.random.seed(42)
    df = create_example_data(2)
    results = ydbt(df.cell_1, df.cell_2)
    expected = pickle.load(open("test_data/ydbt.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    #assert results == expected
    assert actual_truth == expected_truth

def test_wwmcppb():

    np.random.seed(42)
    df = create_example_data(6)
    results = wwmcppb(2, 3, df, trim_mean, tr)
    expected = pickle.load(open("test_data/wwmcppb.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_wwmcpbt():

    np.random.seed(42)
    df = create_example_data(6)
    results = wwmcpbt(2, 3, df, tr)
    expected = pickle.load(open("test_data/wwmcpbt.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_bwamcp():

    np.random.seed(42)
    df = create_example_data(6)
    results = bwamcp(2, 3, df)
    expected = pickle.load(open("test_data/bwamcp.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_bwbmcp():

    np.random.seed(42)
    df = create_example_data(6)
    results = bwbmcp(2, 3, df)
    expected = pickle.load(open("test_data/bwbmcp.pkl", "rb"))

    print(results)
    print(expected)
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_bwmcp():

    np.random.seed(42)
    df = create_example_data(6)
    results = bwmcp(2, 3, df)
    expected = pickle.load(open("test_data/bwmcp.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_bwimcp():

    np.random.seed(42)
    df = create_example_data(6)
    results = bwimcp(2, 3, df)
    expected = pickle.load(open("test_data/bwimcp.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_bwmcppb():

    np.random.seed(42)
    df = create_example_data(6)
    results = bwmcppb(2, 3, df, trim_mean, tr)
    expected = pickle.load(open("test_data/bwmcppb.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_spmcpa():

    np.random.seed(42)
    df = create_example_data(6)
    results = spmcpa(2, 3, df, trim_mean, tr)
    expected = pickle.load(open("test_data/spmcpa.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_spmcpb():

    np.random.seed(42)
    df = create_example_data(6)
    results = spmcpb(2, 3, df, trim_mean, tr)
    expected = pickle.load(open("test_data/spmcpb.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_spmcpi():

    np.random.seed(42)
    df = create_example_data(6)
    results = spmcpi(2, 3, df, trim_mean, tr)
    expected = pickle.load(open("test_data/spmcpi.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_corb():

    np.random.seed(42)
    df = create_example_data(2)
    results = corb(wincor, df.cell_1, df.cell_2, alpha, nboot, tr)
    expected = pickle.load(open("test_data/corb.pkl", "rb"))
    expected_truth = build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    #assert results  == expected
    assert actual_truth == expected_truth

def test_pball():

    np.random.seed(42)
    df = create_example_data(3)
    results = pball(df)
    expected = pickle.load(open("test_data/pball.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_pbcor():

    np.random.seed(42)
    df = create_example_data(2)
    results = pbcor(df.cell_1, df.cell_2)
    expected = pickle.load(open("test_data/pbcor.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    #assert results == expected
    assert actual_truth == expected_truth

def test_winall():

    np.random.seed(42)
    df = create_example_data(3)
    results = winall(df)
    expected = pickle.load(open("test_data/winall.pkl", "rb"))
    expected_truth=build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    assert actual_truth == expected_truth

def test_wincor():

    np.random.seed(42)
    df = create_example_data(2)
    results = wincor(df.cell_1, df.cell_2)
    expected = pickle.load(open("test_data/wincor.pkl", "rb"))
    print(results)
    print(expected)
    expected_truth = build_truth_list(expected)
    actual_truth = check_dict_items_equality(expected, results)

    #assert results  == expected
    assert actual_truth == expected_truth



