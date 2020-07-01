#!/usr/bin/env python

import numpy as np
import scipy.linalg





def get_bins(bin_arr, min_bin_size, max_bin_count):
    """
    Select equal width bins along bin_arr such that there are >min_bin_size rows
    per bin but at most max_bin_count bins.

    Parameters:
    bin_arr : array-like
        Array over which to generate bins. Must have length N > min_bin_size.
    min_bin_size : int
        Target minimum number of rows per bin.
    max_bin_count : int
        Maximum number of bins to consider.

    Returns:
    bins : array-ndarray
        Array of bin boundaries.
    """
    # raise error if there aren't enough rows to make one bin
    if len(bin_arr)<min_bin_size:
        raise ValueError("Not enough rows in bin_arr to make one bin.")
    else:
        # make sure input is in appropriate formats
        bin_arr = np.asarray(bin_arr, dtype=np.float64)
        min_bin_size = int(min_bin_size)
        max_bin_count = int(max_bin_count)

    # bin count to start with
    bin_count = max_bin_count
    # variable in which to store current smallest bin size
    smallest_bin_size = 0.
    # while the bins are too small...
    while smallest_bin_size<min_bin_size:
        # generate new bins
        bins = np.linspace(bin_arr.min()-1E-9, bin_arr.max()+1E-9, bin_count)
        # count rows per bin
        bin_sizes, _ = np.histogram(bin_arr, bins=bins)
        # update current smallest bin size
        smallest_bin_size = np.min(bin_sizes)
        #decrement bin count
        bin_count -= 1
    # return final bin array
    return bins






def Chebfit(x, y,
            val, err,
            ndeg=2,
            domain=[-1, 1],
            fitsubset=None):
    """
    Fit Chebyshev polynomial.
    Authors: Sergey Koposov, Leigh Smith

    Parameters:
    x, y : array_like
        1D arrays of shape (n,) containing floating point coordinates.
    val : array_like
        1D array of shape (n,) containing floating point target values.
    err : array_like
        Weight the fit with the provided 1d array of errors.
    ndeg : int, optional
        Number of degrees to use for the Chebyshev polynomial (default 2).
    domain : (2,) array_like, optional
        The domain to use for the Chebyshev polynomial (default [-1,1]).
    fitsubset : boolean array, optional
        Use only this subset of rows for the LSQ fit.

    Returns:
    coeff : ndarray
        An array containing fitted coefficients.
    val_pred : ndarray
        An array containing predicted values.
    """

    # build design matrix
    M1 = []
    for i in range(ndeg + 1):
        for j in range(i + 1):
            M1.append(
                np.polynomial.Chebyshev(
                    (np.arange(ndeg + 1) == j).astype(int),
                    domain=domain)(x) * np.polynomial.Chebyshev(
                        (np.arange(ndeg + 1) == (i - j)).astype(int),
                        domain=domain)(y))
    M1 = np.array(M1).T

    # fit subset is all if none set
    if fitsubset is None:
        fitsubset = slice(None)

    # divide through by errors
    err = np.asarray(err)
    _val = val/err
    Mv = M1/err[:,None]

    # run the lsq fit, and subsequent prediction
    coeff = scipy.linalg.lstsq(Mv[fitsubset], _val[fitsubset])[0]
    val_pred = np.dot(Mv, coeff)*err
    return coeff, val_pred







def Qn(arr, c_n=2.219144):
    """
    Rousseeuw and Croux Qn scale estimator, Gaussian efficiency 82%.

    Parameters:
    arr : array-like
        Array of size N over which to compute Qn.
    c_n : float, optional
        Multiplier to use. Optimal value is dependant on N. Where N is large the
        optimal value is 2.219144 (the default).

    Returns:
    Q_n : float
        Rousseeuw and Croux Qn scale estimator for arr.
    """

    if len(arr)==0:
        return None
    elif len(arr)==1:
        return 0
    else:
        if len(arr)<10 and c_n==2.219144:
            warnings.warn("Array is small, default Cn might not be appropriate.")
        arr = np.asarray(arr)

    # create array mask where upper triangle are ones, zeros elsewhere
    # i.e. 1 where i<j else 0
    msk = np.tri(arr.size, k=-1, dtype=np.bool).T
    _i, _j = map(lambda a: a[msk].flatten(), np.mgrid[0:arr.size, 0:arr.size])

    # differences
    diff = np.abs(arr[_i] - arr[_j])

    # Qn is c_n * Q_25(diff)
    # c_n is 2.219144 for normally distributed data of moderate length and above
    Q_n = c_n * np.quantile(diff, 0.25, interpolation='nearest')

    return Q_n
