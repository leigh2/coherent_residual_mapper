#!/usr/bin/env python3

import numpy as np
import argparse
from os.path import isfile
import h5py
from utils import get_bins, Chebfit, Qn
from scipy.optimize import minimize


# command line setup and reading
def cmdargs():
    parser = argparse.ArgumentParser(
        description="Solve for spatially coherent residuals global calibration.")
    parser.add_argument("ref_file_path", type=str,
        help="Path to hdf5 archive of reference sources")
    parser.add_argument("out_file_path", type=str,
        help="Path to solution output")
    parser.add_argument("-s", "--spread", choices=["MAD", "STD", "Qn"],
        default="MAD",
        help="The spread measure to use for error estimation (default MAD)")
    parser.add_argument("-f", "--figures", action="store_true", default=False,
        help="Produce diagnostic figures")
    parser.add_argument("-v", "--verbose", action="store_true", default=False,
        help="Verbose output")
    parser.add_argument("--overwrite", action="store_true", default=False,
        help="Overwrite output file if exists")
    return parser.parse_args()




# columns to be read
read_cols = ["ext",         # detector number
             "ucal_mag",    # calibrated mag
             "mean_mag",    # mean mag of all observations
             "eimag",       # instrumental mag error
             "std_mag",     # standard deviation on mags from all observations
             "n_phot",      # number of observations
             "x", "y",      # x and y detector coordinates
             ]
# check which files from the possible combinations exist and then read them
def read_data(infile):
    with h5py.File(infile, "r") as f:
        gp = f["refs"]

        # select rows to read
        read = ~np.isnan(gp["mean_mag"][:])         # at least one good detection
        bad = (   gp["imap_edge"][read]             # not near array edge
               | (gp["n_phot"][read]<2)             # at least two good detections in this band
               | (gp["ucal_min_rcnt"][read]<100)    # more than 100 refs for all ubercal coeffs
               | (gp["objtype"][read]!=1)  )        # perfect dophot detection
        read[read] = ~bad

        # read the required columns and rows
        data = [gp[c][read] for c in read_cols]

    # return numpy record array
    return np.rec.fromarrays(data, names=read_cols)





# the main code block
def run(infile, outfile,
        figs=False, verbose=False, overwrite=False, spread='MAD'):


    if figs:
        # some additional imports if figures are to be produced
        from scipy.stats import binned_statistic_2d as bs2d
        import matplotlib.pyplot as plt


    # choose appropriate spread measure for error estimation
    if spread=="MAD":
        # median absolute deviation
        _std = lambda arr: np.median(np.abs(arr-np.median(arr)))*1.4826
    elif spread=="STD":
        # standard deviation
        _std = np.std
    elif spread=="Qn":
        # Rousseeuw and Croux Qn
        _std = Qn
    else:
        # user defined spread measure
        _std = spread

    if verbose:
        print("filename: {:s}".format(infile))

    ###################################
    ### load, clean and format data ###
    ###################################

    # read files
    data = read_data(infile)
    # slices for each detector
    data = data[np.argsort(data["ext"])]
    chips, idx, cnt = np.unique(data["ext"], return_index=True, return_counts=True)
    slices = [None]*chips.size
    for _ch, _id, _ct in zip(chips, idx, cnt):
        slices[_ch-1] = slice(_id, _id+_ct)
    # delta magnitude, Δmag error, mag error
    Δmag = data["ucal_mag"]-data["mean_mag"]
    σ_Δmag = np.hypot(data["eimag"], data["std_mag"]/np.sqrt(data["n_phot"]))


    ##################################################
    ### identify the best number of degrees to use ###
    ##################################################


    if figs:
        fig = plt.figure(figsize=(11,8))

    # ndegs to test for
    ndegs = [3,5,7,10,13,17,21,25]

    # somewhere to store best ndeg and other parameters and arrays
    best_ndegs = [None]*chips.size
    chip_res = [None]*chips.size
    ndegs_tested = [None]*chips.size
    rt_mn_sq_sig_chips = [None]*chips.size

    if verbose:
        print("\nFinding best ndeg")
        _space_count = int(((len(ndegs)*10)-13)/2)
        print("     |  ndeg  Δ√<(residual/error)²>")
        print("chip | "+" | ".join(["   {:2d}  ".format(nd) for nd in ndegs]))
        print("-----"+"|---------"*len(ndegs))



    for chip in chips:
        if verbose:
            _str = []
            print("{:4d}".format(chip), end=" ")
        # slice for this chip
        _s = slices[chip-1]

        # list to store residuals
        res_store = []
        # list to append root mean square residual over error
        rt_mn_sq_sig = []
        # starting delta fraction
        delta_frac = np.inf
        # starting index
        _i = 0
        # run the loop
        while (delta_frac>0.01 or _i<=3) and _i!=len(ndegs):
            ndeg = ndegs[_i]

            # get cross validated residuals
            res = cv_residuals(data["x"][_s], data["y"][_s], Δmag[_s], σ_Δmag[_s],
                               ncv=5, ndeg=ndeg)
            # store residuals
            res_store += [res]
            # store root mean square residual over error
            rt_mn_sq_sig += [np.sqrt(np.mean((res/σ_Δmag[_s])**2))]

            # calculate fractional change
            if _i>0:
                delta_frac = (rt_mn_sq_sig[_i-1]-rt_mn_sq_sig[_i])/rt_mn_sq_sig[_i]
            else:
                delta_frac = np.inf

            if verbose:
                _str += [char for char in "| {:7.4f} ".format(delta_frac)]

            _i += 1

        # find the best ndeg
        rt_mn_sq_sig = np.asarray(rt_mn_sq_sig)
        idx_best = np.argmin(rt_mn_sq_sig)

        if verbose:
            _str[1+(idx_best*10)] = "⇒"
            _str[9+(idx_best*10)] = "⇐"
            print("".join(_str) + "|         "*(len(ndegs)-len(rt_mn_sq_sig)))
        # store the results for this chip
        best_ndegs[chip-1] = ndegs[idx_best]
        chip_res[chip-1] = res_store[idx_best].copy()
        ndegs_tested[chip-1] = ndegs[:rt_mn_sq_sig.size]
        rt_mn_sq_sig_chips[chip-1] = rt_mn_sq_sig.copy()


        if figs:
            _sp_size = int(np.ceil(np.sqrt(chips.size)))
            plt.subplot(_sp_size,_sp_size,chip)
            plt.plot(ndegs[:len(rt_mn_sq_sig)], rt_mn_sq_sig)
            plt.xlim(0,30)
            plt.xlabel("Polynomial degrees")
            plt.ylabel("√<(residual/error)²>")


    if figs:
        plt.tight_layout()


    ###############################
    ### main distortion fitting ###
    ###############################

    if verbose:
        print("\nMain distortion and error fitting")
        print("chip | ndeg | bins | cal_err | err_scale")
        print("-----|------|------|---------|----------")

    if figs:
        fig = plt.figure(figsize=(11,8))

    # somewhere to store solutions and other parameters and arrays
    coeffs = [None]*chips.size
    cal_errs = [None]*chips.size
    err_scales = [None]*chips.size
    ref_counts = [None]*chips.size
    Δmag_pred = Δmag*np.nan

    for chip in chips:
        # slice for this chip
        _s = slices[chip-1]

        # clip 5 sigma residuals and rerun
        fit = np.abs(chip_res[chip-1])/σ_Δmag[_s] < 5
        ref_counts[chip-1] = np.count_nonzero(fit)
        # run the fitter
        coeffs[chip-1], Δmag_pred[_s] = run_Chebfit(
            data["x"][_s], data["y"][_s], Δmag[_s], σ_Δmag[_s],
            ndeg=best_ndegs[chip-1], fitsubset=fit
        )



        ### estimate error scale factor and calibration error for this chip ###
        # select appropriate bins inside which to measure spread
        bins = get_bins(data["ucal_mag"][_s], 100, 100)
        # get bin numbers
        bn = np.digitize(data["ucal_mag"][_s], bins=bins)
        # slice for each bin
        bin_slices = [bn==_b for _b in np.unique(bn)]

        # build the cost function
        resid = (Δmag[_s].astype(np.float64) - Δmag_pred[_s].astype(np.float64))
        def minfunc(params):
            cal_err, err_scale = params
            # the total error is the instrumental error multiplied by the error
            # scale and added in quadrature to the calibration error
            tot_err = np.hypot(err_scale*data["eimag"][_s].astype(np.float64), cal_err)
            # measure the residual over error spread in each evaluation bin
            stat = np.array([_std(resid/tot_err) for bs in bin_slices])
            # the spread in residual over error should ideally be 1 in each bin
            # sum(ln(stat)^2) is therefore an appropriate cost function
            return np.sum(np.log(stat)**2)

        # run the minimizer
        res = minimize(minfunc, [0.01, 1.0], bounds=[(1E-6,1.0), (0.1,10.0)])
        # store the result
        cal_errs[chip-1], err_scales[chip-1] = res.x

        if verbose:
            print(" {:3d} | {:4d} | {:4d} | {:7.5f} | {:8.6f} ".format(chip, best_ndegs[chip-1], bins.size-1, *res.x))

        if figs:
            _sp_size = int(np.ceil(np.sqrt(chips.size)))
            plt.subplot(_sp_size,_sp_size,chip)
            plt.scatter(data["ucal_mag"][_s], data["eimag"][_s], s=3, alpha=0.5)
            plt.scatter(data["ucal_mag"][_s], np.hypot(res.x[1]*data["eimag"][_s], res.x[0]), s=3, alpha=0.5)
            plt.xlabel("mag")
            plt.ylabel("mag error")

    if figs:
        plt.tight_layout()



    if figs:
        # bins
        step = 100. # pixels
        bins = np.arange(0., 2110., step)

        fig = plt.figure(figsize=(11,8))
        lims = [None]*chips.size
        for chip in chips:
            _d = data[slices[chip-1]]
            _Δmag = Δmag[slices[chip-1]]

            st,_,_,_ = bs2d(_d["x"], _d["y"], _Δmag, statistic='median', bins=[bins]*2)

            _sp_size = int(np.ceil(np.sqrt(chips.size)))
            plt.subplot(_sp_size,_sp_size,chip)
            vl,vu = np.nanpercentile(st, [1,99])
            vmax = np.max([-vl, vu])
            lims[chip-1] = (-vmax,vmax)
            plt.imshow(st.T, cmap='jet', vmin=-vmax, vmax=vmax)
            plt.colorbar()
            plt.gca().get_xaxis().set_visible(False)
            plt.gca().get_yaxis().set_visible(False)

        fig.suptitle('Original', x=0.5, y=1.0)
        plt.tight_layout()


        fig = plt.figure(figsize=(11,8))
        for chip in chips:
            _d = data[slices[chip-1]]
            _Δmag = Δmag_pred[slices[chip-1]]

            st,_,_,_ = bs2d(_d["x"], _d["y"], _Δmag, statistic='median', bins=[bins]*2)

            _sp_size = int(np.ceil(np.sqrt(chips.size)))
            plt.subplot(_sp_size,_sp_size,chip)
            vl,vu = lims[chip-1]
            plt.imshow(st.T, cmap='jet', vmin=vl, vmax=vu)
            plt.colorbar()
            plt.gca().get_xaxis().set_visible(False)
            plt.gca().get_yaxis().set_visible(False)

        fig.suptitle('Model', x=0.5, y=1.0)
        plt.tight_layout()


        fig = plt.figure(figsize=(11,8))
        for chip in chips:
            _d = data[slices[chip-1]]
            _Δmag = (Δmag-Δmag_pred)[slices[chip-1]]

            st,_,_,_ = bs2d(_d["x"], _d["y"], _Δmag, statistic='median', bins=[bins]*2)

            _sp_size = int(np.ceil(np.sqrt(chips.size)))
            plt.subplot(_sp_size,_sp_size,chip)
            vl,vu = lims[chip-1]
            plt.imshow(st.T, cmap='jet', vmin=vl, vmax=vu)
            plt.colorbar()
            plt.gca().get_xaxis().set_visible(False)
            plt.gca().get_yaxis().set_visible(False)

        fig.suptitle('Residuals', x=0.5, y=1.0)
        plt.tight_layout()



    ####################
    ### store result ###
    ####################

    if not isfile(outfile) or overwrite:
        with h5py.File(outfile, "w") as f:
            # store the input filename
            f.attrs["infile"] = infile
            # store Chebyshev domain
            f.attrs["domain"] = [0.,2110.]

            # for each chip
            for chip in chips:
                gp = f.create_group("ext_{:d}".format(chip))
                # store the number of Chebyshev polynomial degrees used
                gp.attrs["ndeg"] = best_ndegs[chip-1]
                # store the number of reference sources
                gp.attrs["ref_count"] = ref_counts[chip-1]
                # store the fitted calibration error
                gp.attrs["calibration_error"] = cal_errs[chip-1]
                # store the fitted instrumental error scale factor
                gp.attrs["error_scale_factor"] = err_scales[chip-1]
                # ndegs tested and their mean residual significance are stored
                gp.create_dataset("ndegs_tested", data=ndegs_tested[chip-1])
                gp.create_dataset("root_mean_sq_sig", data=rt_mn_sq_sig_chips[chip-1])
                # store the solution for each chip
                gp.create_dataset("soln", data=coeffs[chip-1])


    if figs:
        plt.show()


    # done
    return True







def cv_residuals(x, y, Δmag, err, ncv=5, ndeg=6):
    # measure cross validated residuals
    # array in which to store residuals
    res = Δmag*np.nan
    # array of indices
    idx = np.arange(Δmag.size)
    # array of random integers to use for cross validation selection
    xids = np.random.randint(0, ncv, size=res.size)
    # get cross validated residuals
    for j in range(ncv):
        # select fit and validation set
        fit = xids != j
        val = ~fit
        # run the fitting routine
        c, Δmag_pred = Chebfit(
            x, y, Δmag, err,
            fitsubset=fit, ndeg=ndeg,
            domain=[0.,2110.]
        )
        # send these residuals to output
        res[idx[val]] = Δmag[val]-Δmag_pred[val]
    return res





def run_Chebfit(x, y, Δmag, err, ndeg=6, fitsubset=None):
    # run the Chebyshev polynomial fitter

    # array for output
    Δmag_pred = Δmag*np.nan

    #run the fit
    c, pred = Chebfit(
        x, y, Δmag, err,
        fitsubset=fitsubset, ndeg=ndeg, domain=[0.,2110.]
    )

    return c, pred










if __name__=="__main__":
    # get cmdline args
    args = cmdargs()


    # run main code block
    run(args.ref_file_path,
        args.out_file_path,
        figs=args.figures,
        verbose=args.verbose,
        overwrite=args.overwrite,
        spread=args.spread)
