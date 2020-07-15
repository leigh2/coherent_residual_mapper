# High Frequency Atmospheric Distortion (HFAD) Fitter
Solve for residuals to a photometric calibration. For example those caused by
high frequency atmospheric distortion.

## Authors:
Leigh C. Smith,
Sergey E. Koposov

```
usage: hfad_fitter.py [-h] [-s {MAD,STD,Qn}] [-f] [-v] [--overwrite]
                      ref_file_path out_file_path

Solve for spatially coherent residuals global calibration.

positional arguments:
  ref_file_path         Path to hdf5 archive of reference sources
  out_file_path         Path to hfad solution output

optional arguments:
  -h, --help            show this help message and exit
  -s {MAD,STD,Qn}, --spread {MAD,STD,Qn}
                        The spread measure to use for error estimation
                        (default MAD)
  -f, --figures         Produce diagnostic figures
  -v, --verbose         Verbose output
  --overwrite           Overwrite output file if exists
```

## Example usage

```
python hfad_fitter.py v20120320_00503_st_refs.hdf5 test.hdf5 -fv --overwrite
```

v20120320_00503_st_refs.hdf5 is an HDF5 archive containing instantaneous and average magnitudes (m and <m> respectively hereafter) of a number of well measured reference stars for the VIRCAM observation v20120320_00503_st. A map of median m-<m> inside spatial bins of each of the 16 VIRCAM detectors is shown below.
![Original residual map](/figs/original.png)
test  
