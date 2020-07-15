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

### The problem

v20120320_00503_st_refs.hdf5 is an HDF5 archive containing instantaneous and average magnitudes of a number of well measured reference stars for the VIRCAM observation v20120320_00503_st. A map of median residual (i.e. instantaneous minus average magnitude) inside spatial bins of each of the 16 VIRCAM detectors is shown below.

![Original residual map](/figs/original.png)

Clearly there is a lot of spatially coherent structure present in this image. We can reduce the scatter in the observations of individual stars if we can remove this structure.

### The solution

hfad_fitter.py maps the structure of the spatially coherent residuals for each detector by fitting Chebyshev polynomials of increasing degrees until the fractional improvement of sqrt(chisq) drops below 1%. Residuals are measured using 5-fold cross validation for robustness.

Verbose output for ideal Chebyshev polynomial degree finding
```
Finding best ndeg
     |  ndeg  Δ√<(residual/error)²>
chip |     3   |     5   |     7   |    10   |    13   |    17   |    21   |    25  
-----|---------|---------|---------|---------|---------|---------|---------|---------
   1 |     inf |  0.0045 |  0.0006 |⇒ 0.0040⇐|         |         |         |         
   2 |     inf |  0.0054 |  0.0086 |⇒ 0.0066⇐|         |         |         |         
   3 |     inf |  0.0174 |  0.0130 |⇒ 0.0046⇐|         |         |         |         
   4 |     inf |  0.0217 |  0.0078 |⇒ 0.0050⇐|         |         |         |         
   5 |     inf |  0.0190 |  0.0084 |  0.0107 |⇒ 0.0078⇐|         |         |         
   6 |     inf |  0.0129 |  0.0288 |  0.0114 |⇒ 0.0095⇐|         |         |         
   7 |     inf |  0.0223 |  0.0241 |⇒ 0.0096⇐|         |         |         |         
   8 |     inf |  0.0090 |  0.0055 |⇒ 0.0058⇐|         |         |         |         
   9 |     inf |  0.0251 |  0.0128 |  0.0184 |⇒ 0.0054⇐|         |         |         
  10 |     inf |  0.0076 |  0.0032 |  0.0105 |⇒ 0.0013⇐|         |         |         
  11 |     inf |  0.0236 |  0.0363 |  0.0285 |  0.0151 |⇒ 0.0034⇐|         |         
  12 |     inf |  0.0169 |  0.0075 |⇒ 0.0067⇐|         |         |         |         
  13 |     inf |  0.0176 |  0.0043 |⇒ 0.0087⇐|         |         |         |         
  14 |     inf |  0.0251 |  0.0065 |⇒ 0.0057⇐|         |         |         |         
  15 |     inf |  0.0141 |  0.0104 |  0.0146 |⇒ 0.0080⇐|         |         |         
  16 |     inf |  0.0240 |  0.0127 |⇒ 0.0069⇐|         |         |         |         
```

And in figure form:

![Err vs. ndeg](/figs/error_vs_ndeg.png)

Once the ideal number of Chebyshev polynomial degrees is identified (i.e. high enough to capture as much structure as possible but without overfitting the data or being too expensive to compute) it is refit to all the data. The fitted model is shown in the below figure.

![Model](/figs/model.png)

Once this model is subtracted from the data we obtain the new residual maps shown below -- a significant improvement.

![New residual map](/figs/residuals.png)

An additional step is taken to calibrate the raw magnitude errors by fitting for a raw magnitude error scaling factor and a calibration error added in quadrature such that the distribution of the new residuals over their calibrated errors is approximately unit Gaussian in a series of magnitude bins. In the figure below the raw magnitude errors are the blue points, the calibrated magnitude errors are in orange. Generally this calibration results in an increase to the final magnitude error since the raw errors tend to be underestimated.

![Original vs. calibrated errors](/figs/old_new_errs.png)

Verbose output for the error calibration

```
chip | ndeg | bins | cal_err | err_scale
-----|------|------|---------|----------
   1 |   10 |   11 | 0.02633 | 1.127961
   2 |   10 |   10 | 0.02878 | 1.206450
   3 |   10 |   10 | 0.02610 | 1.119989
   4 |   10 |    7 | 0.01951 | 1.019898
   5 |   13 |   12 | 0.02511 | 1.171598
   6 |   13 |    9 | 0.02425 | 1.060683
   7 |   10 |    8 | 0.02288 | 1.055600
   8 |   10 |   13 | 0.01942 | 1.038915
   9 |   13 |    9 | 0.02203 | 1.060033
  10 |   13 |   14 | 0.02306 | 1.131490
  11 |   17 |   12 | 0.02401 | 1.084157
  12 |   10 |    9 | 0.02221 | 1.044167
  13 |   10 |   10 | 0.02238 | 1.072370
  14 |   10 |   11 | 0.02245 | 1.091288
  15 |   13 |   11 | 0.01824 | 1.028085
  16 |   10 |    9 | 0.01638 | 1.015260
```
