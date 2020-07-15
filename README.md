# High Frequency Atmospheric Distortion (HFAD) Fitter
Solve for residuals to a photometric calibration. For example those caused by
high frequency atmospheric distortion.

Authors:
Leigh C. Smith
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
