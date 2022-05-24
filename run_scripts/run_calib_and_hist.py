from core_fct.fct_calib import exec_calib, post_calib


##################################################
##################################################

## OPTIONS
name = 'v1'
method = 'FullRankADVI'
folder_out = 'internal_data/pyMC_calib/'


##################################################
##################################################

## execute Bayesian calibration (usually takes several hours)
exec_calib(folder_out=folder_out, method=method, name=name)

## execute post-processing
post_calib(folder_calib=folder_out + method + '__' + name, make_prior=True, write_csv=True)

