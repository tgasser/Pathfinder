from core_fct.fct_calib import exec_calib, post_calib


##################################################
##################################################

## OPTIONS
name = 'v1'
method = 'FullRankADVI'


##################################################
##################################################

## execute Bayesian calibration (usually takes several hours)
exec_calib(method=method, name=name, redo_prior_calib=True)

## execute post-processing
post_calib(folder_calib='internal_data/pyMC_calib/' + method + '__' + name, save_prior=True, save_Var2=False, write_csv=True, redo_prior_calib=False)

