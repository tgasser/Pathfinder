"""
Copyright IIASA (International Institute for Applied Systems Analysis), 2019-2022
Contributor(s): Thomas Gasser (gasser@iiasa.ac.at)

This software is a computer program whose purpose is simulate a large number of global climate (and environmental) change scenarios.

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
"""


import os
import numpy as np
import xarray as xr
import scipy.stats as st

from scipy.integrate import quad
from scipy.optimize import fsolve


##################################################
##   CALCULATE LATEST TREND
##################################################

def get_trend(var, axis, trend_over):
    
    ## take specified period
    var_tmp = var.isel({axis:slice(-trend_over-1, None)}).values.reshape(trend_over+1, -1)
    
    ## calculate as narray
    d_var = np.empty_like(var_tmp[-1])
    for n in range(len(var_tmp[-1])):
        d_var[n] = st.theilslopes(var_tmp[:,n])[0] 
    d_var = d_var.reshape(var_tmp[-1].shape)

    ## remap to xarray
    d_var = xr.zeros_like(var[-1].drop(axis)) + d_var 

    ## deduce current value from mean and trend
    var = var.isel({axis:slice(-trend_over-1, None)}).mean(axis) + 0.5*trend_over*d_var

    return var, d_var


##################################################
##   PROBABILITY DISTRIBUTIONS
##################################################

## draw from a (truncated) standard normal distribution
f_draw = lambda size, std_max=np.inf, **kwargs: xr.DataArray(st.truncnorm.rvs(-abs(std_max),+abs(std_max), size=size, **kwargs), coords={'config':np.arange(size)}, dims=['config'])


## deduce distribution parameters from known mean and std, for logit-normal distribution
def logitnorm_distrib_param(mean, std, raise_error=False):
    ## error function
    def err(par):
        exp, _ = quad(lambda x, mu, sigma: 1/(1.-x) * 1./np.sqrt(2*np.pi*sigma**2.) * np.exp(-0.5*(np.log(x/(1.-x))-mu)**2./sigma**2.), 0, 1, args=tuple(par), limit=100)
        var, _ = quad(lambda x, mu, sigma: x/(1.-x) * 1./np.sqrt(2*np.pi*sigma**2.) * np.exp(-0.5*(np.log(x/(1.-x))-mu)**2./sigma**2.), 0, 1, args=tuple(par), limit=100)
        return np.array([exp-mean, np.sqrt(var-exp**2)-std])**2
    ## minimize error function
    par, _, fsolve_flag, _ = fsolve(err, [np.log(mean/(1.-mean)), np.sqrt(std/mean)], full_output=True)
    ## unpack, check and format
    mu, sigma = par
    if fsolve_flag==1: 
        sigma = np.abs(sigma)
    else: 
        if raise_error: raise RuntimeError('solver did not converge')
        else: mu, sigma = np.nan, np.nan
    return mu, sigma


## deduce distribution parameters from known mean and std, for log-normal distribution
def lognorm_distrib_param(mean, std):
    mu = np.log(mean / np.sqrt(1. + std**2./mean**2.))
    sigma = np.sqrt(np.log(1. + std**2./mean**2.))
    return mu, sigma


## combine distribution parameters and drawing from normal standard distrib
def apply_distrib(param, draw, distrib):
    if distrib == 'norm':
        return param[0] + param[1]*draw
    elif distrib == 'lognorm':
        par0, par1 = lognorm_distrib_param(param[0], param[1])
        return np.exp(par0 + par1*draw)
    elif distrib == 'logitnorm':
        par0, par1 = logitnorm_distrib_param(param[0], param[1])
        return 1./(1 + np.exp(-(par0 + par1*draw)))
    else:
        raise ValueError('wrong distribution')


##################################################
##   ROLLING MEAN & STD
##################################################

## rolling mean function
def rolling_mean(da, where=True, time_axis='year', time_len=5):
    return da.where(where).rolling({time_axis: time_len}, center=True).mean().dropna(time_axis).values


## rolling std function
def rolling_std(da, where=True, time_axis='year', time_len=5, min_cv=1E-6):
    min_std = min_cv * rolling_mean(da, where=where, time_axis=time_axis, time_len=time_len)
    std = da.where(where).rolling({time_axis: time_len}, center=True).std().dropna(time_axis).values    
    return np.maximum(std, min_std)

