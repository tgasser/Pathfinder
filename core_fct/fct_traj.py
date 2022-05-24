"""
Copyright IIASA (International Institute for Applied Systems Analysis), 2019-2021
Contributor(s): Thomas Gasser (gasser@iiasa.ac.at)

This software is a computer program whose purpose is simulate a large number of global climate (and environmental) change scenarios.

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
"""


import warnings
import numpy as np
import xarray as xr

from random import choices
from scipy.optimize import fsolve
from scipy.stats import theilslopes
from collections import OrderedDict

from core_fct.fct_ancillary import get_trend


## /!\ WARNING: /!\
## old code; need some cleaning and updating! but it does what it's supposed to.


##################################################
##   1. TEMPERATURE PATHWAYS
##################################################

## generate temperature change pathways
def gen_T_path(var_in, tF, T_obj={'Tl':1.5, 'Tp':2.0}, n_traj=0, series='', param_traj={}, trend_over=20):
    print('generating: T pathways')

    ## get time axis
    tH = var_in['T'].year[-1].values
    tt = xr.DataArray(np.arange(tH, tF+1, dtype=float)-tH, coords={'year':np.arange(tH, tF+1, dtype=int)}, dims=['year'])

    ## input temperature data
    ## unpack objectives
    Tl, Tp = [T_obj[var] for var in ['Tl', 'Tp']]
    ## get historical values
    T_hist = var_in['T']

    ## initial values and trends
    T_0, dT_0 = get_trend(T_hist, axis='year', trend_over=trend_over)

    ## priority to given parameters
    ## otherwise draw series randomly
    if len(param_traj) == 0:
        if series in ['1a', '1b', '1c', '2a', '2b', '2c', '3a', '3b']:
            series = [series for _ in range(n_traj)]
        else: series = choices(['1a', '1b', '1c', '2a', '2b', '2c', '3a', '3b'], k=n_traj)
    else: n_traj = 1; series = [series]

    ## number of trajectories in each series
    n_series = OrderedDict([(s, series.count(s)) for s in ['1a', '1b', '1c', '2a', '2b', '2c', '3a', '3b']])
    n_cumsum = [0] + list(np.cumsum(list(n_series.values())))

    ## dummy traj variables to format xarray   
    traj = {s: xr.DataArray(np.ones(n_series[s]), 
        coords={'traj': np.arange(n_cumsum[n], n_cumsum[n+1])}, 
        dims=['traj']) 
        for n,s in enumerate(['1a', '1b', '1c', '2a', '2b', '2c', '3a', '3b'])}

    ## SERIES 1
    ## parameter-sparse global temperature time profiles
    ## (Huntingford et al., 2017; doi:10.5194/esd-8-617-2017)
    traj1 = xr.concat([traj[s] for s in ['1a', '1b', '1c']], dim='traj')

    ## read values when provided
    if len(param_traj) > 0:
        if series[0][:1] == '1':
            Tll = traj1 * Tl
            if set(['mu0', 'mu1']).issubset(set(param_traj.keys())):
                mu0, mu1 = [traj1 * param_traj[var] for var in ['mu0', 'mu1']]
            else: raise ValueError('wrong parameters for series 1')
        else: Tll = mu0 = mu1 = np.nan * traj1

    ## draw them randomly otherwise
    else:
        ## asymptote is long-term target (series 1a)
        Tll_a = 0*T_0 + traj['1a'] * Tl
        ## asymptote above long-term target (series 1b)
        Tll_b = 0*T_0 + traj['1b'] * np.random.uniform(Tl, Tp, size=len(traj['1b']))
        ## asymptote below long-term target (series 1c)
        Tll_c = T_0 + (Tl - T_0) * (traj['1c'] * np.random.uniform(0, 1, size=len(traj['1c'])))
        ## common parameters
        Tll = xr.concat([Tll_a, Tll_b, Tll_c], dim='traj')
        del Tll_a, Tll_b, Tll_c
        mu0 = (traj1 * np.random.triangular(-0.05, 0, 0.10, size=len(traj1)))
        mu1 = (traj1 * np.random.triangular(0, 0.05, 0.10, size=len(traj1)))

    ## compute pathways
    T_1 = 0*tt + T_0 + (Tll - T_0) * (1 - np.exp(-mu0 * tt - mu1**2 * tt**2)) + (dT_0 - mu0 * (Tll - T_0)) * tt * np.exp(-mu0 * tt - mu1**2 * tt**2)
    dT_1 = 0*tt + ((Tll - T_0) * (mu0 + 2*mu1**2 * tt) + (dT_0 - mu0 * (Tll - T_0)) * (1 - mu0 * tt - 2*mu1**2 * tt**2)) * np.exp(-mu0 * tt - mu1**2 * tt**2)
    ddT_1 = 0*tt + ((Tll - T_0) * (2*mu1**2 - mu0**2 - 4*mu0 * mu1**2 * tt - 4*mu1**4 * tt**2) - (dT_0 - mu0 * (Tll - T_0)) * (2*mu0 + (6*mu1**2 - mu0**2) * tt - 4*mu0 * mu1**2 * tt**2 - 4*mu1**4 * tt**3)) * np.exp(-mu0 * tt - mu1**2 * tt**2)
    
    ## clean memory
    del Tll, mu0, mu1

    ## SERIES 2
    ## hand-crafted overshooting trajectories with continuity of T and dT at peak
    traj2 = xr.concat([traj[s] for s in ['2a', '2b', '2c']], dim='traj')

    ## read values when provided
    if len(param_traj) > 0:
        if series[0][:1] == '2':
            Tll = traj2 * Tl
            Tpp = traj2 * Tp
            if set(['tp', 'alpha', 'beta']).issubset(set(param_traj.keys())):
                tp, alpha, beta = [traj2 * param_traj[var] for var in ['tp', 'alpha', 'beta']]
            # /!\ missing parameters!
            else: raise ValueError('wrong parameters for series 2')
        else: Tll = Tpp = tp = to = alpha = beta = tau = sigma = eta = np.nan * traj2

    ## draw them randomly otherwise
    else:
        ## common parameters
        Tll = traj2 * Tl
        Tpp = traj2 * np.random.uniform(Tll, Tp)
        tp = traj2 * np.random.uniform(10, 100, size=len(traj2))
        alpha = traj2 * np.random.uniform(1+1E-3, 20, size=len(traj2))
        beta = 2. 
        ## solve for remaining parameter
        tp_ = (tp + 0*T_0).values.reshape(-1)
        alpha_ = (alpha + 0*T_0).values.reshape(-1)
        tmp_ = (alpha * beta * (Tpp-T_0)/dT_0).values.reshape(-1)
        to = np.empty_like(tmp_)
        with warnings.catch_warnings(): # ignore warnings
            warnings.filterwarnings('ignore')
            for n in range(len(to)):
                to[n] = fsolve(lambda t: (-t)*((1-tp_[n]/t)**alpha_[n]-1) - tmp_[n], -tp_[n])
        to = xr.full_like((tp + 0*T_0), 1) * to.reshape((tp + 0*T_0).shape)
        del tp_, alpha_, tmp_
        ## second half follows normal distrib. (series 2a)
        tau = traj['2a'] * np.sqrt(Tpp-Tll) / np.sqrt(Tpp-T_0) * (tp-to) * (1 - (-to/(tp-to))**alpha)**(beta/2.) * np.sqrt(2./(alpha**2 * beta*(beta-1)))
        ## second half follows log-normal distrib. (series 2b)
        sigma = traj['2b'] * np.sqrt(Tpp-Tll) / np.sqrt(Tpp-T_0) * (1-to/tp) * (1 - (-to/(tp-to))**alpha)**(beta/2.) * np.sqrt(1./(alpha**2 * beta*(beta-1)))
        ## second half follows Gompertz distrib. (series 2c)
        eta = traj['2c'] * np.exp(np.sqrt(Tpp-T_0) / np.sqrt(Tpp-Tll) * (1-to/tp)**-1 * (1 - (-to/(tp-to))**alpha)**(-beta/2.) * np.sqrt(alpha**2 * beta*(beta-1)))

    ## compute pathways
    T_2 = (tt <= tp) * (T_0 + (Tpp - T_0) * (1 - (1 - (to / (to - tp))**alpha)**-beta * (1 - ((to - tt) / (to - tp))**alpha)**beta))
    dT_2 = (tt <= tp) * (alpha * (Tpp - T_0) * beta * (1 - (to / (to - tp))**alpha)**-beta * ((to - tt) / (to - tp))**alpha * (1 - ((to - tt) / (to - tp))**alpha)**(beta - 1)) / (tt - to)
    ddT_2 = (tt <= tp) * (-alpha * (Tpp - T_0) * beta * (1 - (to / (to - tp))**alpha)**-beta * ((to - tt) / (to - tp))**alpha * (1 - ((to - tt) / (to - tp))**alpha)**(beta - 2) * (alpha * (beta * ((to - tt) / (to - tp))**alpha - 1) - ((to - tt) / (to - tp))**alpha + 1)) / (tt - to)**2

    ## (series 2a)
    T_2.loc[{'traj':traj['2a'].traj}] += (tt > tp) * (Tll + (Tpp - Tll) * np.exp(-(tt - tp)**2 / tau**2))
    dT_2.loc[{'traj':traj['2a'].traj}] +=  (tt > tp) * (Tpp - Tll) * (-2 / tau**2) * (tt - tp) * np.exp(-(tt - tp)**2 / tau**2)
    ddT_2.loc[{'traj':traj['2a'].traj}] +=  (tt > tp) * (Tpp - Tll) * (-2 / tau**4) * (tau**2 - 2*(tt - tp)**2) * np.exp(-(tt - tp)**2 / tau**2)
    ## (series 2b)
    T_2.loc[{'traj':traj['2b'].traj, 'year':tt.year[1:]}] += (tt[1:] > tp) * (Tll + (Tpp - Tll) * np.exp(-0.5 * np.log(tt[1:] / tp)**2 / sigma**2))
    dT_2.loc[{'traj':traj['2b'].traj, 'year':tt.year[1:]}] += (tt[1:] > tp) * (Tpp - Tll) * -np.log(tt[1:] / tp) / sigma**2 / tt[1:] * np.exp(-0.5 * np.log(tt[1:] / tp)**2 / sigma**2)
    ddT_2.loc[{'traj':traj['2b'].traj, 'year':tt.year[1:]}] += (tt[1:] > tp) * (Tpp - Tll) * (sigma**2 * np.log(tt[1:] / tp) + np.log(tt[1:] / tp)**2 - sigma**2) / sigma**4 / tt[1:]**2 * np.exp(-0.5 * np.log(tt[1:] / tp)**2 / sigma**2)
    ## (series 2c)
    T_2.loc[{'traj':traj['2c'].traj}] += (tt > tp) * (Tll + (Tpp - Tll) * np.exp(1 - eta**(1 - tt / tp)) * eta**(1 - tt / tp))
    dT_2.loc[{'traj':traj['2c'].traj}] += (tt > tp) * (Tpp - Tll) * -np.log(eta) / tp * np.exp(1 - eta**(1 - tt / tp)) * eta**(1 - 2*tt / tp) * (eta**(tt / tp) - eta)
    ddT_2.loc[{'traj':traj['2c'].traj}] += (tt > tp) * (Tpp - Tll) * np.log(eta)**2 / tp**2 * np.exp(1 - eta**(1 - tt / tp)) * eta**(1 - 3*tt / tp) * (eta**(2*tt / tp) - 3*eta**(1 + tt / tp) + eta**2)

    ## clean memory
    del Tll, Tpp, tp, alpha, beta, tau, sigma, eta

    ## SERIES 3
    ## trajectories based on dampened oscillator differential equations
    traj3 = xr.concat([traj[s] for s in ['3a', '3b']], dim='traj')

    ## read values when provided
    if len(param_traj) > 0:
        if series[0][:1] == '3':
            Tll = traj3 * Tl
            if set(['kappa', 'omega']).issubset(set(param_traj.keys())):
                kappa, omega = [traj1 * param_traj[var] for var in ['kappa', 'omega']]
            else: raise ValueError('wrong parameters for series 3')
            assert omega <= kappa
        else: Tll = kappa = omega = np.nan * traj3

    ## draw them randomly otherwise
    else:
        ## common parameters
        Tll = traj3 * Tl
        kappa = (traj3 * np.random.triangular(0, 0.05, 0.15, size=len(traj3)))
        ## critically dampened oscillator (series 3a)
        omega = kappa.copy()
        ## over-dampened oscillator (series 3b)
        omega.loc[{'traj':traj['3b'].traj}] = traj['3b'] * np.random.uniform(0, kappa.loc[{'traj':traj['3b'].traj}])

    ## compute pathways
    ## (series 3a)
    T_3 = 0*tt + Tll + np.exp(-kappa * tt) * ((T_0 - Tll) + (kappa * (T_0 - Tll) + dT_0) * tt)
    dT_3 = 0*tt + np.exp(-kappa * tt) * (dT_0 - kappa * (kappa * (T_0 - Tll) + dT_0) * tt)
    ddT_3 = 0*tt + kappa * np.exp(-kappa * tt) * (kappa * (T_0 - Tll) * (kappa * tt - 1) + dT_0 * (kappa * tt - 2))
    ## (series 3b)
    T_3.loc[{'traj':traj['3b'].traj}] = (0*tt + 0*traj['3b']) + Tll + np.exp(-kappa * tt) * ((T_0 - Tll) * np.cosh(np.sqrt(kappa**2 - omega**2) * tt) + (kappa * (T_0 - Tll) + dT_0) / np.sqrt(kappa**2 - omega**2) * np.sinh(np.sqrt(kappa**2 - omega**2) * tt))
    dT_3.loc[{'traj':traj['3b'].traj}] = (0*tt + 0*traj['3b']) + np.exp(-kappa * tt) * (dT_0 * np.cosh(np.sqrt(kappa**2 - omega**2) * tt) + -(omega**2 + kappa * dT_0) / np.sqrt(kappa**2 - omega**2) * np.sinh(np.sqrt(kappa**2 - omega**2)*tt))
    ddT_3.loc[{'traj':traj['3b'].traj}] = (0*tt + 0*traj['3b']) + np.exp(-kappa * tt) * (-(omega**2 * (T_0 - Tll) + 2*kappa * dT_0) * np.cosh(np.sqrt(kappa**2 - omega**2) * tt) + omega**2 * (kappa * (T_0 - Tll) - dT_0) * np.sinh(np.sqrt(kappa**2 - omega**2)*tt))

    ## clean memory
    del Tll, kappa, omega

    ## FINALIZE
    ## combine series
    T = xr.concat([T_1, T_2, T_3], dim='traj')
    dT = xr.concat([dT_1, dT_2, dT_3], dim='traj')
    ddT = xr.concat([ddT_1, ddT_2, ddT_3], dim='traj')

    ## remove pathways that don't fit criteria
    keep = np.logical_and((T.max('year') <= Tp) & (T.min('year') >= 0), np.isclose(dT.sel(year=tH), dT_0*xr.concat(traj.values(), dim='traj')))
    T = T.where(keep, drop=True)
    dT = dT.where(keep, drop=True)
    ddT = ddT.where(keep, drop=True)

    ## switch to float32
    T = T.astype(np.float32)
    dT = dT.astype(np.float32)
    ddT = ddT.astype(np.float32)

    ## add units
    T.attrs['units'] = 'K'
    dT.attrs['units'] = 'K yr-1'
    ddT.attrs['units'] = 'K yr-2'

    ## return as dataset with added series variable
    series.sort()
    series = xr.DataArray(series, coords={'traj':np.arange(len(series))}, dims=['traj'])
    return xr.Dataset({'T':T, 'd_T':dT, 'dd_T':ddT, 'series':series})


##################################################
##   2. ATMOSPHERIC CO2 PATHWAYS
##################################################

## generate atm. CO2 pathways
def gen_CO2_path(var_in, param_clim, tF, T_obj={'Tl':1.5}, n_traj=0, param_traj={}, trend_over=6):
    print('generating: CO2 pathways')

    ## get time axis
    tH = var_in['CO2'].year[-1].values
    tt = xr.DataArray(np.arange(tH, tF+1, dtype=float)-tH, coords={'year':np.arange(tH, tF+1, dtype=int)}, dims=['year'])

    ## get input data
    ## unpack objectives
    Tl = T_obj['Tl']
    ## get historical values
    CO2_hist = var_in['CO2']
    ## get climate system parameters
    CO2pi, T2x = [param_clim[var] for var in ['CO2pi', 'T2x']]
    
    ## initial values and trends
    CO2_0, dCO2_0 = get_trend(CO2_hist, axis='year', trend_over=trend_over)

    ## TRAJ.
    ## trajectories based on dampened oscillator differential equations
    traj = xr.DataArray(np.ones(n_traj), coords={'traj':np.arange(n_traj)}, dims=['traj'])
    series = []

    ## read values when provided
    if len(param_traj) > 0:
        if set(['kX', 'muX']).issubset(set(param_traj.keys())):
            kX, muX = [param_traj[var] for var in ['kX', 'muX']]
            k_x = ['-1'] * (-1 <= kX < -0.1) + ['0'] * (-0.1 <= kX <= 0.1) + ['+1'] * (0.1 < kX <= 1)
        if set(['betaX']).issubset(set(param_traj.keys())):
            betaX = param_traj['betaX']
            series.append('0')
        if set(['muXX']).issubset(set(param_traj.keys())):
            muXX = param_traj['muXX']
            series.append('1')
        if len(series) == 0: raise ValueError('wrong parameters for CO2 series')

    ## draw them randomly otherwise
    else:
        k_x = np.sort(choices(['-1', '0', '+1'], k=n_traj))
        kX = traj * np.array([np.random.uniform(-1*(K=='-1') + -0.1*(K=='0') + 0.1*(K=='+1'), -0.1*(K=='-1') + 0.1*(K=='0') + 1*(K=='+1')) for K in k_x])
        muXX = traj * np.random.triangular(-0.05, 0, 0.10, size=n_traj)
        muX = traj * np.random.triangular(0.0, 0.05, 0.10, size=n_traj)
        betaX = traj * np.random.uniform(2+1E-3, 10, size=n_traj)

    ## generate future CO2 pathway
    CO2_l = CO2pi * np.exp((1 - kX) * Tl * np.log(2) / T2x)
    ## old series
    #CO2 = 0*tt + CO2_0 + (CO2pi * np.exp((1 - kX) * Tl * np.log(2) / T2x) - CO2_0) * (1 - np.exp(-muX * tt)) + (dCO2_0 - muX * (CO2pi * np.exp((1 - kX) * Tll * np.log(2) / T2x) - CO2_0)) * tt * np.exp(-muX * tt)
    #dCO2 = (0*tt + 0*CO2_0) + (muX**2 * tt * (CO2pi * np.exp((1 - kX) * Tll * np.log(2) / T2x) - CO2_0) - dCO2_0 * muX * tt + dCO2_0) * np.exp(-muX * tt)
    ## series 1 (same as with temperature)
    CO2_1 = 0*tt + CO2_0 + (CO2_l - CO2_0) * (1 - np.exp(-muXX * tt - muX**2 * tt**2)) + (dCO2_0 - muXX * (CO2_l - CO2_0)) * tt * np.exp(-muXX * tt - muX**2 * tt**2)
    dCO2_1 = 0*tt + ((CO2_l - CO2_0) * (muXX + 2*muX**2 * tt) + (dCO2_0 - muXX * (CO2_l - CO2_0)) * (1 - muXX * tt - 2*muX**2 * tt**2)) * np.exp(-muXX * tt - muX**2 * tt**2)
    ## series 0 (log-logistic function)
    CO2_ = 0*tt + CO2_0 + ((CO2_l - CO2_0) * (muX * tt)**betaX + dCO2_0 * tt) / (1 + (muX * tt)**betaX)
    dCO2_ = (0*tt + 0*CO2_0) + (betaX * (CO2_l - CO2_0) * muX**betaX * tt**(betaX - 1) + dCO2_0 * (1 + (muX * tt)**betaX - betaX * (muX * tt)**betaX)) / (1 + (muX * tt)**betaX)**2
    #CO2 = 0*tt + CO2_0 + ((CO2_l - CO2_0) * (muX * tt)**betaX + dCO2_0 * tt) / (1 + (muX * tt)**betaX)
    
    ## draw or choose series
    if len(param_traj) > 0 and series == ['1']:
        CO2, dCO2 = CO2_1, dCO2_1
        series = xr.DataArray(['1'] * n_traj, coords={'traj':np.arange(n_traj)}, dims=['traj'])
    elif len(param_traj) > 0 and series == ['0']:
        CO2, dCO2 = CO2_, dCO2_
        series = xr.DataArray(['0'] * n_traj, coords={'traj':np.arange(n_traj)}, dims=['traj'])
    else:
        series = xr.DataArray(choices(['1', '0'], k=n_traj), coords={'traj':np.arange(n_traj)}, dims=['traj'])
        CO2 = xr.concat([CO2_1.where(series == '1').dropna('traj', how='all'), CO2_.where(series == '0').dropna('traj', how='all')], dim='traj')
        dCO2 = xr.concat([dCO2_1.where(series == '1').dropna('traj', how='all'), dCO2_.where(series == '0').dropna('traj', how='all')], dim='traj')

    ## remove pathways that don't fit criteria
    keep = (CO2.max('year') <= 2000) & (CO2.min('year') >= CO2pi)
    CO2 = CO2.where(keep, drop=True)
    dCO2 = dCO2.where(keep, drop=True)

    ## switch to float32
    CO2 = CO2.astype(np.float32)
    dCO2 = dCO2.astype(np.float32)

    ## add units
    CO2.attrs['units'] = 'ppm'
    dCO2.attrs['units'] = 'ppm yr-1'

    ## return as dataset with added k_x variable
    k_x = xr.DataArray(k_x, coords={'traj':np.arange(len(k_x))}, dims=['traj'])
    return xr.Dataset({'CO2':CO2, 'd_CO2':dCO2, 'sign':k_x, 'seriesX':series})

