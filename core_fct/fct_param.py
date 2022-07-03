"""
Copyright IIASA (International Institute for Applied Systems Analysis), 2019-2022
Contributor(s): Thomas Gasser (gasser@iiasa.ac.at)

This software is a computer program whose purpose is simulate a large number of global climate (and environmental) change scenarios.

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
"""


import numpy as np
import pandas as pd
import xarray as xr
import scipy.special as sp
import matplotlib.pyplot as plt

from scipy.optimize import fsolve, fmin

from core_fct.fct_load import load_gmst, load_CMIP6_climate, load_CMIP5_ocean, load_CMIP6_ocean, load_TRENDYv7_land, load_CMIP5_land, load_CMIP6_land, load_Edwards_2021_slr


##################################################
## CALIBRATE PARAMETERS
##################################################

##=============
## ECS distrib.
##=============

## get paramteres for distribution of ECS
## AR6: (Forster et al., 2021; doi:10.1017/9781009157896.009) (Executive summary)
## AR5: (Collins et al., 2013; doi:10.1017/CBO9781107415324.024) (Box 12.2)
## AR4: (Meehl et al., 2007) (Box 10.2)
## AR3: (IPCC, 2001; Technical Summary) (Section F.3)
## AR2: (IPCC, 1995; Technical Summary) (Section D.2)
## AR1: (Mitchell et al., 1990) (Section 5.2.1)
def get_ecs_distrib(data='AR6', display_fit=False):

    ## distribution for T2x (and for the feedback factor ff)
    ## (Roe & Baker, 2007; doi:10.1126/science.1144735)
    ## note: ff is logit-normal as restricted to [0,1], thus T2x is shifted log-normal
    pdf = lambda T, mu, sigma, T0: 1 / sigma / np.sqrt(2 * np.pi) / (T - T0) * np.exp(- 0.5 * (np.log(T / T0 - 1) - mu)**2 / sigma**2)
    cdf = lambda T, mu, sigma, T0: 0.5 * (1 - sp.erf((mu - np.log(T / T0 - 1)) / sigma / np.sqrt(2)))
    qtf = lambda p, mu, sigma, T0: np.exp(mu + np.log(T0) + np.sqrt(2) * sigma * sp.erfinv(2 * p - 1)) + T0

    ## first error function (to be used with fsolve)
    def err(par, points, T0=None, get_par=False):
        pts = dict(points)
        mu, sigma = par[0], abs(par[1])
        if 0.50 in pts.keys(): T0 = pts[0.50] / (1 + np.exp(mu)); del pts[0.50]
        elif len(par)>2: T0 = abs(par[2])
        else: T0 = float(T0)
        if get_par: return mu, sigma, T0
        else: return [qtf(pct, mu, sigma, T0) - T for pct, T in pts.items()]
 
    ## second error function (to be used with fmin)
    def err2(par, points, T0=None):
        return np.sum(np.array(err(par, points, T0))**2)

    ## points in the ECS distribution
    if data == 'AR6':
        points = [(0.50, 3.0), (0.83, 4.0), (0.17, 2.5), (0.95, 5.0), (0.05, 2.0), (0.01, 1.5)]
    elif data == 'AR5':
        points = [(0.83, 4.5), (0.17, 1.5), (0.90, 6.0), (0.05, 1.0)]
    elif data == 'AR4':
        points = [(0.50, 3.0), (0.83, 4.5), (0.17, 2.0), (0.10, 1.5)]
    elif data in ['AR3', 'AR2', 'AR1']:
        points = [(0.50, 2.5), (0.83, 4.5), (0.17, 1.5)]

    ## run fit
    if len(points)<=3:
        par, _, flag, _ = fsolve(err, (0., 0.2), args=(points, 1.2), maxfev=10000, full_output=True)
        mu, sigma, T2x0 = err(par, points, 1.2, get_par=True) if flag==1 else 3*[np.nan]
    else:
        par, _, _, _, flag = fmin(err2, (0., 0.2, 1.2), args=(points,), maxiter=1000, maxfun=1000, full_output=True, disp=False)
        mu, sigma, T2x0 = err(par, points, None, get_par=True) if flag==0 else 3*[np.nan]

    ## plot fit
    if display_fit:
        plt.figure()
        plt.plot([val for pct, val in points if type(pct) != str], [pct for pct, val in points if type(pct) != str], ls='none', marker='+', ms=8, mew=2)
        plt.plot(np.arange(0, 10, 0.01), cdf(np.arange(0, 10, 0.01), mu, sigma, T2x0), 'k-')
        plt.plot(np.arange(0, 10, 0.01), pdf(np.arange(0, 10, 0.01), mu, sigma, T2x0) / np.nanmax(pdf(np.arange(0, 10, 0.01), mu, sigma, T2x0)), ls='--', color='0.5')
        plt.title('ECS distribution | T2x0 = {0:.2f},\n mu = {1:.3f}, sigma = {2:.3f}'.format(T2x0, mu, sigma), fontsize='small')

    ## wrap in dataset with units
    Par = xr.Dataset()
    Par['T2x0'], Par['mu_ff'], Par['sigma_ff'] = T2x0, mu, sigma
    Par['T2x0'].attrs['units'] = 'K'
    Par['mu_ff'].attrs['units'] = Par['sigma_ff'].attrs['units'] = '1'

    ## return
    return Par


##========
## Climate
##========

## get climate response parameters
## CMIP5_a (Geoffroy et al., 2013a; doi:10.1175/JCLI-D-12-00195.1) (Tables 3 and 4)
## CMIP5 (Geoffroy et al., 2013b; doi:10.1175/JCLI-D-12-00196.1) (Table 1)
## CMIP6 (Smith et al., 2021; doi:???) (Section 7.SM.2.1) (https://github.com/chrisroadmap/ar6/blob/main/data_output/cmip6_twolayer_tuning_params.csv)
def get_param_clim(data='CMIP6', get_irf_param=False):
    
    ## parameter correspondence
    def box_2_irf(l0, Cs, Cd, g, e):
        Cd, g = e * Cd, e * g
        b = (l0 + g) / Cs + g / Cd
        bb = (l0 + g) / Cs - g / Cd
        delta = b**2 - 4 * l0 * g / Cs / Cd
        phiF = 0.5 * Cs / g * (bb - np.sqrt(delta))
        phiS = 0.5 * Cs / g * (bb + np.sqrt(delta))
        tF = 0.5 * Cs * Cd / l0 / g * (b - np.sqrt(delta))
        tS = 0.5 * Cs * Cd / l0 / g * (b + np.sqrt(delta))
        aF = tF * l0 / Cs * phiS / (phiS - phiF)
        aS = tS * l0 / Cs * phiF / (phiF - phiS)
        return aS, tS, aF, tF, phiS, phiF

    ## initialise local dataset
    Par = xr.Dataset()
    if data not in ['CMIP5_a', 'CMIP5', 'CMIP6']:
        raise ValueError

    ## for CMIP5 data
    elif data in ['CMIP5_a', 'CMIP5']:
        Par.coords['model_clim'] = ['bcc-csm1-1', 'BNU-ESM','CanESM2', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3.6.0', 'FGOALS-s2', 'GFDL-ESM2M', 
        'GISS-E2-R', 'HadGEM2-ES', 'INM-CM4', 'IPSL-CM5A-LR', 'MIROC5', 'MPI-ESM-LR', 'MRI-CGCM3', 'NorESM1-M']
        for var in ['F4x', 'lambda', 'eheat', 'T4x', 'THs', 'THd', 'th']:
            Par[var] = np.nan * xr.zeros_like(Par.model_clim, dtype=float)
    
        ## from first paper, i.e. without the ocean heat uptake feedback
        if data == 'CMIP5_a': 
            Par['F4x'][:] = np.array([6.7, 7.4, 7.6, 7.2, 7.3, 5.1, 7.5, 6.6, 7.3, 5.9, 6.2, 6.4, 8.5, 8.2, 6.6, 6.2])
            Par['lambda'][:] = np.array([1.21, 0.93, 1.03, 1.24, 1.11, 0.61, 0.88, 1.34, 1.70, 0.65, 1.51, 0.79, 1.58, 1.14, 1.26, 1.11])
            Par['eheat'][:] = xr.ones_like(Par['lambda'])
            Par['T4x'][:] = np.array([5.6, 8.0, 7.4, 5.8, 6.5, 8.3, 8.5, 4.9, 4.3, 9.1, 4.1, 8.1, 5.4, 7.3, 5.2, 5.6])
            Par['THs'][:] = np.array([7.6, 7.4, 7.3, 6.1, 8.4, 6.0, 7.0, 8.1, 4.7, 6.5, 8.6, 7.7, 8.3, 7.3, 8.5, 8.0])
            Par['THd'][:] = np.array([53., 90., 71., 69., 99., 69., 127., 105., 126., 82., 317., 95., 145., 71., 64., 105.])
            Par['th'][:] = np.array([0.67, 0.53, 0.59, 0.93, 0.50, 0.88, 0.76, 0.90, 1.16, 0.55, 0.65, 0.59, 0.76, 0.72, 0.66, 0.88])

        ## better use the second set of parameters since they include the ocean heat uptake feedback
        elif data == 'CMIP5':
            Par['F4x'][:] = np.array([7.4, 7.3, 8.2, 8.5, 7.1, 7.0, 8.0, 7.1, 9.1, 6.8, 6.0, 6.7, 8.9, 9.4, 7.1, 7.4])
            Par['lambda'][:] = np.array([1.28, 0.92, 1.06, 1.40, 1.12, 0.68, 0.87, 1.38, 2.03, 0.61, 1.56, 0.79, 1.58, 1.21, 1.31, 1.15])
            Par['eheat'][:] = np.array([1.27, 0.98, 1.28, 1.36, 0.92, 1.82, 1.21, 1.21, 1.44, 1.54, 0.83, 1.14, 1.19, 1.42, 1.25, 1.57])
            Par['T4x'][:] = np.array([5.8, 7.9, 7.8, 6.0, 6.4, 10.2, 9.1, 5.1, 4.5, 11.1, 3.9, 8.5, 5.6, 7.8, 5.4, 6.5])
            Par['THs'][:] = np.array([8.4, 7.3, 8.0, 7.6, 8.3, 8.5, 7.5, 8.8, 6.1, 7.5, 8.5, 8.1, 8.7, 8.5, 9.3, 9.7])
            Par['THd'][:] = np.array([56., 89., 77., 72., 95., 76., 138., 112., 134., 98., 271., 100., 158., 78., 68., 121.])
            Par['th'][:] = np.array([0.59, 0.54, 0.54, 0.81, 0.51, 0.71, 0.72, 0.84, 1.06, 0.49, 0.67, 0.57, 0.73, 0.62, 0.59, 0.76])

    ## for CMIP6 data
    elif data == 'CMIP6': 
        TMP = xr.Dataset.from_dataframe(pd.read_csv('input_data/AR6/cmip6_twolayer_tuning_params.csv', index_col=0))
        TMP = TMP.sel(index=[model for model in TMP.index.values if model not in ['Mean', 'StDev']]).rename({'index': 'model_clim'})
        Par['F4x'] = TMP['q4x']
        Par['lambda'] = TMP['q4x'] / TMP['t4x']
        Par['eheat'] = TMP['eff']
        Par['T4x'] = TMP['t4x']
        Par['THs'] = TMP['cmix']
        Par['THd'] = TMP['cdeep']
        Par['th'] = TMP['gamma_2l']

    ## units
    Par['T4x'].attrs['units'] = 'K'
    Par['F4x'].attrs['units'] = 'W m-2'
    Par['THs'].attrs['units'] = Par['THd'].attrs['units'] = 'W yr m-2 K-1'
    Par['th'].attrs['units'] = Par['lambda'].attrs['units'] = 'W m-2 K-1'
    Par['eheat'].attrs['units'] = '1'

    ## also get irf parameters if requested
    if get_irf_param:
        for par, val in zip(['aS', 'tS', 'aF', 'tF', 'phiS', 'phiF'], box_2_irf(*[Par[par] for par in ['lambda', 'THs', 'THd', 'th', 'eheat']])):
            Par[par] = val
        Par['tF'].attrs['units'] = Par['tS'].attrs['units'] = 'yr'
        Par['aF'].attrs['units'] = Par['aS'].attrs['units'] = Par['phiS'].attrs['units'] = Par['phiF'].attrs['units'] = '1'
    else:
        Par = Par.drop(['lambda'])

    ## return
    return Par


##=============
## Ocean carbon
##=============

## get structural parameters
## (Strassmann et al., 2018; doi:10.5194/gmd-11-1887-2018) (Tables A2 and A3)
## & sensitivities of ocean carbon cycle
## (Arora et al., 2013; doi:10.1175/JCLI-D-12-00494.1) (from CMIP5 models' outputs)
## (Arora et al., 2019; doi:10.5194/bg-17-4173-2020) (from CMIP6 models' outputs)
def get_param_ocean(struct='Bern2.5D', get_original=False, data='CMIP6', display_fit=False):
    folder_out = 'internal_data/prior_calib/'

    ## initialise local dataset
    Par = xr.Dataset()
    Par.coords['model_struct'] = ['HILDA', 'Bern2.5D','Princeton']
    for var in ['Ho', 'Ao', 'vgx', 'To'] + [var+'_'+str(n) for var in ['aoc', 'toc'] for n in range(1, 1+5)]:
        Par[var] = np.nan * xr.zeros_like(Par['model_struct'], dtype=float)

    ## primary parameters
    Par['Ho'][:] = [75., 50.0, 50.9]
    Par['Ao'][:] = np.array([3.62, 3.5375, 3.55]) * 1E14
    Par['vgx'][:] = np.array([1/9.06, 1/7.46, 1/7.66]) * 2.123
    Par['To'][:] = [18.17, 18.30, 17.70]
    Par['gdic'] = 0.023

    ## secondary parameters
    Par['adic'] = 1E15 / Par['Ho'] / Par['Ao'] / 1026.5 / 12.0107E-6

    ## response function parameters
    ## allocation coefficients (forced to sum to 1)
    Par['aoc_1'][:] = [1 - 0.13733 - 0.051541 - 0.035033 - 0.022936, 
                       1 - 0.10292 - 0.0392835 - 0.01298 - 0.013691, 
                       1 - 0.061618 - 0.037265 - 0.019565 - 0.014818]
    Par['aoc_2'][:] = [0.13733, 0.10292, 0.061618]
    Par['aoc_3'][:] = [0.051541, 0.0392835, 0.037265]
    Par['aoc_4'][:] = [0.035033, 0.012986, 0.019565]
    Par['aoc_5'][:] = [0.022936, 0.013691, 0.014818]
    ## turnover times (merged for fastest pools)
    Par['toc_1'][:] = [(0.27830 * 0.45254 + 0.24014 * 0.03855 + 0.23337 * 2.1990) / (0.27830 + 0.24014 + 0.23337), 
                       (0.27022 * 0.07027 + 0.45937 * 0.57621 + 0.094671 * 2.6900) / (0.27022 + 0.45937 + 0.094671), 
                       (2.2745 * 1.1976 + -2.7093 * 1.5521 + 1.2817 * 2.0090) / (2.2745 + -2.7093 + 1.2817)]
    Par['toc_2'][:] = [12.038, 13.617, 16.676]
    Par['toc_3'][:] = [59.584, 86.797, 65.102]
    Par['toc_4'][:] = [237.31, 337.30, 347.58]
    Par['toc_5'][:] = [1E9, 1E9, 1E9]

    ## units
    Par['Ho'].attrs['units'] = 'm'
    Par['Ao'].attrs['units'] = 'm2'
    Par['vgx'].attrs['units'] = 'PgC ppm-1 yr-1'
    Par['To'].attrs['units'] = 'degC'
    Par['adic'].attrs['units'] = 'umol kg-1 PgC-1'
    Par['gdic'].attrs['units'] = 'K-1'
    for n in range(1, 1+5): Par['aoc_'+str(n)].attrs['units'] = '1'
    for n in range(1, 1+5): Par['toc_'+str(n)].attrs['units'] = 'yr'

    ## return structural parameters
    if get_original: return Par
    ## or choose one model and calibrate
    else: Par = Par.sel(model_struct=struct, drop=True)

    ## =======

    ## load CMIP data 
    ## (and list of simulations)
    if data == 'CMIP5': 
        Var = load_CMIP5_ocean()
        simu_list = ['1pctCO2', 'esmFixClim1', 'esmFdbk1']
    elif data == 'CMIP6': 
        Var = load_CMIP6_ocean()
        simu_list = ['1pctCO2', '1pctCO2-bgc', '1pctCO2-rad']
    else: 
        raise ValueError

    ## additional variable needed for calibration
    ## surface layer carbon (using the IRF from the chosen structure)
    toc = xr.DataArray([Par['toc_'+str(n)].values for n in range(1, 1+5)], coords={'box':range(1, 1+5)}, dims=['box'])
    aoc = xr.DataArray([Par['aoc_'+str(n)].values for n in range(1, 1+5)], coords={'box':range(1, 1+5)}, dims=['box'])
    F = aoc * (Var.fgco2 - Var.fgco2.sel(simu='piControl', drop=True))
    IRF = np.exp(-(Var.year - Var.year[0] + 0.5) / toc)
    Co = np.array([[np.array([np.convolve(F.sel(model_ocean=model, simu=simu, box=box), IRF.sel(box=box))[:len(Var.year)] for box in IRF.box]).sum(0) for model in Var.model_ocean] for simu in Var.simu])
    Var['Co'] = (('year', 'model_ocean', 'simu'), Co.transpose())

    ## pCO2 function required for calibration
    ## /!\ WARNING: must be same as in main model
    f_p = lambda dic, To: (1.5568-1.3993E-2*To)*dic + (7.4706-0.20207*To)*1E-3*dic**2. - (1.2748-0.12015*To)*1E-5*dic**3. + (2.4491-0.12639*To)*1E-7*dic**4. - (1.5768-0.15326*To)*1E-10*dic**5.

    ## initialise new parameters
    for var in ['To', 'vgx', 'ggx', 'bdic', 'gdic']:
        Par[var] = xr.zeros_like(Var.model_ocean, dtype=float)

    ## calculation of To
    Par['To'][:] = Var['tos'].sel(simu='piControl').mean('year') - 273.15

    ## calibration of vgx and ggx
    if display_fit: plt.figure()
    for m, model in enumerate(Var.model_ocean):
        
        ## get data
        dpco2 = Var.dpco2.sel(simu=simu_list, model_ocean=model).dropna('year', how='all')
        dpco2_0 = Var.dpco2.sel(simu='piControl', model_ocean=model).where(dpco2.notnull())
        fgco2 = Var.fgco2.sel(simu=simu_list, model_ocean=model).dropna('year', how='all')
        fgco2_0 = Var.fgco2.sel(simu='piControl', model_ocean=model).where(fgco2.notnull())
        tas = Var.tas.sel(simu=simu_list, model_ocean=model).dropna('year', how='all')
        tas_0 = Var.tas.sel(simu='piControl', model_ocean=model).where(tas.notnull())
        year_ok = Var.year.where((dpco2_0.notnull() & fgco2_0.notnull() & tas_0.notnull()).sum('simu') >= 2).dropna('year') - Var.year[0] + 0.5
        
        ## ad hoc restriction of calibration range
        #if model == 'NorESM2-LM': year_ok = year_ok[:70]

        ## fit data
        if len(year_ok) > 0:
            err1 = lambda par: (( -abs(par[0]) * (1 + par[1] * (tas - tas_0)) * (dpco2 - dpco2_0) - (fgco2 - fgco2_0) )**2.).where(year_ok).sum().values
            par, _, _, _, fmin_flag = fmin(err1, [0.1, 0.], disp=False, full_output=True)
            Par['vgx'][m] = abs(par[0]) if fmin_flag == 0 else np.nan
            Par['ggx'][m] = par[1] if fmin_flag == 0 else np.nan
        else:
            Par['vgx'][m], Par['ggx'][m] = np.nan, np.nan
        
        ## plot fit
        if display_fit:
            N = int(np.ceil(max(np.roots((1, 1, -len(Var.model_ocean))))))
            plt.subplot(N, N+1, m+1)
            plt.plot(fgco2 - fgco2_0, marker='+', ls='none')
            plt.plot(-Par['vgx'][m] * (1 + Par['ggx'][m] * (tas - tas_0)) * (dpco2 - dpco2_0), color='k')
            plt.title('dpco2 => fgco2 | ' + str(model.values) + '\n' + 
                'vgx = {0:.4f}, ggx = {1:.4f}'.format(Par['vgx'][m].values, Par['ggx'][m].values), fontsize='small')

    ## calibration of bdic and gdic
    if display_fit: plt.figure()
    for m, model in enumerate(Var.model_ocean):
        
        ## get data
        To = Par['To'].sel(model_ocean=model)
        Co = Var.Co.sel(simu=simu_list, model_ocean=model).dropna('year', how='all')
        Co_0 = Var.Co.sel(simu='piControl', model_ocean=model).where(Co.notnull())
        spco2 = Var.spco2.sel(simu=simu_list, model_ocean=model).dropna('year', how='all')
        spco2_0 = Var.spco2.sel(simu='piControl', model_ocean=model).where(spco2.notnull())
        tas = Var.tas.sel(simu=simu_list, model_ocean=model).dropna('year', how='all')
        tas_0 = Var.tas.sel(simu='piControl', model_ocean=model).where(tas.notnull())
        year_ok = Var.year.where((Co_0.notnull() & spco2_0.notnull() & tas_0.notnull()).sum('simu') >= 2).dropna('year') - Var.year[0] + 0.5

        ## fit data
        if len(year_ok) > 0:
            err2 = lambda par: (( (f_p(Par['adic'] / abs(par[0]) * (Co - Co_0), To) + spco2_0) * np.exp(par[1] * (tas - tas_0)) - spco2 )**2.).where(year_ok).sum().values
            par, _, _, _, fmin_flag = fmin(err2, [1., 0.], disp=False, full_output=True)
            Par['bdic'][m] = abs(par[0]) if fmin_flag == 0 else np.nan
            Par['gdic'][m] = par[1] if fmin_flag == 0 else np.nan
        else:
            Par['bdic'][m], Par['gdic'][m] = np.nan, np.nan
        
        ## plot fit
        if display_fit:
            N = int(np.ceil(max(np.roots((1, 1, -len(Var.model_ocean))))))
            plt.subplot(N, N+1, m+1)
            plt.plot(spco2, marker='+', ls='none')
            plt.plot((f_p(Par['adic'] / Par['bdic'][m] * (Co - Co_0), To) + spco2_0) * np.exp(Par['gdic'][m] * (tas - tas_0)), color='k')
            plt.title('Co, tas => spco2 | ' + str(model.values) + '\n' + 
                'bdic = {0:.4f}, gdic = {1:.4f}'.format(Par['bdic'][m].values, Par['gdic'][m].values), fontsize='small')

    ## units
    Par['To'].attrs['units'] = 'degC'
    Par['vgx'].attrs['units'] = 'PgC ppm-1 yr-1'
    Par['ggx'].attrs['units'] = Par['gdic'].attrs['units'] = 'K-1'
    Par['bdic'].attrs['units'] = '1'

    ## save and return
    Par = Par.drop(['Ho', 'Ao'])
    Par.to_netcdf(folder_out + 'param_ocean.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Par})
    return Par


##===============
## PI Land carbon
##===============

## get preindustrial land carbon cycle
## (LeQuéré et al., 2018; doi:10.5194/essd-10-2141-2018) (from TRENDYv7 models' ouputs)
def get_param_landPI(data='TRENDYv7'):

    ## initialise local dataset
    Par = xr.Dataset()

    ## load TRENDY data
    if data == 'TRENDYv7': Var = load_TRENDYv7_land()
    else: raise ValueError

    ## take average over control
    Var = Var.mean('year')

    ## intermediate variables (using surrogate when unavailable)
    Fmort = (Var.fVegLitter.fillna(0) + Var.fVegSoil.fillna(0)).where(Var.fVegLitter.notnull() | Var.fVegSoil.notnull(), Var.npp - Var.fFire.fillna(0) - Var.fLateral)
    Fstab = Var.fLitterSoil.where(Var.fLitterSoil.notnull(), 0.5 * Var.rh)

    ## calculate paramaters
    Par['npp0'] = Var.npp
    Par['vfire'] = Var.fFire / Var.cVeg
    Par['vharv'] = Var.fLateral / Var.cVeg
    Par['vmort'] = Fmort / Var.cVeg
    Par['vstab'] = Fstab / Var.cLitter
    Par['vrh1'] = (Var.rh - Fstab) / Var.cLitter
    Par['vrh2'] = Fstab / Var.cSoil

    ## units
    for var in Par:
        if var in ['npp0']: Par[var].attrs['units'] = 'PgC yr-1'
        else: Par[var].attrs['units'] = 'yr-1'

    ## return
    return Par


##============
## Land carbon
##============

## get sensitivities of land carbon cycle
## (Arora et al., 2013; doi:10.1175/JCLI-D-12-00494.1) (from CMIP5 outputs)
## (Arora et al., 2019; doi:10.5194/bg-17-4173-2020) (from CMIP6 outputs)
def get_param_land(data='CMIP6', log_npp=False, display_fit=False):
    folder_out = 'internal_data/prior_calib/'

    ## load CMIP data
    ## (and list of simulations)
    if data == 'CMIP5': 
        Var = load_CMIP5_land()
        simu_list = ['1pctCO2', 'esmFixClim1', 'esmFdbk1']
    elif data == 'CMIP6': 
        Var = load_CMIP6_land()
        simu_list = ['1pctCO2', '1pctCO2-bgc', '1pctCO2-rad']
    else: 
        raise ValueError

    ## initialise local dataset
    Par = xr.Dataset()
    for var in ['bnpp', 'anpp', 'gnpp', 'bfire', 'gfire', 'brh', 'grh']:
        Par[var] = xr.zeros_like(Var.model_land, dtype=float)

    ## calibration of bnpp and gnpp
    if display_fit: plt.figure()
    for m, model in enumerate(Var.model_land):
        
        ## get data
        npp = Var.npp.sel(simu=simu_list, model_land=model).dropna('year', how='all')
        npp_0 = Var.npp.sel(simu='piControl', model_land=model).where(npp.notnull())
        co2 = Var.co2.sel(simu=simu_list, model_land=model).dropna('year', how='all')
        co2_0 = Var.co2.sel(simu='piControl', model_land=model).where(co2.notnull())
        tas = Var.tas.sel(simu=simu_list, model_land=model).dropna('year', how='all')
        tas_0 = Var.tas.sel(simu='piControl', model_land=model).where(tas.notnull())
        year_ok = Var.year.where((npp_0.notnull() & co2_0.notnull() & tas_0.notnull()).sum('simu') >= 2).dropna('year') - Var.year[0] + 0.5
        
        ## fit data
        if len(year_ok) > 0:
            if log_npp:
                err1 = lambda par: ((( 1 + abs(par[0]) / 1E-9 * (1 - (co2 / co2_0)**-1E-9)) * (1 + par[2] * (tas - tas_0))  - (npp / npp_0) )**2.).sum().values
            else:
                err1 = lambda par: ((( 1 + abs(par[0]) / abs(par[1]) * (1 - (co2 / co2_0)**-abs(par[1]))) * (1 + par[2] * (tas - tas_0))  - (npp / npp_0) )**2.).sum().values
            par, _, _, _, fmin_flag = fmin(err1, [0.6, 1., 0.], disp=False, full_output=True)
            Par['bnpp'][m] = abs(par[0]) if fmin_flag == 0 else np.nan
            Par['anpp'][m] = abs(par[1]) * (1-log_npp) + 1E-9 * log_npp if fmin_flag == 0 else np.nan
            Par['gnpp'][m] = par[2] if fmin_flag == 0 else np.nan
        else:
            Par['bnpp'][m], Par['anpp'][m], Par['gnpp'][m] = np.nan, np.nan, np.nan
        
        ## plot fit
        if display_fit:
            N = int(np.ceil(max(np.roots((1, 1, -len(Var.model_land))))))
            plt.subplot(N, N+1, m+1)
            plt.plot(npp / npp_0, marker='+', ls='none')
            plt.plot(( 1 + Par['bnpp'][m] / Par['anpp'][m] * (1 - (co2 / co2_0)**-Par['anpp'][m])) * (1 + Par['gnpp'][m] * (tas - tas_0)), color='k')
            plt.title('co2, tas => npp | ' + str(model.values) + '\n' + 
                'bnpp = {0:.4f}, anpp = {1:.4f}, gnpp = {2:.4f}'.format(Par['bnpp'][m].values, Par['anpp'][m].values, Par['gnpp'][m].values), fontsize='small')

    ## calibration of bef and gef
    if display_fit: plt.figure()
    for m, model in enumerate(Var.model_land):
        
        ## get data
        fFire = Var.fFire.sel(simu=simu_list, model_land=model).dropna('year', how='all')
        fFire_0 = Var.fFire.sel(simu='piControl', model_land=model).where(fFire.notnull())
        cVeg = Var.cVeg.sel(simu=simu_list, model_land=model).dropna('year', how='all')
        cVeg_0 = Var.cVeg.sel(simu='piControl', model_land=model).where(cVeg.notnull())
        co2 = Var.co2.sel(simu=simu_list, model_land=model).dropna('year', how='all')
        co2_0 = Var.co2.sel(simu='piControl', model_land=model).where(co2.notnull())
        tas = Var.tas.sel(simu=simu_list, model_land=model).dropna('year', how='all')
        tas_0 = Var.tas.sel(simu='piControl', model_land=model).where(tas.notnull())
        year_ok = Var.year.where((fFire_0.notnull() & cVeg_0.notnull() & co2_0.notnull() & tas_0.notnull()).sum('simu') >= 2).dropna('year') - Var.year[0] + 0.5

        ## fit data
        if len(year_ok) > 0:
            err2 = lambda par: (( (1 + par[0] * (co2 / co2_0 - 1)) * (1 + par[1] * (tas - tas_0))  - (fFire / fFire_0 * cVeg_0 / cVeg) )**2.).sum().values
            par, _, _, _, fmin_flag = fmin(err2, [0., 0.], disp=False, full_output=True)
            Par['bfire'][m] = par[0] if fmin_flag == 0 else np.nan
            Par['gfire'][m] = par[1] if fmin_flag == 0 else np.nan
        else:
            Par['bfire'][m], Par['gfire'][m] = np.nan, np.nan
        
        ## plot fit
        if display_fit:
            N = int(np.ceil(max(np.roots((1, 1, -len(Var.model_land))))))
            plt.subplot(N, N+1, m+1)
            plt.plot(fFire / fFire_0 * cVeg_0 / cVeg, marker='+', ls='none')
            plt.plot((1 + Par['bef'][m] * (co2 / co2_0 - 1)) * (1 + Par['gef'][m] * (tas - tas_0)), color='k')
            plt.title('co2, tas => fFire | ' + str(model.values) + '\n' + 
                'bfire = {0:.4f}, gfire = {1:.4f}'.format(Par['bef'][m].values, Par['gef'][m].values), fontsize='small')

    ## calibration of brh and grh
    if display_fit: plt.figure()
    for m, model in enumerate(Var.model_land):
        
        ## get data
        rh = Var.rh.sel(simu=simu_list, model_land=model).dropna('year', how='all')
        rh_0 = Var.rh.sel(simu='piControl', model_land=model).where(rh.notnull())
        cLitter = Var.cLitter.sel(simu=simu_list, model_land=model).dropna('year', how='all')
        cLitter_0 = Var.cLitter.sel(simu='piControl', model_land=model).where(cLitter.notnull())
        cSoil = Var.cSoil.sel(simu=simu_list, model_land=model).dropna('year', how='all')
        cSoil_0 = Var.cSoil.sel(simu='piControl', model_land=model).where(cSoil.notnull())
        co2 = Var.co2.sel(simu=simu_list, model_land=model).dropna('year', how='all')
        co2_0 = Var.co2.sel(simu='piControl', model_land=model).where(co2.notnull())
        tas = Var.tas.sel(simu=simu_list, model_land=model).dropna('year', how='all')
        tas_0 = Var.tas.sel(simu='piControl', model_land=model).where(tas.notnull())
        year_ok = Var.year.where((rh_0.notnull() & cLitter_0.notnull() & cSoil_0.notnull() & co2_0.notnull() & tas_0.notnull()).sum('simu') >= 2).dropna('year') - Var.year[0] + 0.5

        ## fit data
        if len(year_ok) > 0:
            err3 = lambda par: (( (1 + abs(par[0]) * (cLitter / (cSoil + cLitter) * (cSoil_0 + cLitter_0) / cLitter_0 - 1)) * np.exp(abs(par[1]) * (tas - tas_0))  - (rh / rh_0 * (cSoil_0 + cLitter_0) / (cSoil + cLitter)) )**2.).sum().values
            par, _, _, _, fmin_flag = fmin(err3, [0., 0.], disp=False, full_output=True)
            Par['brh'][m] = abs(par[0]) if fmin_flag == 0 else np.nan
            Par['grh'][m] = abs(par[1]) if fmin_flag == 0 else np.nan
        else:
            Par['brh'][m], Par['grh'][m] = np.nan, np.nan
        
        ## plot fit
        if display_fit:
            N = int(np.ceil(max(np.roots((1, 1, -len(Var.model_land))))))
            plt.subplot(N, N+1, m+1)
            plt.plot(rh / rh_0 * (cSoil_0 + cLitter_0) / (cSoil + cLitter), marker='+', ls='none')
            plt.plot((1 + Par['brh'][m] * (cLitter / (cSoil + cLitter) * (cSoil_0 + cLitter_0) / cLitter_0 - 1)) * np.exp(Par['grh'][m] * (tas - tas_0)), color='k')
            plt.title('co2, tas => rh | ' + str(model.values) + '\n' + 
                'brh = {0:.4f}, grh = {1:.4f}'.format(Par['brh'][m].values, Par['grh'][m].values), fontsize='small')

    ## units
    Par['bnpp'].attrs['units'] = Par['anpp'].attrs['units'] = Par['bfire'].attrs['units'] = Par['brh'].attrs['units'] = '1'
    Par['gnpp'].attrs['units'] = Par['gfire'].attrs['units'] = Par['grh'].attrs['units'] = 'K-1'

    ## save and return
    Par.to_netcdf(folder_out + 'param_land.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Par})
    return Par


##===============
## Sea level rise
##===============

## get structural parameters
## (Mengel et al., 2016; doi:10.1073/pnas.1500515113)
## & sensitivities of sea level rise
## (Edwards et al., 2021; doi:10.1038/s41586-021-03302-y) (from compiled outputs)
def get_param_slr(fixed_param, display_fit=False):
    folder_out='internal_data/prior_calib/'

    ## get fixed parameters
    tgla_bg, tgis_bg, tais_bg = [fixed_param[var] for var in ['tgla', 'tgis', 'tais']]
    Lgla_bg, Lais_bg = [fixed_param[var] for var in ['Lgla', 'Lais']]

    ## load Edwards et al. data
    GLA = load_Edwards_2021_slr(ice='gla')
    GIS = load_Edwards_2021_slr(ice='gis')
    AIS = load_Edwards_2021_slr(ice='ais')

    ## restrict config/model to available data and more than one experiment
    GLA_model = xr.concat([mod for mod in GLA.model if (GLA.slr.sel(model=mod).isel(year=-1).notnull().sum('region') >= len(GLA.region)-1).sum() > 1], dim='model')
    GIS_config = xr.concat([cfg for cfg in GIS.config if GIS.slr.sel(config=cfg).isel(year=-1).sel(region='ALL', drop=True).notnull().sum() > 1], dim='config')
    AIS_config = xr.concat([cfg for cfg in AIS.config if (AIS.slr.sel(config=cfg).isel(year=-1).notnull().sum('region') == len(AIS.region)).sum() > 1], dim='config')

    ## initialise local dataset
    Par = xr.Dataset()
    for var in ['lgla0', 'Lgla', 'Ggla1', 'Ggla3', 'tgla', 'ggla']:
        Par[var] = np.nan * xr.zeros_like(GLA_model, dtype=float).rename({'model': 'cfg_gla'})
    for var in ['lgis0', 'Lgis1', 'Lgis3', 'tgis']:
        Par[var] = np.nan * xr.zeros_like(GIS_config, dtype=float).rename({'config': 'cfg_gis'}).drop('model')
    for var in ['lais0', 'Lais_smb', 'Lais', 'tais', 'aais']:
        Par[var] = np.nan * xr.zeros_like(AIS_config, dtype=float).rename({'config': 'cfg_ais'}).drop('model')

    ## calibration of GLA
    if display_fit: plt.figure()
    for m, model in enumerate(GLA_model):
        
        ## get data
        H = GLA.slr.sum('region', min_count=1).sel(model=model).dropna('year', how='all').dropna('exp')
        H_0 = 67.2 - 3 * 0.62 # AR6 WG1 Table 9.5
        d_H = H.differentiate('year')
        d_H0 = 0.58 # AR6 WG1 Table 9.5
        D_tas = GLA.tas.sel(year=H.year, exp=H.exp)
       
        ## fit data
        err1 = lambda par: (( d_H0 + (Lgla_bg * (1 - np.exp(-par[2] * D_tas - par[3] * D_tas**3)) - (H + H_0)) / tgla_bg * np.exp(par[5] * D_tas)  - d_H )**2.).sum().values
        par, _, _, _, fmin_flag = fmin(err1, [d_H0, Lgla_bg, 0.2, 0.0, tgla_bg, 0.0], disp=False, full_output=True, maxiter=2000, maxfun=2000)
        Par['lgla0'][m] = d_H0 if fmin_flag == 0 else np.nan
        Par['Lgla'][m] = Lgla_bg if fmin_flag == 0 else np.nan
        Par['Ggla1'][m] = par[2] if fmin_flag == 0 else np.nan
        Par['Ggla3'][m] = par[3] if fmin_flag == 0 else np.nan
        Par['tgla'][m] = tgla_bg if fmin_flag == 0 else np.nan
        Par['ggla'][m] = par[5] if fmin_flag == 0 else np.nan
        
        ## plot fit
        if display_fit:
            N = int(np.ceil(max(np.roots((1, 1, -len(GLA_model))))))
            plt.subplot(N, N+1, m+1)
            plt.plot(d_H, marker='+', ls='none')
            plt.plot(Par['lgla0'][m] + (Par['Lgla'][m] * (1 - np.exp(-Par['Ggla1'][m] * D_tas - Par['Ggla3'][m] * D_tas**3)) - (H + H_0)) / Par['tgla'][m] * np.exp(Par['ggla'][m] * D_tas), color='k')
            plt.title('tas => d_H | ' + str(model.values) + '\n' + 
                'Ggla1 = {0:.1f}, Ggla3 = {1:.1f}, ggla = {2:.1f}'.format(*[Par[var][m].values for var in ['Ggla1', 'Ggla3', 'ggla']]), fontsize='small')

    ## calibration of GIS
    if display_fit: plt.figure()
    for m, model in enumerate(GIS_config):
        
        ## get data
        H = GIS.slr.sum('region', min_count=1).sel(config=model).dropna('year', how='all').dropna('exp')
        H_0 = 40.4 - 3 * 0.63 # AR6 WG1 Table 9.5
        d_H = H.differentiate('year')
        d_H0 = 0.33 # AR6 WG1 Table 9.5
        D_tas = GIS.tas.sel(year=H.year, exp=H.exp)
       
        ## fit data
        err2 = lambda par: (( d_H0 + (par[1] * D_tas + par[2] * D_tas**3 - (H + H_0)) / tgis_bg  - d_H )**2.).sum().values
        par, _, _, _, fmin_flag = fmin(err2, [d_H0, 3., 2., tgis_bg], disp=False, full_output=True, maxiter=1000, maxfun=1000)
        Par['lgis0'][m] = d_H0 if fmin_flag == 0 else np.nan
        Par['Lgis1'][m] = par[1] if fmin_flag == 0 else np.nan
        Par['Lgis3'][m] = par[2] if fmin_flag == 0 else np.nan
        Par['tgis'][m] = tgis_bg if fmin_flag == 0 else np.nan
        
        ## plot fit
        if display_fit:
            N = int(np.ceil(max(np.roots((1, 1, -len(GIS_config))))))
            plt.subplot(N, N+1, m+1)
            plt.plot(d_H, marker='+', ls='none')
            plt.plot(Par['lgis0'][m] + (Par['Lgis1'][m] * D_tas + Par['Lgis3'][m] * D_tas**3 - (H + H_0)) / Par['tgis'][m], color='k')
            plt.title('tas => d_H | ' + str(model.values) + '\n' + 
                'Lgis1 = {0:.1f}, Lgis3 = {1:.1f}'.format(*[Par[var][m].values for var in ['Lgis1', 'Lgis3']]), fontsize='small')

    ## calibration of AIS
    if display_fit: plt.figure()
    for m, model in enumerate(AIS_config):
        
        ## get data
        H = AIS.slr.sum('region', min_count=1).sel(config=model).dropna('year', how='all').dropna('exp')
        H_0 = 6.7 - 3 * 0.37 # AR6 WG1 Table 9.5
        d_H = H.differentiate('year')
        d_H0 = 0.00 # AR6 WG1 Table 9.5
        D_tas = AIS.tas.sel(year=H.year, exp=H.exp)
        sum_tas_0 = load_gmst().T.sel(year=slice(1901, 2014)).sum('year').mean('data').values
        sum_tas = D_tas.cumsum('year') + sum_tas_0
       
        ## fit data
        err3 = lambda par: (( -abs(par[1]) * D_tas + d_H0 + (Lais_bg * D_tas - (H + H_0 - -abs(par[1]) * sum_tas)) / tais_bg * (1 + par[4] * (H + H_0 - -abs(par[1]) * sum_tas)) - d_H )**2.).sum().values
        par, _, _, _, fmin_flag = fmin(err3, [d_H0, 0.1, Lais_bg, tais_bg, 0.], disp=False, full_output=True, maxiter=3000, maxfun=3000)
        Par['lais0'][m] = d_H0 if fmin_flag == 0 else np.nan
        Par['Lais_smb'][m] = abs(par[1]) if fmin_flag == 0 else np.nan
        Par['Lais'][m] = Lais_bg if fmin_flag == 0 else np.nan
        Par['tais'][m] = tais_bg if fmin_flag == 0 else np.nan
        Par['aais'][m] = par[4] if fmin_flag == 0 else np.nan

        ## plot fit
        if display_fit:
            N = int(np.ceil(max(np.roots((1, 1, -len(AIS_config))))))
            plt.subplot(N, N+1, m+1)
            plt.plot(d_H, marker='+', ls='none')
            plt.plot(-Par['Lais_smb'][m] * D_tas + Par['lais0'][m] + (Par['Lais'][m] * D_tas - (H + H_0 - -Par['Lais_smb'][m] * sum_tas)) / Par['tais'][m] * (1 + Par['aais'][m] * (H + H_0 - -Par['Lais_smb'][m] * sum_tas)), color='k')
            plt.title('tas => d_H | ' + str(model.values) + '\n' + 
                'Lais_smb = {0:.1f}, aais = {1:.1f}'.format(*[Par[var][m].values for var in ['Lais_smb', 'aais']]), fontsize='small')

    ## remove config axis is constant
    for var in Par:
        if all(np.isclose(Par[var], Par[var].mean())):
            Par[var] = Par[var].mean()

    ## units
    Par['lgla0'].attrs['units'] = Par['lgis0'].attrs['units'] = Par['lais0'].attrs['units'] = 'mm yr-1'
    Par['Lgla'].attrs['units'] = 'mm'
    Par['Ggla1'].attrs['units'] = Par['ggla'].attrs['units'] = 'K-1'
    Par['Ggla3'].attrs['units'] = 'K-3'
    Par['tgla'].attrs['units'] = Par['tgis'].attrs['units'] = Par['tais'].attrs['units'] = 'yr'
    Par['Lgis1'].attrs['units'] = Par['Lais'].attrs['units'] = 'mm K-1'
    Par['Lgis3'].attrs['units'] = 'mm K-3'
    Par['Lais_smb'].attrs['units'] = 'mm yr-1 K-1'
    Par['aais'].attrs['units'] = 'mm-1'

    ## save and return
    Par.to_netcdf(folder_out + 'param_slr.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Par})
    return Par

