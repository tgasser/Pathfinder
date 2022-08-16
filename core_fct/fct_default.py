"""
Copyright IIASA (International Institute for Applied Systems Analysis), 2019-2022
Contributor(s): Thomas Gasser (gasser@iiasa.ac.at), Thomas Bossy (thomas.bossy@lsce.ipsl.fr)

This software is a computer program whose purpose is simulate a large number of global climate (and environmental) change scenarios.

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
"""


import numpy as np
import xarray as xr

from core_fct.fct_load import load_GCB, load_IPCC_AR5_WG1, load_IPCC_AR6_WG1, load_NOAA_ESRL, load_gmst
from core_fct.fct_param import get_ecs_distrib, get_param_clim, get_param_ocean, get_param_landPI, get_param_land, get_param_slr


##################################################
## 1. DEFAULT PARAMETERS
##################################################

## get default parameters
## (a lot of hard-coded stuff)
def get_dflt_param(redo_prior_calib=True, ipcc='AR6', passiveC_on=True, toc_adjust=True, calib_ecs=True):
    assert ipcc in ['AR6', 'AR5']

    ## INITIALISATION
    ## list of uncertain parameters
    param_unc = ['phi', 'T2x', 'THs', 'THd', 'th', 'eheat'] # climate
    param_unc += ['aOHC', 'Lthx', 'Lgis1', 'Lgis3', 'Ggla1', 'Ggla3', 'ggla', 'Lais_smb', 'aais'] + ['lgla0', 'lgis0', 'lais0']# slr
    param_unc += (toc_adjust) * ['k_toc'] + ['vgx', 'ggx', 'To', 'bdic', 'gdic'] # ocean
    param_unc += ['npp0', 'vfire', 'vharv', 'vmort', 'vstab', 'vrh1', 'vrh2', 'bnpp', 'anpp', 'gnpp', 'bfire', 'gfire', 'brh', 'grh'] + (passiveC_on) * ['apass'] # land
    param_unc += ['ga', 'ka', 'k_tth', 'Cfr0'] # pf
    param_unc += ['CO2pi'] # atmo

    ## list of fixed parameters
    param_fix = (calib_ecs) * ['T2x0'] # climate
    param_fix = ['Lgla', 'Lais'] + [] + ['tgla', 'tgis', 'tais']  # slr
    param_fix += ['adic', 'aoc_1', 'aoc_2', 'aoc_3', 'aoc_4', 'aoc_5', 'toc_1', 'toc_2', 'toc_3', 'toc_4', 'toc_5'] + (not toc_adjust) * ['k_toc'] # ocean
    param_fix += ['vrh3'] + (not passiveC_on) * ['apass'] # land
    param_fix += ['aLST', 'grt1', 'grt2', 'krt', 'amin', 'vthaw', 'vfroz', 'ath_1', 'ath_2', 'ath_3', 'tth_1', 'tth_2', 'tth_3'] # pf
    param_fix += ['aCO2', 'k_pH'] # atmo

    ## creating final dataset
    Par0 = xr.Dataset()
    Par0.coords['stat'] = ['mean', 'std']
    for var in param_unc:
        Par0[var] = np.nan * xr.zeros_like(Par0['stat'], dtype=float)
    for var in param_fix:
        Par0[var] = np.nan


    ##========
    ## Climate
    ##========

    ## CO2 radiative forcing parameter
    ## AR6: expression requires N2O (Table 7.SM.1), so ignored for now and use AR5
    ## AR5 (bg): (Myhre et al., 2013; doi:10.1017/CBO9781107415324.018) (Table 8.SM.1)
    ## AR5 (unc): taken from main RF CO2 estimates in main chapter
    Par0['phi'][:] = [5.35, 5.35 * 0.1]
    Par0['phi'].attrs['units'] = 'W m-2'
    Par0['phi'].attrs['bounds'] = (0., np.inf)


    ## climate response parameters
    ## load from CMIP models
    ## load and calculate
    Par_tmp = get_param_clim(data='CMIP6' * (ipcc=='AR6') + 'CMIP5' * (ipcc=='AR5'))
    Par_tmp['T2x'] = Par_tmp['T4x'] / 2.
    Par_tmp['T2x'].attrs['units'] = 'K'
    ## take mean and unbiased std
    for var in ['T2x', 'THs', 'THd', 'th', 'eheat']:
        Par0[var][:] = [Par_tmp[var].mean(), Par_tmp[var].std() * np.sqrt((len(Par_tmp[var]) - 1.) / ((len(Par_tmp[var]) - 1.5)))]
        Par0[var].attrs['units'] = Par_tmp[var].units
        Par0[var].attrs['bounds'] = (0., np.inf)


    ## minimum ECS in the distribution (somewhat akin to T2x from Plank feedback)
    ## note: not really a parameter of the model, but required to calibrate T2x
    Par0['T2x0'] = get_ecs_distrib(data=ipcc)['T2x0']
    Par0['T2x0'].attrs['units'] = 'K'


    ##==========
    ## Sea level
    ##==========

    ## fraction of energy warming up the ocean
    ## AR6 (bg): (Forster et al., 2021; doi:10.1017/9781009157896.011) (Table 7.1)
    ## AR6 (unc): derived from same table assuming logit-norm distribution
    ## AR5 (bg): (Rhein et al., 2013; doi:10.1017/CBO9781107415324.010) (Box 3.1)
    ## AR5 (unc): assumed
    if ipcc == 'AR6': 
        Par0['aOHC'][:] = [0.910, 0.015]
    elif ipcc == 'AR5': 
        Par0['aOHC'][:] = [0.93, 0.02]
    Par0['aOHC'].attrs['units'] = '1'
    Par0['aOHC'].attrs['bounds'] = (0., 1.)


    ## linear thermosteric SLR
    ## AR6: (Fox-Kemper et al., 2021; doi:10.1017/9781009157896.011) (Section 9.2.4.1)
    ## AR5: (Kuhlbrodt & Gregory, 2012; doi:10.1029/2012GL052952) (CMIP5 value)
    if ipcc == 'AR6': 
        Par0['Lthx'][:] = np.array([113, 13]) * (3600*24*365.25) * 510.1E12 / 1E24
    elif ipcc == 'AR5': 
        Par0['Lthx'][:] = np.array([110, 10]) * (3600*24*365.25) * 510.1E12 / 1E24
    Par0['Lthx'].attrs['units'] = 'mm m2 W-1 yr-1'
    Par0['Lthx'].attrs['bounds'] = (0., np.inf)


    ## characteristic times for ice components
    ## (Mengel et al., 2016; doi:10.1073/pnas.1500515113) (Table S1)
    ## note: derived assuming log-norm distribution and 90% range provided
    Par0['tgla'] = 190. # [190., 62.] # from (98, 295)
    Par0['tgis'] = 481. # [481., 292.] # from (99.7, 927)
    Par0['tais'] = 2093. # [2093., 482.] # from (1350, 2910)
    Par0['tgla'].attrs['units'] = Par0['tgis'].attrs['units'] = Par0['tais'].attrs['units'] = 'yr'
    #Par0['tgla'].attrs['bounds'] = Par0['tgis'].attrs['bounds'] = Par0['tais'].attrs['bounds'] = (0., np.inf)


    ## holocene trends in ice components
    ## (Fox-Kemper et al., 2021; doi:10.1017/9781009157896.011) (Table 9.5)
    ## note: taken as earliest period available; likely overestimated but range increased (from 90% to 1-sigma)
    Par0['lgla0'][:] = [0.58, 0.5*(0.82-0.34)]
    Par0['lgis0'][:] = [0.33, 0.5*(0.47-0.18)]
    Par0['lais0'][:] = [0.00, 0.5*(0.11+0.10)]
    Par0['lgla0'].attrs['units'] = Par0['lgis0'].attrs['units'] = Par0['lais0'].attrs['units'] = 'mm'
    Par0['lgla0'].attrs['bounds'] = Par0['lgis0'].attrs['bounds'] = Par0['lais0'].attrs['bounds'] = (-np.inf, np.inf)


    ## maximum contribution from glaciers
    ## (Fox-Kemper et al., 2021; doi:10.1017/9781009157896.011) (Section 9.6.3.2 & Table 9.5)
    Par0['Lgla'] = 320. + np.round(67.2 - 7.5)
    Par0['Lgla'].attrs['units'] = 'mm'


    ## equilibrium sensitivity of AIS
    ## (Church et al., 2013; doi:10.1017/CBO9781107415324.026) (Figure 13.14)
    Par0['Lais'] = 1200.
    Par0['Lais'].attrs['units'] = 'mm K-1'


    ## SLR sensitivity parameters
    ## load from Edwards et al. models
    if redo_prior_calib:
        print('recalibrating SLR')
        Par_tmp = get_param_slr(fixed_param={var: Par0[var] if Par0[var].size==1 else Par0[var][0] for var in ['tgla', 'tgis', 'tais', 'Lgla', 'Lais']})
    else:
        with xr.open_dataset('internal_data/prior_calib/param_slr.nc') as TMP: Par_tmp = TMP.load()
    ## take mean and unbiased std
    for var in ['Ggla1', 'Ggla3', 'ggla', 'Lgis1', 'Lgis3', 'Lais_smb', 'aais']:
        Par0[var][:] = [Par_tmp[var].mean(), Par_tmp[var].std() * np.sqrt((len(Par_tmp[var]) - 1.) / ((len(Par_tmp[var]) - 1.5)))]
        Par0[var].attrs['units'] = Par_tmp[var].units
    Par0['Ggla1'].attrs['bounds'] = Par0['aais'].attrs['bounds'] = (-np.inf, np.inf)
    Par0['Ggla3'].attrs['bounds'] = Par0['ggla'].attrs['bounds'] = Par0['Lgis1'].attrs['bounds'] = Par0['Lgis3'].attrs['bounds'] = Par0['Lais_smb'].attrs['bounds'] = (0., np.inf)


    ##=============
    ## Ocean carbon
    ##=============

    ## ocean carbon structural parameters
    ## (Strassmann et al., 2018; doi:10.5194/gmd-11-1887-2018) (Tables A2 and A3; Princeton preferred as more stable without solver)
    ## load
    if redo_prior_calib:
        print('recalibrating ocean C')
        Par_tmp = get_param_ocean(struct='Princeton', data='CMIP6' * (ipcc=='AR6') + 'CMIP5' * (ipcc=='AR5'))
    else:
        with xr.open_dataset('internal_data/prior_calib/param_ocean.nc') as TMP: Par_tmp = TMP.load()
    ## assign
    for var in ['adic', 'aoc_1', 'aoc_2', 'aoc_3', 'aoc_4', 'aoc_5', 'toc_1', 'toc_2', 'toc_3', 'toc_4', 'toc_5']:
        Par0[var] = Par_tmp[var]
    

    ## adjustment of ocean timescales
    ## note: arbitrary uncertainty if turned on
    if toc_adjust:
        Par0['k_toc'][:] = [1., 0.2]
        Par0['k_toc'].attrs['bounds'] = (0., np.inf)
    else:
        Par0['k_toc'] = 1.
    Par0['k_toc'].attrs['units'] = '1'


    ## ocean sensitivity parameters
    ## already loaded from CMIP models
    ## take mean and unbiased std
    for var in ['vgx', 'ggx', 'To', 'bdic', 'gdic']:
        Par0[var][:] = [Par_tmp[var].mean(), Par_tmp[var].std() * np.sqrt((len(Par_tmp[var]) - 1.) / ((len(Par_tmp[var]) - 1.5)))]
        Par0[var].attrs['units'] = Par_tmp[var].units
    for var in ['vgx', 'To', 'bdic']: Par0[var].attrs['bounds'] = (0., np.inf)
    for var in ['ggx', 'gdic']: Par0[var].attrs['bounds'] = (-np.inf, np.inf)


    ##============
    ## Land carbon
    ##============

    ## land preindustrial state
    ## load from TRENDY models
    Par_tmp = get_param_landPI(data='TRENDYv7')
    ## take mean and unbiased std
    for var in ['npp0', 'vfire', 'vharv', 'vmort', 'vstab', 'vrh1', 'vrh2']:
        Par0[var][:] = [Par_tmp[var].mean(), Par_tmp[var].std() * np.sqrt((len(Par_tmp[var]) - 1.) / ((len(Par_tmp[var]) - 1.5)))]
        Par0[var].attrs['units'] = Par_tmp[var].units
        Par0[var].attrs['bounds'] = (0., np.inf)
    ## rename respiration rate of the last soil box
    Par0 = Par0.rename({'vrh2':'vrh23'})

    ## turnover rate of passive carbon
    ## (He et al., 2016; doi:10.1126/science.aad4273) (Table 1; central value only)
    Par0['vrh3'] = (1185 * 10.2)**-1
    Par0['vrh3'].attrs['units'] = 'yr-1'

    ## fraction of passive carbon
    ## (He et al., 2016; doi:10.1126/science.aad4273) (Table S5)
    if passiveC_on:
        Par0['apass'][0] = np.mean([666. / (65 + 666), 873. / (492 + 873), 817. / (737 + 817)])
        Par0['apass'][1] = np.std([666. / (65 + 666), 873. / (492 + 873), 817. / (737 + 817)]) * np.sqrt((3 - 1) / (3 - 1.5))
        Par0['apass'].attrs['bounds'] = (0., 1.)
    else:
        Par0['apass'] = 0.
    Par0['apass'].attrs['units'] = '1'

    ## land sensitivity parameters
    ## load from CMIP models
    if redo_prior_calib:
        print('recalibrating land C')
        Par_tmp = get_param_land(log_npp=False, data='CMIP6' * (ipcc=='AR6') + 'CMIP5' * (ipcc=='AR5'))
    else:
        with xr.open_dataset('internal_data/prior_calib/param_land.nc') as TMP: Par_tmp = TMP.load()
    ## take mean and unbiased std
    for var in ['bnpp', 'anpp', 'gnpp', 'bfire', 'gfire', 'brh', 'grh']:
        Par0[var][:] = [Par_tmp[var].mean(), Par_tmp[var].std() * np.sqrt((len(Par_tmp[var]) - 1.) / ((len(Par_tmp[var]) - 1.5)))]
        Par0[var].attrs['units'] = Par_tmp[var].units
    for var in ['bnpp', 'anpp', 'brh', 'grh']: Par0[var].attrs['bounds'] = (0., np.inf)
    for var in ['gnpp', 'bfire', 'gfire']: Par0[var].attrs['bounds'] = (-np.inf, np.inf)


    ##==================
    ## Permafrost carbon
    ##==================

    ## fixed and best-guess parameters
    ## (Gasser et al., 2018; ; doi:10.1038/s41561-018-0227-0)
    ## note: recalibrated on global aggregate and multi-model average, with one extra model (UVic)
    Par0['Cfr0'][0] = 545.7
    Par0['aLST'] = 1.872
    Par0['grt1'] = 0.1223
    Par0['grt2'] = 0.002946
    Par0['krt'] = 1.341
    Par0['amin'] = 0.9811
    Par0['ka'][0] = 2.655
    Par0['ga'][0] = 0.1312
    Par0['vthaw'] = 0.1411
    Par0['vfroz'] = 0.01089
    Par0['ath_1'] = 0.05086
    Par0['ath_2'] = 0.1167
    Par0['ath_3'] = 0.8325
    Par0['tth_1'] = 18.23
    Par0['tth_2'] = 251.5
    Par0['tth_3'] = 3494.
    Par0['k_tth'][0] = 1.

    ## units
    for var in ['aLST', 'krt', 'amin', 'ka', 'ath_1', 'ath_2', 'ath_3', 'k_tth']: Par0[var].attrs['units'] = '1'
    for var in ['grt1', 'ga']: Par0[var].attrs['units'] = 'K-1'
    for var in ['grt2']: Par0[var].attrs['units'] = 'K-2'
    for var in ['vthaw', 'vfroz']: Par0[var].attrs['units'] = 'yr-1'
    for var in ['tth_1', 'tth_2', 'tth_3']: Par0[var].attrs['units'] = 'yr'
    for var in ['Cfr0']: Par0[var].attrs['units'] = 'PgC'


    ## uncertainty in varying parameters
    ## derived from the same ensemble of 5 models
    Par0['Cfr0'][1] = 110.97 * np.sqrt((5 - 1) / (5 - 1.5))
    Par0['ka'][1] = 2.108 * np.sqrt((5 - 1) / (5 - 1.5))
    Par0['ga'][1] = 0.03446 * np.sqrt((5 - 1) / (5 - 1.5))
    Par0['k_tth'][1] = 1.059 * np.sqrt((5 - 1) / (5 - 1.5))
    for var in ['ka', 'ga', 'k_tth', 'Cfr0']: Par0[var].attrs['bounds'] = (0., np.inf)


    ##===========
    ## Atmosphere
    ##===========

    ## atmospheric conversion factor
    ## same as recent GCBs (e.g. Le Quere et al., 2018; doi:10.5194/essd-10-405-2018) (Table 1)
    Par0['aCO2'] = 2.124
    Par0['aCO2'].attrs['units'] = 'PgC ppm-1'


    ## preindustrial CO2 concentration
    ## AR6 (bg): (Gulev et al., 2021; doi:10.1017/9781009157896.004) (Section 2.2.3.2.1)
    ## AR6 (unc): same, but taking 0-1850 period for uncertainty range
    ## AR5: (Ciais et al., 2013; doi:10.1017/CBO9781107415324.015)
    if ipcc == 'AR6':
        Par0['CO2pi'][:] = [278.3, 0.5 * (285 - 274) / 1.645]
    elif ipcc == 'AR5':
        Par0['CO2pi'][:] = [278., 5 / 1.645]
    Par0['CO2pi'].attrs['units'] = 'ppm'
    Par0['CO2pi'].attrs['bounds'] = (0., np.inf)


    ## surface ocean acidification scaling factor
    ## assumed known
    Par0['k_pH'] = 1.
    Par0['k_pH'].attrs['units'] = '1'


    ## RETURN
    return Par0


##################################################
## 2. DEFAULT DRIVERS
##################################################

## get default drivers
## (for all versions of the model)
def get_dflt_forcing(ipcc='AR6'):

    ## INITIALISATION
    ## lists of drivers
    forcing_unc = ['Eco2', 'CO2', 'ERFx', 'T']
    forcing_fix = ['d_CO2', 'd_ERFx', 'd_T', 'dd_T']

    ## get last historical year
    TMP = load_gmst()
    TMP2 = load_NOAA_ESRL()
    indH = min(TMP.year.values[-1], TMP2.year.values[-1])

    ## creating final dataset
    For0 = xr.Dataset()
    For0.coords['year'] = np.arange(1750, indH+1)
    For0.coords['stat'] = ['mean', 'std']
    for var in forcing_unc:
        For0[var] = np.nan * xr.zeros_like(For0['year'], dtype=float) + xr.zeros_like(For0['stat'], dtype=float)
    for var in forcing_fix:
        For0[var] = np.nan * xr.zeros_like(For0['year'], dtype=float)


    ##==============
    ## CO2 emissions
    ##==============

    ## best guess: GCB historical data, extrapolated backward to zero for land-use emissions
    ## note: see fct_load for reference
    TMP = load_GCB().sel(data='hist', drop=True)
    TMP.Eluc.loc[:] = np.interp(TMP.year.values, [1750] + list(TMP.Eluc.dropna('year').year.values), [0] + list(TMP.Eluc.dropna('year').values))
    TMP.Eff.loc[1750] = 0.
    For0.Eco2[:len(TMP.Eff), 0] = TMP.Eff + TMP.Eluc

    ## uncertainty: GCB uncertainties, 50% before budget period for land-use emissions
    ## note: see fct_load for reference
    For0.Eco2[:len(TMP.Eff), 1] = np.sqrt((0.05 * TMP.Eff)**2 + (0.5 * TMP.Eluc)**2)
    For0.Eco2[1959-1750:len(TMP.Eff), 1] = np.sqrt((0.05 * TMP.Eff[1959-1750:])**2 + 0.7**2)

    ## units
    For0['Eco2'].attrs['units'] = 'PgC yr-1'


    ##================
    ## Atmospheric CO2
    ##================

    ## best guess: combination of IPCC and latest NOAA/ESRL global measurement estimates
    if ipcc == 'AR6':
        TMP = load_IPCC_AR6_WG1()
    elif ipcc == 'AR5':
        TMP = load_IPCC_AR5_WG1()
    TMP2 = load_NOAA_ESRL().sel(site='global', drop=True).dropna('year')
    For0.CO2[:, 0] = TMP.CO2.combine_first(TMP2.CO2)

    ## uncertainty: ~0.1 ppm over measurement period, extrapolated backward to value assessed by IPCC
    ## AR6: (Gulev et al., 2021; doi:10.1017/9781009157896.004)
    ## AR5: (Ciais et al., 2013; doi:10.1017/CBO9781107415324.015)
    if ipcc == 'AR6':
        For0.CO2[:, 1] = np.interp(For0.year.values, [1750] + list(TMP2.year.values), [2.9 / 1.645] + [0.1 for _ in TMP2.year])
    elif ipcc == 'AR5':
        For0.CO2[:, 1] = np.interp(For0.year.values, [1750] + list(TMP2.year.values), [2.] + [0.1 for _ in TMP2.year])

    ## derivative
    For0.d_CO2[:len(For0.CO2[:,0].dropna('year'))] = np.gradient(For0.CO2[:,0].dropna('year'), edge_order=2)

    ## units
    For0['CO2'].attrs['units'] = 'ppm'
    For0['d_CO2'].attrs['units'] = 'ppm yr-1'


    ##===========
    ## Non-CO2 RF
    ##===========

    ## best guess: IPCC estimates
    ## AR6: no correction needed
    ## AR5: volcanic forcing corrected for average volcano and warming efficacy
    ## (Prather et al., 2013; doi:10.1017/CBO9781107415324.030)
    ## (Gregory et al., 2016; doi:10.1007/s00382-016-3055-1)
    if ipcc == 'AR6':
        TMP = load_IPCC_AR6_WG1()
    elif ipcc == 'AR5':
        TMP = load_IPCC_AR5_WG1()
        TMP.RF.loc[:, 'volc'] -= TMP.RF.loc[:, 'volc'].mean('year')
        TMP.RF.loc[:, 'volc'] *= 0.6
        TMP = TMP.rename({'RF':'ERF'})
    For0.ERFx[:len(TMP.ERF), 0] = TMP.ERF[:, 1:].sum('forcing')

    ## uncertainty: deduced and rounded from IPCC, volcanoes arbitrary
    ## AR6: (Forster et al., 2021; doi:10.1017/9781009157896.009) (Table 7.8)
    ## AR5: (Myhre et al., 2013; doi:10.1017/CBO9781107415324.018) (Tables 8.2, 8.3, 8.6, and text)
    if ipcc == 'AR6':
        unc_ipcc = {'CO2': 0.5*(2.41-1.90)/2.16, 'CH4': 0.5*(0.65-0.43)/0.54, 'N2O': 0.5*(0.24-0.18)/0.21,
            'other_wmghg': 0.5*(0.49-0.33)/0.41, 'O3': 0.5*(0.71-0.24)/0.27, 'H2O_stratospheric': 0.5*(0.10-0.00)/0.05, 
            'aerosol-radiation_interactions': 0.5*(0.04--0.47)/-0.22, 'aerosol-cloud_interactions': 0.5*(-0.25--1.45)/-0.84, 
            'bc_on_snow': 0.5*(0.18-0.00)/0.08, 'land_use': 0.5*(-0.30--0.10)/-0.20, 'contrails': 0.5*(0.10-0.02)/0.06, 
            'solar': 0.5*(0.08--0.06)/0.01, 'volcanic': 0}
    elif ipcc == 'AR5':
        unc_ipcc = {'CO2': 0.1, 'GHG_nonCO2': 0.1, 'O3_trop': 0.20/0.40, 'O3_strat': 0.10/0.05, 
            'AER_tot': np.sqrt(0.5**2 + 0.6**2)/0.9, 'LCC_alb': 0.10/0.15, 'H2O_strat': 0.05/0.07, 'BC_snow': 0.035/0.04, 
            'contrails': 0.025/0.01, 'solar': 0.05/0.05, 'volc': 0.}
    For0.ERFx[:len(TMP.ERF), 1] = np.sqrt(sum([(TMP.ERF.sel(forcing=forc).values * unc_ipcc[forc] / 1.645)**2 for forc in TMP.forcing[1:].values]))

    ## derivative
    For0.d_ERFx[:len(For0.ERFx[:,0].dropna('year'))] = np.gradient(For0.ERFx[:,0].dropna('year'), edge_order=2)

    ## units
    For0['ERFx'].attrs['units'] = 'W m-2'
    For0['d_ERFx'].attrs['units'] = 'W m-2 yr-1'


    ##=====
    ## GMST
    ##=====

    ## best guess: average of all GMST datasets (HadCRUT4: limited coverage; Cowtan_and_Way: not updated to 2021)
    TMP = load_gmst(ref_period=[1850, 1900], ignore_data=['HadCRUT4', 'Cowtan_and_Way'])
    For0.T[1850-1750:, 0] = TMP.T.mean('data')

    ## uncertainty: debiased standard deviation between same datasets
    For0.T[1850-1750:, 1] = TMP.T.std('data') * np.sqrt((TMP.T.notnull().sum('data') - 1.) / (TMP.T.notnull().sum('data') - 1.5))
    For0['T'] = For0.T.fillna(0.)

    ## derivatives
    For0.d_T[:len(For0.T[:,0].dropna('year'))] = np.gradient(For0.T[:,0].dropna('year'), edge_order=2)
    For0.dd_T[:len(For0.d_T.dropna('year'))] = np.gradient(For0.d_T.dropna('year'), edge_order=2)

    ## units
    For0['T'].attrs['units'] = 'K'
    For0['d_T'].attrs['units'] = 'K yr-1'
    For0['dd_T'].attrs['units'] = 'K yr-2'


    ## RETURN
    return For0


##################################################
## 3. DEFAULT CONSTRAINTS
##################################################

## get default constraints
def get_dflt_constr(ipcc='AR6', corr_peat=True, corr_unch_glacier=True):

    ## INITIALISATION
    ## lists of constraints
    constr_all = ['Eco2', 'd_CO2', 'Fland', 'Focean', 'Cv', 'Cs', 'NPP', 'CO2', 'ERFx', 'T', 'd_T', 'd_OHC', 'd_Hthx', 'd_Hgla', 'd_Hgis', 'd_Hais', 'Htot', 'Hlia', 'logit_ff']

    ## load default forcings
    For0 = get_dflt_forcing(ipcc=ipcc)

    ## creating final dataset
    Con0 = xr.Dataset()
    Con0.coords['stat'] = ['mean', 'std']
    for var in constr_all:
        Con0[var] = np.nan * xr.zeros_like(Con0['stat'], dtype=float)
        Con0[var].attrs['is_sum'] = False
        Con0[var].attrs['is_mean'] = False
        Con0[var].attrs['is_diff'] = False


    ##==============
    ## Carbon budget
    ##==============

    ## average fluxes over latest decade
    ## (Friedlingstein et al., 2020; doi:10.5194/essd-12-3269-2020) (Table 6)
    ## note: added FOS and LUC together
    Con0.Eco2[:] = [9.5 + 1.1, np.sqrt(0.5**2 + 0.7**2)]
    Con0.d_CO2[:] = np.array([5.1, 0.02]) / 2.124
    Con0.Eco2.attrs['period'] = Con0.d_CO2.attrs['period'] = (2011, 2020)
    Con0.Eco2.attrs['is_mean'] = Con0.d_CO2.attrs['is_mean'] = True
    
    ## cumulative fluxes
    ## (Friedlingstein et al., 2020; doi:10.5194/essd-12-3269-2020) (Table 8)
    ## note: LASC is excluded here...
    Con0.Focean[:] = [115, 25]
    Con0.Fland[:] = [135, 25]
    Con0.Focean.attrs['period'] = Con0.Fland.attrs['period'] = (1960, 2020)
    Con0.Focean.attrs['is_sum'] = Con0.Fland.attrs['is_sum'] = True
    
    ## units
    Con0.Eco2.attrs['units'] = 'PgC yr-1'
    Con0.d_CO2.attrs['units'] = 'ppm yr-1'
    Con0.Focean.attrs['units'] = 'PgC'
    Con0.Fland.attrs['units'] = 'PgC'


    ##=============
    ## Land C cycle
    ##=============

    ## preindustrial land carbon pools, subtracting peatland for soils
    ## AR6 (bg): (Canadell et al., 2021; doi:10.1017/9781009157896.007) (Figure 5.12)
    ## AR6 (unc): assumed same relative uncertainty as AR5
    ## AR5: (Ciais et al., 2013; doi:10.1017/CBO9781107415324.015) (Figure 6.1)
    if ipcc == 'AR6':
        Con0.Cv[:] = [450, 0.11 * 450]
        Con0.Cs[:] = [1700, 0.14 * 1700]
    elif ipcc == 'AR5':
        Con0.Cv[:] = [550, 100 / 1.645]
        Con0.Cs[:] = [1950, 450 / 1.645]
    Con0.Cv.attrs['period'] = Con0.Cs.attrs['period'] = (1750, 1750)
    Con0.Cv.attrs['is_mean'] = Con0.Cs.attrs['is_mean'] = True
    Con0.Cv.attrs['units'] = Con0.Cs.attrs['units'] = 'PgC'


    ## correction for peatland as unaccounted for in the model
    ## (Yu et al., 2010;  doi:10.1029/2010GL043584) (Table 1)
    if corr_peat:
        Con0.Cs[0] = Con0.Cs[0] - (547 + 50 + 15)
        Con0.Cs[1] = np.sqrt(Con0.Cs[1]**2 + 74**2 + 5.5**2 + 2.5**2)


    ## present-day net primary production
    ## (Ciais et al., 2013; doi:10.1017/CBO9781107415324.015) (text; see also next ref.)
    ## (Zhao et al., 2005; doi:10.1016/j.rse.2004.12.011)
    Con0.NPP[:] = [55, 5]
    Con0.NPP.attrs['period'] = (1998, 2002)
    Con0.NPP.attrs['is_mean'] = True
    Con0.NPP.attrs['units'] = 'PgC yr-1'


    ##================
    ## Atmospheric CO2
    ##================

    ## average concentration from NOAA/ESRL, over latest decade
    ## note: uncertainty assumed 100% correlated
    TMP = For0.CO2.dropna('year')
    Con0.CO2[:] = TMP[-10:].mean('year')
    Con0.CO2.attrs['period'] = (TMP.year.values[-10], TMP.year.values[-1])
    Con0.CO2.attrs['is_mean'] = True
    Con0.CO2.attrs['units'] = TMP.units


    ##===========
    ## Non-CO2 RF
    ##===========

    ## from IPCC latest 10 years
    ## note: uncertainty assumed 100% correlated
    TMP = For0.ERFx.dropna('year')
    Con0.ERFx[:] = TMP[-10:].mean('year')
    Con0.ERFx.attrs['period'] = (TMP.year.values[-10], TMP.year.values[-1])
    Con0.ERFx.attrs['is_mean'] = True
    Con0.ERFx.attrs['units'] = TMP.units


    ##=====
    ## GMST
    ##=====
    
    ## average value: from our GMST datasets, over latest two decades
    TMP = load_gmst().dropna('year')
    TMP2 = TMP.T[-20:].mean('year')
    Con0.T[:] = [np.mean(TMP2), np.std(TMP2) * np.sqrt((len(TMP.data) - 1.) / (len(TMP.data) - 1.5))]
    Con0.T.attrs['period'] = (TMP.year.values[-20], TMP.year.values[-1])
    Con0.T.attrs['is_mean'] = True
    Con0.T.attrs['units'] = 'K'


    ## derivative: same, using second-order accuracy gradient
    TMP2 = [np.gradient(TMP.T.sel(data=data), edge_order=2)[-20-1:-1].mean() for data in TMP.data.values]
    Con0.d_T[:] = [np.mean(TMP2), np.std(TMP2) * np.sqrt((len(TMP.data) - 1.) / (len(TMP.data) - 1.5))]
    Con0.d_T.attrs['period'] = (TMP.year.values[-20-1], TMP.year.values[-1-1])
    Con0.d_T.attrs['is_mean'] = True
    Con0.d_T.attrs['units'] = 'K yr-1'


    ##==================
    ## Ocean heat uptake
    ##==================

    ## from IPCC assessments
    ## AR6: (Gulev et al., 2021; doi:10.1017/9781009157896.004) (Table 2.7)
    ## AR5: (Rhein et al., 2013; doi:10.1017/CBO9781107415324.010) (Box 3.1)
    ## note: converted from ZJ yr-1 to W m-2
    if ipcc == 'AR6':
        Con0.d_OHC[:] = np.array([11.57, 0.5*(15.94-7.20) / 1.645]) * 1E21 / (3600*24*365.25) / 510.1E12
        Con0.d_OHC.attrs['period'] = (2006, 2018)
    elif ipcc == 'AR5':
        Con0.d_OHC[:] = np.array([163, 0.5*(201-127) / 1.645]) / (2010-1993) * 0.93 * 1E21 / (3600*24*365.25) / 510.1E12
        Con0.d_OHC.attrs['period'] = (1993, 2010)
    Con0.d_OHC.attrs['is_mean'] = True
    Con0.d_OHC.attrs['units'] =  'W m-2'


    ##===============
    ## Sea level rise
    ##===============

    ## speed of change, from IPCC assessments
    ## AR6: (Fox-Kemper et al., 2021; doi:10.1017/9781009157896.011) (Table 9.5)
    ## SROCC: (Oppenheimer et al., 2019; doi:10.1017/9781009157964.006) (Table 4.1)
    ## AR5: (Church et al., 2013; doi:10.1017/CBO9781107415324.026) (Table 13.1)
    if ipcc == 'AR6':
        #Con0.d_Htot[:] = [3.69, 0.5*(4.17-3.21) / 1.645]
        Con0.d_Hthx[:] = [1.39, 0.5*(2.05-0.74) / 1.645]
        Con0.d_Hgla[:] = [0.62, 0.5*(0.68-0.57) / 1.645]
        Con0.d_Hgis[:] = [0.63, 0.5*(0.74-0.51) / 1.645]
        Con0.d_Hais[:] = [0.37, 0.5*(0.50-0.24) / 1.645]
        Con0.d_Hthx.attrs['period'] = Con0.d_Hgla.attrs['period'] = Con0.d_Hgis.attrs['period'] = Con0.d_Hais.attrs['period'] = (2006, 2018)
    elif ipcc == 'SROCC':
        #Con0.d_Htot[:] = [3.58, 0.5*(4.06-3.10) / 1.645]
        Con0.d_Hthx[:] = [1.40, 0.5*(1.72-1.08) / 1.645]
        Con0.d_Hgla[:] = [0.61, 0.5*(0.69-0.53) / 1.645]
        Con0.d_Hgis[:] = [0.77, 0.5*(0.82-0.72) / 1.645] 
        Con0.d_Hais[:] = [0.43, 0.5*(0.52-0.34) / 1.645] 
        Con0.d_Hthx.attrs['period'] = Con0.d_Hgla.attrs['period'] = Con0.d_Hgis.attrs['period'] = Con0.d_Hais.attrs['period'] = (2006, 2015)
    elif ipcc == 'AR5':
        #Con0.d_Htot[:] = [3.2, 0.5*(3.6-2.8) / 1.645]
        Con0.d_Hthx[:] = [1.1, 0.5*(1.4-0.8) / 1.645]
        Con0.d_Hgla[:] = [0.76, 0.5*(1.13-0.39) / 1.645]
        Con0.d_Hgis[:] = [0.33, 0.5*(0.25-0.41) / 1.645]
        Con0.d_Hais[:] = [0.27, 0.5*(0.38-0.16) / 1.645] 
        Con0.d_Hthx.attrs['period'] = Con0.d_Hgla.attrs['period'] = Con0.d_Hgis.attrs['period'] = Con0.d_Hais.attrs['period'] = (1993, 2010)
    Con0.d_Hthx.attrs['is_mean'] = Con0.d_Hgla.attrs['is_mean'] = Con0.d_Hgis.attrs['is_mean'] = Con0.d_Hais.attrs['is_mean'] = True
    Con0.d_Hthx.attrs['units'] = Con0.d_Hgla.attrs['units'] = Con0.d_Hgis.attrs['units'] = Con0.d_Hais.attrs['units'] = 'mm yr-1'


    ## sea level, from IPCC assessments (tide gauges only)
    ## idem
    if ipcc == 'AR6':
        Con0.Htot[:] = [120.1, 0.5*(170.8-69.3) / 1.645]
    elif ipcc == 'SROCC':
        Con0.Htot[:] = np.array([1.38, 0.5*(1.95-0.81) / 1.645]) * (1990-1901)
    elif ipcc == 'AR5':
        Con0.Htot[:] = np.array([1.5, 0.5*(1.7-1.3) / 1.645]) * (1990-1901)
    Con0.Htot.attrs['period'] = (1901, 1990)
    Con0.Htot.attrs['is_diff'] = True
    Con0.Htot.attrs['units'] = 'mm'


    ## correction for uncharted glaciers
    ## (Parkes and Marzeion, 2018; doi: 10.1038/s41586-018-0687-9)
    if corr_unch_glacier:
        Con0.Htot[0] = Con0.Htot[0] - 0.5*(0.17+0.53) * (1990-1901)
        Con0.Htot[1] = np.sqrt(Con0.Htot[1]**2 + (0.5*(0.53-0.17) * (1990-1901) / 1.645)**2)


    ## sea level rise from relaxation after LIA
    ## (Slangen et al., 2016; doi:10.1038/nclimate2991)
    Con0.Hlia[:] = [30., 13.]
    Con0.Hlia.attrs['period'] = (1750, 1750)
    Con0.Hlia.attrs['is_mean'] = True
    Con0.Hlia.attrs['units'] = 'mm'


    ##====================
    ## Climate sensitivity
    ##====================

    ## distribution of logit feedback factor
    ## from IPCC AR and Roe & Baker formulation
    Con0.logit_ff[:] = [get_ecs_distrib(data=ipcc)[var] for var in ['mu_ff', 'sigma_ff']]
    Con0.logit_ff.attrs['period'] = (1750, 1750)
    Con0.logit_ff.attrs['is_mean'] = True
    Con0.logit_ff.attrs['units'] = '1'


    ## RETURN
    return Con0

