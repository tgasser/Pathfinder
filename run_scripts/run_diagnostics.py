import numpy as np
import xarray as xr

from core_fct.fct_load import load_SSPs

from core_fct.mod_Cdriven import PF_Cdriven as PF
from core_fct.mod_OBSdriven import PF_OBSdriven as PF_obs
from core_fct.mod_Cdriven_bgc import PF_Cdriven_bgc as PF_bgc
from core_fct.mod_Cdriven_rad import PF_Cdriven_rad as PF_rad


##################################################
##################################################

## OPTIONS
name = 'v1'
folder_calib = 'internal_data/pyMC_calib/'

## length or last year of diagnostic runs
len_2x = 2000
len_1p = 70
yr_ssp = 2100


##################################################
##################################################

## load historical run from calib
with xr.open_dataset(folder_calib + 'Par_' + name + '.nc') as TMP: Par = TMP.load()
with xr.open_dataset(folder_calib + 'Var_' + name + '.nc') as TMP: Var = TMP.load()
Var = xr.merge([Var, PF_obs.get_Var2(Var, Par, Var)])

## create idealized drivers
## time
year = Var.year.interp(year=np.arange(1750, 1750 + max(150, min(len_2x, 10000)) + 1)).year
## 2x CO2
For_2x = xr.Dataset()
For_2x['CO2'] = (1. + (year > 1750)) * Par.CO2pi
For_2x['d_CO2'] = (year == 1750) * Par.CO2pi
For_2x['ERFx'] = 0. * year
## 1% CO2
For_1p = xr.Dataset()
For_1p['CO2'] = 1.01 ** (year.isel(year=slice(max(70, len_1p) + 1)) - 1750) * Par.CO2pi
For_1p['d_CO2'] = np.log(1.01) * For_1p.CO2
For_1p['ERFx'] = 0. * year.isel(year=slice(max(70, len_1p) + 1))

## load SSP drivers from IPCC AR6
For_ssp = load_SSPs().sel(year=slice(max(2100, yr_ssp)))
For_ssp['d_CO2'] = For_ssp.CO2.differentiate('year')
## create initial state for SSPs
Ini_ssp = Var.sel(year=slice(int(For_ssp.year[0])-2, int(For_ssp.year[0])+2)).mean('year')


##################################################
##################################################

## run diagnostic experiments
## note: 'nt' increased for high-CO2 experiments
print('RUN: 2xCO2')
Out_2x = PF.run_xarray(Par=Par, For=For_2x, get_Var2=True)
print('RUN: 1%CO2')
Out_1p = PF.run_xarray(Par=Par, For=For_1p, get_Var2=True, nt=16)
Out_1p_bgc = PF_bgc.run_xarray(Par=Par, For=For_1p, get_Var2=True, nt=16)
Out_1p_rad = PF_rad.run_xarray(Par=Par, For=For_1p, get_Var2=True)
print('RUN: SSPs')
Out_ssp1 = PF.run_xarray(Par=Par, For=For_ssp.sel(scen=[scen for scen in For_ssp.scen.values if scen not in ['ssp585', 'ssp370', 'ssp370-lowNTCF']]), Ini=Ini_ssp, get_Var2=True)
Out_ssp2 = PF.run_xarray(Par=Par, For=For_ssp.sel(scen=['ssp370', 'ssp370-lowNTCF']), Ini=Ini_ssp, get_Var2=True, nt=16)
Out_ssp3 = PF.run_xarray(Par=Par, For=For_ssp.sel(scen=['ssp585']), Ini=Ini_ssp, get_Var2=True, nt=24)
Out_ssp = xr.concat([Out_ssp1, Out_ssp2, Out_ssp3], dim='scen')


##################################################
##################################################

## calculate diagnostics
## climate
ECS = Out_2x.T.isel(year=-1)
TCR = Out_1p.T.sel(year=int(Out_1p.year[0]) + 70)
TCRE = 1E3 * Out_1p.T.sel(year=int(Out_1p.year[0]) + 70) / Out_1p.Eco2.cumsum('year').sel(year=int(Out_1p.year[0]) + 70)
## carbon
beta_ocean = (Out_1p_bgc.Focean.cumsum('year') / (For_1p.CO2 - Par.CO2pi)).sel(year=int(Out_1p.year[0]) + 70)
gamma_ocean = (((Out_1p.Focean - Out_1p_bgc.Focean).cumsum('year')) / Out_1p.T).sel(year=int(Out_1p.year[0]) + 70)
beta_land = ((Out_1p_bgc.Fland - Out_1p_bgc.Epf).cumsum('year') / (For_1p.CO2 - Par.CO2pi)).sel(year=int(Out_1p.year[0]) + 70)
gamma_land = (((Out_1p.Fland - Out_1p.Epf - Out_1p_bgc.Fland + Out_1p_bgc.Epf).cumsum('year')) / Out_1p.T).sel(year=int(Out_1p.year[0]) + 70)

## key scenario values
## climate
T_2030 = (Out_ssp.T.sel(year=slice(2021, 2040)).mean('year') - Var.T.sel(year=slice(1850, 1900)).mean('year'))
T_2050 = (Out_ssp.T.sel(year=slice(2041, 2060)).mean('year') - Var.T.sel(year=slice(1850, 1900)).mean('year'))
T_2090 = (Out_ssp.T.sel(year=slice(2081, 2100)).mean('year') - Var.T.sel(year=slice(1850, 1900)).mean('year'))
## carbon
C_ocean = Out_ssp.Focean.sel(year=slice(2015, 2100)).sum('year')
C_land = Out_ssp.Fland.sel(year=slice(2015, 2100)).sum('year')
## slr
H_2050 = (Out_ssp.Htot.sel(year=2050) - Var.Htot.sel(year=slice(1995, 2014)).mean('year')) / 1000
H_2100 = (Out_ssp.Htot.sel(year=2100) - Var.Htot.sel(year=slice(1995, 2014)).mean('year')) / 1000
dH_2050 = Out_ssp.d_Htot.sel(year=slice(2041, 2060)).mean('year')
dH_2090 = Out_ssp.d_Htot.sel(year=slice(2081, 2100)).mean('year')


##################################################
##################################################

## scenarios to compare to
scen_ref = ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']

## quick display functions
print_sym0 = lambda da: '{:.0f} +- {:.0f}'.format(da.mean().values, da.std().values)
print_sym1 = lambda da: '{:.1f} +- {:.1f}'.format(da.mean().values, da.std().values)
print_sym2 = lambda da: '{:.2f} +- {:.2f}'.format(da.mean().values, da.std().values)
print_asym0 = lambda da: '{:.0f} ({:.0f}, {:.0f})'.format(da.median().values, da.quantile(0.17).values, da.quantile(0.83).values)
print_asym1 = lambda da: '{:.1f} ({:.1f}, {:.1f})'.format(da.median().values, da.quantile(0.17).values, da.quantile(0.83).values)
print_asym2 = lambda da: '{:.2f} ({:.2f}, {:.2f})'.format(da.median().values, da.quantile(0.17).values, da.quantile(0.83).values)

## ECS & TCR
## CMIP6: (Meehl et al., 2020; doi:10.1126/sciadv.aba1981) (Table 2)
## CMIP6 data: mean +- std
## AR6: (Forster et al., 2021; doi:10.1017/9781009157896.009) (Executive Summary)
## AR6 data: median (likely range)
print(6*20*'-')
print((6*'{:>20}').format('', 'Pathfinder', 'CMIP6', '|', 'Pathfinder', 'AR6'))
print((6*'{:>20}').format('ECS', print_sym1(ECS), '3.7 +- 1.0', '|', print_asym1(ECS), '3.0 (2.0, 4.5)'))
print((6*'{:>20}').format('TCR', print_sym1(TCR), '2.0 +- 0.4', '|', print_asym1(TCR), '1.8 (1.4, 2.2)'))

## TCRE
## CMIP6: (Arora et al., 2020; doi:10.5194/bg-17-4173-2020) (Table 7)
## CMIP6 data: mean +- std
## AR6: (Canadell et al., 2021; doi:10.1017/9781009157896.007) (Table 5.7)
## AR6 data: median (likely range)
print((6*'{:>20}').format('TCRE', print_sym2(TCRE), '1.77 +- 0.37', '|', print_asym2(TCRE), '1.65 (1.0, 2.3)'))
print(6*20*'-')

## beta/gamma
## (Arora et al., 2020; doi:10.5194/bg-17-4173-2020) (Table A1)
## data:  mean +- std
print((3*'{:>20}').format('', 'Pathfinder', 'CMIP6'))
print((3*'{:>20}').format('beta_ocean', print_sym2(beta_ocean), '0.91 +- 0.09'))
print((3*'{:>20}').format('gamma_ocean', print_sym1(gamma_ocean), '-8.6 +- 2.9'))
print((3*'{:>20}').format('beta_land', print_sym2(beta_land), '1.22 +- 0.40'))
print((3*'{:>20}').format('gamma_land', print_sym1(gamma_land), '-34.1 -+ 38.4'))
print(6*20*'-')

## climate
## (Lee et al., 2021; doi:10.1017/9781009157896.006) (Table 4.5)
## data: central (very likely range)
print(6*20*'-')
print((6*'{:>20}').format('', *scen_ref))
print((6*'{:>20}').format('T_2030 Pathfinder', *[print_asym1(T_2030.sel(scen=scen)) for scen in scen_ref]))
print((6*'{:>20}').format('T_2030 AR6', *['1.5 (1.2, 1.7)', '1.5 (1.2, 1.8)', '1.5 (1.2, 1.8)', '1.5 (1.2, 1.8)', '1.6 (1.3, 1.9)']))
print(6*20*'-')
print((6*'{:>20}').format('T_2050 Pathfinder', *[print_asym1(T_2050.sel(scen=scen)) for scen in scen_ref]))
print((6*'{:>20}').format('T_2050 AR6', *['1.6 (1.2, 2.0)', '1.7 (1.3, 2.2)', '2.0 (1.6, 2.5)', '2.1 (1.7, 2.6)', '2.4 (1.9, 3.0)']))
print(6*20*'-')
print((6*'{:>20}').format('T_2090 Pathfinder', *[print_asym1(T_2090.sel(scen=scen)) for scen in scen_ref]))
print((6*'{:>20}').format('T_2090 AR6', *['1.4 (1.0, 1.8)', '1.8 (1.3, 2.4)', '2.7 (2.1, 3.5)', '3.6 (2.8, 4.6)', '4.4 (3.3, 5.7)']))

## carbon
## (Liddicoat et al., 2021; doi:10.1175/JCLI-D-19-0991.1) (Tables S5 & S6)
## note: removed LUC emissions from ScenarioMIP database without uncertainty (CE_luc = -13, -27, -14, 73, 34)
## data: mean +- std
print(6*20*'-')
print((6*'{:>20}').format('', *scen_ref))
print((6*'{:>20}').format('C_ocean Pathfinder', *[print_sym0(C_ocean.sel(scen=scen)) for scen in scen_ref]))
print((6*'{:>20}').format('C_ocean CMIP6/AR6', *['111 +- 11', '162 +- 8', '252 +- 11', '338 +- 15', '398 +- 17']))
print(6*20*'-')
print((6*'{:>20}').format('C_land Pathfinder', *[print_sym0(C_land.sel(scen=scen)) for scen in scen_ref]))
print((6*'{:>20}').format('C_land CMIP6/AR6', *['73 +- 33', '120 +- 50', '178 +- 76', '269 +- 124', '311 +- 162']))

## slr
## (Fox-Kemper et al., 2021; doi:10.1017/9781009157896.011) (Table 9.9)
## data: median (likely range)
print(6*20*'-')
print((6*'{:>20}').format('', *scen_ref))
print((6*'{:>20}').format('H_2050 Pathfinder', *[print_asym2(H_2050.sel(scen=scen)) for scen in scen_ref]))
print((6*'{:>20}').format('H_2050 AR6', *['0.18 (0.15, 0.23)', '0.19 (0.16, 0.25)', '0.20 (0.17, 0.26)', '0.22 (0.18, 0.27)', '0.23 (0.20, 0.29)']))
print(6*20*'-')
print((6*'{:>20}').format('H_2100 Pathfinder', *[print_asym2(H_2100.sel(scen=scen)) for scen in scen_ref]))
print((6*'{:>20}').format('H_2100 AR6', *['0.38 (0.28, 0.55)', '0.44 (0.32, 0.62)', '0.56 (0.44, 0.76)', '0.68 (0.55, 0.90)', '0.77 (0.63, 1.01)']))
print(6*20*'-')
print((6*'{:>20}').format('dH_2050 Pathfinder', *[print_asym1(dH_2050.sel(scen=scen)) for scen in scen_ref]))
print((6*'{:>20}').format('dH_2050 AR6', *['4.1 (2.8, 6.0)', '4.8 (3.5, 6.8)', '5.8 (4.4, 8.0)', '6.4 (5.0, 8.7)', '7.2 (5.6, 9.7)']))
print(6*20*'-')
print((6*'{:>20}').format('dH_2090 Pathfinder', *[print_asym1(dH_2090.sel(scen=scen)) for scen in scen_ref]))
print((6*'{:>20}').format('dH_2090 AR6', *['4.2 (2.4, 6.6)', '5.2 (3.2, 8.0)', '7.7 (5.2, 11.6)', '10.4 (7.4, 14.8)', '12.1 (8.6, 17.6)']))
print(6*20*'-')

