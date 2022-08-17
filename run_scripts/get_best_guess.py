import numpy as np
import xarray as xr

from core_fct.mod_Edriven import PF_Edriven as PFe
from core_fct.mod_OBSdriven import PF_OBSdriven as PF


##################################################
##################################################

## OPTIONS
name = 'v1'
folder_calib = 'internal_data/pyMC_calib/'
years_csv = slice(2015-2, 2015+2)


##################################################
##################################################

## load historical run from calib
with xr.open_dataset(folder_calib + 'Par_' + name + '.nc') as TMP: Par = TMP.load()
with xr.open_dataset(folder_calib + 'Var_' + name + '.nc') as TMP: Var = TMP.load()

## get years for initial state
year_ini = float(Var.year.sel(year=years_csv).mean('year').values)
year_str = str(int(year_ini)) if year_ini.is_integer() else '{:.1f}'.format(year_ini)

## make best guess parameters and input forcings
Par_bg = Par.mean('config')
For_bg = Var.mean('config').drop([var for var in Var if var not in PF.For_name])

## run Pathfinder to get best guess output
Out_bg = PF.run_xarray(Par_bg, For_bg, get_Var2=True)

## group variables of interest together
Var_bg = xr.merge([For_bg, Out_bg]).drop([var for var in Out_bg if var not in PF.Var_name + PFe.For_name])

## get exact initial state in 1750
for var in PF.For_name: 
    if var in PFe.Var_name:
        Var_bg[var].loc[{'year': int(Var_bg.year[0])}] = PFe.Ini_dflt([np.float32(Par_bg[var].values) for var in PFe.Par_name])[PFe.Var_name.index(var)]
    elif var[:2] == 'd_': # assuming nil derivative
        Var_bg[var].loc[{'year': int(Var_bg.year[0])}] = 0.
for var in PFe.For_name:
    Var_bg[var].loc[{'year': int(Var_bg.year[0])}] = 0.

## get initial state over selected years
Ini_bg = Var_bg.sel(year=years_csv).mean('year')

## save best guess
Par_bg.assign_coords(config='').expand_dims('config', -1).astype(np.float32).to_dataframe().to_csv(folder_calib + '/Par_' + name + '_bg.csv')
Var_bg.astype(np.float32).to_dataframe().to_csv(folder_calib + '/Var_' + name + '_bg.csv')
Ini_bg.assign_coords(config='').expand_dims('config', -1).astype(np.float32).to_dataframe().to_csv(folder_calib + '/Ini_' + year_str + '_' + name + '_bg.csv')

