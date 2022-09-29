# Pathfinder

*Manual*

## Run a simulation

Basic instructions to run a simulation with `xarray` from the main folder:

1. open python in a terminal, import numpy and xarray: 
```python
import numpy as np
import xarray as xr
```

2. import the model from one of the `mod_` files in `core_fct`, for instance the emission-driven mode:
```python
from core_fct.mod_Edriven import PF_Edriven as PF
```

3. create or load parameters (`Par`), for instance the saved posterior parameters:
```python
Par = xr.open_dataset('internal_data/pyMC_calib/Par_v1.nc').load()
```

4. create or load an initial state (`Ini`), for instance the end of the posterior historical simulation averaged over 11 years: (*Note: this is optional, without a specific initial state, the model will assume it starts from the end of the preindustrial era in 1750.*) 
```python
Ini = xr.open_dataset('internal_data/pyMC_calib/Var_v1.nc')
yr_ini = int(Ini.year.isel(year=slice(-11, None)).mean('year'))
Ini = Ini.isel(year=slice(-11, None)).mean('year')
```

5. create or load some forcings (`For`), for instance one experiment that keeps forcings constant and another that relaxes them:
```python
For = xr.merge([Ini[var] for var in ['Eco2', 'ERFx']])
years = xr.DataArray(yr_ini + np.arange(100), coords={'year': yr_ini + np.arange(100)}, dims=['year'])
For = 0*years + xr.DataArray([1., 0.], coords={'exp': ['const', 'relax']}, dims=['exp']) * For
```

6. run the model by calling the `run_xarray` method (turn on `get_Var2` to save more output variables):
```python
Out = PF.run_xarray(Par, For, Ini, get_Var2=True)
```

7. display results with whatever package you prefer, for instance matplotlib:
```python
import matplotlib.pyplot as plt
plt.figure()
for n, var in enumerate(['T', 'CO2']):
    plt.subplot(1, 2, n+1)
    for exp in Out.exp.values:
        plt.plot(Out.year, Out[var].sel(exp=exp).mean('config'), label=exp)
    plt.legend(loc=0)
plt.show()
```

Check the `run_scripts/run_diagnostics.py` file for more complex examples!


<!--------->
<!--------->
## Code structure

The main folder contains the following subfolders:
| Folder | Content |
| --- | --- |
| `core_fct`| core classes (`cls_`), functions (`fct_`) and objects (`mod_`) that constitute Pathfinder |
| `input_data` | input data used either for calibration or to run Pathfinder |
| `internal_data` | internally generated data, notably during calibration |
| `results` | output data, normally empty in the open-source version |
| `run_scripts` | scripts used to execute calibration and diagnostic (and possibly more) |

The `core_fct` subfolder contains the following files:
| File | Content |
| --- | --- |
| `cls_calib`| `pymc3` subclasses for use during calibration |
| `cls_model`| specific class used to wrap models and solve them |
| `fct_ancillary` | compilation of ancillary functions |
| `fct_calib` | main functions used for the Bayesian calibration |
| `fct_default` | functions setting the default (i.e. prior) values of parameters, forcings and constraints |
| `fct_load` | functions loading and formatting `input_data` |
| `fct_param` | functions doing the OLS calibration of the prior parameters |
| `fct_traj` | functions to generate future trajectories of `T` and `CO2` |
| `mod_Cdriven_bgc` | Pathfinder model in concentration-driven mode with only the biogeochemical effect of CO2 |
| `mod_Cdriven_rad` | Pathfinder model in concentration-driven mode with only the radiative effect of CO2 |
| `mod_Cdriven` | Pathfinder model in concentration-driven mode |
| `mod_Edriven` | Pathfinder model in emission-driven mode |
| `mod_OBSdriven` | Pathfinder model in observation-driven mode (used for calibration) |
| `mod_Tdriven` | Pathfinder model in temperature-driven mode |

The `input_data` subfolder contains the following data folders: 
(*Note that data was kept as close to the original format as possible. Sources are detailed in the core functions that call the files.*)
| Folder | Content |
| --- | --- |
| `AR5`| Data from the IPCC 5th report |
| `AR6`| Data from the IPCC 6th report |
| `CMIP5`| Outputs of the CMIP5 models compiled for Pathfinder |
| `CMIP6`| Outputs of the CMIP6 models compiled for Pathfinder |
| `Edwards_2021`| Data from Edwards et al. (2021) for the SLR module |
| `GCB`| Data from the Global Carbon Budget |
| `obs_CO2`| Observations of atmospheric CO2 |
| `obs_T`| Observations of global mean surface temperature |
| `RCPs`| RCP scenarios |
| `SSPs`| SSP scenarios |
| `TRENDYv7`| Outputs of the TRENDYv7 models compiled for Pathfinder |

The `internal_data` subfolder contains the following folders:
| Folder | Content |
| --- | --- |
| `prior_calib`| Results of the OLS calibration for the prior parameters (for faster loading, can be re-calculated on the fly) |
| `pyMC_calib`| Results of the Bayesian calibration, including historical simulations |

The `run_scripts` subfolder contains the following files:
| File | Content |
| --- | --- |
| `get_best_guess`| script to extract best-guess configuration after calibration |
| `plot_calib_check`| script to display results of the Bayesian calibration |
| `run_calib_and_hist`| script to run the Bayesian calibration |
| `run_diagnostics`| script to run idealized and scenario experiments for model evaluation |


<!--------->
<!--------->
## Notations

<!----------------------->
### Forcings and Variables 

| In manual| In code | Description | Units | Prog? | Dims |
| --- | --- | --- | --- | --- | --- |
| $$ R_c $$ | `RFco2` | CO2 (effective) radiative forcing | W m<sup>-2</sup> |||
| $$ R_x $$ | `ERFx` | Non-CO2 effective radiative forcing | W m<sup>-2</sup>  |||
| $$ R $$ | `ERF` | Effective radiative forcing | W m<sup>-2</sup> |||
| $$ T $$ | `T` | Global surface temperature anomaly | K | *yes* ||
| $$ T_d $$ | `Td` | Deep ocean temperature anomaly | K | *yes* ||
| $$ \mathrm{logit} (\mathrm{ff}) $$ | `logit_ff` | Logit of the climate feedback factor (for calib.) | 1 |||
|||||||
| $$ U_\mathrm{ohc} $$ | `OHC` | Ocean heat content (anomaly) | W yr m<sup>-2</sup> |||
| $$ H_\mathrm{thx} $$ | `Hthx` | Thermosteric  sea level rise | mm |||
| $$ H_\mathrm{gla} $$ | `Hgla` | Glaciers' contribution to sea level rise | mm | *yes* ||
| $$ H_\mathrm{gis} $$ | `Hgis` | Grenland ice sheet's contribution to sea level rise | mm | *yes* ||
| $$ H_\mathrm{ais,smb} $$ | `Hais_smb` | Surface mass balance component of `Hais` | mm |||
| $$ H_\mathrm{ais} $$ | `Hais` | Antartica ice sheet's contribution to sea level rise | mm | *yes* ||
| $$ H_\mathrm{tot} $$ | `Htot` | Total sea level rise | mm |||
| $$ H_\mathrm{lia} $$ | `Hlia` | Sea level rise from relaxation after LIA between 1900 and 2005 (for calib.) | mm |||
|||||||
| $$ C_{o,j} $$ | `Co_j` | Change in surface ocean carbon subpools | PgC | *yes* | $$ j\in[\![1,5]\!] $$ |
| $$ C_o $$ | `Co` | Change in surface ocean carbon pool | PgC |||
| $$ C_d $$ | `Cd` | Change in deep ocean carbon pool | *yes* ||
| $$ c_\mathrm{dic} $$ | `dic` | Change in surface DIC | µmol kg<sup>-1</sup> |||
| $$ p_\mathrm{dic} $$ | `pdic` | Subcomponent of `pCO2` | ppm |||
| $$ p_\mathrm{CO2} $$ | `pCO2` | CO2 partial pressure at the ocean surface | ppm |||
| $$ F_\mathrm{ocean} $$ | `Focean` | Ocean carbon sink | PgC yr<sup>-1</sup> |||
|||||||
| $$ r_\mathrm{npp} $$ | `r_npp` | Relative change in NPP | 1 |||
| $$ r_\mathrm{fire} $$ | `r_fire` | Relative change in wildfire intensity | 1 |||
| $$ r_\mathrm{rh} $$ | `r_rh` | Relative change in heterotrophic respiration rate | 1 |||
| $$ F_\mathrm{npp} $$ | `NPP` | Net primary productivity | PgC yr<sup>-1</sup> |||
| $$ E_\mathrm{fire} $$ | `Efire` | Emissions from wildfire | PgC yr<sup>-1</sup> |||
| $$ E_\mathrm{harv} $$ | `Eharv` | Emissions from harvest and grazing | PgC yr<sup>-1</sup> |||
| $$ F_\mathrm{mort} $$ | `Fmort` | Mortality flux | PgC yr<sup>-1</sup> |||
| $$ E_\mathrm{rh1} $$ | `RH1` | Litter heterotrophic respiration | PgC yr<sup>-1</sup> |||
| $$ F_\mathrm{stab} $$ | `Fstab` | Stabilization flux | PgC yr<sup>-1</sup> |||
| $$ E_\mathrm{rh2} $$ | `RH2` | Active soil heterotrophic respiration | PgC yr<sup>-1</sup> |||
| $$ F_\mathrm{pass} $$ | `Fpass` | Passivization flux | PgC yr<sup>-1</sup> |||
| $$ E_\mathrm{rh3} $$ | `RH3` | Passive soil heterotrophic respiration | PgC yr<sup>-1</sup> |||
| $$ F_\mathrm{land} $$ | `Fland` | Land carbon sink | PgC yr<sup>-1</sup> |||
| $$ E_\mathrm{rh} $$ | `RH` | Heterotrophic respiration | PgC yr<sup>-1</sup> |||
| $$ C_v $$ | `Cv` | Vegetation carbon pool | PgC | *yes* ||
| $$ C_{s1} $$ | `Cs1` | Litter carbon pool | PgC | *yes* ||
| $$ C_{s2} $$ | `Cs2` | Active soil carbon pool | PgC | *yes* ||
| $$ C_{s3} $$ | `Cs3` | Passive soil carbon pool | PgC | *yes* ||
| $$ C_s $$ | `Cs` | Total soil carbon pool | PgC |||
|||||||
| $$ r_\mathrm{rt} $$ | `r_rt` | Relative change in permafrost respiration rate | 1 |||
| $$ \bar{a} $$ | `abar` | Theoretical thawed fraction | 1 |||
| $$ a $$ | `a` | Actual thawed fraction | 1 | *yes* ||
| $$ E_\mathrm{pf} $$ | `Epf` | Emissions from permafrost | PgC yr<sup>-1</sup> |||
| $$ C_{\mathrm{th},j} $$ | `Cth_j` | Thawed permafrost carbon subpools | PgC | *yes* | $$ j\in[\![1,3]\!] $$ |
| $$ C_\mathrm{fr} $$ | `Cfr` | Frozen permafrost carbon pool | PgC |||
|||||||
| $$ E_\mathrm{CO2} $$ | `Eco2` | Anthropogenic CO2 emissions | PgC yr<sup>-1</sup> |||
| $$ C $$ | `CO2` | Atmospheric CO2 concentration | ppm | *yes* ||
| $$ \mathrm{pH} $$ | `pH` | Surface ocean pH | 1 |||

<!----------->
### Parameters 

| In manual | In code | Description | Units | Dims |
| --- | --- | --- | --- | --- |
| $$ \phi $$ | `phi` | Radiative parameter of CO2 | W m<sup>-2</sup> ||
| $$ T_\mathrm{2\times} $$ | `T2x` | Equilibrium climate sensitivity | K ||
| $$ \Theta_s $$ | `THs` | Heat capacity of the surface | W yr m<sup>-2</sup> K<sup>-1</sup> ||
| $$ \Theta_d $$ | `THd` | Heat capacity of the deep ocean | W yr m<sup>-2</sup> K<sup>-1</sup> ||
| $$ \theta $$ | `th` | Heat exchange coefficient | W m<sup>-2</sup> K<sup>-1</sup> ||
| $$ \epsilon_\mathrm{heat} $$ | `eheat` | Deep ocean heat uptake efficacy | 1 ||
| $$ T_\mathrm{2\times}^* $$ | `T2x0` | Minimal value of the ECS distribution (for calib.) | K ||
||||||
| $$ \alpha_\mathrm{ohc} $$ | `aOHC` | Fraction of energy warming the ocean | 1 ||
| $$ \Lambda_\mathrm{thx} $$ | `Lthx` | Proportionality factor of thermosteric SLR | mm m<sup>2</sup> W<sup>-1</sup> yr<sup>-1</sup> ||
| $$ \lambda_\mathrm{gla} $$ | `lgla0` | Initial imbalance in SLR from Glaciers | mm yr<sup>-1</sup> ||
| $$ \Lambda_\mathrm{gla} $$ | `Lgla` | Maximum contribution to SLR from Glaciers | mm ||
| $$ \Gamma_\mathrm{gla1} $$ | `Ggla1` | Linear sensitivity of steady-state Glaciers SLR to climate | K<sup>-1</sup> ||
| $$ \Gamma_\mathrm{gla3} $$ | `Ggla3` | Cubic sensitivity of steady-state Glaciers SLR to climate | K<sup>-3</sup> ||
| $$ \tau_\mathrm{gla} $$ | `tgla` | Timescale of Glaciers' contribution to SLR | yr ||
| $$ \gamma_\mathrm{gla} $$ | `ggla` | Sensitivity of Glaciers' timescale to climate | K<sup>-1</sup> ||
| $$ \lambda_\mathrm{gis} $$ | `lgis0` | Initial imbalance in SLR from GIS | mm yr<sup>-1</sup> ||
| $$ \Lambda_\mathrm{gis1} $$ | `Lgis1` | Linear sensitivity of steady-state GIS SLR to climate| mm K<sup>-1</sup> ||
| $$ \Lambda_\mathrm{gis3} $$ | `Lgis3` | Cubic sensitivity of steady-state GIS SLR to climate | mm K<sup>-3</sup> ||
| $$ \tau_\mathrm{gis} $$ | `tgis` | Timescale of GIS contribution to SLR | yr ||
| $$ \Lambda_\mathrm{ais,smb} $$ | `Lais_smb` | Sensitivity of AIS SMB increase due to climate | mm yr<sup>-1</sup> K<sup>-1</sup> ||
| $$ \lambda_\mathrm{ais} $$ | `lais` | Initial imbalance in SLR from AIS | mm yr<sup>-1</sup> ||
| $$ \Lambda_\mathrm{ais} $$ | `Lais` | Sensitivity of steady-state AIS SLR to climate | mm K<sup>-1</sup> ||
| $$ \tau_\mathrm{ais} $$ | `tais` | Timescale of AIS contribution to SLR | yr ||
| $$ \alpha_\mathrm{ais} $$ | `aais` | Sensitivity of AIS timescale to AIS SLR | mm<sup>-1</sup> ||
||||||
| $$ \alpha_\mathrm{dic} $$ | `adic` | Conversion factor for DIC | µmol kg<sup>-1</sup> PgC<sup>-1</sup> ||
| $$ \beta_\mathrm{dic} $$ | `bdic` | Inverse-scaling factor for DIC | 1 ||
| $$ \gamma_\mathrm{dic} $$ | `gdic` | Sensitivity of pCO2 to climate | K<sup>-1</sup> ||
| $$ T_o $$ | `To` | Preindustrial surface ocean temperature | °C ||
| $$ \nu_\mathrm{gx} $$ | `vgx` | Surface ocean gas exchange rate | yr<sup>-1</sup> ||
| $$ \gamma_\mathrm{gx} $$ | `ggx` | Sensitivity of gas exchange to climate | K<sup>-1</sup> ||
| $$ \alpha_{o,j} $$ | `aoc_j` | Surface ocean subpools fractions | 1 | $$ j\in[\![1,5]\!] $$ |
| $$ \tau_{o,j} $$ | `toc_j` | Timescales of surface ocean subpools | yr | $$ j\in[\![1,5]\!] $$ |
| $$ \kappa_{\tau_o} $$ | `k_toc` | Scaling factor for timescales of surface ocean subpools | 1 ||
||||||
| $$ \beta_\mathrm{npp} $$ | `bnpp` | Sensitivity of NPP to CO2 (= fertilization effect) | 1 ||
| $$ \alpha_\mathrm{npp} $$ | `anpp` | Shape parameter for fertilization effect | 1 ||
| $$ \gamma_\mathrm{npp} $$ | `gnpp` | Sensitivity of NPP to climate | K<sup>-1</sup> ||
| $$ \beta_\mathrm{fire} $$ | `bfire` | Sensitivity of wildfire intensity to CO2 | 1 ||
| $$ \gamma_\mathrm{fire} $$ | `gfire` | Sensitivity of wildfire intensity to climate | K<sup>-1</sup> ||
| $$ \beta_\mathrm{rh} $$ | `brh` | Sensitivity of heterotrophic respiration to fresh organic matter | 1 ||
| $$ \gamma_\mathrm{rh} $$ | `grh` | Sensitivity of heterotrophic respiration to climate | K<sup>-1</sup> ||
| $$ F_\mathrm{npp,0} $$ | `npp0` | Preindustrial NPP | PgC yr<sup>-1</sup> ||
| $$ \nu_\mathrm{fire} $$ | `vfire` | Wildfire intensity | yr<sup>-1</sup> ||
| $$ \nu_\mathrm{harv} $$ | `vharv` | Harvest and grazing rate | yr<sup>-1</sup> ||
| $$ \nu_\mathrm{mort} $$ | `vmort` | Mortality rate | yr<sup>-1</sup> ||
| $$ \nu_\mathrm{stab} $$ | `vstab` | Stabilization rate | yr<sup>-1</sup> ||
| $$ \nu_\mathrm{rh1} $$ | `vrh1` | Litter heterotrophic respiration rate | yr<sup>-1</sup> ||
| $$ \nu_\mathrm{rh23} $$ | `vrh23` | Soil (active and passive) respiration rate | yr<sup>-1</sup> ||
| $$ \nu_\mathrm{rh3} $$ | `vrh3` | Passive soil respiration rate | yr<sup>-1</sup> ||
| $$ \alpha_\mathrm{pass} $$ | `apass` | Fraction of passive soil | 1 ||
||||||
| $$ \alpha_\mathrm{lst} $$ | `aLST` | Climate scaling factor over permafrost regions | 1 ||
| $$ \gamma_\mathrm{rt1} $$ | `grt1` | Sensitivity of (boreal) heterotrophic respiration to climate | K<sup>-1</sup> ||
| $$ \gamma_\mathrm{rt2} $$ | `grt2` | Sensitivity of (boreal) heterotrophic respiration to climate (quadratic) | K<sup>-2</sup> ||
| $$ \kappa_\mathrm{rt} $$ | `krt` | Scaling factor for sensitivity of permafrost respiration to climate | 1 ||
| $$ a_\mathrm{min} $$ | `amin` | Minimal thawed fraction | 1 ||
| $$ \kappa_a $$ | `ka` | Shape parameter for theoretical thawed fraction | 1 ||
| $$ \gamma_a $$ | `ga` | Sensitivity of theoretical thawed fraction to climate | K<sup>-1</sup> ||
| $$ \nu_\mathrm{thaw} $$ | `vthaw` | Thawing rate | yr<sup>-1</sup> ||
| $$ \nu_\mathrm{froz} $$ | `vfroz` | Freezing rate | yr<sup>-1</sup> ||
| $$ \alpha_{\mathrm{th},j} $$ | `ath_j` | Thawed permafrost carbon subpools fractions | 1 | $$ j\in[\![1,3]\!] $$ |
| $$ \tau_{\mathrm{th},j} $$ | `tth_j` | Timescales of thawed permafrost carbon subpools | yr | $$ j\in[\![1,3]\!] $$ |
| $$ \kappa_{\tau_\mathrm{th}} $$ | `k_tth` | Scaling factor for timescales of surface ocean subpools | 1 ||
| $$ C_\mathrm{fr,0} $$ | `Cfr0` | Preindustrial frozen permafrost carbon pool | PgC ||
||||||
| $$ \alpha_C $$ | `aCO2` | Conversion factor for atmospheric CO2 | PgC ppm<sup>-1</sup> ||
| $$ C_\mathrm{pi} $$ | `CO2pi` | Preindustrial CO2 concentration | ppm ||
| $$ \kappa_\mathrm{pH} $$ | `k_pH` | Scaling factor for surface ocean pH | 1 ||
||||||
| $$ \tilde{\sigma}_C $$ | `std_CO2` | Relative standard deviation of the historical `CO2` time series (for calib.) | 1 ||
| $$ \epsilon_C $$ | `ampl_CO2` | Noise amplitude of the historical `CO2` time series (for calib.) | ppm ||
| $$ \rho_C $$ | `corr_CO2` | Autocorrelation of the historical `CO2` time series (for calib.) | 1 ||
| $$ \tilde{\sigma}_T $$ | `std_T` | Relative standard deviation of the historical `T` time series (for calib.) | 1 ||
| $$ \epsilon_T $$ | `ampl_T` | Noise amplitude of the historical `T` time series (for calib.) | K ||
| $$ \rho_T $$ | `corr_T` | Autocorrelation of the historical `T` time series (for calib.) | 1 ||

<!--------->
<!--------->
## Equations

<!----------->
### 1. Climate

###### diagnostic

* $$ R_c = \phi \: \ln \! \left( \frac{C}{C_\mathrm{pi}} \right) $$

* $$ R = R_c + R_x $$

###### prognostic

* $$ \Theta_s \: \frac{\mathrm{d} T}{\mathrm{d} t} = R - \frac{\phi \: \ln (2)}{T_{2\times}} \: T - \epsilon_\mathrm{heat} \: \theta \: (T - T_d) $$

* $$ \Theta_d \: \frac{\mathrm{d} T_d}{\mathrm{d} t} = \theta \: (T - T_d) $$

###### diagnostic (2nd; for calib.)

* $$ \mathrm{logit} (\mathrm{ff}) = \ln \! \left( \frac{T_{2\times}}{T_{2\times}^*} - 1 \right) $$

<!------------->
### 2. Sea level

###### diagnostic

* $$ U_\mathrm{ohc} = \alpha_\mathrm{ohc} \: (\Theta_s \: T + \Theta_d \: T_d) $$

* $$ \frac{\mathrm{d} U_\mathrm{ohc}}{\mathrm{d} t} = \alpha_\mathrm{ohc} \: (\Theta_s \: \frac{\mathrm{d} T}{\mathrm{d} t} + \Theta_d \: \frac{\mathrm{d} T_d}{\mathrm{d} t}) $$

* $$ H_\mathrm{thx} = \Lambda_\mathrm{thx} \: U_\mathrm{ohc} $$

* $$ \frac{\mathrm{d} H_\mathrm{thx}}{\mathrm{d} t} = \Lambda_\mathrm{thx} \: \frac{\mathrm{d} U_\mathrm{ohc}}{\mathrm{d} t} $$

###### prognostic

* $$ \frac{\mathrm{d} H_\mathrm{gla}}{\mathrm{d} t} = \lambda_\mathrm{gla} + \frac{\exp (\gamma_\mathrm{gla} \: T)}{\tau_\mathrm{gla}} \: (\Lambda_\mathrm{gla} \: (1 - \exp (-\Gamma_\mathrm{gla1} \: T - \Gamma_\mathrm{gla3} \: T^3)) - H_\mathrm{gla}) $$

* $$ \frac{\mathrm{d} H_\mathrm{gis}}{\mathrm{d} t} = \lambda_\mathrm{gis} + \frac{1}{\tau_\mathrm{gis}} \: (\Lambda_\mathrm{gis1} \: T + \Lambda_\mathrm{gis3} \: T^3 - H_\mathrm{gis}) $$

* $$ \frac{\mathrm{d} H_\mathrm{ais,smb}}{\mathrm{d} t} = -\Lambda_\mathrm{ais_smb} \: T $$

* $$ \frac{\mathrm{d} H_\mathrm{ais}}{\mathrm{d} t} = \frac{\mathrm{d} H_\mathrm{ais,smb}}{\mathrm{d} t} + \lambda_\mathrm{ais} + \frac{1 + \alpha_\mathrm{ais} \: (H_\mathrm{ais} - H_\mathrm{ais,smb})}{\tau_\mathrm{ais}} \: (\Lambda_\mathrm{ais} \: T - (H_\mathrm{ais} - H_\mathrm{ais,smb})) $$

###### diagnostic (2nd)

* $$ H_\mathrm{tot} = H_\mathrm{thx} + H_\mathrm{gla} + H_\mathrm{gis} + H_\mathrm{ais} $$

* $$ \frac{\mathrm{d} H_\mathrm{tot}}{\mathrm{d} t} = \frac{\mathrm{d} H_\mathrm{thx}}{\mathrm{d} t} + \frac{\mathrm{d} H_\mathrm{gla}}{\mathrm{d} t} + \frac{\mathrm{d} H_\mathrm{gis}}{\mathrm{d} t} + \frac{\mathrm{d} H_\mathrm{ais}}{\mathrm{d} t} $$

###### diagnostic (3rd; for calib.)

* $$ H_\mathrm{lia} = \sum_{\mathrm{ice} \in \{\mathrm{gla}, \mathrm{gis}, \mathrm{ais}\}} \lambda_\mathrm{ice} \: \tau_\mathrm{ice} \: (\exp (-150 / \tau_\mathrm{ice}) - \exp (-205 / \tau_\mathrm{ice})) $$

<!---------------->
### 3. Ocean carbon

###### diagnostic

* $$ C_o = \sum_j C_{o,j} $$

* $$ c_\mathrm{dic} = \frac{\alpha_\mathrm{dic}}{\beta_\mathrm{dic}} \: C_o $$

* $$ \begin{aligned} p_\mathrm{dic} & = (1.5568 - 0.013993 \: T_o) \: c_\mathrm{dic} \\ & + (7.4706 - 0.20207 \: T_o) \: 10^{-3} \: {c_\mathrm{dic}}^2 \\ & - (1.2748 - 0.12015 \: T_o) \: 10^{-5} \: {c_\mathrm{dic}}^3 \\ & + (2.4491 - 0.12639 \: T_o) \: 10^{-7} \: {c_\mathrm{dic}}^4 \\ & - (1.5768 - 0.15326 \: T_o) \: 10^{-10} \: {c_\mathrm{dic}}^5 \end{aligned} $$

* $$ p_\mathrm{CO2} = (p_\mathrm{dic} + C_\mathrm{pi}) \: \exp (\gamma_\mathrm{dic} \: T) $$

* $$ F_\mathrm{ocean} = \nu_\mathrm{gx} \: (1 + \gamma_\mathrm{gx} \: T) \: (C - p_\mathrm{CO2}) $$

###### prognostic

* $$ \frac{\mathrm{d} C_{o,j}}{\mathrm{d} t} = -\frac{C_{o,j}}{\kappa_{\tau_o} \: \tau_{o,j}} + \alpha_{o,j} \: F_\mathrm{ocean} $$

* $$ \frac{\mathrm{d} C_d}{\mathrm{d} t} = \sum_j \frac{C_{o,j}}{\kappa_{\tau_o} \: \tau_{o,j}} $$

<!--------------->
### 4. Land carbon

###### diagnostic

* $$ r_\mathrm{npp} = \left( 1 + \frac{\beta_\mathrm{npp}}{\alpha_\mathrm{npp}} \: \left( 1 - \left( \frac{C}{C_\mathrm{pi}} \right) ^{-\alpha_\mathrm{npp}} \right) \right) (1 + \gamma_\mathrm{npp} \: T) $$

* $$ r_\mathrm{fire} = \left( 1 + \beta_\mathrm{fire} \left( \frac{C}{C_\mathrm{pi}} - 1 \right) \right) (1 + \gamma_\mathrm{fire} \: T) $$

* $$ r_\mathrm{rh} = \left( 1 + \beta_\mathrm{rh} \left( \frac{C_{s1}}{C_{s1} + C_{s2} + C_{s3}} \left( 1 + \frac{\nu_{stab}}{\nu_{rh23}} \right) - 1 \right) \right) \exp (\gamma_\mathrm{rh} \: T) $$

* $$ F_\mathrm{npp} = F_\mathrm{npp,0} \: r_\mathrm{npp} $$

* $$ E_\mathrm{fire} = \nu_\mathrm{fire} \: r_\mathrm{fire} \: C_v $$

* $$ E_\mathrm{harv} = \nu_\mathrm{harv} \: C_v $$

* $$ F_\mathrm{mort} = \nu_\mathrm{mort} \: C_v $$

* $$ E_\mathrm{rh1} = \nu_\mathrm{rh1} \: r_\mathrm{rh} \: C_{s1} $$

* $$ F_\mathrm{stab} = \nu_\mathrm{stab} \: r_\mathrm{rh} \: C_{s1} $$

* $$ E_\mathrm{rh2} = \frac{\nu_\mathrm{rh23} - \nu_\mathrm{rh3} \: \alpha_\mathrm{pass}}{1 - \alpha_\mathrm{pass}} \: r_\mathrm{rh} \: C_{s2} $$

* $$ F_\mathrm{pass} = \nu_\mathrm{rh3} \: \frac{\alpha_\mathrm{pass}}{1 - \alpha_\mathrm{pass}} \: r_\mathrm{rh} \: C_{s2} $$

* $$ E_\mathrm{rh3} = \nu_\mathrm{rh3} \: r_\mathrm{rh} \: C_{s3} $$

* $$ F_\mathrm{land} = F_\mathrm{npp} - E_\mathrm{fire} - E_\mathrm{harv} - E_\mathrm{rh1} - E_\mathrm{rh2} - E_\mathrm{rh3} $$

###### prognostic

* $$ \frac{\mathrm{d} C_v}{\mathrm{d} t} = F_\mathrm{npp} - E_\mathrm{fire} - E_\mathrm{harv} - F_\mathrm{mort} $$

* $$ \frac{\mathrm{d} C_{s1}}{\mathrm{d} t} = F_\mathrm{mort} - F_\mathrm{stab} - E_\mathrm{rh1} $$

* $$ \frac{\mathrm{d} C_{s2}}{\mathrm{d} t} = F_\mathrm{stab} - F_\mathrm{pass} - E_\mathrm{rh2} $$

* $$ \frac{\mathrm{d} C_{s3}}{\mathrm{d} t} = F_\mathrm{pass} - E_\mathrm{rh3} $$

###### diagnostic (2nd)

* $$ E_\mathrm{rh} = E_\mathrm{rh1} + E_\mathrm{rh2} + E_\mathrm{rh3} $$

* $$ C_s = C_{s1} + C_{s2} + C_{s3} $$

<!--------------------->
### 5. Permafrost carbon

###### diagnostic

* $$ r_\mathrm{rt} = \exp \! \left( \kappa_\mathrm{rt} \: \gamma_\mathrm{rt1} \: \alpha_\mathrm{lst} \: T - \kappa_\mathrm{rt} \: \gamma_\mathrm{rt2} \: (\alpha_\mathrm{lst} \: T)^{2} \right) $$

* $$  \bar{a} = -a_\mathrm{min} + \frac{(1 + a_\mathrm{min})} { \left( 1 + \left( \left( 1 + \frac{1}{a_\mathrm{min}} \right) ^{\kappa_a} - 1 \right) \exp (-\gamma_a \: \kappa_a \: \alpha_\mathrm{lst} \: T) \right) ^\frac{1}{\kappa_a} } $$

* $$ E_\mathrm{pf} = \sum_j \frac{C_{\mathrm{th},j}}{\kappa_{\tau_\mathrm{th}} \: \tau_{\mathrm{th},j}} \: r_\mathrm{rt} $$

###### prognostic

* $$ \frac{\mathrm{d} a}{\mathrm{d} t} = 0.5 \: (\nu_\mathrm{thaw} + \nu_\mathrm{froz}) \: (\bar{a} - a) + 0.5 \: |(\nu_\mathrm{thaw} - \nu_\mathrm{froz}) \: (\bar{a} - a)| $$

* $$ \frac{\mathrm{d} C_{\mathrm{th},j}}{\mathrm{d} t} = \alpha_{\mathrm{th},j} \: \frac{\mathrm{d} a}{\mathrm{d} t} \: C_\mathrm{fr,0} - \frac{C_{\mathrm{th},j}}{\kappa_{\tau_\mathrm{th}} \: \tau_{\mathrm{th},j}} \: r_\mathrm{rt} $$

###### diagnostic (2nd)

* $$ C_\mathrm{fr} = (1 - a) \: C_\mathrm{fr,0} $$

* $$ \frac{\mathrm{d} C_\mathrm{fr}}{\mathrm{d} t} = -\frac{\mathrm{d} a}{\mathrm{d} t} \: C_\mathrm{fr,0} $$

<!------------------->
### 6. Atmospheric CO2

###### diagnostic

* $$ \mathrm{pH} = \kappa_\mathrm{pH} \: (8.5541 - 0.00173 \: C + 1.3264 \: 10^{-6} \: C^2 - 4.4943 \: 10^{-10} \: C^3) $$

###### prognostic

* $$ \alpha_C \: \frac{\mathrm{d} C}{\mathrm{d} t} = E_\mathrm{CO2} + E_\mathrm{pf} - F_\mathrm{land} - F_\mathrm{ocean} $$

