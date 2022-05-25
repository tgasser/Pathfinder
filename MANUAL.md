

<!--------->
<!--------->
## Notations

<!----------------------->
### Forcings and Variables 

| In manual| In code | Description | Units | Prog? | Dims |
| --- | --- | --- | --- | --- | --- |
| <img src="https://latex.codecogs.com/gif.latex?R_c" /> | `RFco2` | CO2 (effective) radiative forcing | W m<sup>-2</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?R_x" /> | `ERFx` | Non-CO2 effective radiative forcing | W m<sup>-2</sup>  |||
| <img src="https://latex.codecogs.com/gif.latex?R" /> | `ERF` | Effective radiative forcing | W m<sup>-2</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?T" /> | `T` | Global surface temperature anomaly | K | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?T_d" /> | `Td` | Deep ocean temperature anomaly | K | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?\mathrm{logit} (\mathrm{ff})" /> | `logit_ff` | Logit of the climate feedback factor (for calib.) | 1 |||
|||||||
| <img src="https://latex.codecogs.com/gif.latex?U_\mathrm{ohc}" /> | `OHC` | Ocean heat content (anomaly) | W yr m<sup>-2</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?H_\mathrm{lin}" /> | `Hlin` | Linear part of thermosteric sea level rise | mm |||
| <img src="https://latex.codecogs.com/gif.latex?H_\mathrm{thx}" /> | `Hthx` | Total thermosteric sea level rise | mm | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?H_\mathrm{ice}" /> | `Hice` | Ice contributions to sea level rise | mm |||
| <img src="https://latex.codecogs.com/gif.latex?H_\mathrm{tot}" /> | `Htot` | Total sea level rise | mm | *yes* ||
|||||||
| <img src="https://latex.codecogs.com/gif.latex?C_{o,j}" /> | `Co_j` | Change in surface ocean carbon subpools | PgC | *yes* | <img src="https://latex.codecogs.com/gif.latex?j\in[\![1,5]\!]" /> |
| <img src="https://latex.codecogs.com/gif.latex?C_o" /> | `Co` | Change in surface ocean carbon pool | PgC |||
| <img src="https://latex.codecogs.com/gif.latex?C_d" /> | `Cd` | Change in deep ocean carbon pool | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?c_\mathrm{dic}" /> | `dic` | Change in surface DIC | µmol kg<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?p_\mathrm{dic}" /> | `pdic` | Subcomponent of `pCO2` | ppm |||
| <img src="https://latex.codecogs.com/gif.latex?p_\mathrm{CO2}" /> | `pCO2` | CO2 partial pressure at the ocean surface | ppm |||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{ocean}" /> | `Focean` | Ocean carbon sink | PgC yr<sup>-1</sup> |||
|||||||
| <img src="https://latex.codecogs.com/gif.latex?r_\mathrm{npp}" /> | `r_npp` | Relative change in NPP | 1 |||
| <img src="https://latex.codecogs.com/gif.latex?r_\mathrm{fire}" /> | `r_fire` | Relative change in wildfire intensity | 1 |||
| <img src="https://latex.codecogs.com/gif.latex?r_\mathrm{rh}" /> | `r_rh` | Relative change in heterotrophic respiration rate | 1 |||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{npp}" /> | `NPP` | Net primary productivity | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{fire}" /> | `Efire` | Emissions from wildfire | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{harv}" /> | `Eharv` | Emissions from harvest and grazing | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{mort}" /> | `Fmort` | Mortality flux | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{rh1}" /> | `RH1` | Litter heterotrophic respiration | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{stab}" /> | `Fstab` | Stabilization flux | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{rh2}" /> | `RH2` | Active soil heterotrophic respiration | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{pass}" /> | `Fpass` | Passivization flux | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{rh3}" /> | `RH3` | Passive soil heterotrophic respiration | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{land}" /> | `Fland` | Land carbon sink | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{rh}" /> | `RH` | Heterotrophic respiration | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?C_v" /> | `Cv` | Vegetation carbon pool | PgC | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?C_{s1}" /> | `Cs1` | Litter carbon pool | PgC | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?C_{s2}" /> | `Cs2` | Active soil carbon pool | PgC | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?C_{s3}" /> | `Cs3` | Passive soil carbon pool | PgC | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?C_s" /> | `Cs` | Total soil carbon pool | PgC |||
|||||||
| <img src="https://latex.codecogs.com/gif.latex?r_\mathrm{rt}" /> | `r_rt` | Relative change in permafrost respiration rate | 1 |||
| <img src="https://latex.codecogs.com/gif.latex?\bar{a}" /> | `abar` | Theoretical thawed fraction | 1 |||
| <img src="https://latex.codecogs.com/gif.latex?a" /> | `a` | Actual thawed fraction | 1 | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{pf}" /> | `Epf` | Emissions from permafrost | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?C_{\mathrm{th},j}" /> | `Cth_j` | Thawed permafrost carbon subpools | PgC | *yes* | <img src="https://latex.codecogs.com/gif.latex?j\in[\![1,3]\!]" /> |
| <img src="https://latex.codecogs.com/gif.latex?C_\mathrm{fr}" /> | `Cfr` | Frozen permafrost carbon pool | PgC |||
|||||||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{CO2}" /> | `Eco2` | Anthropogenic CO2 emissions | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?C" /> | `CO2` | Atmospheric CO2 concentration | ppm | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?\mathrm{pH}" /> | `pH` | Surface ocean pH | 1 |||

<!----------->
### Parameters 

| In manual | In code | Description | Units | Dims |
| --- | --- | --- | --- | --- |
| <img src="https://latex.codecogs.com/gif.latex?\phi" /> | `phi` | Radiative parameter of CO2 | W m<sup>-2</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?T_\mathrm{2\times}" /> | `T2x` | Equilibrium climate sensitivity | K ||
| <img src="https://latex.codecogs.com/gif.latex?\Theta_s" /> | `THs` | Heat capacity of the surface | W yr m<sup>-2</sup> K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\Theta_d" /> | `THd` | Heat capacity of the deep ocean | W yr m<sup>-2</sup> K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\theta" /> | `th` | Heat exchange coefficient | W m<sup>-2</sup> K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\epsilon_\mathrm{heat}" /> | `eheat` | Deep ocean heat uptake efficacy | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?T_\mathrm{2\times}^*" /> | `T2x0` | Minimal value of the ECS distribution (for calib.) | K ||
||||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{ohc}" /> | `aOHC` | Fraction of energy warming the ocean | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\Lambda_\mathrm{lin}" /> | `Llin` | Linear factor for thermosteric SLR | mm m<sup>2</sup> W<sup>-1</sup> yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\Lambda_\mathrm{thx}" /> | `Lthx` | Equilibrium thermosteric SLR | mm K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\Lambda_\mathrm{tot1}" /> | `Ltot1` | Linear equilibrium ice-related SLR | mm K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\Lambda_\mathrm{tot2}" /> | `Ltot2` | Quadratic equilibrium ice-related SLR | mm K<sup>-2</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{thx}" /> | `tthx` | Timescale of thermosteric SLR | yr ||
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{ice}" /> | `tice` | Timescale of ice-related SLR | yr ||
||||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{dic}" /> | `adic` | Conversion factor for DIC | µmol kg<sup>-1</sup> PgC<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\beta_\mathrm{dic}" /> | `bdic` | Inverse-scaling factor for DIC | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{dic}" /> | `gdic` | Sensitivity of pCO2 to climate | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?T_o" /> | `To` | Preindustrial surface ocean temperature | °C ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{gx}" /> | `vgx` | Surface ocean gas exchange rate | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{gx}" /> | `ggx` | Sensitivity of gas exchange to climate | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_{o,j}" /> | `aoc_j` | Surface ocean subpools fractions | 1 | <img src="https://latex.codecogs.com/gif.latex?j\in[\![1,5]\!]" /> |
| <img src="https://latex.codecogs.com/gif.latex?\tau_{o,j}" /> | `toc_j` | Timescales of surface ocean subpools | yr | <img src="https://latex.codecogs.com/gif.latex?j\in[\![1,5]\!]" /> |
| <img src="https://latex.codecogs.com/gif.latex?\kappa_{\tau_o}" /> | `k_toc` | Scaling factor for timescales of surface ocean subpools | 1 ||
||||||
| <img src="https://latex.codecogs.com/gif.latex?\beta_\mathrm{npp}" /> | `bnpp` | Sensitivity of NPP to CO2 (= fertilization effect) | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{npp}" /> | `anpp` | Shape parameter for fertilization effect | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{npp}" /> | `gnpp` | Sensitivity of NPP to climate | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\beta_\mathrm{fire}" /> | `bfire` | Sensitivity of wildfire intensity to CO2 | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{fire}" /> | `gfire` | Sensitivity of wildfire intensity to climate | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\beta_\mathrm{rh}" /> | `brh` | Sensitivity of heterotrophic respiration to fresh organic matter | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{rh}" /> | `grh` | Sensitivity of heterotrophic respiration to climate | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{npp,0}" /> | `npp0` | Preindustrial NPP | PgC yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{fire}" /> | `vfire` | Wildfire intensity | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{harv}" /> | `vharv` | Harvest and grazing rate | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{mort}" /> | `vmort` | Mortality rate | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{stab}" /> | `vstab` | Stabilization rate | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{rh1}" /> | `vrh1` | Litter heterotrophic respiration rate | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{rh23}" /> | `vrh23` | Soil (active and passive) respiration rate | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{rh3}" /> | `vrh3` | Passive soil respiration rate | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{pass}" /> | `apass` | Fraction of passive soil | 1 ||
||||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{lst}" /> | `aLST` | Climate scaling factor over permafrost regions | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{rt1}" /> | `grt1` | Sensitivity of (boreal) heterotrophic respiration to climate | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{rt2}" /> | `grt2` | Sensitivity of (boreal) heterotrophic respiration to climate (quadratic) | K<sup>-2</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{rt}" /> | `krt` | Scaling factor for sensitivity of permafrost respiration to climate | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?a_\mathrm{min}" /> | `amin` | Minimal thawed fraction | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\kappa_a" /> | `ka` | Shape parameter for theoretical thawed fraction | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_a" /> | `ga` | Sensitivity of theoretical thawed fraction to climate | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{thaw}" /> | `vthaw` | Thawing rate | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{froz}" /> | `vfroz` | Freezing rate | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_{\mathrm{th},j}" /> | `ath_j` | Thawed permafrost carbon subpools fractions | 1 | <img src="https://latex.codecogs.com/gif.latex?j\in[\![1,3]\!]" /> |
| <img src="https://latex.codecogs.com/gif.latex?\tau_{\mathrm{th},j}" /> | `tth_j` | Timescales of thawed permafrost carbon subpools | yr | <img src="https://latex.codecogs.com/gif.latex?j\in[\![1,3]\!]" /> |
| <img src="https://latex.codecogs.com/gif.latex?\kappa_{\tau_\mathrm{th}}" /> | `k_tth` | Scaling factor for timescales of surface ocean subpools | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?C_\mathrm{fr,0}" /> | `Cfr0` | Preindustrial frozen permafrost carbon pool | PgC ||
||||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_C" /> | `aCO2` | Conversion factor for atmospheric CO2 | PgC ppm<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?C_\mathrm{pi}" /> | `CO2pi` | Preindustrial CO2 concentration | ppm ||
| <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{pH}" /> | `k_pH` | Scaling factor for surface ocean pH | 1 ||
||||||
| <img src="https://latex.codecogs.com/gif.latex?\tilde{\sigma}_C" /> | `std_CO2` | Relative standard deviation of the historical `CO2` time series (for calib.) | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\epsilon_C" /> | `ampl_CO2` | Noise amplitude of the historical `CO2` time series (for calib.) | ppm ||
| <img src="https://latex.codecogs.com/gif.latex?\rho_C" /> | `corr_CO2` | Autocorrelation of the historical `CO2` time series (for calib.) | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\tilde{\sigma}_T" /> | `std_T` | Relative standard deviation of the historical `T` time series (for calib.) | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\epsilon_T" /> | `ampl_T` | Noise amplitude of the historical `T` time series (for calib.) | K ||
| <img src="https://latex.codecogs.com/gif.latex?\rho_T" /> | `corr_T` | Autocorrelation of the historical `T` time series (for calib.) | 1 ||

<!--------->
<!--------->
## Equations

<!----------->
### 1. Climate

###### diagnostic

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
R_c = \phi \: \ln \! \left( \frac{C}{C_\mathrm{pi}} \right)
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
R = R_c + R_x
" />

###### prognostic

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\Theta_s \: \frac{\mathrm{d} T}{\mathrm{d} t} = 
R - \frac{\phi \: \ln (2)}{T_{2\times}} \: T - \epsilon_\mathrm{heat} \: \theta \: (T - T_d)
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\Theta_d \: \frac{\mathrm{d} T_d}{\mathrm{d} t} = 
\theta \: (T - T_d)
" />

###### diagnostic (2nd; for calib.)

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
\mathrm{logit} (\mathrm{ff}) = 
\ln \! \left( \frac{T_{2\times}}{T_{2\times}^*} - 1 \right)
" />

<!------------->
### 2. Sea level

###### diagnostic

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
U_\mathrm{ohc} = 
\alpha_\mathrm{ohc} \: (\Theta_s \: T + \Theta_d \: T_d)
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} U_\mathrm{ohc}}{\mathrm{d} t} = 
\alpha_\mathrm{ohc} \: (\Theta_s \: \frac{\mathrm{d} T}{\mathrm{d} t} + \Theta_d \: \frac{\mathrm{d} T_d}{\mathrm{d} t})
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
H_\mathrm{lin} = 
\Lambda_\mathrm{lin} \: U_\mathrm{ohc}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} H_\mathrm{lin}}{\mathrm{d} t} = 
\Lambda_\mathrm{lin} \: \frac{\mathrm{d} U_\mathrm{ohc}}{\mathrm{d} t}
" />

###### prognostic

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} H_\mathrm{thx}}{\mathrm{d} t} = 
\frac{\mathrm{d} H_\mathrm{lin}}{\mathrm{d} t} - \frac{H_\mathrm{thx} - H_\mathrm{lin}}{\tau_\mathrm{thx}} + \frac{\Lambda_\mathrm{thx} - \Lambda_\mathrm{lin} \: \alpha_\mathrm{ohc} \: (\Theta_s + \Theta_d)}{\tau_\mathrm{thx}} \: T_d
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} H_\mathrm{tot}}{\mathrm{d} t} = 
\frac{\mathrm{d} H_\mathrm{thx}}{\mathrm{d} t} - \frac{H_\mathrm{tot} - H_\mathrm{thx}}{\tau_\mathrm{ice}} + \frac{\Lambda_\mathrm{tot1} + \Lambda_\mathrm{tot2} \: T - \Lambda_\mathrm{thx}}{\tau_\mathrm{ice}} \: T
" />

###### diagnostic (2nd)

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
H_\mathrm{ice} = 
H_\mathrm{tot} - H_\mathrm{thx}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} H_\mathrm{ice}}{\mathrm{d} t} = 
\frac{\mathrm{d} H_\mathrm{tot}}{\mathrm{d} t} - \frac{\mathrm{d} H_\mathrm{thx}}{\mathrm{d} t}
" />

<!---------------->
### 3. Ocean carbon

###### diagnostic

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
C_o = 
\sum_j C_{o,j}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
c_\mathrm{dic} = 
\frac{\alpha_\mathrm{dic}}{\beta_\mathrm{dic}} \: C_o
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
\begin{align*}
p_\mathrm{dic} = 
\: & (1.5568 - 0.013993 \: T_o) \: c_\mathrm{dic}\\ + 
\: & (7.4706 - 0.20207 \: T_o) \: 10^{-3} \: {c_\mathrm{dic}}^2\\ - 
\: & (1.2748 - 0.12015 \: T_o) \: 10^{-5} \: {c_\mathrm{dic}}^3\\ + 
\: & (2.4491 - 0.12639 \: T_o) \: 10^{-7} \: {c_\mathrm{dic}}^4\\ - 
\: & (1.5768 - 0.15326 \: T_o) \: 10^{-10} \: {c_\mathrm{dic}}^5
\end{align}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
p_\mathrm{CO2} = 
(p_\mathrm{dic} + C_\mathrm{pi}) \: \exp (\gamma_\mathrm{dic} \: T)
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
F_\mathrm{ocean} = 
\nu_\mathrm{gx} \: (1 + \gamma_\mathrm{gx} \: T) \: (C - p_\mathrm{CO2})
" />

###### prognostic

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_{o,j}}{\mathrm{d} t} = 
-\frac{C_{o,j}}{\kappa_{\tau_o} \: \tau_{o,j}} + \alpha_{o,j} \: F_\mathrm{ocean}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_d}{\mathrm{d} t} = 
\sum_j \frac{C_{o,j}}{\kappa_{\tau_o} \: \tau_{o,j}} 
" />

<!--------------->
### 4. Land carbon

###### diagnostic

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
r_\mathrm{npp} = 
\left( 1 + \frac{\beta_\mathrm{npp}}{\alpha_\mathrm{npp}} \: \left( 1 - \left( \frac{C}{C_\mathrm{pi}} \right) ^{-\alpha_\mathrm{npp}} \right) \right) (1 + \gamma_\mathrm{npp} \: T)
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
r_\mathrm{fire} = 
\left( 1 + \beta_\mathrm{fire} \left( \frac{C}{C_\mathrm{pi}} - 1 \right) \right) (1 + \gamma_\mathrm{fire} \: T)
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
r_\mathrm{rh} = 
\left( 1 + \beta_\mathrm{rh} \left( \frac{C_{s1}}{C_{s1} + C_{s2} + C_{s3}} \left( 1 + \frac{\nu_{stab}}{\nu_{rh23}} \right) - 1 \right) \right) \exp (\gamma_\mathrm{rh} \: T)
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
F_\mathrm{npp} = 
F_\mathrm{npp,0} \: r_\mathrm{npp}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
E_\mathrm{fire} = 
\nu_\mathrm{fire} \: r_\mathrm{fire} \: C_v 
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
E_\mathrm{harv} = 
\nu_\mathrm{harv} \: C_v 
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
F_\mathrm{mort} = 
\nu_\mathrm{mort} \: C_v
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
E_\mathrm{rh1} = 
\nu_\mathrm{rh1} \: r_\mathrm{rh} \: C_{s1}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
F_\mathrm{stab} = 
\nu_\mathrm{stab} \: r_\mathrm{rh} \: C_{s1}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
E_\mathrm{rh2} = 
\frac{\nu_\mathrm{rh23} - \nu_\mathrm{rh3} \: \alpha_\mathrm{pass}}{1 - \alpha_\mathrm{pass}} \: r_\mathrm{rh} \: C_{s2}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
F_\mathrm{pass} = 
\nu_\mathrm{rh3} \: \frac{\alpha_\mathrm{pass}}{1 - \alpha_\mathrm{pass}} \: r_\mathrm{rh} \: C_{s2}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
E_\mathrm{rh3} = 
\nu_\mathrm{rh3} \: r_\mathrm{rh} \: C_{s3}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex?
F_\mathrm{land} = 
F_\mathrm{npp} - E_\mathrm{fire} - E_\mathrm{harv} - E_\mathrm{rh1} - E_\mathrm{rh2} - E_\mathrm{rh3}
" />

###### prognostic

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_v}{\mathrm{d} t} = 
F_\mathrm{npp} - E_\mathrm{fire} - E_\mathrm{harv} - F_\mathrm{mort}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_{s1}}{\mathrm{d} t} = 
F_\mathrm{mort} - F_\mathrm{stab} - E_\mathrm{rh1}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_{s2}}{\mathrm{d} t} = 
F_\mathrm{stab} - F_\mathrm{pass} - E_\mathrm{rh2}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_{s3}}{\mathrm{d} t} = 
F_\mathrm{pass} - E_\mathrm{rh3}
" />

###### diagnostic (2nd)

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
E_\mathrm{rh} = 
E_\mathrm{rh1} + E_\mathrm{rh2} + E_\mathrm{rh3}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
C_s = 
C_\mathrm{s1} + C_\mathrm{s2} + C_\mathrm{s3}
" />

<!--------------------->
### 5. Permafrost carbon

###### diagnostic

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
r_\mathrm{rt} = 
\exp \! \left( \kappa_\mathrm{rt} \: \gamma_\mathrm{rt1} \: \alpha_\mathrm{lst} \: T - \kappa_\mathrm{rt} \: \gamma_\mathrm{rt2} \: (\alpha_\mathrm{lst} \: T)^{2} \right)
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\bar{a} = 
-a_\mathrm{min} + \frac{(1 + a_\mathrm{min})} { \left( 1 + \left( \left( 1 + \frac{1}{a_\mathrm{min}} \right) ^{\kappa_a} - 1 \right) \exp (-\gamma_a \: \kappa_a \: \alpha_\mathrm{lst} \: T) \right) ^\frac{1}{\kappa_a} }
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
E_\mathrm{pf} = 
\sum_j \frac{C_{\mathrm{th},j}}{\kappa_{\tau_\mathrm{th}} \: \tau_{\mathrm{th},j}} \: r_\mathrm{rt}
" />

###### prognostic

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} a}{\mathrm{d} t} = 
0.5 \: (\nu_\mathrm{thaw} + \nu_\mathrm{froz}) \: (\bar{a} - a) + 0.5 \: |(\nu_\mathrm{thaw} - \nu_\mathrm{froz}) \: (\bar{a} - a)|
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_{\mathrm{th},j}}{\mathrm{d} t} = 
\alpha_{\mathrm{th},j} \: \frac{\mathrm{d} a}{\mathrm{d} t} \: C_\mathrm{fr,0} - \frac{C_{\mathrm{th},j}}{\kappa_{\tau_\mathrm{th}} \: \tau_{\mathrm{th},j}} \: r_\mathrm{rt}
" />

###### diagnostic (2nd)

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
C_\mathrm{fr} = 
(1 - a) \: C_\mathrm{fr,0}
" />

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_\mathrm{fr}}{\mathrm{d} t} = 
-\frac{\mathrm{d} a}{\mathrm{d} t} \: C_\mathrm{fr,0}
" />

<!------------------->
### 6. Atmospheric CO2

###### diagnostic

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\mathrm{pH} = 
\kappa_\mathrm{pH} \: (8.5541 - 0.00173 \: C + 1.3264 \: 10^{-6} \: C^2 - 4.4943 \: 10^{-10} \: C^3)
" />

###### prognostic

* <img style="vertical-align:middle" src="https://latex.codecogs.com/gif.latex? 
\alpha_C \: \frac{\mathrm{d} C}{\mathrm{d} t} = 
E_\mathrm{CO2} + E_\mathrm{pf} - F_\mathrm{land} - F_\mathrm{ocean}
" />

