
<!------------------->
## Differential System

<!-------->
### 1. Climate

###### diagnostic

<img src="https://latex.codecogs.com/gif.latex? 
R_c = \phi \: \ln \! \left( \frac{C}{C_\mathrm{pi}} \right)
" />

<img src="https://latex.codecogs.com/gif.latex? 
R = R_c + R_x
" />

###### prognostic

<img src="https://latex.codecogs.com/gif.latex? 
\Theta_s \: \frac{\mathrm{d} T}{\mathrm{d} t} = 
R - \frac{\phi \: \ln (2)}{T_{2\times}} \: T - \epsilon_\mathrm{heat} \: \theta \: (T - T_d)
" />

<img src="https://latex.codecogs.com/gif.latex? 
\Theta_d \: \frac{\mathrm{d} T_d}{\mathrm{d} t} = 
\theta \: (T - T_d)
" />

<!---------->
### 2. Sea level

###### diagnostic

<img src="https://latex.codecogs.com/gif.latex? 
U_\mathrm{ohc} = 
\alpha_\mathrm{ohc} \: (\Theta_s \: T + \Theta_d \: T_d)
" />

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} U_\mathrm{ohc}}{\mathrm{d} t} = 
\alpha_\mathrm{ohc} \: (\Theta_s \: \frac{\mathrm{d} T}{\mathrm{d} t} + \Theta_d \: \frac{\mathrm{d} T_d}{\mathrm{d} t})
" />

<img src="https://latex.codecogs.com/gif.latex? 
H_\mathrm{lin} = 
\Lambda_\mathrm{lin} \: U_\mathrm{ohc}
" />

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} H_\mathrm{lin}}{\mathrm{d} t} = 
\Lambda_\mathrm{lin} \: \frac{\mathrm{d} U_\mathrm{ohc}}{\mathrm{d} t}
" />

###### prognostic

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} H_\mathrm{thx}}{\mathrm{d} t} = 
\frac{\mathrm{d} H_\mathrm{lin}}{\mathrm{d} t} - \frac{H_\mathrm{thx} - H_\mathrm{lin}}{\tau_\mathrm{thx}} + \frac{\Lambda_\mathrm{thx} - \Lambda_\mathrm{lin} \: \alpha_\mathrm{ohc} \: (\Theta_s + \Theta_d)}{\tau_\mathrm{thx}} \: T_d
" />

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} H_\mathrm{tot}}{\mathrm{d} t} = 
\frac{\mathrm{d} H_\mathrm{thx}}{\mathrm{d} t} - \frac{H_\mathrm{tot} - H_\mathrm{thx}}{\tau_\mathrm{ice}} + \frac{\Lambda_\mathrm{tot1} + \Lambda_\mathrm{tot2} \: T - \Lambda_\mathrm{thx}}{\tau_\mathrm{ice}} \: T
" />

###### diagnostic (2nd)

<img src="https://latex.codecogs.com/gif.latex? 
H_\mathrm{ice} = 
H_\mathrm{tot} - H_\mathrm{thx}
" />

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} H_\mathrm{ice}}{\mathrm{d} t} = 
\frac{\mathrm{d} H_\mathrm{tot}}{\mathrm{d} t} - \frac{\mathrm{d} H_\mathrm{thx}}{\mathrm{d} t}
" />

<!------------->
### 3. Ocean carbon

###### diagnostic

<img src="https://latex.codecogs.com/gif.latex? 
C_o = 
\sum_j C_{o,j}
" />

<img src="https://latex.codecogs.com/gif.latex? 
c_\mathrm{dic} = 
\frac{\alpha_\mathrm{dic}}{\beta_\mathrm{dic}} \: C_o
" />

<img src="https://latex.codecogs.com/gif.latex?
\begin{align*}
p_\mathrm{dic} = 
\: & (1.5568 - 0.013993 \: T_o) \: c_\mathrm{dic}\\ + 
\: & (7.4706 - 0.20207 \: T_o) \: 10^{-3} \: {c_\mathrm{dic}}^2\\ - 
\: & (1.2748 - 0.12015 \: T_o) \: 10^{-5} \: {c_\mathrm{dic}}^3\\ + 
\: & (2.4491 - 0.12639 \: T_o) \: 10^{-7} \: {c_\mathrm{dic}}^4\\ - 
\: & (1.5768 - 0.15326 \: T_o) \: 10^{-10} \: {c_\mathrm{dic}}^5
\end{align}
" />

<img src="https://latex.codecogs.com/gif.latex?
p_\mathrm{CO2} = 
(p_\mathrm{dic} + C_\mathrm{pi}) \: \exp (\gamma_\mathrm{dic} \: T)
" />

<img src="https://latex.codecogs.com/gif.latex?
F_\mathrm{ocean} = 
\nu_\mathrm{gx} \: (1 + \gamma_\mathrm{gx} \: T) \: (C - p_\mathrm{CO2})
" />

###### prognostic

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_{o,j}}{\mathrm{d} t} = 
-\frac{C_{o,j}}{\kappa_{\tau_o} \: \tau_{o,j}} + \alpha_{o,j} \: F_\mathrm{ocean}
" />

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_d}{\mathrm{d} t} = 
\sum_j \frac{C_{o,j}}{\kappa_{\tau_o} \: \tau_{o,j}} 
" />

<!------------>
### 4. Land carbon

###### diagnostic

<img src="https://latex.codecogs.com/gif.latex?
r_\mathrm{npp} = 
\left( 1 + \frac{\beta_\mathrm{npp}}{\alpha_\mathrm{npp}} \: \left( 1 - \left( \frac{C}{C_\mathrm{pi}} \right) ^{-\alpha_\mathrm{npp}} \right) \right) (1 + \gamma_\mathrm{npp} \: T)
" />

<img src="https://latex.codecogs.com/gif.latex?
r_\mathrm{ef} = 
\left( 1 + \beta_\mathrm{ef} \left( \frac{C}{C_\mathrm{pi}} - 1 \right) \right) (1 + \gamma_\mathrm{ef} \: T)
" />

<img src="https://latex.codecogs.com/gif.latex?
r_\mathrm{rh} = 
\left( 1 + \beta_\mathrm{rh} \left( \frac{C_{s1}}{C_{s1} + C_{s2} + C_{s3}} \left( 1 + \frac{\nu_{met}}{\nu_{cs2}} \right) - 1 \right) \right) \exp (\gamma_\mathrm{rh} \: T)
" />

<img src="https://latex.codecogs.com/gif.latex?
F_\mathrm{npp} = 
F_\mathrm{npp,0} \: r_\mathrm{npp}
" />

<img src="https://latex.codecogs.com/gif.latex?
E_\mathrm{fire} = 
\nu_\mathrm{fire} \: r_\mathrm{ef} \: C_v 
" />

<img src="https://latex.codecogs.com/gif.latex?
E_\mathrm{harv} = 
\nu_\mathrm{harv} \: C_v 
" />

<img src="https://latex.codecogs.com/gif.latex?
F_\mathrm{mort} = 
\nu_\mathrm{mort} \: C_v
" />

<img src="https://latex.codecogs.com/gif.latex?
E_\mathrm{rh1} = 
\nu_\mathrm{rh1} \: r_\mathrm{rh} \: C_{s1}
" />

<img src="https://latex.codecogs.com/gif.latex?
F_\mathrm{met} = 
\nu_\mathrm{met} \: r_\mathrm{rh} \: C_{s1}
" />

<img src="https://latex.codecogs.com/gif.latex?
E_\mathrm{rh2} = 
\frac{\nu_\mathrm{cs2} - \nu_\mathrm{rh3} \: \alpha_\mathrm{pass}}{1 - \alpha_\mathrm{pass}} \: r_\mathrm{rh} \: C_{s2}
" />

<img src="https://latex.codecogs.com/gif.latex?
F_\mathrm{pass} = 
\nu_\mathrm{rh3} \: \frac{\alpha_\mathrm{pass}}{1 - \alpha_\mathrm{pass}} \: r_\mathrm{rh} \: C_{s2}
" />

<img src="https://latex.codecogs.com/gif.latex?
E_\mathrm{rh3} = 
\nu_\mathrm{rh3} \: r_\mathrm{rh} \: C_{s3}
" />

<img src="https://latex.codecogs.com/gif.latex?
F_\mathrm{land} = 
F_\mathrm{npp} - E_\mathrm{fire} - E_\mathrm{harv} - E_\mathrm{rh1} - E_\mathrm{rh2} - E_\mathrm{rh3}
" />

###### prognostic

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_v}{\mathrm{d} t} = 
F_\mathrm{npp} - E_\mathrm{fire} - E_\mathrm{harv} - F_\mathrm{mort}
" />

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_{s1}}{\mathrm{d} t} = 
F_\mathrm{mort} - F_\mathrm{met} - E_\mathrm{rh1}
" />

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_{s2}}{\mathrm{d} t} = 
F_\mathrm{met} - F_\mathrm{pass} - E_\mathrm{rh2}
" />

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_{s3}}{\mathrm{d} t} = 
F_\mathrm{pass} - E_\mathrm{rh3}
" />

###### diagnostic (2nd)

<img src="https://latex.codecogs.com/gif.latex? 
E_\mathrm{rh} = 
E_\mathrm{rh1} + E_\mathrm{rh2} + E_\mathrm{rh3}
" />

<img src="https://latex.codecogs.com/gif.latex? 
C_s = 
C_\mathrm{s1} + C_\mathrm{s2} + C_\mathrm{s3}
" />

<!------------------>
### 5. Permafrost carbon

###### diagnostic

<img src="https://latex.codecogs.com/gif.latex? 
r_\mathrm{rt} = 
\exp \! \left( \kappa_\mathrm{rt} \: \gamma_\mathrm{rt1} \: \alpha_\mathrm{lst} \: T - \kappa_\mathrm{rt} \: \gamma_\mathrm{rt2} \: (\alpha_\mathrm{lst} \: T)^{2} \right)
" />

<img src="https://latex.codecogs.com/gif.latex? 
\bar{a} = 
-a_\mathrm{min} + \frac{(1 + a_\mathrm{min})} { \left( 1 + \left( \left( 1 + \frac{1}{a_\mathrm{min}} \right) ^{\kappa_a} - 1 \right) \exp (-\gamma_a \: \kappa_a \: \alpha_\mathrm{lst} \: T) \right) ^\frac{1}{\kappa_a} }
" />

<img src="https://latex.codecogs.com/gif.latex? 
E_\mathrm{pf} = 
\sum_j \frac{C_{\mathrm{th},j}}{\kappa_{\tau_\mathrm{th}} \: \tau_{\mathrm{th},j}} \: r_\mathrm{rt}
" />

###### prognostic

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} a}{\mathrm{d} t} = 
0.5 \: (\nu_\mathrm{thaw} + \nu_\mathrm{froz}) \: (\bar{a} - a) + 0.5 \: |(\nu_\mathrm{thaw} - \nu_\mathrm{froz}) \: (\bar{a} - a)|
" />

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_{\mathrm{th},j}}{\mathrm{d} t} = 
\alpha_{\mathrm{th},j} \: \frac{\mathrm{d} a}{\mathrm{d} t} \: C_\mathrm{fr,0} - \frac{C_{\mathrm{th},j}}{\kappa_{\tau_\mathrm{th}} \: \tau_{\mathrm{th},j}} \: r_\mathrm{rt}
" />

###### diagnostic (2nd)

<img src="https://latex.codecogs.com/gif.latex? 
C_\mathrm{fr} = 
(1 - a) \: C_\mathrm{fr,0}
" />

<img src="https://latex.codecogs.com/gif.latex? 
\frac{\mathrm{d} C_\mathrm{fr}}{\mathrm{d} t} = 
-\frac{\mathrm{d} a}{\mathrm{d} t} \: C_\mathrm{fr,0}
" />

<!---------------->
### 6. Atmospheric CO2

###### diagnostic

<img src="https://latex.codecogs.com/gif.latex? 
\mathrm{pH} [C] = 
\kappa_\mathrm{pH} \: (8.5541 - 0.00173 \: C + 1.3264 \: 10^{-6} \: C^2 - 4.4943 \: 10^{-10} \: C^3)
" />

###### prognostic

<img src="https://latex.codecogs.com/gif.latex? 
\alpha_C \: \frac{\mathrm{d} C}{\mathrm{d} t} = 
E_\mathrm{CO2} + E_\mathrm{pf} - F_\mathrm{land} - F_\mathrm{ocean}
" />


<!--------->
<!--------->
## Notations

<!-------->
### Drivers

| In manual | In code | Units |
| --- | --- | --- |
| <img src="https://latex.codecogs.com/gif.latex?R_x" /> | `ERFx` | W m<sup>-2</sup> |
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{CO2}" /> | `Eco2` | PgC yr<sup>-1</sup> |

<!---------->
### Variables 

| In manual | In code | Units | Prog? | Dims |
| --- | --- | --- | --- | --- |
| <img src="https://latex.codecogs.com/gif.latex?R_c" /> | `RFco2` | W m<sup>-2</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?R" /> | `ERF` | W m<sup>-2</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?T" /> | `T` | K | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?T_d" /> | `Td` | K | *yes* ||
||||||
| <img src="https://latex.codecogs.com/gif.latex?U_\mathrm{ohc}" /> | `OHC` | W yr m<sup>-2</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?H_\mathrm{lin}" /> | `Hlin` | mm |||
| <img src="https://latex.codecogs.com/gif.latex?H_\mathrm{thx}" /> | `Hthx` | mm | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?H_\mathrm{ice}" /> | `Hice` | mm |||
| <img src="https://latex.codecogs.com/gif.latex?H_\mathrm{tot}" /> | `Htot` | mm | *yes* ||
||||||
| <img src="https://latex.codecogs.com/gif.latex?C_{o,j}" /> | `Co_j` | PgC | *yes* | <img src="https://latex.codecogs.com/gif.latex?j\in[\![1,5]\!]" /> |
| <img src="https://latex.codecogs.com/gif.latex?C_o" /> | `Co` | PgC |||
| <img src="https://latex.codecogs.com/gif.latex?C_d" /> | `Cd` | PgC | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?c_\mathrm{dic}" /> | `dic` | µmol kg<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?p_\mathrm{dic}" /> | `pdic` | ppm |||
| <img src="https://latex.codecogs.com/gif.latex?p_\mathrm{CO2}" /> | `pCO2` | ppm |||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{ocean}" /> | `Focean` | PgC yr<sup>-1</sup> |||
||||||
| <img src="https://latex.codecogs.com/gif.latex?r_\mathrm{npp}" /> | `r_npp` | 1 |||
| <img src="https://latex.codecogs.com/gif.latex?r_\mathrm{ef}" /> | `r_ef` | 1 |||
| <img src="https://latex.codecogs.com/gif.latex?r_\mathrm{rh}" /> | `r_rh` | 1 |||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{npp}" /> | `NPP` | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{fire}" /> | `Efire` | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{harv}" /> | `Eharv` | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{mort}" /> | `Fmort` | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{rh1}" /> | `RH1` | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{met}" /> | `Fmet` | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{rh2}" /> | `RH2` | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{pass}" /> | `Fpass` | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{rh3}" /> | `RH3` | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{land}" /> | `Fland` | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{rh}" /> | `RH` | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?C_v" /> | `Cv` | PgC | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?C_{s1}" /> | `Cs1` | PgC | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?C_{s2}" /> | `Cs2` | PgC | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?C_{s3}" /> | `Cs3` | PgC | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?C_s" /> | `Cs` | PgC |||
||||||
| <img src="https://latex.codecogs.com/gif.latex?r_\mathrm{rt}" /> | `r_rt` | 1 |||
| <img src="https://latex.codecogs.com/gif.latex?\bar{a}" /> | `abar` | 1 |||
| <img src="https://latex.codecogs.com/gif.latex?a" /> | `a` | 1 | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?E_\mathrm{pf}" /> | `Epf` | PgC yr<sup>-1</sup> |||
| <img src="https://latex.codecogs.com/gif.latex?C_{\mathrm{th},j}" /> | `Cth_j` | PgC | *yes* | <img src="https://latex.codecogs.com/gif.latex?j\in[\![1,3]\!]" /> |
| <img src="https://latex.codecogs.com/gif.latex?C_\mathrm{fr}" /> | `Cfr` | PgC |||
||||||
| <img src="https://latex.codecogs.com/gif.latex?C" /> | `CO2` | ppm | *yes* ||
| <img src="https://latex.codecogs.com/gif.latex?\mathrm{pH}" /> | `pH` | 1 |||

<!---------->
### Parameters 

| In manual | In code | Units | Dims |
| --- | --- | --- | --- |
| <img src="https://latex.codecogs.com/gif.latex?\phi" /> | `phi` | W m<sup>-2</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?T_\mathrm{2\times}" /> | `T2x` | K ||
| <img src="https://latex.codecogs.com/gif.latex?\Theta_s" /> | `THs` | W yr m<sup>-2</sup> K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\Theta_d" /> | `THd` | W yr m<sup>-2</sup> K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\theta" /> | `th` | W m<sup>-2</sup> K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\epsilon_\mathrm{heat}" /> | `eheat` | 1 ||
|||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{ohc}" /> | `aOHC` | 1  ||
| <img src="https://latex.codecogs.com/gif.latex?\Lambda_\mathrm{lin}" /> | `Llin` | mm m<sup>2</sup> W<sup>-1</sup> yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\Lambda_\mathrm{thx}" /> | `Lthx` | mm K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\Lambda_\mathrm{tot1}" /> | `Ltot1` | mm K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\Lambda_\mathrm{tot2}" /> | `Ltot2` | mm K<sup>-2</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{thx}" /> | `tthx` | yr ||
| <img src="https://latex.codecogs.com/gif.latex?\tau_\mathrm{ice}" /> | `tice` | yr ||
|||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{dic}" /> | `adic` | µmol kg<sup>-1</sup> PgC<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\beta_\mathrm{dic}" /> | `bdic` | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{dic}" /> | `gdic` | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?T_o" /> | `To` | °C ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{gx}" /> | `vgx` | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{gx}" /> | `ggx` | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_{o,j}" /> | `aoc_j` | 1 | <img src="https://latex.codecogs.com/gif.latex?j\in[\![1,5]\!]" /> |
| <img src="https://latex.codecogs.com/gif.latex?\tau_{o,j}" /> | `toc_j` | yr | <img src="https://latex.codecogs.com/gif.latex?j\in[\![1,5]\!]" /> |
| <img src="https://latex.codecogs.com/gif.latex?\kappa_{\tau_o}" /> | `k_toc` | 1 ||
|||||
| <img src="https://latex.codecogs.com/gif.latex?\beta_\mathrm{npp}" /> | `bnpp` | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{npp}" /> | `anpp` | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{npp}" /> | `gnpp` | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\beta_\mathrm{ef}" /> | `bef` | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{ef}" /> | `gef` | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\beta_\mathrm{rh}" /> | `brh` | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{rh}" /> | `grh` | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?F_\mathrm{npp,0}" /> | `npp0` | PgC yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{fire}" /> | `vfire` | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{harv}" /> | `vharv` | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{mort}" /> | `vmort` | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{met}" /> | `vmet` | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{rh1}" /> | `vrh1` | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{cs2}" /> | `vcs2` | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{rh3}" /> | `vrh3` | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{pass}" /> | `apass` | 1 ||
|||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_\mathrm{lst}" /> | `aLST` | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{rt1}" /> | `grt1` | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_\mathrm{rt2}" /> | `grt2` | K<sup>-2</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{rt}" /> | `krt` | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?a_\mathrm{min}" /> | `amin` | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\kappa_a" /> | `ka` | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?\gamma_a" /> | `ga` | K<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{thaw}" /> | `vthaw` | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{froz}" /> | `vfroz` | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\nu_\mathrm{froz}" /> | `vfroz` | yr<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_{\mathrm{th},j}" /> | `ath_j` | 1 | <img src="https://latex.codecogs.com/gif.latex?j\in[\![1,3]\!]" /> |
| <img src="https://latex.codecogs.com/gif.latex?\tau_{\mathrm{th},j}" /> | `tth_j` | yr | <img src="https://latex.codecogs.com/gif.latex?j\in[\![1,3]\!]" /> |
| <img src="https://latex.codecogs.com/gif.latex?\kappa_{\tau_\mathrm{th}}" /> | `k_tth` | 1 ||
| <img src="https://latex.codecogs.com/gif.latex?C_\mathrm{fr,0}" /> | `Cfr0` | PgC ||
|||||
| <img src="https://latex.codecogs.com/gif.latex?\alpha_C" /> | `aCO2` | PgC ppm<sup>-1</sup> ||
| <img src="https://latex.codecogs.com/gif.latex?C_\mathrm{pi}" /> | `CO2pi` | ppm ||
| <img src="https://latex.codecogs.com/gif.latex?\kappa_\mathrm{pH}" /> | `k_pH` | 1 ||

