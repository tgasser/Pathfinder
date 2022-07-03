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
import theano.tensor as tt

from core_fct.cls_model import WrapModel


##################################################
## EMISSION DRIVEN MODEL
##################################################

## prognostic variables
Var_name = ['Td', 
            'Hgla', 'Hgis', 'Hais_smb', 'Hais',
            'Co_1', 'Co_2', 'Co_3', 'Co_4', 'Co_5', 'Cd', 
            'Cv', 'Cs1', 'Cs2', 'Cs3', 
            'a', 'Cth_1', 'Cth_2', 'Cth_3']

## diagnostic variables (non-exhaustive)
Var2_name = ['RFco2', 'ERF', 'ERFx', 'logit_ff', 
             'OHC', 'd_OHC', 'Hthx', 'd_Hthx', 'Htot', 'd_Htot', 'Hlia', 
             'Co', 'dic', 'pCO2', 'Focean', 
             'NPP', 'Efire', 'RH', 'Cs', 'Fland', 
             'abar', 'Epf', 'Cfr', 
             'pH', 'Eco2']

## parameters
Par_name = ['phi', 'T2x', 'THs', 'THd', 'th', 'eheat', 'T2x0', 
            'aOHC', 'Lthx', 'lgla0', 'Lgla', 'Ggla1', 'Ggla3', 'tgla', 'ggla', 'lgis0', 'Lgis1', 'Lgis3', 'tgis', 'Lais_smb', 'lais0', 'Lais', 'tais', 'aais', 
            'adic', 'aoc_1', 'aoc_2', 'aoc_3', 'aoc_4', 'aoc_5', 'toc_1', 'toc_2', 'toc_3', 'toc_4', 'toc_5', 'k_toc', 'vgx', 'ggx', 'To', 'bdic', 'gdic', 
            'npp0', 'vfire', 'vharv', 'vmort', 'vstab', 'vrh1', 'vrh23', 'vrh3', 'apass', 'bnpp', 'anpp', 'gnpp', 'bfire', 'gfire', 'brh', 'grh', 
            'aLST', 'grt1', 'grt2', 'krt', 'amin', 'ka', 'ga', 'vthaw', 'vfroz', 'ath_1', 'ath_2', 'ath_3', 'tth_1', 'tth_2', 'tth_3', 'k_tth', 'Cfr0', 
            'aCO2', 'CO2pi', 'k_pH']
    
## forcings
For_name = ['T', 'd_T', 'CO2', 'd_CO2']


##################################################
##################################################

## default initial values
def Ini_dflt(Par):

    ## unpack parameters
    [_, _, _, _, _, _, _,
    _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, 
    _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
    npp0, vfire, vharv, vmort, vstab, vrh1, vrh23, _, apass, _, _, _, _, _, _, _, 
    _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
    _, _, _] = [Par[n] for n in range(len(Par_name))]

    ## calculate variables
    Cv0 = npp0 / (vfire + vharv + vmort)
    Cl0 = Cv0 * vmort / (vrh1 + vstab)
    Cs0 = Cl0 * vstab / vrh23 * (1 - apass)
    Cp0 = Cl0 * vstab / vrh23 * apass

    ## pack initial values
    Ini = [0, 
    0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 
    Cv0, Cl0, Cs0, Cp0, 
    0, 0, 0, 0]

    ## return
    return Ini


## speed of linear system
def v_linear(Par):

    ## unpack parameters
    [_, _, _, THd, th, _, _, 
    _, _, _, _, _, _, tgla, _, _, _, _, tgis, _, _, _, tais, _, 
    adic, aoc_1, aoc_2, aoc_3, aoc_4, aoc_5, toc_1, toc_2, toc_3, toc_4, toc_5, k_toc, vgx, _, To, bdic, _,
    _, vfire, vharv, vmort, vstab, vrh1, vrh23, vrh3, apass, _, _, _, _, _, _, _, 
    _, _, _, _, _, _, _, vthaw, vfroz, _, _, _, tth_1, tth_2, tth_3, k_tth, _,
    _, _, _] = [Par[n] for n in range(len(Par_name))]

    ## calculate linear speeds
    v_Td = th / THd
    v_Hgla = 1. / tgla
    v_Hgis = 1. / tgis
    v_Hais_smb = 1E-9
    v_Hais = 1. / tais
    v_Co_1 = 1. / toc_1 / k_toc  + aoc_1 * vgx * (1.5568 - 1.3993E-2 * To) * adic / bdic
    v_Co_2 = 1. / toc_2 / k_toc  + aoc_2 * vgx * (1.5568 - 1.3993E-2 * To) * adic / bdic
    v_Co_3 = 1. / toc_3 / k_toc  + aoc_3 * vgx * (1.5568 - 1.3993E-2 * To) * adic / bdic
    v_Co_4 = 1. / toc_4 / k_toc  + aoc_4 * vgx * (1.5568 - 1.3993E-2 * To) * adic / bdic
    v_Co_5 = 1. / toc_5 / k_toc  + aoc_5 * vgx * (1.5568 - 1.3993E-2 * To) * adic / bdic
    v_Cd = 1E-9
    v_Cv = vfire + vharv + vmort 
    v_Cs1 = vstab + vrh1
    v_Cs2 = vrh23 / (1 - apass)
    v_Cs3 = vrh3
    v_a = 0.5 * (vthaw + vfroz)
    v_Cth_1 = 1. / tth_1 / k_tth
    v_Cth_2 = 1. / tth_2 / k_tth
    v_Cth_3 = 1. / tth_3 / k_tth

    ## pack speeds
    v = [v_Td,
    v_Hgla, v_Hgis, v_Hais_smb, v_Hais,
    v_Co_1, v_Co_2, v_Co_3, v_Co_4, v_Co_5, v_Cd,
    v_Cv, v_Cs1, v_Cs2, v_Cs3, 
    v_a, v_Cth_1, v_Cth_2, v_Cth_3]

    ## return
    return v


##################################################
##################################################

## differential system
def d_Var(t, Var, Par, For=None, autonomous=False, tensor=False, expost=False):

    ## unpack variables
    [Td,
    Hgla, Hgis, Hais_smb, Hais,
    Co_1, Co_2, Co_3, Co_4, Co_5, _,
    Cv, Cs1, Cs2, Cs3, 
    a, Cth_1, Cth_2, Cth_3] = [Var[n] for n in range(len(Var_name))]

    ## unpack parameters
    [phi, T2x, THs, THd, th, eheat, T2x0,
    aOHC, Lthx, lgla0, Lgla, Ggla1, Ggla3, tgla, ggla, lgis0, Lgis1, Lgis3, tgis, Lais_smb, lais0, Lais, tais, aais, 
    adic, aoc_1, aoc_2, aoc_3, aoc_4, aoc_5, toc_1, toc_2, toc_3, toc_4, toc_5, k_toc, vgx, ggx, To, bdic, gdic,
    npp0, vfire, vharv, vmort, vstab, vrh1, vrh23, vrh3, apass, bnpp, anpp, gnpp, bfire, gfire, brh, grh, 
    aLST, grt1, grt2, krt, amin, ka, ga, vthaw, vfroz, ath_1, ath_2, ath_3, tth_1, tth_2, tth_3, k_tth, Cfr0,
    aCO2, CO2pi, k_pH] = [Par[n] for n in range(len(Par_name))]

    ## unpack forcings
    if For is None: T, d_T, CO2, d_CO2 = 0, 0, 0, 0
    else: T, d_T, CO2, d_CO2 = [For[n] for n in range(len(For_name))]

    ## get forcings over t-th year (if not autonomous)
    if not autonomous and For is not None:
        if tensor: T, d_T, CO2, d_CO2 = [tt.as_tensor_variable(Z)[tt.cast(t+1-1E-9, 'int64')] for Z in [T, d_T, CO2, d_CO2]]
        else: T, d_T, CO2, d_CO2 = [Z[np.array(t+1-1E-9, dtype=int)] for Z in [T, d_T, CO2, d_CO2]]

    ## get correct mathematical functions
    if tensor: exp, log, abs_ = tt.exp, tt.log, tt.abs_
    else: exp, log, abs_ = np.exp, np.log, np.abs

    ## =======

    ## 1. CLIMATE
    ## diagnostic
    RFco2 = phi * log(CO2 / CO2pi)
    ERF = THs * d_T + phi * log(2) / T2x * T + eheat * th * (T - Td)
    ERFx = ERF - RFco2
    ## prognostic
    d_Td = 1. / THd * th * (T - Td)
    ## diagnostic 2 (for calib)
    logit_ff = 0*T + log(T2x / T2x0 - 1.)

    ## 2. SEA LEVEL
    ## diagnostic
    OHC = aOHC * (THs * T + THd * Td)
    d_OHC = aOHC * (THs * d_T + THd * d_Td)
    Hthx = Lthx * OHC
    d_Hthx = Lthx * d_OHC
    ## prognostic
    d_Hgla = lgla0 + (Lgla * (1. - exp(-Ggla1 * T - Ggla3 * T**3)) - Hgla) / tgla * exp(ggla * T)
    d_Hgis = lgis0 + (Lgis1 * T + Lgis3 * T**3 - Hgis) / tgis
    d_Hais_smb = -Lais_smb * T
    d_Hais = d_Hais_smb + lais0 + (Lais * T - (Hais - Hais_smb)) / tais * (1 + aais * (Hais - Hais_smb))
    ## diagnostic 2
    Htot = Hthx + Hgla + Hgis + Hais
    d_Htot = d_Hthx + d_Hgla + d_Hgis + d_Hais
    ## diagnostic 3 (for calib)
    Hlia = 0*Htot + lgla0 * tgla * (exp(-150. / tgla) - exp(-255. / tgla)) + lgis0 * tgis * (exp(-150. / tgis) - exp(-255. / tgis)) + lais0 * tais * (exp(-150. / tais) - exp(-255. / tais))

    ## 3. OCEAN CARBON
    ## diagnostic
    Co = Co_1 + Co_2 + Co_3 + Co_4 + Co_5
    dic = adic / bdic * Co
    pdic = ((1.5568 - 1.3993E-2 * To) * dic 
        + (7.4706 - 0.20207 * To) * 1E-3 * dic**2. 
        - (1.2748 - 0.12015 * To) * 1E-5 * dic**3. 
        + (2.4491 - 0.12639 * To) * 1E-7 * dic**4. 
        - (1.5768 - 0.15326 * To) * 1E-10 * dic**5.)
    pCO2 = (pdic + CO2pi) * exp(gdic * T)
    Focean = vgx * (1. + ggx * T) * (CO2 - pCO2)
    ## prognostic
    d_Co_1 = aoc_1 * Focean - Co_1 / toc_1 / k_toc
    d_Co_2 = aoc_2 * Focean - Co_2 / toc_2 / k_toc
    d_Co_3 = aoc_3 * Focean - Co_3 / toc_3 / k_toc
    d_Co_4 = aoc_4 * Focean - Co_4 / toc_4 / k_toc
    d_Co_5 = aoc_5 * Focean - Co_5 / toc_5 / k_toc
    d_Cd = (Co_1 / toc_1 + Co_2 / toc_2 + Co_3 / toc_3 + Co_4 / toc_4 + Co_5 / toc_5) / k_toc

    ## 4. LAND CARBON
    ## convenience
    r_npp = (1. + bnpp / anpp * (1. - (CO2 / CO2pi)**-anpp)) * (1. + gnpp * T)
    r_fire = (1. + bfire * (CO2 / CO2pi - 1.)) * (1. + gfire * T)
    r_rh = (1. + brh * (Cs1 / (Cs1 + Cs2 + Cs3) * (1. + vstab / vrh23) - 1.)) * exp(grh * T)
    ## diagnostic
    NPP = npp0 * r_npp
    Efire = vfire * r_fire * Cv 
    Eharv = vharv * Cv 
    Fmort = vmort * Cv
    RH1 = vrh1 * r_rh * Cs1
    Fstab = vstab * r_rh * Cs1
    RH2 = (vrh23 - vrh3 * apass) / (1. - apass) * r_rh * Cs2
    Fpass = vrh3 * apass / (1. - apass) * r_rh * Cs2
    RH3 = vrh3 * r_rh * Cs3
    Fland = NPP - Efire - Eharv - RH1 - RH2 - RH3
    ## prognostic
    d_Cv = NPP - Efire - Eharv - Fmort
    d_Cs1 = Fmort - Fstab - RH1
    d_Cs2 = Fstab - Fpass - RH2
    d_Cs3 = Fpass - RH3
    ## diagnostic 2
    RH = RH1 + RH2 + RH3
    Cs = Cs1 + Cs2 + Cs3

    ## 5. PERMAFROST
    ## convenience
    r_rt = exp(krt * grt1 * aLST * T - krt * grt2 * (aLST * T)**2.)
    ## diagnostic
    abar = -amin + (1. + amin) / (1. + ((1. + 1. / amin)**ka - 1.) * exp(-ga * ka * aLST * T))**(1. / ka)
    Epf = (Cth_1 / tth_1 + Cth_2 / tth_2 + Cth_3 / tth_3) * r_rt / k_tth
    ## prognostic
    d_a = 0.5 * (vthaw + vfroz) * (abar - a) + 0.5 * abs_((vthaw - vfroz) * (abar - a))
    d_Cth_1 = ath_1 * d_a * Cfr0 - Cth_1 / tth_1 * r_rt / k_tth
    d_Cth_2 = ath_2 * d_a * Cfr0 - Cth_2 / tth_2 * r_rt / k_tth
    d_Cth_3 = ath_3 * d_a * Cfr0 - Cth_3 / tth_3 * r_rt / k_tth
    ## diagnostic 2
    Cfr = (1. - a) * Cfr0
    #d_Cfr = -d_a * Cfr0

    ## 6. ATMOSPHERIC CO2
    ## diagnostic
    pH = k_pH * (8.5541 - 0.00173 * CO2 + 1.3264E-6 * CO2**2 - 4.4943E-10 * CO2**3)
    Eco2 = aCO2 * d_CO2 - Epf + Fland + Focean

    ## =======

    ## pack derivatives
    d_Var = [d_Td,
    d_Hgla, d_Hgis, d_Hais_smb, d_Hais,
    d_Co_1, d_Co_2, d_Co_3, d_Co_4, d_Co_5, d_Cd,
    d_Cv, d_Cs1, d_Cs2, d_Cs3, 
    d_a, d_Cth_1, d_Cth_2, d_Cth_3]

    ## pack secondary variables
    if expost:
        Var2 = [RFco2, ERF, ERFx, logit_ff, 
        OHC, d_OHC, Hthx, d_Hthx, Htot, d_Htot, Hlia, 
        Co, dic, pCO2, Focean,
        NPP, Efire, RH, Cs, Fland,
        abar, Epf, Cfr, 
        pH, Eco2]

    ## return
    if expost: return Var2 + d_Var
    else: return d_Var


##################################################
##################################################

## wrapping
PF_OBSdriven = WrapModel(Var_name, Var2_name, Par_name, For_name, Ini_dflt, v_linear, d_Var)

