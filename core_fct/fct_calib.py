"""
Copyright IIASA (International Institute for Applied Systems Analysis), 2019-2022
Contributor(s): Thomas Gasser (gasser@iiasa.ac.at)

This software is a computer program whose purpose is simulate a large number of global climate (and environmental) change scenarios.

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
"""


import theano
import datetime
import pymc3 as pm
import numpy as np
import xarray as xr
import theano.tensor as tt

from time import perf_counter

from core_fct.cls_calib import my_AR1
from core_fct.mod_OBSdriven import PF_OBSdriven as PF
from core_fct.mod_Edriven import PF_Edriven as PFe
from core_fct.fct_ancillary import lognorm_distrib_param, logitnorm_distrib_param
from core_fct.fct_default import get_dflt_param, get_dflt_forcing, get_dflt_constr


##################################################
## ANCILLARY FUNCTIONS
##################################################

## function to get 2nd-order accurate differential
## same as np.gradient(..., edge_order=2)
gradient = lambda X, n: [-1.5*X[0] + 2*X[1] - 0.5*X[2]] + [0.5*X[t+1] - 0.5*X[t-1] for t in range(1, n-1)] + [1.5*X[-1] + -2*X[-2] + 0.5*X[-3]]

## functions to iterate over the differential system
## implicit-explicit scheme with four steps in total
## note: better check consistency with PF call in post_calib function
def imex_4steps(F, X, p, v):
    Xi = (X + 0.25*(tt.stack(PF(None, X, p, F, autonomous=True, tensor=True, expost=False)) + v*X)) / (1 + v*0.25)
    Xii = (Xi + 0.25*(tt.stack(PF(None, Xi, p, F, autonomous=True, tensor=True, expost=False)) + v*Xi)) / (1 + v*0.25)
    Xiii = (Xii + 0.25*(tt.stack(PF(None, Xii, p, F, autonomous=True, tensor=True, expost=False)) + v*Xii)) / (1 + v*0.25)
    return (Xiii + 0.25*(tt.stack(PF(None, Xiii, p, F, autonomous=True, tensor=True, expost=False)) + v*Xiii)) / (1 + v*0.25)


##################################################
## MAIN CALIBRATION FUNCTIONS
##################################################

## function returning pyMC3 model 
## note: quite a bunch of options are hard-coded in there
def build_model(redo_prior_calib=True, std_max=5, ignore_constr=[], uniform_prior=False, with_noise=True, f_scheme=imex_4steps):

    ## load parameters, drivers and constraints
    t_counter = perf_counter()
    print('** loading inputs **')
    Par0 = get_dflt_param(redo_prior_calib=redo_prior_calib)
    For0 = get_dflt_forcing()
    Con0 = get_dflt_constr()
    print('(done in {0:.0f} seconds)'.format(perf_counter() - t_counter))

    ## define model
    with pm.Model() as model:
        t_counter = perf_counter()
        print('** building model **')

        ##===========
        ## Parameters
        ##===========

        ## create all parameters
        param = []
        for par in PF.Par_name:

            ## calibrated parameters
            if 'stat' in Par0[par].dims:
                mean, std = float(Par0[par][0]), float(Par0[par][1])
                bounds = [max(Par0[par].bounds[0], mean - std_max * std), min(Par0[par].bounds[1], mean + std_max * std)]
                if par == 'T2x' and 'T2x0' in Par0: bounds[0] = max(bounds[0], float(Par0.T2x0))
                                    
                ## uniform distrib
                if uniform_prior:
                    param.append( pm.Uniform(par, lower=bounds[0], upper=bounds[1]) )
                
                ## log-normal distrib
                elif Par0[par].bounds == (0., np.inf):
                    mu, sigma = lognorm_distrib_param(mean, std)
                    param.append( pm.Bound(pm.Lognormal, lower=bounds[0], upper=bounds[1])(par, mu=mu, sigma=sigma) )
                
                ## logit-normal distrib
                elif Par0[par].bounds == (0., 1.):
                    mu, sigma = logitnorm_distrib_param(float(Par0[par][0]), float(Par0[par][1]), raise_error=True)
                    param.append( pm.Bound(pm.LogitNormal, lower=bounds[0], upper=bounds[1])(par, mu=mu, sigma=sigma) )
                
                ## normal distrib
                else:
                    mu, sigma = mean, std
                    param.append( pm.Bound(pm.Normal, lower=bounds[0], upper=bounds[1])(par, mu=mu, sigma=sigma) )
            
            ## fixed parameters
            else:
                param.append( float(Par0[par]) )


        ##========
        ## Drivers
        ##========

        ## years
        t0 = int(For0.year[0])
        n_year = len(For0.year)

        ## temperature
        T_mean = tt.as_tensor_variable(For0.T.sel(stat='mean').values)
        T_std = tt.as_tensor_variable(For0.T.sel(stat='std').values)
        if with_noise:
            if uniform_prior:
                ampl_T = pm.Uniform('ampl_T', lower=0, upper=std_max*0.05)
            else:
                ampl_T = pm.Bound(pm.HalfNormal, lower=0, upper=std_max*0.05)('ampl_T', sigma=0.05)
            corr_T = pm.Uniform('corr_T', lower=0, upper=1)
            ar_T = my_AR1('ar_T', k=corr_T, shape=n_year)
        else:
            ar_T, ampl_T = 0., 0.
        std_T = pm.Bound(pm.Normal, lower=-std_max*1, upper=+std_max*1)('std_T', mu=0, sigma=1)

        ## CO2
        CO2_mean = tt.as_tensor_variable(For0.CO2.sel(stat='mean').values)
        CO2_std = tt.as_tensor_variable(For0.CO2.sel(stat='std').values)
        if with_noise:
            if uniform_prior:
                ampl_CO2 = pm.Uniform('ampl_CO2', lower=0, upper=std_max*0.5)
            else:
                ampl_CO2 = pm.Bound(pm.HalfNormal, lower=0, upper=std_max*0.5)('ampl_CO2', sigma=0.5)
            corr_CO2 = pm.Uniform('corr_CO2', lower=0, upper=1)
            ar_CO2 = my_AR1('ar_CO2', k=corr_CO2, shape=n_year)
        else:
            ar_CO2, ampl_CO2 = 0., 0.
        std_CO2 = pm.Bound(pm.Normal, lower=-std_max*1, upper=+std_max*1)('std_CO2', mu=0, sigma=1)

        ## pack drivers in dict
        forcing = {}
        forcing['T'] = pm.Deterministic('T', T_mean + std_T * T_std + ampl_T * ar_T)
        forcing['d_T'] = pm.Deterministic('d_T', tt.stack(gradient(forcing['T'], n_year)))
        forcing['CO2'] = pm.Deterministic('CO2', CO2_mean + std_CO2 * CO2_std + ampl_CO2 * ar_CO2)
        forcing['d_CO2'] = pm.Deterministic('d_CO2', tt.stack(gradient(forcing['CO2'], n_year)))


        ##=============
        ## Diff. system
        ##=============

        ## forcing (as tt.matrix)
        forc = tt.stack([forcing[var] for var in PF.For_name])
        ## initial state (as tt.vector)
        init = tt.stack(PF.Ini_dflt(param))
        ## linear system speeds
        v_lin = PF.v_linear(param)

        ## prognostic variables(as tt.matrix)
        var_prog = theano.scan(f_scheme, sequences=[tt.transpose(forc)], outputs_info=[init], non_sequences=[param, v_lin])[0]

        ## diagnostic variables (as list of tt.vector)
        var_diag = PF(None, Var=tt.transpose(var_prog), Par=param, For=forc, autonomous=True, tensor=True, expost=True)


        ##============
        ## Constraints
        ##============

        ## add constraints to model
        for var in [var for var in Con0 if var not in ignore_constr]:
            
            ## indices
            inds = [Con0[var].period[0] - t0, Con0[var].period[1] - t0 + 1]

            ## select variable
            var_ = var.split('__')[0]
            if var_ in PF.For_name:
                var_tmp = forcing[var_]
            elif var_ in PF.Var_name:
                var_tmp = var_prog[:, PF.Var_name.index(var_)]
            elif var_ in PF.Var2_name or var in ['d_'+name for name in PF.Var_name]:
                var_tmp = var_diag[(PF.Var2_name + ['d_'+name for name in PF.Var_name]).index(var_)]
            else:
                raise RuntimeError(var_ + ' in constraints but not in output variables or forcings')

            ## create constraint
            if Con0[var].is_sum:
                pm.Normal(var+'_obs', mu=var_tmp[slice(*inds)].sum(), sigma=float(Con0[var][1]), observed=float(Con0[var][0]))
            elif Con0[var].is_mean:
                pm.Normal(var+'_obs', mu=var_tmp[slice(*inds)].mean(), sigma=float(Con0[var][1]), observed=float(Con0[var][0]))
            elif Con0[var].is_diff:
                pm.Normal(var+'_obs', mu=var_tmp[slice(*inds)][-1]-var_tmp[slice(*inds)][0], sigma=float(Con0[var][1]), observed=float(Con0[var][0]))


    ## RETURN
    print('(done in {0:.0f} seconds)'.format(perf_counter() - t_counter))
    return model


##################################################
##################################################

## function to execute Bayesian calibration
## note: most options are hard-coded
def exec_calib(method='FullRankADVI', name=None, n_sample=2000, **model_args):
    folder_out = 'internal_data/pyMC_calib/'

    ## get model
    model = build_model(**model_args)

    ## calibrate!
    print('** Bayesian magic **')
    with model:

        ## MCMC
        if method in ['NUTS', 'Metropolis']:
            if method == 'NUTS': 
                method_mcmc = pm.NUTS(target_accept=0.99)
                trace = pm.sample(n_sample, tune=3000, step=method_mcmc, cores=1, chains=2)
            elif method == 'Metropolis': 
                method_mcmc = pm.Metropolis()
                trace = pm.sample(n_sample, tune=4000, step=method_mcmc, cores=1, chains=4)
            sample = trace

        ## VI
        elif method in ['ADVI', 'FullRankADVI', 'SVGD']:
            if method == 'ADVI': 
                method_vi = pm.ADVI()
                n_iter = 100000
            elif method == 'FullRankADVI': 
                method_vi = pm.FullRankADVI()
                n_iter = 100000
            elif method == 'SVGD': 
                method_vi = pm.SVGD(n_particles=1000)
                n_iter = 5000
            tracker = pm.callbacks.Tracker(mean=method_vi.approx.mean.eval, std=method_vi.approx.std.eval)
            approx = method_vi.fit(n_iter, callbacks=[tracker])
            sample = approx.sample(n_sample)
        
        ## else raise error
        else: raise ValueError('wrong method; must be NUTS, Metropolis, ADVI, FullRankADVI, or SVGD')

    ## save folder
    date = datetime.datetime.today().strftime("%Y_%m_%d__%H_%M_%S")
    folder_calib = folder_out + method + '__' + name if name is not None else date
    
    ## trace/samples
    pm.backends.save_trace(sample, folder_calib, overwrite=True)
    with open(folder_calib + '/sample.npz', 'wb') as f:
        np.savez(f, **{var:sample[var] for var in sample.varnames}, allow_pickle=False)
    
    ## convergence info
    if method in ['ADVI', 'FullRankADVI']:
        np.save(folder_calib + '/elbo.npy', approx.hist, allow_pickle=False)
        np.save(folder_calib + '/tracker_mean.npy', np.array(tracker['mean']), allow_pickle=False)
        np.save(folder_calib + '/tracker_std.npy', np.array(tracker['std']), allow_pickle=False)

    ## RETURN
    return sample


##################################################
## POST-PROCESSING FUNCTIONS
##################################################

## function to post-process one set of calibration
def post_calib(folder_calib, name=None, save_prior=False, save_Var2=False, write_csv=False, years_csv=slice(2015-2, 2015+2), **model_args):
    
    ## get ignored constraints (if any)
    ignore_constr = model_args['ignore_constr'] if 'ignore_constr' in model_args else []

    ## AR1 parameters (names as in build_model, with units)
    ar1_pars = {'std_T': '1', 'ampl_T': 'K', 'corr_T': '1', 'std_CO2': '1', 'ampl_CO2': 'ppm', 'corr_CO2': '1'}

    ##==========
    ## Posterior
    ##==========

    ## load default
    Par0 = get_dflt_param(redo_prior_calib=model_args['redo_prior_calib'] if 'redo_prior_calib' in model_args else True)
    For0 = get_dflt_forcing()
    Con0 = get_dflt_constr()

    ## load calibration
    with np.load(open(folder_calib + '/sample.npz', 'rb')) as TMP:
        Par_post = xr.Dataset({var:(('config'), TMP[var]) for var in TMP.keys() if var in PF.Par_name + list(ar1_pars.keys())})
        For_post = xr.Dataset({var:(('year', 'config'), TMP[var].T) for var in TMP.keys() if var in PF.For_name})

    ## format
    Par_post.coords['config'] = np.arange(Par_post.dims['config'])
    For_post.coords['config'] = np.arange(For_post.dims['config'])
    For_post.coords['year'] = 1750 + np.arange(For_post.dims['year'])

    ## add fixed parameters
    Par_post = xr.merge([Par_post, Par0.drop([var for var in Par_post if var in Par0] + ['stat'])])

    ## re-run historical period
    Out_post = PF.run_xarray(Par_post, For_post, scheme='imex', nt=4, get_Var2=True)

    ## get constraints values
    Con_post = xr.Dataset()
    for var in [var for var in Con0 if var not in ignore_constr]:
        if var.split('__')[0] in PF.For_name: TMP = For_post[var.split('__')[0]].sel(year=slice(*Con0[var].period))
        else: TMP = Out_post[var.split('__')[0]].sel(year=slice(*Con0[var].period))
        if Con0[var].is_sum: Con_post[var] = TMP.sum('year')
        elif Con0[var].is_mean: Con_post[var] = TMP.mean('year')
        elif Con0[var].is_diff: Con_post[var] = TMP[-1] - TMP[0]

    ## add units
    for var in Par_post: Par_post[var].attrs['units'] = Par0[var].units if var in Par0 else ar1_pars[var]
    for var in For_post: For_post[var].attrs['units'] = For0[var].units
    for var in Con_post: Con_post[var].attrs['units'] = Con0[var].units

    ##======
    ## Prior
    ##======

    ## make prior if requested
    Par_prior = For_prior = Out_prior = None
    if save_prior:

        ## reload pyMC3 model
        model = build_model(**model_args)

        ## load default drivers
        T_mean = For0.T.sel(stat='mean').values
        T_std = For0.T.sel(stat='std').values
        CO2_mean = For0.CO2.sel(stat='mean').values
        CO2_std = For0.CO2.sel(stat='std').values

        ## sample prior parameters
        with model:
            sample = pm.sample_prior_predictive(len(Par_post.config), var_names=[var for var in PF.Par_name if 'stat' in Par0[var].dims] + list(ar1_pars.keys()))
            
        ## AR1 timeseries (because sampling not possible...)
        ar_T = [np.random.normal(0, 1, size=len(Par_post.config))]
        for _ in range(1, len(For0.year)): 
            ar_T.append( sample['corr_T'] * ar_T[-1] + np.random.normal(0, np.sqrt(1 - sample['corr_T']**2)) )
        ar_CO2 = [np.random.normal(0, 1, size=len(Par_post.config))]
        for _ in range(1, len(For0.year)): 
            ar_CO2.append( sample['corr_CO2'] * ar_CO2[-1] + np.random.normal(0, np.sqrt(1 - sample['corr_CO2']**2)) )

        ## add deterministic variables to sample
        sample['T'] = T_mean[:,np.newaxis] + sample['std_T'][np.newaxis,:] * T_std[:,np.newaxis] +  sample['ampl_T'][np.newaxis,:] * np.array(ar_T)
        sample['d_T'] = np.gradient(sample['T'], edge_order=2, axis=0)
        sample['CO2'] = CO2_mean[:,np.newaxis] + sample['std_CO2'][np.newaxis,:] * CO2_std[:,np.newaxis] + sample['ampl_CO2'][np.newaxis,:] * np.array(ar_CO2)
        sample['d_CO2'] = np.gradient(sample['CO2'], edge_order=2, axis=0)

        ## put in xarray
        Par_prior = xr.Dataset({var:(('config'), sample[var]) for var in sample.keys() if var in PF.Par_name + list(ar1_pars.keys())})
        For_prior = xr.Dataset({var:(('year', 'config'), sample[var]) for var in sample.keys() if var in PF.For_name})
    
        ## format
        Par_prior.coords['config'] = np.arange(Par_prior.dims['config'])
        For_prior.coords['config'] = np.arange(For_prior.dims['config'])
        For_prior.coords['year'] = 1750 + np.arange(For_prior.dims['year'])

        ## as previously
        Par_prior = xr.merge([Par_prior, Par0.drop([var for var in Par_prior if var in Par0] + ['stat'])])
        Out_prior = PF.run_xarray(Par_prior, For_prior, scheme='imex', nt=4, get_Var2=True)

        ## get constraints values
        Con_prior = xr.Dataset()
        for var in [var for var in Con0 if var not in ignore_constr]:
            if var.split('__')[0] in PF.For_name: TMP = For_prior[var.split('__')[0]].sel(year=slice(*Con0[var].period))
            else: TMP = Out_prior[var.split('__')[0]].sel(year=slice(*Con0[var].period))
            if Con0[var].is_sum: Con_prior[var] = TMP.sum('year')
            elif Con0[var].is_mean: Con_prior[var] = TMP.mean('year')
            elif Con0[var].is_diff: Con_prior[var] = TMP[-1] - TMP[0]

        ## add units
        for var in Par_prior: Par_prior[var].attrs['units'] = Par0[var].units if var in Par0 else ar1_pars[var]
        for var in For_prior: For_prior[var].attrs['units'] = For0[var].units
        for var in Con_prior: Con_prior[var].attrs['units'] = Con0[var].units

    ##=====
    ## Save
    ##=====
    
    ## decide name
    name = name if name is not None else folder_calib.split('/')[-1].split('__')[-1]

    ## save everything as netcdf
    ## save posterior
    Par2_post = Par_post.drop(var for var in Par_post if var not in ar1_pars.keys())
    Par2_post.to_netcdf(folder_calib + '/../Par2_' + name + '.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Par2_post})
    Par1_post = Par_post.drop(ar1_pars.keys())
    Par1_post.to_netcdf(folder_calib + '/../Par_' + name + '.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Par1_post})
    Var_post = xr.merge([For_post, Out_post.drop([var for var in Out_post if var not in PF.Var_name + PFe.For_name])])
    Var_post.to_netcdf(folder_calib + '/../Var_' + name + '.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Var_post})
    if save_Var2:
        Var2_post = Out_post.drop([var for var in Out_post if var in PF.Var_name + PFe.For_name])
        Var2_post.to_netcdf(folder_calib + '/../Var2_' + name + '.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Var2_post})
    Con_post.to_netcdf(folder_calib + '/../Con_' + name + '.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Con_post})

    ## save prior
    if save_prior:
        Par2_prior = Par_prior.drop(var for var in Par_prior if var not in ar1_pars.keys())
        Par2_prior.to_netcdf(folder_calib + '/../Par2_' + name + '_prior.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Par2_prior})
        Par1_prior = Par_prior.drop(ar1_pars.keys())
        Par1_prior.to_netcdf(folder_calib + '/../Par_' + name + '_prior.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Par1_prior})
        Var_prior = xr.merge([For_prior, Out_prior.drop([var for var in Out_prior if var not in PF.Var_name + PFe.For_name])])
        Var_prior.to_netcdf(folder_calib + '/../Var_' + name + '_prior.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Var_prior})
        if save_Var2:
            Var2_prior = Out_prior.drop([var for var in Out_prior if var in PF.Var_name + PFe.For_name])
            Var2_prior.to_netcdf(folder_calib + '/../Var2_' + name + '_prior.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Var2_prior})
        Con_prior.to_netcdf(folder_calib + '/../Con_' + name + '_prior.nc', encoding={var:{'zlib':True, 'dtype':np.float32} for var in Con_prior})

    ## also write only partial info in csv if requested
    if write_csv:

        ## initial year info
        year_ini = float(Out_post.year.sel(year=years_csv).mean('year').values)
        year_str = str(int(year_ini)) if year_ini.is_integer() else '{:.1f}'.format(year_ini)
        
        ## save posterior
        Par1_post.astype(np.float32).to_dataframe().to_csv(folder_calib + '/../Par_' + name + '.csv')
        Ini_post = Var_post.sel(year=years_csv).mean('year')
        Ini_post.astype(np.float32).to_dataframe().to_csv(folder_calib + '/../Ini_' + year_str + '_' + name + '.csv')
        
        ## save prior
        if save_prior:
            Par1_prior.astype(np.float32).to_dataframe().to_csv(folder_calib + '/../Par_' + name + '_prior.csv')
            Ini_prior = Var_prior.sel(year=years_csv).mean('year')
            Ini_prior.astype(np.float32).to_dataframe().to_csv(folder_calib + '/../Ini_' + year_str + '_' + name + '_prior.csv')

    ## add info on constraints
    Con_stat =[Con0,
        xr.Dataset({var: xr.DataArray(Con0[var].attrs['is_sum'], coords={'stat':['is_sum']}, dims=['stat']) for var in Con0}), 
        xr.Dataset({var: xr.DataArray(Con0[var].attrs['is_mean'], coords={'stat':['is_mean']}, dims=['stat']) for var in Con0}), 
        xr.Dataset({var: xr.DataArray(Con0[var].attrs['is_diff'], coords={'stat':['is_diff']}, dims=['stat']) for var in Con0}), 
        xr.Dataset({var: xr.DataArray(str(Con0[var].attrs['period']), coords={'stat':['period']}, dims=['stat']) for var in Con0}), 
        xr.Dataset({var: xr.DataArray(Con0[var].attrs['units'], coords={'stat':['units']}, dims=['stat']) for var in Con0})]
    xr.concat(Con_stat, dim='stat').drop([var for var in Con0 if var in ignore_constr]).to_dataframe().to_csv(folder_calib + '/../Con0_' + name + '.csv')

