"""
Copyright IIASA (International Institute for Applied Systems Analysis), 2019-2021
Contributor(s): Thomas Gasser (gasser@iiasa.ac.at)

This software is a computer program whose purpose is simulate a large number of global climate (and environmental) change scenarios.

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
"""


import numpy as np
import xarray as xr

from scipy.integrate import solve_ivp


## wrapper class
class WrapModel:
    
    ## initialization
    def __init__(self, Var_name, Var2_name, Par_name, For_name, Ini_dflt, v_linear, d_Var):
        ## safety checks
        assert len(Ini_dflt(len(Par_name)*[.1])) == len(Var_name)
        assert len(v_linear(len(Par_name)*[.1])) == len(Var_name)
        assert len(d_Var(None, len(Var_name)*[.1], len(Par_name)*[.1], len(For_name)*[.1], autonomous=True, tensor=False, expost=False)) == len(Var_name)
        assert len(d_Var(None, len(Var_name)*[.1], len(Par_name)*[.1], len(For_name)*[.1], autonomous=True, tensor=False, expost=True)) == len(Var2_name) + len(Var_name)
        ## initialize list of names, and functions
        self.Var_name, self.Var2_name, self.Par_name, self.For_name = Var_name, Var2_name, Par_name, For_name
        self.Ini_dflt, self.v_linear, self.d_Var = Ini_dflt, v_linear, d_Var


    ## call shortcut
    def __call__(self, *args, **kwargs):
        return self.d_Var(*args, **kwargs)


    ## run with solver (only one extra dim okay)
    def run_solver(self, param, forc, init=None, method='LSODA'):
        
        ## turn into ndarray
        param = {var: np.array(param[var]).squeeze() for var in self.Par_name}
        forc = {var: np.array(forc[var]).squeeze() for var in self.For_name}
        if init is not None:
            init = {var: np.array(init[var]).squeeze() for var in self.Var_name}

        ## check dimensions and lengths
        if not all([param[var].ndim <= 1 for var in param.keys()]):
            raise ValueError('all parameters should have ndim <= 1')
        if not all([1 <= forc[var].ndim <= 2 for var in forc.keys()]):
            raise ValueError('all forcings should have 1 <= ndim <= 2')
        if init is not None:
            if not all([init[var].ndim <= 1 for var in init.keys()]):
                raise ValueError('all initial states should have ndim <= 1')

        ## get minimal lengths
        n_year = min([forc[var].shape[0] for var in forc.keys()])
        n_config_param = min([param[var].shape[0] if param[var].ndim == 1 else np.inf for var in param.keys()])
        n_config_forc = min([forc[var].shape[1] if forc[var].ndim == 2 else np.inf for var in forc.keys()])
        n_config_init = min([init[var].shape[0] if init[var].ndim == 1 else np.inf for var in init.keys()]) if init is not None else np.inf
        if not all(np.isinf([n_config_param, n_config_forc, n_config_init])):
            n_config = min([n_config_param, n_config_forc, n_config_init])
        else:
            n_config = 1

        ## create correct input
        param_ok = np.array([param[var][:n_config] if param[var].ndim == 1 else param[var] * np.ones(n_config) for var in param.keys()])
        forc_ok = np.array([forc[var][:n_year, :n_config] if forc[var].ndim == 2 else forc[var][:n_year, np.newaxis] * np.ones(n_config)[np.newaxis, :] for var in forc.keys()])
        if init is None:
            init_ok = np.array([self.Ini_dflt(param_ok[:, n]) for n in range(n_config)]).transpose()
        else:
            init_ok = np.array([init[var][:n_config] if init[var].ndim == 1 else init[var] * np.ones(n_config) for var in init.keys()])

        ## create output
        var_prog = np.zeros((len(self.Var_name), n_year, n_config))
        var_diag = np.zeros((len(self.Var2_name) + len(self.Var_name), n_year, n_config))

        ## loop on config
        print('*solving*')
        for n in range(n_config):
            print('config: ' + str(n+1) + '/' + str(n_config), end='\r')
       
            ## call solver
            sol_solver = solve_ivp(lambda t, y: self.d_Var(t, y, param_ok[..., n], forc_ok[..., n], autonomous=False, tensor=False, expost=False), 
                t_span=(0, n_year-1), y0=init_ok[..., n], t_eval=np.arange(n_year), method=method)

            ## get results if convergence
            if sol_solver.success: 
                var_prog[..., n] = sol_solver.y
                var_diag[..., n] = self.d_Var(None, var_prog[..., n], param_ok[..., n], forc_ok[..., n], autonomous=True, tensor=False, expost=True)
            else: 
                print('solver failed! (config = ' + str(n) + ')')
                var_prog[..., n] *= np.nan
                var_diag[..., n] *= np.nan

        ## pack in dictionnaries
        var_prog = {var: var_prog[n, ...].squeeze() for n, var in enumerate(self.Var_name)}
        var_diag = {'d_' * (var in self.Var_name) + var: var_diag[n, ...].squeeze() for n, var in enumerate(self.Var2_name + self.Var_name)}

        ## return
        print('')
        return {'t': np.arange(n_year), **var_prog, **var_diag}
    

    ## run with solving scheme for big xarrays
    def run_xarray(self, Par, For, Ini=None, time_axis='year', scheme='imex', nt=4):

        ## get time
        time0 = For[time_axis][0]
        time = For[time_axis]

        ## get input in proper order
        Par_ok = [Par[var] for var in self.Par_name]
        For_ok = [For[var] for var in self.For_name]
        Ini_ok = self.Ini_dflt(Par_ok) if Ini is None else [Ini[var] for var in self.Var_name]
        v_lin = self.v_linear(Par_ok)

        ## function solving next step
        if scheme == 'ExpInt': next_step = lambda X, dX, v, dt: X * np.exp(-v*dt) + np.expm1(-v*dt) / -v * (dX + v*X)
        elif scheme == 'imex': next_step = lambda X, dX, v, dt: (X + dt * (dX + v*X)) / (1 + v*dt)
        elif scheme == 'ex': next_step = lambda X, dX, v, dt: X + dt * dX 

        ## initialise
        Var_t = Ini_ok.copy()
        Var_out = [xr.Dataset({var: val for var, val in zip(self.Var_name, Var_t)}) .assign_coords(**{time_axis:time[0]}) .expand_dims(time_axis, 0)]

        ## solve system
        print('*solving*')
        for t in (time[1:] - time0):
            print(time_axis + ': ' + str(time[t].values), end='\r')
            
            ## get time increment
            dt = float(time[t] - time[t-1])

            ## subloop & function call
            for _ in range(nt):
                Var_t = [next_step(X, dX, v, dt/nt) for X, dX, v in zip(Var_t, self.d_Var(t, Var_t, Par_ok, For_ok, autonomous=False, tensor=False, expost=False), v_lin)]
        
            ## save outputs
            Var_out.append(xr.Dataset({var: val for var, val in zip(self.Var_name, Var_t)}) .assign_coords(**{time_axis:time[t]}) .expand_dims(time_axis, 0))

        ## merge outputs
        Var_out = xr.concat(Var_out, dim=time_axis)
        Var_out = Var_out.astype(np.float32)

        ## one last call for derivatives
        Var_out2 = self.d_Var(time - time0, [Var_out[var] for var in self.Var_name], Par_ok, For_ok, autonomous=False, tensor=False, expost=True)
        Var_out2 = xr.Dataset({'d_' * (var in self.Var_name) + var: 0.*time + val for var, val in zip(self.Var2_name + self.Var_name, Var_out2)})
        Var_out2 = Var_out2.astype(np.float32)

        ## return
        print('')
        return xr.merge([Var_out, Var_out2])

