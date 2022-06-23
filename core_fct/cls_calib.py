"""
Copyright IIASA (International Institute for Applied Systems Analysis), 2019-2022
Contributor(s): Thomas Gasser (gasser@iiasa.ac.at)

This software is a computer program whose purpose is simulate a large number of global climate (and environmental) change scenarios.

This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free software.  You can use, modify and/ or redistribute the software under the terms of the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL "http://www.cecill.info". 

As a counterpart to the access to the source code and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty and the software's author, the holder of the economic rights, and the successive licensors have only limited liability. 

In this respect, the user's attention is drawn to the risks associated with loading, using, modifying and/or developing or reproducing the software by the user in light of its specific status of free software, that may mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data to be ensured and,  more generally, to use and operate it in the same conditions as regards security. 

The fact that you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
"""


import pymc3 as pm
import numpy as np
import theano.tensor as tt

from scipy.integrate import solve_ivp
from pymc3.distributions import Continuous
from pymc3.distributions.continuous import Normal


##################################################
##################################################

## subclass defining a normalised AR1
## note: the 'unit' condition seems to be useless
class my_AR1(Continuous):

    ## initiliased with normalised innovation
    def __init__(self, k, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.k = tt.as_tensor_variable(k)
        self.tau = 1 / (1 - k ** 2)
        self.mode = tt.as_tensor_variable(0.)

    ## log-probability based on regular AR1
    def logp(self, x):

        x_im1 = x[:-1]
        x_i = x[1:]
        
        innov_like = Normal.dist(0., tau=self.tau).logp(x_i - self.k * x_im1)
        unit_like = Normal.dist(0., 1.).logp(x)

        return tt.sum(innov_like) + tt.sum(unit_like)


##################################################
##################################################

## create subclass of the original Theano Op
## change: uses the newest scipy ode integrator
## note: not used in the end, but kept here just in case
class DiffEqOp(pm.ode.DifferentialEquation):

    ## new initialisation
    ## change: ignore t0; option to choose method; option to raise error when fail
    def __init__(self, func, times, n_states, n_theta, method=None, raise_when_fail=True):
        super().__init__(func, times, n_states=n_states, n_theta=n_theta)
        self._augmented_times = self._augmented_times[1:]
        self.method = method
        self.raise_when_fail = raise_when_fail

    ## new simulate function
    ## change: use solve_ivp; size of time vector; catch integration failure
    def _simulate(self, y0, theta):

        ## initial condition
        s0 = np.concatenate([y0, self._sens_ic])

        ## perform the integration
        out = solve_ivp(fun=lambda t, Y, p: self._system(Y, t, p), t_span=(self._augmented_times[0], self._augmented_times[-1]), y0=s0, 
            t_eval=self._augmented_times, args=(np.concatenate([y0, theta]),), **dict((self.method != None) * [('method', self.method)]))

        ## take output if successful integration
        if out.success: 
            sol = out.y.transpose()
        ## or raise error (default behaviour)
        elif self.raise_when_fail: 
            raise RuntimeError('ODE solver failed')
        ## or move on with zeroed output (not recommended!)
        else:
            print('*solver failed*')
            sol = np.zeros((len(self._augmented_times), len(s0)))

        ## main solution and reshaped sensitivities
        y = sol[:, :self.n_states]
        sens = sol[:, self.n_states:].reshape(self.n_times, self.n_states, self.n_p)

        ## return
        return y, sens


## example of use of the above Theano Op
'''
from core_fct.cls_calib import DiffEqOp
''''''
ode_model = DiffEqOp(lambda X, t, c: PF(t, X, c[:len(PF.Par_name)], 
        For=[c[len(PF.Par_name)+n*len(For0.year):len(PF.Par_name)+(n+1)*len(For0.year)] for n in range(len(PF.For_name))], 
        autonomous=False, tensor=True, expost=False), 
    times=For0.year.values - For0.year.values[0], 
    n_states=len(PF.Var_name), 
    n_theta=len(PF.Par_name) + len(For0.year) * len(PF.For_name))
''''''
with model: 
    var_prog = ode_model(PF.Ini_dflt(param), param + [forcing[var][t] for var in PF.For_name for t in range(n_year)])
'''

