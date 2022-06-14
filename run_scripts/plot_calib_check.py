import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


##################################################
##################################################

## OPTIONS
name = 'v1'
name_extra = []
folder_calib = 'internal_data/pyMC_calib/'

## general figure options
plt.style.use('seaborn-colorblind')


##################################################
## LOAD: calib. results
##################################################

## prior
with xr.open_dataset(folder_calib + 'Par_' + name + '_prior.nc') as TMP: Par_prior = TMP.load()
with xr.open_dataset(folder_calib + 'Par2_' + name + '_prior.nc') as TMP: Par2_prior = TMP.load()
with xr.open_dataset(folder_calib + 'Var_' + name + '_prior.nc') as TMP: Var_prior = TMP.load()
with xr.open_dataset(folder_calib + 'Con_' + name + '_prior.nc') as TMP: Con_prior = TMP.load()
par2_list = [par for par in Par2_prior]
Par_prior = xr.merge([Par_prior, Par2_prior]); del Par2_prior

## posterior
with xr.open_dataset(folder_calib + 'Par_' + name + '.nc') as TMP: Par_post = TMP.load()
with xr.open_dataset(folder_calib + 'Par2_' + name + '.nc') as TMP: Par2_post = TMP.load()
with xr.open_dataset(folder_calib + 'Var_' + name + '.nc') as TMP: Var_post = TMP.load()
with xr.open_dataset(folder_calib + 'Con_' + name + '.nc') as TMP: Con_post = TMP.load()
Par_post = xr.merge([Par_post, Par2_post]); del Par2_post

## default constraints (= obs.)
Con_obs = xr.Dataset.from_dataframe(pd.read_csv(folder_calib + 'Con0_' + name + '.csv', index_col=0))

## extra posterior (prior ignored!)
Par_extra, Par2_extra, Var_extra, Con_extra = [], [], [], []
for name_ in name_extra:
    with xr.open_dataset(folder_calib + 'Par_' + name_ + '.nc') as TMP: Par_extra.append(TMP.load())
    with xr.open_dataset(folder_calib + 'Par2_' + name_ + '.nc') as TMP: Par2_extra.append(TMP.load())
    with xr.open_dataset(folder_calib + 'Var_' + name_ + '.nc') as TMP: Var_extra.append(TMP.load())
    with xr.open_dataset(folder_calib + 'Con_' + name_ + '.nc') as TMP: Con_extra.append(TMP.load())
Par_extra = [xr.merge([Par, Par2]) for Par, Par2 in zip(Par_extra, Par2_extra)]; del Par2_extra


#########################################
## FIGURE: distrib. param.
#########################################

## figure info
par_list = [par for par in Par_post if 'config' in Par_post[par].dims]
Par_list = [Par_prior, Par_post] + Par_extra
Par_name = ['prior', 'post.'] + ['extra'+str(n+1) for n in range(len(name_extra))]
color_list = ['0.2', 'C0'] + ['C'+str(n+1) for n in range(len(name_extra))]

## figure options
n_col = 5

## figure
fig = plt.figure()
fig.set_figwidth(7.2, forward=True)
fig.set_figheight(8.7, forward=True)

## subplots
for n_par, par in enumerate(par_list):
    ax = plt.subplot(int(np.ceil(len(par_list)/ n_col)), n_col, n_par + 1)

    ## distributions
    for n, Par in enumerate(Par_list):
        plt.hist(Par[par].values, density=True, bins=50, histtype='step', color=color_list[n], alpha=0.8, label=Par_name[n])

    ## axis
    if par in ['rho_T', 'rho_CO2']:
        plt.title(par + ' (1)', fontsize='x-small', va='top')
    else:
        plt.title(par + ' (' + Par[par].units.replace('-1', r'$^{-1}$').replace('-2', r'$^{-2}$').replace('deg', r'$^\circ$') + ')', fontsize='x-small', va='top')
    plt.xticks(fontsize='x-small')
    plt.minorticks_on()
    plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

    ## legend
    if n_par==0:
        plt.legend(loc=0, ncol=1, fontsize='xx-small', frameon=False, borderpad=0., labelspacing=0.4, handlelength=0.6, handletextpad=0.4)

##
plt.subplots_adjust(left=0.02, bottom=0.03, right=0.98, top=0.98, wspace=0.10, hspace=0.60)


#########################################
## FIGURE: distrib. obs.
#########################################

## figure info
Con_list = [Con_prior, Con_post] + Con_extra
Con_name = ['prior', 'post.'] + ['extra'+str(n+1) for n in range(len(name_extra))]
color_list = ['0.2', 'C0'] + ['C'+str(n+1) for n in range(len(name_extra))]

## figure options
n_col = 4

## figure
fig = plt.figure()
fig.set_figwidth(7.2, forward=True)
fig.set_figheight(5.2, forward=True)

## subplots
for n_var, var in enumerate(Con_prior):
    ax = plt.subplot(int(np.ceil(len(Con_prior) / n_col)), n_col, n_var + 1)

    ## observation
    mu, sigma = [float(Con_obs[var].loc[stat].values) for stat in ['mean', 'std']]
    xx = np.linspace(mu - 4 * sigma, mu + 4 * sigma, 1000)
    yy = 1/np.sqrt(2*np.pi) / sigma * np.exp(-0.5 * (xx - mu)**2 / sigma**2)
    plt.plot(xx, yy, lw=2, ls='--', color='k', alpha=0.8, label='obs.')

    ## distributions
    for n, Con in enumerate(Con_list):
        if var in Con:
            plt.hist(Con[var].values, density=True, bins=50, histtype='step', color=color_list[n], alpha=0.8, label=Con_name[n])
    
    ## axis
    units = str(Con_obs[var].loc['units'].values).replace('-1', r'$^{-1}$').replace('-2', r'$^{-2}$')
    is_mean, is_sum, is_diff = [eval(str(Con_obs[var].loc[stat].values)) for stat in ['is_mean', 'is_sum', 'is_diff']]
    taken = is_mean * 'mean' + is_sum * 'sum' + is_diff * 'diff' + ': {0}-{1}'.format(*[str(eval(str(Con_obs[var].loc['period'].values))[n]) for n in range(2)])
    plt.title('{0} ({1})\n[{2}]'.format(var, units, taken), fontsize='x-small', va='top')
    plt.xticks(fontsize='x-small')
    plt.minorticks_on()
    plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

    ## legend
    if n_var==0:
        plt.legend(loc=0, ncol=1, fontsize='xx-small', frameon=False, borderpad=0., labelspacing=0.4, handlelength=0.6, handletextpad=0.4)

##
plt.subplots_adjust(left=0.02, bottom=0.05, right=0.98, top=0.93, wspace=0.15, hspace=0.60)


#########################################
## FIGURE: correl. param.
#########################################

## figure options
par_list = [par for par in Par_post if 'config' in Par_post[par].dims and par not in par2_list]
Par_list = [Par_prior, Par_post] + Par_extra
Par_name = ['prior', 'post.'] + ['extra'+str(n) for n in range(len(name_extra))]

## nb of parameters per group
grp_nb = [6, 7, 6, 15, 4, 1]
grp_name = ['climate', 'sea level', 'ocean C', 'land C', 'pf. C', '']

## figure
fig = plt.figure()
fig.set_figwidth(7.2, forward=True)
fig.set_figheight(4.2, forward=True)

## subplots
for n, Par in enumerate(Par_list):
    
    ## matrix
    ax = plt.subplot(1, len(Par_list), n+1)
    mapp = ax.matshow(np.corrcoef(np.array([Par[par].values for par in par_list])), cmap='RdBu_r', vmin=-1, vmax=1) # RdBu_r / bwr / seismic / coolwarm
    
    ## color bar
    cb = plt.colorbar(mapp, fraction=0.10, pad=0.03, shrink=0.9, ticks=np.arange(-1, 1+0.01, 0.2), orientation="horizontal")
    cb.ax.tick_params(labelsize='x-small')

    ## axis
    plt.xticks(range(len(par_list)), par_list, fontsize='xx-small', rotation=90)
    plt.yticks(range(len(par_list)), par_list, fontsize='xx-small')

    ## legend
    if n+1==len(Par_list):
        grp_val = np.append(0, np.cumsum(grp_nb)) / len(par_list)
        for nx, x0, x1 in zip(np.arange(len(grp_val)-1)[:-1], grp_val[:-1][:-1], grp_val[1:][:-1]):
            plt.annotate('', xy=(x0, 1.12 + 0.0 * (nx % 2)), xycoords='axes fraction', xytext=(x1, 1.12 + 0.0 * (nx % 2)), arrowprops=dict(arrowstyle="-", linewidth=0.5, color='k'))
            plt.figtext((x0+x1)/2, 1.12 + 0.0 * (nx % 2) + 0.01, grp_name[nx], transform=ax.transAxes, fontsize='xx-small', ha='center', va='bottom')

##
plt.subplots_adjust(left=0.06, bottom=0.01, right=0.98, top=0.92, wspace=0.15, hspace=0.15)

