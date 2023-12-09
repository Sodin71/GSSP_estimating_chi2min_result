# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 22:02:42 2023

@author: JeongMJ (mjj0055@hanmail.net, chungbuk national university)

This program is for estimating results from the grid-search result file of GSSP.
ver. 0.0

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import time


def Setting_Design_figure(xrange='auto',yrange = 'auto',ytickloc='left',xtickexi='on', yinv = 'off'):
    if xrange!='auto': # for adj. axes 1~3
        plt.gca().set_xlim(xrange)
    if yrange!='auto':plt.gca().set_ylim(yrange)
    if ytickloc=='right':
        plt.gca().yaxis.tick_right()
        plt.gca().yaxis.set_label_position("right")
    if xtickexi=='off':plt.gca().xaxis.set_ticklabels([])
    if yinv == 'on': plt.gca().invert_yaxis()
    
    plt.minorticks_on()
    plt.gca().xaxis.set_ticks_position('both')
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator(2))
    plt.gca().yaxis.set_ticks_position('both')
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator(2))
    plt.gca().tick_params(which='major', direction = 'in',length=8)
    plt.gca().tick_params(which='minor', direction = 'in',length=5)
    
    plt.gcf().tight_layout()

def get_min_chi2(data, sel_para_label, chi2_label):
    # finding the minimum sigma for each fixed value
    
    val_list = np.sort(np.unique(data[sel_para_label]))
    results = []
    
    for val in val_list:
        chi2 = np.min(data[chi2_label][data[sel_para_label] == val])
        results.append((val,chi2))
        
    results = np.array(results, dtype= [(sel_para_label,np.float64),(chi2_label,np.float64)])
    return results

    
def get_val_err_with_polyfitting(x, y, ylim, dimension = 3, savename = 'result', x_ref = [], y_ref = [], xlabel = 'value', ylabel = r'\chi^2'):
    # estimating the value having minimum sigma and its error range corresponding to 1-sigma
    
    # polynomial fitting
    p = np.polyfit(x, y, dimension)
    fx = np.poly1d(p)
    
    # estimate the X value in minimum chi2
    dfx1 = np.polyder(fx)
    root = np.roots(dfx1)
    root_real = root[np.isreal(root)]
    solX = np.real(root_real[(root_real < np.max(x)) & (root_real > np.min(x))])
    
    # x range, smaller than 1-sigma, to estimate the error
    gx = fx - ylim
    root_ylim = np.roots(gx)
    root_ylim_real = root_ylim[np.isreal(root_ylim)]
    possible_X_range = np.real(root_ylim_real[(root_ylim_real < np.max(x)) & (root_ylim_real > np.min(x))])
    diff_X = possible_X_range - solX
    
    if len(diff_X) == 2: error = [np.min(diff_X), np.max(diff_X)]  
    elif len(diff_X) ==1:
        if diff_X[0] >= 0: error = [np.nan, diff_X[0]]
        else: error = [diff_X[0], np.nan]
    else: 
        error = [np.nan, np.nan] 
        print('       WARNING!!: cannot estimate the errors')
    
    # create a result figure
    x_c_add = (np.max(x)-np.min(x))*0.1
    x_c = np.linspace(np.min(x)-x_c_add, np.max(x)+x_c_add,200)
    y_c = fx(x_c)
    plt.figure()
    plt.plot(x_ref,y_ref,'.', c = '0.8', ms = 1)
    plt.plot(x,y,'o')
    plt.plot(x_c, y_c,'-k')
    plt.axvline(solX, c ='red', ls = '-', label = 'minimum point')
    plt.axvline(possible_X_range[0], c ='gray', ls = '-')
    plt.axvline(possible_X_range[1], c ='gray', ls = '-')
    plt.axhline(ylim, c ='gray', ls = ':', label = r'1-$\sigma$')
    plt.title(r'%s = $%.4f_{%.4f}^{+%.4f}$'%(xlabel, solX, error[0], error[1]))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ylim([np.min(y)-np.std(y)*0.2, np.max(y)+np.std(y)*0.2])
    plt.savefig(savename+'.png')
    Setting_Design_figure()
    plt.show()
    
    
    # save the results
    
    f = open('%s.txt'%savename,'w')
    f.write('# %s\n'%time.strftime("%Y-%m-%d %H:%M:%S"))
    f.write('# dimension  = %.0f\n'%dimension)
    f.write('# %s = %.4f [%.4f, +%.4f]\n\n'%(xlabel, solX, error[0], error[1]))
    f.write('# 1. minimu chi2 value\n')
    f.write('# %s   chi2\n'%xlabel)
    for k in range(len(x)):
        f.write('%.5f %.5f\n'%(x[k], y[k]))
    f.write('\n')
    f.write('# 2. polynomial fitting result\n')
    f.write('# %s   chi2\n'%xlabel)
    for k in range(len(x_c)):
        f.write('%.10f %.5f\n'%(x_c[k], y_c[k]))
    
    f.close()
    
    return solX, error



def main(input_file, work_parameters = []):
    
    chi2_lab = 'reduced_chi2'
    
    # load data from [input_file]
    names = ['MH','Teff','logg','micro_turbulence','vsini','chi2_inter','contin_factor','reduced_chi2','chi2_1sigma']
    data = np.genfromtxt(input_file, names = names)   
    
    # get chi2_1sigma value
    chi2_1sig = data['chi2_1sigma'][0]
    
    # run for each parameter
    result = {} # make directory for saving results
    for param in work_parameters:
        
        # find minimum chi2
        MINchi_val = get_min_chi2(data, param, chi2_lab)
        if len(MINchi_val[param]) == 1:            
            print('  - %s = %.4f (fixed)'%(param, MINchi_val[param]))
            print('       WARNING!!: %s was fixed > result figure was not created'%(param))
            result[param] = {'val': MINchi_val[param], 'error': [0.,0.]}
            continue
        
        # get estimated value and its error
        param_val, param_err = get_val_err_with_polyfitting(MINchi_val[param], MINchi_val[chi2_lab], chi2_1sig, 
                                                      savename = param, x_ref = data[param], y_ref = data[chi2_lab], 
                                                      xlabel = param, ylabel = r'$\chi^2$')
        print('  - %s = %.4f [%.4f, +%.4f]'%(param, param_val, param_err[0], param_err[1]))
        # save the results in dictionary
        result[param] = {'val': param_val, 'error': param_err}
    

if __name__ == "__main__":
    
    # 1. put result file name
    filename = 'Chi2_table.dat'
    
    # 2. put the parameter. It should be the list! select: ['MH','Teff','logg','micro_turbulence','vsini']
    work_parameters = ['MH','Teff','logg','micro_turbulence','vsini'] 
    
    # run ...  Polynomial fitting to estimate parameters is performed with dimension 3 (default).    
    # If you want to change the dimension, add dimension = any value as follow:
    # main(filename, work_parameters, dimension = 5)
    
    main(filename, work_parameters)