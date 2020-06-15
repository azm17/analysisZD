# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:16:55 2020

@author: Azumi Mamiya
"""

import matplotlib.pyplot as plt
import numpy as np

payoff = {'T':1.5, 'R':1.0, 'P':0.0, 'S':-0.5}
payoff = {'T':2.0, 'R':1.2, 'P':0.1, 'S':-0.5}

T,R,P,S = payoff['T'],payoff['R'],payoff['P'],payoff['S']

def my_plot(eta, w, name):# generate the figure and setting of the figure
    plt.ylim(-1.0, 2.0)
    plt.xlim(0, 20)
    plt.xlabel('$\chi$',
               fontsize = 18)
    plt.ylabel('$\kappa$',
               fontsize = 18)
    
    plt.rcParams["xtick.direction"] = "in"  
    plt.rcParams["ytick.direction"] = "in"  
    plt.plot(eta, w, '-', 
             # color = '#1f77b4', 
             markersize = 1,
             label = name)
    #plt.savefig('../data/figure/nondiscount_hao/{}.png'.format(n))

def kappa1(chi, delta, p0):
    return R - (1 / delta - 1) * (1 - p0) * (T - R + (T-S) / (chi - 1))

def kappa2(chi, delta, p0):
    return P + (1 - delta) * p0 * (T - P + (T-S) / (chi - 1))

def kappa3(chi, delta, p0):
    return R - (1 - delta) * (1 - p0) * (R - S + (T-S) / (chi - 1))

def kappa4(chi, delta, p0):
    return P + (1 / delta - 1) * p0 * (P - S + (T-S) / (chi - 1))

if __name__ == "__main__":
    delta = 0.55
    p0 = 0.5
    
    chi_c1 = 1 + (1 - delta) * (T - S) / (delta * (T - P) - (T - R))
    chi_c2 = 1 + (1 - delta) * (T - S) / (delta * (R - S) - (P - S))
    chi_c = max(chi_c1, chi_c2)
    
    chi_list = [i for i in np.linspace(chi_c, 20, 10) if i <= 20 and i >= 1]
    kappa_list1 = []; kappa_list2 = []
    kappa_list3 = []; kappa_list4 = []
    
    for chi in chi_list:
        kappa_list1.append(kappa1(chi, delta, p0))
        kappa_list2.append(kappa2(chi, delta, p0))
        
        kappa_list3.append(kappa3(chi, delta, p0))
        kappa_list4.append(kappa4(chi, delta, p0))
        
    my_plot(chi_list, kappa_list1, 'kappa1')
    my_plot(chi_list, kappa_list2, 'kappa2')
    
    my_plot(chi_list, kappa_list3, 'kappa3')
    my_plot(chi_list, kappa_list4, 'kappa4')
    
    plt.grid()
    plt.legend()
    # plt.savefig('./equa_kappa_range.pdf')