# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:16:55 2020

@author: Azumi Mamiya
"""

import matplotlib.pyplot as plt
import numpy as np

payoff={'T':1.5, 'R':1.0, 'P':0.0, 'S':-0.5}
TE, RE, PE, SE = payoff['T'],payoff['R'],payoff['P'],payoff['S']

def my_plot(eta, chi,delta):# generate the figure and setting of the figure
    plt.ylim(0, 5)
    plt.xlim(-0.03, 0.3)
    
    plt.yticks([1,2,3,4,5])
    plt.xticks([0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
    
    plt.xlabel('$\epsilon+\\xi$',
               fontsize = 18)
    plt.ylabel('$\chi_c$',
               fontsize = 18)
    
    plt.rcParams["xtick.direction"] = "in"  
    plt.rcParams["ytick.direction"] = "in"  
    plt.plot(eta, chi, '-', 
             #color = '#1f77b4', 
             markersize = 1,
             label = delta)
    plt.legend(title='$\delta$')
    #plt.savefig('../data/figure/nondiscount_hao/{}.png'.format(n))

def func1(eta, delta):
    mu = 1 - eta
    val = 1 + (1 - delta + 2 * delta * eta) * (TE - SE) \
               / (delta * (mu * (RE - PE) - eta * (TE - SE)) - (1 - delta) * (TE - RE))

    return val

def func2(eta, delta):
    mu = 1 - eta
    val = 1 + (1 - delta + 2 * delta * eta) * (TE - SE) \
               / (delta * (mu * (RE - PE) - eta * (TE - SE)) - (1 - delta) * (PE - SE))

    return val

if __name__ == "__main__":
    for delta in [1.0,0.9,0.8]:
        eta_list = []; chi_list = []
        for eta in [i for i in np.linspace(0, 0.3, 100)]:
            mu = 1 - eta
            
            chi = max(func1(eta, delta), func2(eta, delta))
            
            delta1 = (TE - RE) / (mu * (RE - PE) - eta * (TE - SE) + TE -RE)
            delta2 = (PE - SE) / (mu * (RE - PE) - eta * (TE - SE) + PE -SE)
            
            delta_c = max(delta1, delta2)
            if delta < delta_c:
                break
            else:
                eta_list.append(eta)
                chi_list.append(chi)
        
        my_plot(eta_list, chi_list, delta)
    
    
    eta_list2 = [0.1, 0.1, 0.1]
    chi_list2 = [1.571, 1.966, 2.565]
    
    plt.plot(eta_list2, chi_list2 , color='k',marker='o', linewidth=0)
    
    eta_list3 = [0, 0, 0]
    chi_list3 = [1.0, 1.235, 1.571]
    plt.plot(eta_list3, chi_list3 , color='k', marker='o', linewidth=0)


    
    plt.grid()
    # plt.savefig('./equa_kappa_range.pdf')