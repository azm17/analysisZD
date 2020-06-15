# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:16:55 2020

@author: Azumi Mamiya
"""

import matplotlib.pyplot as plt
import numpy as np

payoff={'T':1.5, 'R':1.0, 'P':0.0, 'S':-0.5}
#payoff={'T':2, 'R':1.0, 'P':0.0, 'S':-0.5}

T,R,P,S=payoff['T'],payoff['R'],payoff['P'],payoff['S']

def my_plot(eta, w):# generate the figure and setting of the figure
    plt.ylim(-1.0, 2.0)
    plt.xlim(0,20)
    plt.xlabel('$\chi$',
               fontsize = 18)
    plt.ylabel('$\kappa$',
               fontsize = 18)
    
    plt.rcParams["xtick.direction"] = "in"  
    plt.rcParams["ytick.direction"] = "in"  
    plt.plot(eta, w, '-', 
             color = '#1f77b4', 
             markersize = 1)
    #plt.savefig('../data/figure/nondiscount_hao/{}.png'.format(n))

def kappa1(chi,eta):
    mu = 1 - eta
    RE = R * mu + S * eta
    TE = T * mu + P * eta
    SE = S * mu + R * eta
    PE = P * mu + T * eta
    
    return RE - eta*((RE-SE)+chi*(TE-RE))/((mu-eta)*(chi-1))

def kappa2(chi,eta):
    mu = 1 - eta
    # RE = R * mu + S * eta
    TE = T * mu + P * eta
    SE = S * mu + R * eta
    PE = P * mu + T * eta
    
    return PE + eta * (TE-PE+chi*(PE-SE))/((mu-eta)*(chi-1))

def kappa3(chi,eta):
    mu = 1 - eta
    RE = R * mu + S * eta
    TE = T * mu + P * eta
    SE = S * mu + R * eta
    PE = P * mu + T * eta
    
    return RE + mu*((RE-SE) + chi*(TE-RE))/((mu-eta)*(chi-1))

def kappa4(chi,eta):
    mu = 1 - eta
    # RE = R * mu + S * eta
    TE = T * mu + P * eta
    SE = S * mu + R * eta
    PE = P * mu + T * eta
    
    return PE - mu * (TE-PE + chi*(PE-SE))/((mu-eta)*(chi-1))

if __name__ == "__main__":
    eta = 0.99
    mu = 1 - eta
    RE = R * mu + S * eta
    TE = T * mu + P * eta
    SE = S * mu + R * eta
    PE = P * mu + T * eta
    
    if mu > eta:
        chi_c = 1 + 2 * eta * (TE - SE)/(mu*(RE-PE)-eta*(TE-SE))
    elif eta == 0.5:
        chi_c = None
    else:
        chi_c = 1 + 2 * mu * (TE - SE)/(eta*(RE-PE)-mu*(TE-SE))
    
    x_list = [i for i in np.linspace(chi_c, 20, 500) if i <= 20 and i>=1]
    y_list1 = []; y_list2 = []
    
    print(x_list)
    for x in x_list:
        if eta < 0.5:
            y_list1.append(kappa1(x,eta))
            y_list2.append(kappa2(x,eta))
        elif eta == 0.5:
            y_list1.append(None)
            y_list2.append(None)
        else:
            y_list1.append(kappa3(x,eta))
            y_list2.append(kappa4(x,eta))
        
    my_plot(x_list, y_list1)
    my_plot(x_list, y_list2)
    plt.grid()
    # plt.savefig('./equa_kappa_range.pdf')