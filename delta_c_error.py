# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:16:55 2020

@author: Azumi Mamiya
"""

import matplotlib.pyplot as plt
import numpy as np

payoff={'T':1.5, 'R':1.0, 'P':0.0, 'S':-0.5}
T,R,P,S=payoff['T'],payoff['R'],payoff['P'],payoff['S']

def my_plot(eta, w):# generate the figure and setting of the figure
    plt.ylim(0, 1)
    plt.xlim(0,1/3)
    plt.xlabel('$\epsilon+\\xi$',
               fontsize = 18)
    plt.ylabel('$\delta_c$',
               fontsize = 18)
    
    plt.rcParams["xtick.direction"] = "in"  
    plt.rcParams["ytick.direction"] = "in"  
    plt.plot(eta, w, '-', 
             #color = '#1f77b4', 
             markersize = 1)
    #plt.savefig('../data/figure/nondiscount_hao/{}.png'.format(n))

def func1(eta):
    mu = 1 - eta
    RE = R * mu + S * eta
    TE = T * mu + P * eta
    SE = S * mu + R * eta
    PE = P * mu + T * eta
    
    # val = 1 + (TE-SE)*(1-delta+2*delta*eta)/(delta*(TE-PE)-delta*eta*(RE-SE+TE-PE)-(TE-RE))
    val = (TE-RE) / (mu*(RE-PE)-eta*(TE-SE)+TE-RE)
    #if val < 0:
    #    val=50
    return val

def func2(eta):
    mu = 1 - eta
    RE = R * mu + S * eta
    TE = T * mu + P * eta
    SE = S * mu + R * eta
    PE = P * mu + T * eta
    
    # val = 1 + (TE-SE)*(1-delta+2*delta*eta)/(delta*(mu*(RE-PE)-eta*(TE-SE))-(1-delta)*(PE-SE))
    val = (PE-SE) / (mu*(RE-PE)-eta*(TE-SE)+PE-SE)
    #if val < 0:
    #    val=50
    return val

if __name__ == "__main__":
    x_list = [i for i in np.linspace(0, 0.3, 100)]
    y_list1 = []
    for x in x_list:
        y_list1.append(func2(x))
            
    my_plot(x_list, y_list1)
    plt.grid()
    # plt.savefig('./equa_kappa_range.pdf')