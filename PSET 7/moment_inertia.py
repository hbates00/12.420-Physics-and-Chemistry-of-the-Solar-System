# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 21:29:16 2017

@author: Haley
"""

from matplotlib import pyplot as plt
import numpy as np

plt.axhline(y=1.0, xmin=0.005, xmax=1, linewidth=1.5, color = '#9467bd', label = 'C/MR^2 = 0.4' )
plt.axvline(x=1.00001, ymin=0, ymax = .989, linewidth=1.5, color='#9467bd') 
plt.legend()

gamma = np.array([0.375, 0.35, 0.325, 0.3, 0.275, 0.25, 0.225])

x = np.arange(0, 1, 0.01)


for i in gamma:
    
    A = (2*x**5.0 - 5 * (x**3.0) * i + 5 *i - 2) / (x**3.0 * (2 * x**2.0 - 5 * i))
    
    rhoavg = 1 / (A * x**3.0 + 1 - x**3)
    Arhoavg = rhoavg * A
    
    plt.plot(Arhoavg, rhoavg, label = "C/MR^2 = %s" %  i)
    plt.legend()
    plt.xlabel('Arho/rho_avg')
    plt.ylabel('rho/rho_avg')
    plt.title('Core Density vs Crust Density for Different Moments of Inertia')
    plt.xlim(0.99, 3.0)
    plt.ylim(0.05, 1.01)
    
