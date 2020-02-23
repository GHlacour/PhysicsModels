import numpy as np

N=900
R=8.3144598e-3 # J/mol/K Gas constant

# Non-Interacting
E310=2327.70
E300=2252.56
CVni=(E310-E300)/10/N/R

# Attractive
E310=1869.83
E300=1792.65
CVat=(E310-E300)/10/N/R

# Repulsive
E310=2390.51
E300=2313.21
CVre=(E310-E300)/10/N/R

CVig=0.5
print([CVni,CVat,CVre,CVig])
