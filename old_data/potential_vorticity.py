import math
import numpy as np
from bathymetric_depth import nearest_depth


####### potential vorticity, planetary vorticity is f/H param
def PV(latd, lond, latg, long): # can be ll of data, grid for Cdg and data, data for Cdd

    # coriolis force of data points and grid point
    fd = 2 * 7.29e-5 * np.sin(np.deg2rad(latd))
    fg = 2 * 7.29e-5 * np.sin(np.deg2rad(latg))

    # nearest bathmetric depth for given data or grid point
    Zg = -nearest_depth(latg, long)
    Zd = -nearest_depth(latd, lond)
                                        
    PV = abs(fd/Zd - fg/Zg) / math.sqrt((fd/Zd)**2 + (fg/Zg)**2)
    
    return PV