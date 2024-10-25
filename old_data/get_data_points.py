import math
from distance import dist
import numpy as np

####### Finding potential data points
# find <Od> spatial mean of ssh over all data points within L, Od array and subset of latlon
def get_data_points(hydr_data, latg, long, L1, L2):
    
    subset1 = [ll for ll, val in hydr_data.items() 
                    if dist(ll[0], ll[1], latg, long) <= L1
                    and math.isnan(val['ssh']) is False]
    Od1 = {ll:hydr_data[ll]['ssh'] for ll in subset1}
    print(len(Od1))

    subset2 = [ll for ll, val in hydr_data.items() 
                    if dist(ll[0], ll[1], latg, long) <= L2 
                    and math.isnan(val['ssh']) is False]
    Od2 = {ll:hydr_data[ll]['ssh'] for ll in subset2}
    print(len(Od2))
    
    return subset1, Od1, subset2, Od2


####### Outlier elimination
# for the potential data points, neglect those outside 2 SD
def outlier_elimination(Od1, Od2):
    
    Od1 = {ll:val for ll,val in Od1.items() 
        if abs(val - np.mean(list(Od1.values()))) < 2*np.std(list(Od1.values()))}
    
    Od2 = {ll:val for ll,val in Od2.items() 
        if abs(val - np.mean(list(Od2.values()))) < 2*np.std(list(Od2.values()))}
    
    return Od1, Od2
