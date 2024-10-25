from distance import dist
import random
import numpy as np

# partition only when data points are more than 60 else take full data as it is
def pick_data_points(hydr_data, subset1, latg, long, Od1, Od2):
    if len(subset1) > 60:
        # pick random 1/3 points here   
        subset = random.sample(subset1, 20)
        # print(len(subset))

        # pick 1/3 points from L=600km radius with weight to closest points 
        # sort ll array based on distance closest and select 1/3 points
        dist1 = [(ll, dist(ll[0], ll[1], latg, long)) for ll in list(Od1.keys())]
        sorted1 = sorted(dist1, key=lambda x: x[1])
        sorted1 = [point for point, distance in sorted1 if point not in subset]
        subset = subset + sorted1[0:20]
        # print(len(subset))

        # pick remaining 1/3 points from L=300km radius 
        dist2 = [(ll, dist(ll[0], ll[1], latg, long)) for ll in list(Od2.keys())]
        sorted2 = sorted(dist2, key=lambda x: x[1])
        sorted2 = [point for point, distance in sorted2 if point not in subset]
        subset = subset + sorted2[0:20]
        print(len(subset))
    else:
        subset = subset1
    
    Od = np.array([hydr_data[ll]['ssh'] for ll in subset])
    n = len(Od)
    # print('Od:', Od, n)
    print('n', n)
    print('mean', np.mean(Od))
    
    return subset, Od, n
