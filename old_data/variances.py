import numpy as np
from scipy.spatial import cKDTree

####### signal variance <s2>
def signal_variance(Od, n):
    s2 = np.sum([(dval - np.mean(Od))**2 for dval in Od])/n
 
    print('s2:', s2)
    return s2

####### noise variance <n2>
def noise_variance(hydr_data, subset, Od, n):
    # Finding noise between a data point with its nearest neighbour using k-d tree
    tree = cKDTree(np.array(subset))
    nn = []
    for i, point in enumerate(subset):
        dist, idx = tree.query(point, k=2)  # k=2 because the closest point to a point is itself
        nearest_idx = idx[1]  # idx[0] is the point itself, idx[1] is the nearest neighbor
        nearest_neighbor = subset[nearest_idx]
        nn.append(nearest_neighbor)
        # print(nearest_neighbor)
        
    n2 = np.sum([(Od[idx] - hydr_data[nn[idx]]['ssh'])**2 for idx in range(len(subset))])/(2*n)
    print('n2:', n2)
    return n2