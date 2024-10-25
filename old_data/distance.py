import geopy

####### distance between the target grid point and all other points
def dist(lat1, lon1, lat2, lon2): 
    return geopy.distance.geodesic((lat1, lon1), (lat2, lon2)).km