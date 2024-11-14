import numpy as np
import csv
import concurrent.futures
import itertools
from objmap_opt import *
import pandas as pd
from datetime import datetime

gp = pd.read_csv('grid_points.csv')
    
# Function to apply to each combination
def proc(lat_lon_pair):
    lat, lon = lat_lon_pair
    t = datetime(2011, 1, 1)
    
    (Og, Og_err) = objmap(lat, lon, t)
    print(t.date(), lat, lon, Og, Og_err)
    return t.date(), lat, lon, Og, Og_err  
            
# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(proc, [(row['Latitude'], row['Longitude']) for _, row in list(gp.iterrows())]))

output_file = "results/grd_dh_2011_01_wnn.csv"  # another run with distance sorted poiints selected within L1 L2

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Datetime", "Latitude", "Longitude", "Dynamic_Height", "DH_error"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")

