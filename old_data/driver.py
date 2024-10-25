import numpy as np
import csv
import concurrent.futures
import itertools
from old.objmap_grid import objmap

# run objective mapping in batches for whole grid
# lats = np.linspace(60, 88, 113)
# lons = np.linspace(-179.25, 180, 480)

# for beaufort gyre
lats = np.linspace(70, 80, 41)
lons = np.linspace(-155, -125, 61)

combinations = np.array(list(itertools.product(lats, lons)))
print(combinations)
# print(len(combinations))

# Function to apply to each combination
def process_combination(lat_lon_pair):
    lat, lon = lat_lon_pair
    result = objmap(lat, lon)
    print(lat, lon, result)
    return lat, lon, result  

# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_combination, combinations))

output_file = "grid_dh.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Latitude", "Longitude", "Dynamic_Height"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")

