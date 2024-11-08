import numpy as np
import csv
import concurrent.futures
import itertools
from objmap import *
import pandas as pd
from scipy.spatial import KDTree
from datetime import datetime

# run objective mapping in batches for whole grid
lats = np.linspace(60, 90, 121)
lons = np.linspace(-179.25, 180, 480)
# ts = pd.date_range(start=datetime(2012, 1, 1), end=datetime(2012, 1, 1), freq='MS')
# print(lats, lons,)

combinations = np.array(list(itertools.product(lats, lons)))
print(combinations)

landmask_data = pd.read_csv("data/landsea_04.csv")
landmask_coordinates = list(zip(landmask_data['Latitude'], landmask_data['Longitude']))
kdtree = KDTree(landmask_coordinates)

def is_land_or_ocean(lat, lon):
    _, idx = kdtree.query((lat, lon))
    depth = landmask_data.iloc[idx]['Bottom_Standard_level']
    return 'land' if depth == 1 else 'ocean'

# Function to apply to each combination
def process_combination(lat_lon_pair):
    lat, lon = lat_lon_pair
    t = datetime(2011, 1, 1)
    
    if lon>-20 and lon<15 and lat<82: # remove Fram strait from mapping
        return t, lat, lon, np.nan, np.nan
    else:
        if is_land_or_ocean(lat, lon) == 'ocean':  
            (Og, Og_err) = objmap(lat, lon, t)
            return t, lat, lon, Og, Og_err  
        
    return t, lat, lon, np.nan, np.nan
    
# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_combination, combinations))

output_file = "results/grd_dh_2011_01.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Datetime", "Latitude", "Longitude", "Dynamic_Height", "DH_error"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")

# Function to apply to each combination
def process_combination(lat_lon_pair):
    lat, lon = lat_lon_pair
    t = datetime(2011, 2, 1)
    
    if lon>-20 and lon<15 and lat<82: # remove Fram strait from mapping
        return t, lat, lon, np.nan, np.nan
    else:
        if is_land_or_ocean(lat, lon) == 'ocean':  
            Og, Og_err = objmap(lat, lon, t)
            return t, lat, lon, Og, Og_err  
        
    return t, lat, lon, np.nan, np.nan
    
# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_combination, combinations))

output_file = "results/grd_dh_2011_02.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Datetime", "Latitude", "Longitude", "Dynamic_Height", "DH_error"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")

# Function to apply to each combination
def process_combination(lat_lon_pair):
    lat, lon = lat_lon_pair
    t = datetime(2011, 3, 1)
    
    if lon>-20 and lon<15 and lat<82: # remove Fram strait from mapping
        return t, lat, lon, np.nan, np.nan
    else:
        if is_land_or_ocean(lat, lon) == 'ocean':  
            Og, Og_err = objmap(lat, lon, t)
            return t, lat, lon, Og, Og_err  
        
    return t, lat, lon, np.nan, np.nan
    
# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_combination, combinations))

output_file = "results/grd_dh_2011_03.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Datetime", "Latitude", "Longitude", "Dynamic_Height", "DH_error"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")

# Function to apply to each combination
def process_combination(lat_lon_pair):
    lat, lon = lat_lon_pair
    t = datetime(2011, 4, 1)
    
    if lon>-20 and lon<15 and lat<82: # remove Fram strait from mapping
        return t, lat, lon, np.nan, np.nan
    else:
        if is_land_or_ocean(lat, lon) == 'ocean':  
            Og, Og_err = objmap(lat, lon, t)
            return t, lat, lon, Og, Og_err  
        
    return t, lat, lon, np.nan, np.nan
    
# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_combination, combinations))

output_file = "results/grd_dh_2011_04.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Datetime", "Latitude", "Longitude", "Dynamic_Height", "DH_error"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")

# Function to apply to each combination
def process_combination(lat_lon_pair):
    lat, lon = lat_lon_pair
    t = datetime(2011, 5, 1)
    
    if lon>-20 and lon<15 and lat<82: # remove Fram strait from mapping
        return t, lat, lon, np.nan, np.nan
    else:
        if is_land_or_ocean(lat, lon) == 'ocean':  
            Og, Og_err = objmap(lat, lon, t)
            return t, lat, lon, Og, Og_err  
        
    return t, lat, lon, np.nan, np.nan
    
# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_combination, combinations))

output_file = "results/grd_dh_2011_05.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Datetime", "Latitude", "Longitude", "Dynamic_Height", "DH_error"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")

# Function to apply to each combination
def process_combination(lat_lon_pair):
    lat, lon = lat_lon_pair
    t = datetime(2011, 6, 1)
    
    if lon>-20 and lon<15 and lat<82: # remove Fram strait from mapping
        return t, lat, lon, np.nan, np.nan
    else:
        if is_land_or_ocean(lat, lon) == 'ocean':  
            Og, Og_err = objmap(lat, lon, t)
            return t, lat, lon, Og, Og_err  
        
    return t, lat, lon, np.nan, np.nan
    
# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_combination, combinations))

output_file = "results/grd_dh_2011_06.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Datetime", "Latitude", "Longitude", "Dynamic_Height", "DH_error"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")


# Function to apply to each combination
def process_combination(lat_lon_pair):
    lat, lon = lat_lon_pair
    t = datetime(2011, 7, 1)
    
    if lon>-20 and lon<15 and lat<82: # remove Fram strait from mapping
        return t, lat, lon, np.nan, np.nan
    else:
        if is_land_or_ocean(lat, lon) == 'ocean':  
            Og, Og_err = objmap(lat, lon, t)
            return t, lat, lon, Og, Og_err  
        
    return t, lat, lon, np.nan, np.nan
    
# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_combination, combinations))

output_file = "results/grd_dh_2011_07.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Datetime", "Latitude", "Longitude", "Dynamic_Height", "DH_error"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")

# Function to apply to each combination
def process_combination(lat_lon_pair):
    lat, lon = lat_lon_pair
    t = datetime(2011, 8, 1)
    
    if lon>-20 and lon<15 and lat<82: # remove Fram strait from mapping
        return t, lat, lon, np.nan, np.nan
    else:
        if is_land_or_ocean(lat, lon) == 'ocean':  
            Og, Og_err = objmap(lat, lon, t)
            return t, lat, lon, Og, Og_err  
        
    return t, lat, lon, np.nan, np.nan
    
# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_combination, combinations))

output_file = "results/grd_dh_2011_08.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Datetime", "Latitude", "Longitude", "Dynamic_Height", "DH_error"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")

# Function to apply to each combination
def process_combination(lat_lon_pair):
    lat, lon = lat_lon_pair
    t = datetime(2011, 9, 1)
    
    if lon>-20 and lon<15 and lat<82: # remove Fram strait from mapping
        return t, lat, lon, np.nan, np.nan
    else:
        if is_land_or_ocean(lat, lon) == 'ocean':  
            Og, Og_err = objmap(lat, lon, t)
            return t, lat, lon, Og, Og_err  
        
    return t, lat, lon, np.nan, np.nan
    
# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_combination, combinations))

output_file = "results/grd_dh_2011_09.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Datetime", "Latitude", "Longitude", "Dynamic_Height", "DH_error"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")

# Function to apply to each combination
def process_combination(lat_lon_pair):
    lat, lon = lat_lon_pair
    t = datetime(2011, 10, 1)
    
    if lon>-20 and lon<15 and lat<82: # remove Fram strait from mapping
        return t, lat, lon, np.nan, np.nan
    else:
        if is_land_or_ocean(lat, lon) == 'ocean':  
            Og, Og_err = objmap(lat, lon, t)
            return t, lat, lon, Og, Og_err  
        
    return t, lat, lon, np.nan, np.nan
    
# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_combination, combinations))

output_file = "results/grd_dh_2011_10.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Datetime", "Latitude", "Longitude", "Dynamic_Height", "DH_error"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")

# Function to apply to each combination
def process_combination(lat_lon_pair):
    lat, lon = lat_lon_pair
    t = datetime(2011, 11, 1)
    
    if lon>-20 and lon<15 and lat<82: # remove Fram strait from mapping
        return t, lat, lon, np.nan, np.nan
    else:
        if is_land_or_ocean(lat, lon) == 'ocean':  
            Og, Og_err = objmap(lat, lon, t)
            return t, lat, lon, Og, Og_err  
        
    return t, lat, lon, np.nan, np.nan
    
# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_combination, combinations))

output_file = "results/grd_dh_2011_11.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Datetime", "Latitude", "Longitude", "Dynamic_Height", "DH_error"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")

# Function to apply to each combination
def process_combination(lat_lon_pair):
    lat, lon = lat_lon_pair
    t = datetime(2011, 12, 1)
    
    if lon>-20 and lon<15 and lat<82: # remove Fram strait from mapping
        return t, lat, lon, np.nan, np.nan
    else:
        if is_land_or_ocean(lat, lon) == 'ocean':  
            Og, Og_err = objmap(lat, lon, t)
            return t, lat, lon, Og, Og_err  
        
    return t, lat, lon, np.nan, np.nan
    
# Use ThreadPoolExecutor to parallelize computation and collect results
results = []
with concurrent.futures.ThreadPoolExecutor() as executor:
    results = list(executor.map(process_combination, combinations))

output_file = "results/grd_dh_2011_12.csv"

with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Datetime", "Latitude", "Longitude", "Dynamic_Height", "DH_error"])  
    writer.writerows(results)  

print(f"Results written to {output_file}")
