import numpy as np
import pandas as pd
# from helpers import *
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Filter for the section 2 - TPD FS
# df = df[
#     (round(df['Latitude'], 3) == 82.494) & (round(df['Longitude'], 3) == 17.904) |
#     (round(df['Latitude'], 3) == 82.620) & (round(df['Longitude'], 3) == 14.578) |
#     (round(df['Latitude'], 3) == 82.721) & (round(df['Longitude'], 3) == 11.148) |
#     (round(df['Latitude'], 3) == 82.794) & (round(df['Longitude'], 3) == 7.635) |
#     (round(df['Latitude'], 3) == 82.840) & (round(df['Longitude'], 3) == 4.064) |
#     (round(df['Latitude'], 3) == 82.858) & (round(df['Longitude'], 3) == 0.461) |
#     (round(df['Latitude'], 3) == 82.847) & (round(df['Longitude'], 3) == -3.145) |
#     (round(df['Latitude'], 3) == 82.809) & (round(df['Longitude'], 3) == -6.727) |
#     (round(df['Latitude'], 3) == 82.742) & (round(df['Longitude'], 3) == -10.257)
# ]

# bg section center
# df = df[
#     (round(df['Latitude'], 3) == 72.724) & (round(df['Longitude'], 3) == -131.033) |
#     (round(df['Latitude'], 3) == 73.064) & (round(df['Longitude'], 3) == -132.036) |
#     (round(df['Latitude'], 3) == 73.399) & (round(df['Longitude'], 3) == -133.079) |
#     (round(df['Latitude'], 3) == 73.728) & (round(df['Longitude'], 3) == -134.164) |
#     (round(df['Latitude'], 3) == 74.05) & (round(df['Longitude'], 3) == -135.293) | 
#     (round(df['Latitude'], 3) == 74.366) & (round(df['Longitude'], 3) == -136.468) | 
#     (round(df['Latitude'], 3) == 74.676) & (round(df['Longitude'], 3) == -137.69) |
#     (round(df['Latitude'], 3) == 74.977) & (round(df['Longitude'], 3) == -138.962) |
#     (round(df['Latitude'], 3) == 75.271) & (round(df['Longitude'], 3) == -140.285) | 
#     (round(df['Latitude'], 3) == 75.556) & (round(df['Longitude'], 3) == -141.661) | 
#     (round(df['Latitude'], 3) == 75.833) & (round(df['Longitude'], 3) == -143.092) | 
#     (round(df['Latitude'], 3) == 76.1) & (round(df['Longitude'], 3) == -144.577) | 
#     (round(df['Latitude'], 3) == 76.358) & (round(df['Longitude'], 3) == -146.12) |
#     (round(df['Latitude'], 3) == 76.605) & (round(df['Longitude'], 3) == -147.721) |
#     (round(df['Latitude'], 3) == 76.841) & (round(df['Longitude'], 3) == -149.38) |
#     (round(df['Latitude'], 3) == 77.066) & (round(df['Longitude'], 3) == -151.098) 
# ]


# print(df)

# dx = 50 * 1e3  # Segment width in meters (assume 50 km spacing)
# # conversion_factor = 1e-6  # Convert m³/s to Sv

# annual_transport = {}
# for year, group in merged_df.groupby('Year'):
#     vg = group['velocity'].values  # Meridional velocity (m/s)
#     segment_transport = vg * dx  # Transport for each segment in m³/s
#     # total_transport = np.sum(segment_transport) * conversion_factor # Convert to Sv
#     total_transport = np.sum(segment_transport) # Convert to Sv
#     # annual_transport[year] = total_transport/(14*dx)
#     annual_transport[year] = total_transport/(6*dx)

# annual_transport_df = pd.DataFrame.from_dict(annual_transport, orient='index', columns=['Transport (Sv)'])
# print(annual_transport_df)

# dummy = {year:np.nan for year in range(1980, 2018) if year not in list(annual_transport.keys())}
# annual_transport = annual_transport | dummy
# print(annual_transport)


# fig, ax1 = plt.subplots(figsize=(10, 6))
# ax1.plot(range(1980, 2019), annual_transport_df.values, marker='o', linestyle='-', color='red')
# ax1.set_xlabel('Year', fontsize=12)
# ax1.set_ylabel('Transport [Sv]', fontsize=12)
# plt.title('Meridional Transport along center of BG')
# ax1.tick_params(axis='y')
# # ax1.legend()
# plt.show()
# plt.savefig('bg_1980_2018_transport.png', dpi=300)

# Example DataFrame
# df = pd.read_csv('grd_dh_2012_01.csv')
# df['ug'] = np.full(len(df), np.nan)
# df['vg'] = np.full(len(df), np.nan)

###### checking nearest neigh function
# target_lat = 77.81912879787636
# target_lon = -141.38977105221358
# lat = round(target_lat, 3)
# lon = round(target_lon,3)

# Find the 4 nearest coordinates within the 49-51 km range
# nearest_coords = find_nearest_coords(target_lat, target_lon, df)
# print(nearest_coords)

# Print results
# idx = df.loc[(df['Latitude']==target_lat ) & (df['Longitude']==target_lon)].index[0]
# print(idx, df['X_meters'][idx], df['Y_meters'][idx], df['Surf_DH'][idx])

# nidx = coord[0]
# print(nidx, df['X_meters'][nidx], df['Y_meters'][nidx], coord[3], df['Surf_DH'][nidx])

# dx_plus = df.loc[(df['X_meters']==df['X_meters'][idx]+50000 ) & (df['Y_meters']==df['Y_meters'][idx])].index[0]
# dh_plus = df['Surf_DH'][dx_plus]
# dx_minus = df.loc[(df['X_meters']==df['X_meters'][idx]-50000 ) & (df['Y_meters']==df['Y_meters'][idx])].index[0]
# dh_minus = df['Surf_DH'][dx_minus]

# grad_x = (dh_plus-dh_minus)/(2*50000)

# dy_plus = df.loc[(df['X_meters']==df['X_meters'][idx] ) & (df['Y_meters']==df['Y_meters'][idx]+50000)].index[0]
# dh_plus = df['Surf_DH'][dy_plus]
# dy_minus = df.loc[(df['X_meters']==df['X_meters'][idx] ) & (df['Y_meters']==df['Y_meters'][idx]-50000)].index[0]
# dh_minus = df['Surf_DH'][dy_minus]

# grad_y = (dh_plus-dh_minus)/(2*50000)

# g = get_g(target_lat)
# f = coriolis(target_lat)
# print(grad_x, grad_y)

# ug = -g*grad_y/f
# vg = g*grad_x/f
# print(ug, vg)
    

##### combine different mapping vars to single file
# df = pd.read_csv('grd_dh_2012_01.csv')

# hfw_df = pd.read_csv('grd_dh_2012_01_hfw.csv')
# df['hFW'] = hfw_df['hFW']
# df['hFW_err'] = hfw_df['hFW_err']

# disoh_df = pd.read_csv('grd_dh_2012_01_disoh.csv')
# df['D_Siso'] = disoh_df['D_Siso']
# df['D_Siso_err'] = disoh_df['D_Siso_err']

# # print(df.tail)
# df.to_csv('grd_dh_2012_01.csv', index=False)


###### combine all months in a year to one - 2012
import os
# print(os.listdir('./'))
csv_files = [file for file in os.listdir('./results/gsc/') if file.startswith('grd_gsc_') and file not in ['grid_50km_nplaea.csv']]
print(csv_files)

# Initialize an empty list to hold DataFrames
dataframes = []

# Loop through the CSV files and read each into a DataFrame
for csv_file in csv_files:
    file_path = os.path.join('./results/gsc/', csv_file)
    df = pd.read_csv(file_path)
    dataframes.append(df)

merged_df = pd.concat(dataframes, ignore_index=True)
# merged_df = merged_df.drop(columns=['D_Siso', 'hFW'])
print(merged_df.head)

output_file = "grd_gsc_1980_2018.csv"
merged_df.to_csv(output_file, index=False)


####### distribution of UDAH 2011-2018
# import pandas as pd
# import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap

# # Load data
# df = pd.read_csv('data_points.csv')

# # Convert 'Datetime' column to datetime format
# df['Datetime'] = pd.to_datetime(df['Datetime'])

# # Define time periods
# time_periods = {
#     # "1980-1989": ('1980-01-01', '1989-12-31'),
#     # "1990-1999": ('1990-01-01', '1999-12-31'),
#     # "2000-2009": ('2000-01-01', '2009-12-31'),
#     "2011-2018": ('2011-01-01', '2018-12-31')
# }

# # Create a 2x2 subplot
# fig, ax = plt.subplots(1, 1, figsize=(12, 12))
# # axes = axes.flatten()  # Flatten for iteration

# # Iterate over time periods and plot
# # for ax, (label, (start, end)) in zip(axes, time_periods.items()):
#     # Filter data for the time period
# df_filtered = df[(df['Datetime'] >= '2011-01-01') & (df['Datetime'] <= '2018-12-31')]
# df_filtered['Month'] = df_filtered['Datetime'].dt.month

# # Set up Basemap
# m = Basemap(
#     projection='npstere',
#     boundinglat=70,
#     lon_0=0,
#     resolution='l',
#     round=True,
#     ax=ax
# )
    
# # Draw coastlines and gridlines
# m.drawcoastlines()
# m.drawparallels(range(70, 91, 5), labels=[1, 0, 0, 0])
# m.drawmeridians(range(-180, 180, 60), labels=[0, 0, 0, 0])

# # Convert lat/lon to map projection
# x, y = m(df_filtered['Longitude'].values, df_filtered['Latitude'].values)

# # Scatter plot
# sc = m.scatter(x, y, c=df_filtered['Month'], cmap='viridis', s=30, edgecolor='black', alpha=0.8)

# # Add title to each subplot
# ax.set_title(f"{2011-2018}", fontsize=12)

# # Add a colorbar
# cbar = fig.colorbar(sc, ax=ax, shrink=0.5)
# cbar.set_label("Month")
# cbar.set_ticks(range(1, 13))

# # Add overall title
# fig.suptitle("Month-Wise Distribution of Profiles (2011-2018)", fontsize=16)

# # Adjust layout
# # plt.tight_layout(rect=[0, 0.03, 1, 0.95])

# # Show plot
# plt.show()



