import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from helpers import get_g, coriolis, D_mat
from scipy.stats import pearsonr

# Load data
file_path = 'grd_gsc_1980_2018.csv'  # Replace with your file path
df = pd.read_csv(file_path)

######## section 1 - TPD LR side start
df1 = df[(round(df['Latitude'], 3) == 83.432) & (round(df['Longitude'], 3) == 123.829)]
df2 = df[(round(df['Latitude'], 3) == 83.071) & (round(df['Longitude'], 3) == 156.552)]

###### GSC annual mean
merged_df = pd.merge(df1, df2, on="Datetime", how="inner")
merged_df = merged_df.drop(columns=['Depth_x', 'Depth_y', 'Surf_DH_err_x', 'Surf_DH_err_y', 'Numobs_x', 'Numobs_y', 'Y_meters_x', 'X_meters_x', 'Y_meters_y', 'X_meters_y', 'ug_x', 'vg_x', 'ug_y', 'vg_y'])
dist = D_mat((83.071, 156.552), target=(83.432, 123.829))*1000  # in meters
merged_df['grad_dh'] = (merged_df['Surf_DH_y'] - merged_df['Surf_DH_x'])/dist
# print(merged_df)

g, f = get_g(83), coriolis(83)
merged_df['velocity'] = (g/f) * (merged_df['grad_dh'])
print(merged_df)

merged_df['Year'] = pd.to_datetime(merged_df['Datetime']).dt.year
merged_df = merged_df.sort_values(by='Year')

# Group by 'Year' and calculate the mean of 'velocity'
annual_mean_velocity1 = merged_df.groupby('Year')['velocity'].mean().reset_index()
print(annual_mean_velocity1)

plt.plot(annual_mean_velocity1['Year'], annual_mean_velocity1['velocity']*100, marker='o', color='red', label='TPD start')
# plt.title('Annual Mean Velocity at TPD Start')
# plt.xlabel('Year')
# plt.ylabel('Mean Velocity (cm/s)')
# plt.show()

####### DH annual mean 
# lat_min, lat_max = 81.0, 82.0
# lon_min, lon_max = 120.0, 160.0
# region_data = df[(df['Latitude'] >= lat_min) & (df['Latitude'] <= lat_max) &
#                    (df['Longitude'] >= lon_min) & (df['Longitude'] <= lon_max) ]
# region_data['Year'] = pd.to_datetime(region_data['Datetime']).dt.year

# annual_mean_dh = region_data.groupby('Year')['Surf_DH'].mean()
# print(max(annual_mean_dh))

# # region_data['GSC'] = np.sqrt(region_data['ug']**2 + region_data['vg']**2)
# # annual_mean_gsc = region_data.groupby('Year')['GSC'].mean()

# plt.plot(annual_mean_velocity['Year'], annual_mean_dh*100, marker='o')
# plt.show()


####### section 2 - TPD side
df1 = df[(round(df['Latitude'], 3) == 80.330) & (round(df['Longitude'], 3) == -169.614)]
df2 = df[(round(df['Latitude'], 3) == 85.540) & (round(df['Longitude'], 3) == -157.012)]

###### GSC annual mean
merged_df = pd.merge(df1, df2, on="Datetime", how="inner")
merged_df = merged_df.drop(columns=['Depth_x', 'Depth_y', 'Surf_DH_err_x', 'Surf_DH_err_y', 'Numobs_x', 'Numobs_y', 'Y_meters_x', 'X_meters_x', 'Y_meters_y', 'X_meters_y', 'ug_x', 'vg_x', 'ug_y', 'vg_y'])
dist = D_mat((80.330, -169.614), target=(85.540, -157.012))*1000  # in meters
merged_df['grad_dh'] = (merged_df['Surf_DH_y'] - merged_df['Surf_DH_x'])/dist
# print(merged_df)

g, f = get_g(83), coriolis(83)
merged_df['velocity'] = -(g/f) * (merged_df['grad_dh'])
print(merged_df)

merged_df['Year'] = pd.to_datetime(merged_df['Datetime']).dt.year
merged_df = merged_df.sort_values(by='Year')

# Group by 'Year' and calculate the mean of 'velocity'
annual_mean_velocity2 = merged_df.groupby('Year')['velocity'].mean().reset_index()
print(annual_mean_velocity2)

plt.plot(annual_mean_velocity2['Year'], annual_mean_velocity2['velocity']*100, marker='o', color='blue', label='TPD side (MB)')


####### section 3 - TPD end
df1 = df[(round(df['Latitude'], 3) == 86.568) & (round(df['Longitude'], 3) == -50.317)]
df2 = df[(round(df['Latitude'], 3) == 85.542) & (round(df['Longitude'], 3) == -5.049)]

###### GSC annual mean
merged_df = pd.merge(df1, df2, on="Datetime", how="inner")
merged_df = merged_df.drop(columns=['Depth_x', 'Depth_y', 'Surf_DH_err_x', 'Surf_DH_err_y', 'Numobs_x', 'Numobs_y', 'Y_meters_x', 'X_meters_x', 'Y_meters_y', 'X_meters_y', 'ug_x', 'vg_x', 'ug_y', 'vg_y'])
dist = D_mat((86.568, -50.317), target=(85.542, -5.049))*1000  # in meters
merged_df['grad_dh'] = (merged_df['Surf_DH_y'] - merged_df['Surf_DH_x'])/dist
# print(merged_df)

g, f = get_g(86), coriolis(86)
merged_df['velocity'] = (g/f) * (merged_df['grad_dh'])
print(merged_df)

merged_df['Year'] = pd.to_datetime(merged_df['Datetime']).dt.year
merged_df = merged_df.sort_values(by='Year')

# Group by 'Year' and calculate the mean of 'velocity'
annual_mean_velocity3 = merged_df.groupby('Year')['velocity'].mean().reset_index()
print(annual_mean_velocity3)

plt.plot(annual_mean_velocity3['Year'], annual_mean_velocity3['velocity']*100, marker='o', color='green', label='TPD end')
plt.title('Annual Mean Velocity normal to TPD sections')
plt.xlabel('Year')
plt.ylabel('Mean Velocity (cm/s)')



df = pd.read_csv('ao_1980_2018.csv')
plt.plot(df['Year'], df['AO'], marker='.', color='black', label='AO index')

# print(df['AO'][10:], annual_mean_velocity2['velocity']*100)
correlation, p_value = pearsonr(df['AO'][10:], annual_mean_velocity1['velocity']*100)    # annual=0.776, summer=0.230, winter=0.830
print(f"Pearson Correlation Coefficient: {correlation:.3f}")
print(f"P-value: {p_value:.3e}")

correlation, p_value = pearsonr(df['AO'][10:], annual_mean_velocity2['velocity']*100)    # annual=0.776, summer=0.230, winter=0.830
print(f"Pearson Correlation Coefficient: {correlation:.3f}")
print(f"P-value: {p_value:.3e}")

correlation, p_value = pearsonr(df['AO'], annual_mean_velocity3['velocity']*100)    # annual=0.776, summer=0.230, winter=0.830
print(f"Pearson Correlation Coefficient: {correlation:.3f}")
print(f"P-value: {p_value:.3e}")

# plt.legend()
# plt.show()

###### correlation with lag
lag = 2  # 1-year lag
ao_series = df['AO'][10:][lag:]  # Remove last `lag` values
velocity_series = annual_mean_velocity2['velocity'][:-lag]*100  # Remove first `lag` values

correlation, p_value = pearsonr(ao_series, velocity_series)    
print(f"Pearson Correlation Coefficient: {correlation:.3f}")
print(f"P-value: {p_value:.3e}")


##### computing lag
# Normalize the series (optional for better interpretation)
x = df['AO'][10:]
y = annual_mean_velocity2['velocity']*100
x = (x - np.mean(x)) / np.std(x)
y = (y - np.mean(y)) / np.std(y)

# Compute cross-correlation
lags = np.arange(-len(x) + 1, len(y))  # Possible lags
ccf = np.correlate(x, y, mode='full')

# Find the lag with the maximum cross-correlation
max_lag = lags[np.argmax(ccf)]
print(f"Maximum correlation occurs at lag: {max_lag}")
