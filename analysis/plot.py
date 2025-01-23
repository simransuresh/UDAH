from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

df = pd.read_csv('grd_gsc_1980_2018.csv')
df = df[df['Datetime']>='1990-01-01']
X = df['X_meters'].values
Y = df['Y_meters'].values
lons = df['Longitude'].values
lats = df['Latitude'].values
dhs = df['Surf_DH'].values
ug = df['ug'].values
vg = df['vg'].values

print(df.head)

#### Plot dh
# m = Basemap(projection='nplaea', boundinglat=70, lon_0=0, resolution='l', round=True)

# fig, ax = plt.subplots(figsize=(8, 8))
# m.drawcoastlines(linewidth=1.2)
# m.drawparallels([80], labels=[True], linewidth=0.5, color='gray', fontsize=6)
# m.drawmeridians([-180, -90, 0, 90], labels=[True, True, True, True], linewidth=0.5, color='gray', fontsize=6)

# x, y = m(lons, lats)
# sc = m.scatter(x, y, s=20, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.02
# q = m.quiver(x, y, ug, vg, scale=scale, color='black', width=0.003)

# cbar = plt.colorbar(sc, ax=ax, shrink=0.7)
# cbar.set_label("Dynamic height (m)")
# plt.show()


fig = plt.figure(figsize=(6, 4))
projection = ccrs.NorthPolarStereo()
ax = plt.axes(projection=projection)

# Set the extent to focus on the region of interest
ax.set_extent([-180, 180, 70, 90], crs=ccrs.PlateCarree())

# Add map features
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
ax.add_feature(cfeature.LAND, color='lightgrey')

# Scatter plot for mean DH
sc = ax.scatter(lons, lats, c=dhs, cmap='YlGnBu', s=180, transform=ccrs.PlateCarree(), label='Mean DH')
cb = plt.colorbar(sc, orientation='vertical', pad=0.05, shrink=0.7)
cb.set_label('Mean Dynamic Height (m)', fontsize=12)

# Quiver plot for geostrophic currents (u_g, v_g)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.03
# quiver = ax.quiver(lons, lats, ug, vg, transform=ccrs.PlateCarree(), scale=scale, color='black', width=0.003)

# Add labels and title
# plt.savefig('bg_saga_2015.png', dpi=300)
plt.show()

# plt.savefig('Sep2015_seas.png', dpi=300)


##### plot gsc 
# m = Basemap(projection='nplaea', boundinglat=70, lon_0=0, resolution='l', round=True)

# fig, ax = plt.subplots(figsize=(8, 8))
# m.drawcoastlines(linewidth=1.2)
# m.drawparallels([80], labels=[True], linewidth=0.5, color='gray', fontsize=6)
# m.drawmeridians([-180, -90, 0, 90], labels=[True, True, True, True], linewidth=0.5, color='gray', fontsize=6)

# x, y = m(lons, lats)
# sc = m.scatter(x, y, s=30, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)

# cbar = plt.colorbar(sc, ax=ax, shrink=0.7)
# cbar.set_label("Dynamic height (m)")

# # calculate scale based on ug vg values for quiver
# velocity_magnitude = np.sqrt(ug**2 + vg**2)
# max_velocity = np.nanmax(velocity_magnitude)
# scale = max_velocity / 0.05
# print(max_velocity, scale)

# Q = m.quiver(x, y, ug, vg, scale=scale, width=0.003)
# qk = plt.quiverkey(Q, 0.1, 0.1, 0.1, '10 cm/s', labelpos='W')

# plt.show()
# plt.savefig('2012_08_gsc.png', dpi=300)


# #### plot all fw/diso of all data points
# df = pd.read_csv('data_500m_fw.csv')
# lats = df['Latitude'].values
# lons = df['Longitude'].values
# hfw = df['hFW'].values
# diso = df['D_Siso'].values

# m = Basemap(projection='nplaea', boundinglat=70, lon_0=0, resolution='l', round=True)

# fig, ax = plt.subplots(figsize=(8, 8))
# m.drawcoastlines(linewidth=1.2)
# m.drawparallels([80], labels=[True], linewidth=0.5, color='gray', fontsize=6)
# m.drawmeridians([-180, -90, 0, 90], labels=[True, True, True, True], linewidth=0.5, color='gray', fontsize=6)

# x, y = m(lons, lats)
# print(min(diso), max(diso))

# sc = m.scatter(x, y, s=10, c=diso, cmap='jet', vmin=50, vmax=300)
# sc = m.scatter(x, y, s=10, c=hfw, cmap='jet', vmin=0, vmax=25)

# cbar = plt.colorbar(sc, ax=ax, shrink=0.7)
# cbar.set_label("Depth of 34 isohaline (m)")
# plt.show()
# plt.savefig('diso_2011_2018.png', dpi=300)


####### Plot mapped hfw/disoh
# m = Basemap(projection='nplaea', boundinglat=70, lon_0=0, resolution='l', round=True)

# fig, ax = plt.subplots(figsize=(8, 8))
# m.drawcoastlines(linewidth=1.2)
# m.drawparallels([80], labels=[True], linewidth=0.5, color='gray', fontsize=6)
# m.drawmeridians([-180, -90, 0, 90], labels=[True, True, True, True], linewidth=0.5, color='gray', fontsize=6)

# x, y = m(lons, lats)
# # print(min(hFW), max(hFW))
# print(min(D_Siso), max(D_Siso))

# sc = m.scatter(x, y, s=20, c=D_Siso, cmap='jet', vmin=30, vmax=250)
# # sc = m.scatter(x, y, s=20, c=hFW, cmap='jet', vmin=0, vmax=25)

# cbar = plt.colorbar(sc, ax=ax, shrink=0.7)
# cbar.set_label("Depth of 34 isohaline (m)")
# # plt.show()
# plt.savefig('figures/2012_01_disoh.png', dpi=300)


##### plot monthly gsc 
# fig, axes = plt.subplots(2, 6, figsize=(16, 8))
# plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05, hspace=0.1, wspace=0.05)

# for ax in axes.flatten():
#     m = Basemap(projection='nplaea', ax=ax, boundinglat=70, lon_0=0, resolution='l', round=True)
#     m.drawcoastlines(linewidth=1.2)
#     m.drawmeridians([-180, -90, 0, 90], linewidth=0.5, color='gray', fontsize=6)
#     x, y = m(lons, lats)

# df = pd.read_csv('results/grd_dh_201.csv')
# dhs, ug, vg = df['Surf_DH'], df['ug'], df['vg']
# sc = axes[0,0].scatter(x, y, s=30, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05
# Q = axes[0,0].quiver(x, y, ug, vg, scale=scale, width=0.003)
# print(np.mean(dhs), np.mean(np.sqrt(ug**2 + vg**2)))

# df = pd.read_csv('results/grd_dh_2012_02.csv')
# dhs, ug, vg = df['Surf_DH'], df['ug'], df['vg']
# sc = axes[0,1].scatter(x, y, s=30, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05
# Q = axes[0,1].quiver(x, y, ug, vg, scale=scale, width=0.003)
# print(np.mean(dhs), np.mean(np.sqrt(ug**2 + vg**2)))

# df = pd.read_csv('results/grd_dh_2012_03.csv')
# dhs, ug, vg = df['Surf_DH'], df['ug'], df['vg']
# sc = axes[0,2].scatter(x, y, s=30, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05
# Q = axes[0,2].quiver(x, y, ug, vg, scale=scale, width=0.003)
# print(np.mean(dhs), np.mean(np.sqrt(ug**2 + vg**2)))

# df = pd.read_csv('results/grd_dh_2012_04.csv')
# dhs, ug, vg = df['Surf_DH'], df['ug'], df['vg']
# sc = axes[0,3].scatter(x, y, s=30, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05
# Q = axes[0,3].quiver(x, y, ug, vg, scale=scale, width=0.003)
# print(np.mean(dhs), np.mean(np.sqrt(ug**2 + vg**2)))

# df = pd.read_csv('results/grd_dh_2012_05.csv')
# dhs, ug, vg = df['Surf_DH'], df['ug'], df['vg']
# sc = axes[0,4].scatter(x, y, s=30, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05
# Q = axes[0,4].quiver(x, y, ug, vg, scale=scale, width=0.003)
# print(np.mean(dhs), np.mean(np.sqrt(ug**2 + vg**2)))

# df = pd.read_csv('results/grd_dh_2012_06.csv')
# dhs, ug, vg = df['Surf_DH'], df['ug'], df['vg']
# sc = axes[0,5].scatter(x, y, s=30, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05
# Q = axes[0,5].quiver(x, y, ug, vg, scale=scale, width=0.003)
# print(np.mean(dhs), np.mean(np.sqrt(ug**2 + vg**2)))

# df = pd.read_csv('results/grd_dh_2012_07.csv')
# dhs, ug, vg = df['Surf_DH'], df['ug'], df['vg']
# sc = axes[1,0].scatter(x, y, s=30, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05
# Q = axes[1,0].quiver(x, y, ug, vg, scale=scale, width=0.003)
# print(np.mean(dhs), np.mean(np.sqrt(ug**2 + vg**2)))

# df = pd.read_csv('results/grd_dh_2012_08.csv')
# dhs, ug, vg = df['Surf_DH'], df['ug'], df['vg']
# sc = axes[1,1].scatter(x, y, s=30, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05
# Q = axes[1,1].quiver(x, y, ug, vg, scale=scale, width=0.003)
# print(np.mean(dhs), np.mean(np.sqrt(ug**2 + vg**2)))

# df = pd.read_csv('results/grd_dh_2012_09.csv')
# dhs, ug, vg = df['Surf_DH'], df['ug'], df['vg']
# sc = axes[1,2].scatter(x, y, s=30, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05
# Q = axes[1,2].quiver(x, y, ug, vg, scale=scale, width=0.003)
# print(np.mean(dhs), np.mean(np.sqrt(ug**2 + vg**2)))

# df = pd.read_csv('results/grd_dh_2012_10.csv')
# dhs, ug, vg = df['Surf_DH'], df['ug'], df['vg']
# sc = axes[1,3].scatter(x, y, s=30, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05
# Q = axes[1,3].quiver(x, y, ug, vg, scale=scale, width=0.003)
# print(np.mean(dhs), np.mean(np.sqrt(ug**2 + vg**2)))

# df = pd.read_csv('results/grd_dh_2012_11.csv')
# dhs, ug, vg = df['Surf_DH'], df['ug'], df['vg']
# sc = axes[1,4].scatter(x, y, s=30, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05
# Q = axes[1,4].quiver(x, y, ug, vg, scale=scale, width=0.003)
# print(np.mean(dhs), np.mean(np.sqrt(ug**2 + vg**2)))

# df = pd.read_csv('results/grd_dh_2012_12.csv')
# dhs, ug, vg = df['Surf_DH'], df['ug'], df['vg']
# sc = axes[1,5].scatter(x, y, s=30, c=dhs, cmap='YlGnBu', vmin=0, vmax=0.8)
# scale = np.nanmax(np.sqrt(ug**2 + vg**2)) / 0.05
# Q = axes[1,5].quiver(x, y, ug, vg, scale=scale, width=0.003)
# print(np.mean(dhs), np.mean(np.sqrt(ug**2 + vg**2)))

# # cbar = plt.colorbar(sc, shrink=0.7)
# # cbar.set_label("Dynamic height (m)")

# # qk = plt.quiverkey(Q, 0.1, 0.1, 0.1, '10 cm/s', labelpos='W')
# plt.tight_layout
# plt.show()