import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# saved from SRC: https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/monthly.ao.index.b50.current.ascii.table
df = pd.read_csv('ao_1980_2018.csv')
print(df.head)

# to compute annual mean of AO (mean of Jan-Dec rows)
# df["AO"] = df.iloc[:, 1:].mean(axis=1)
# df.to_csv('ao_1980_2018.csv', index=False)
# p = plt.plot(df['Year'], df['AO'], 'r-o')
# plt.show()

# years = df['Year']
# months_labels = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
# monthly_ao = np.array([ list(df[mon]) for mon in months_labels ])

# years_extended = np.linspace(2011, 2018, num=monthly_ao.size)  # Extended time axis for all months
# monthly_ao_flat = monthly_ao.flatten()  # Flatten the monthly AO data

# Plotting
# plt.figure(figsize=(14, 6))

# # monthly data
# plt.plot(years_extended, monthly_ao_flat, c='blue', label='Monthly')

# # Labels and title
# plt.xlabel("Year", fontsize=12)
# plt.ylabel("AO Index", fontsize=12)
# plt.title("Arctic Oscillation (AO) Index (2011-2018)", fontsize=14)
# plt.grid(True, linestyle="--", alpha=0.6)
# for year in years:
#     plt.axvline(year, color='gray', linestyle='--', alpha=0.5)

# # annual mean plot
# plt.plot(df['Year'], df['AO'], c='red', marker='o', linestyle='-', label='Annual mean')
# plt.legend(fontsize=10)
# # plt.show()

# plt.savefig('ao_2011_2018.png', dpi=300)
