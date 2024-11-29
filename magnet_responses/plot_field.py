import matplotlib.pyplot as plt
import csv
import pandas as pd
import numpy as np

data = []

pd.set_option('display.max_rows', 50000)

with open('QPQ4.dat', newline='') as csvfile:
    field = csv.reader(csvfile, delimiter=' ')
    for row in field:
        if(len(row) > 2 and row[0]!='!'):
            data.append([float(entry) for entry in row])


df = pd.DataFrame(data, columns=['x', 'y', 'z', 'Bx', 'By', 'Bz'])
line = df[(df['x'].between(10.6, 11)) & (df['y'].between(-0.1, 0.1))]
print(line)
plt.plot(line['z'], line['By'])
plt.show()


unique_z = np.abs(df['z'].unique())
centre = np.min(unique_z)

delta = 1e-5

slice = df[df['z'].between(-centre-delta, centre+delta )]

print(slice)





pivot_table = slice.pivot(index='y', columns='x', values='Bx')

plt.figure(figsize=(8, 6))
plt.imshow(pivot_table, origin='lower', cmap='viridis', interpolation='nearest')
plt.colorbar(label='Z Value')
plt.xlabel('X Bin')
plt.ylabel('Y Bin')
plt.title('2D Histogram')
plt.show()


#fig, ax = plt.subplots()
#hist = ax.hist2d(slice['x'], slice['y'], weights=slice['By'])
#fig.colorbar(hist[3], ax=ax)
#plt.show()

