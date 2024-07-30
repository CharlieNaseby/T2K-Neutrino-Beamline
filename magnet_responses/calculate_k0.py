import matplotlib.pyplot as plt
import csv
import pandas as pd
import numpy as np

data = []

pd.set_option('display.max_rows', 50000)

with open('BPD2.dat', newline='') as csvfile:
    field = csv.reader(csvfile, delimiter=' ')
    for row in field:
        if(len(row) > 2 and row[0]!='!'):
            data.append([float(entry) for entry in row])


df = pd.DataFrame(data, columns=['z', 'Bx', 'By', 'Bz'])
df['zdiff'] = df['z'].diff()

zint = np.sum(df['By']*df['zdiff']/100.)

print(zint)
