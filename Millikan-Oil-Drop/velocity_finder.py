# -*- coding: utf-8 -*-
"""
import pandas as pd
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt

data = pd.read_excel(".xlsx", na_values = '#NV').fillna(0).to_numpy(np.int64).flatten()
time = np.linspace(0, (len(data)-1), len(data))

plt.figure(0)
plt.scatter(time, data)
plt.xlabel("Frame")
plt.ylabel("Pixels from Top")

print(sig.argrelextrema(data, np.less))
print(sig.argrelextrema(data, np.greater))

vt = 10*(data[X]-data[Y])/(520*(time[X]-time[Y]))
v2 = 10*(data[A]-data[B])/(520*(time[A]-time[B]))
print("vt: " + str(vt))
print("v2: " + str(v2))
"""
import pandas as pd
import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt

data = pd.read_excel("9.xlsx", na_values = '#NV').fillna(0).to_numpy(np.int64).flatten()
time = np.linspace(0, (len(data)-1), len(data))

plt.figure(0)
plt.scatter(time, data)
plt.xlabel("Frame")
plt.ylabel("Pixels from Top")

print(sig.argrelextrema(data, np.less))
print(sig.argrelextrema(data, np.greater))

vt = 10*(data[88]-data[21])/(520*(time[88]-time[21]))
v2 = 10*(data[265]-data[195])/(520*(time[265]-time[247]))
print("vt: " + str(vt))
print("v2: " + str(v2))

