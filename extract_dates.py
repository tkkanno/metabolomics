import numpy as np
import pandas as pd

from matplotlib import pyplot as plt


#converting collection data into integers for PCA and PLS

dates =pd.read_csv('dates.csv', header = None)
dates= np.array(dates, dtype = str)
dates =dates.reshape(dates.shape[0])
dates = pd.to_datetime(dates, format = '%m/%d/%Y')
d= dates.map(lambda x: x.strftime('%Y-%m-%d'))
e = np.unique(d)
e_ints= np.array(range(e.shape[0]))

f = np.array([e_ints[j] for i in d for j in range(len(e)) if i == e[j]])

np.savetxt('timecodes', f, fmt = '%i')
