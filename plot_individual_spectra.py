import numpy as np
import pandas as pd

def pd_to_np(data):
  ppm = data.columns
  ppm = np.array(ppm, dtype = float)
  data = np.array(data)
  return  ppm, data

def plt_spect(ppm, data):
  #scale data to 1000
  for i in range(len(data)):
    data[i] = (data[i]/data[i].mean())*100
  
  fig,ax = plt.subplots()
  count = 0
  ppmcount = 0
  for i in range(len(data)):
    plt.plot(ppm+ppmcount, data[i] + count)
    count+=400*(1000/data.mean())
    ppmcount += 0.025
  plt.gca().invert_xaxis()
  plt.show()
  
plt_spect(ppm, pd)

