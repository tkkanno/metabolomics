import numpy as np
import matplotlib.pyplot as plt#
import pandas as pd

print "this script will take NMR data along with a csv of metabolites and ppm"
print "It will then label the NMR data to allow visualisation of you metabolite list."

directory = raw_input('Enter the directory of the data: ')

spct = raw_input('Enter the spectra text file: ')
mets = raw_input('Enter name of metabolite csv file. must be in format Compound then ppm: ')
data= pd.read_csv(directory+spct, sep = ' ')
mets= pd.read_csv(directory+mets, sep = '\t')

ppm = data.columns

ppm = ppm[1:]
ppm = np.array(ppm, dtype = float)
data = np.array(data)
data = data[:,1:]
new_ppms = []
peaks_to_label =[]
mets.columns = ['metabolite', 'ppm']
for i in range(len(mets)):
  idx 	= np.argmin(np.abs(ppm - float(mets.ix[i].ppm)))
  old_idx = idx
  idxmin = np.argmin(np.abs(ppm - float(mets.ix[i].ppm) - 0.002)) #signs reversed b/c nmr spectra x-axis is reversed
  idxmax = np.argmin(np.abs(ppm - float(mets.ix[i].ppm) + 0.002))
  sppm = ppm[idxmin:idxmax]
  max_spect = data[:, idxmin:idxmax].mean(axis = 0)
  try:
    if max_spect.max() >= data[:,idx].mean(axis = 0).max():
      b = np.where(float(max_spect.max()) ==max_spect )
    else:
      b = idx
    #take slices of data where peak is highest
    if len(b)==1:
      b = b[0][0]
      idx = b
      new_ppm = float(sppm[idx])
      new_idx = np.argmin(np.abs(ppm-new_ppm))
      new_ppms.append(new_ppm)
    else:
      print """Too many equal values, taking the first peak at %.3f and moving onto the next one.""" %(sppm[b[0][0]])
      print "There are equal peaks at %i lications." %(len(b)-1)
      new_ppm = float(ppm[idx])
    peakY = data[:,new_idx].max()
    peaks_to_label.append(peakY)
    print "Peak changed %.3f, plotting %s at (%.3f, %.3f)" %(np.abs(ppm[old_idx]-new_ppm), mets.ix[i].metabolite, new_ppm, peakY)
  except ValueError:
    print "something went wrong with %s. Onto the next one:" %(mets.ix[i].metabolite)
    next
new_ppms = np.array(new_ppms)
peaks_to_label = np.array(peaks_to_label)
if len(np.unique(new_ppms)) < len(mets):
  print "there are values that overlap each other"
  for i in range(len(new_ppms)):
    for j in range(len(new_ppms)):
      if i!=j and new_ppms[i]==new_ppms[j]:
	print "%.3f is the same in positions %i and %i. These are %s and %s" %(new_ppms[i], i, j, mets.ix[i].metabolite, mets.ix[j].metabolite)

fig,ax = plt.subplots()
for i in range(len(data)):
  plt.plot(ppm, data[i], 'b', alpha = 0.75)
  
plt.gca().invert_xaxis()

for i in range(len(mets)):
  plt.annotate(mets.ix[i].metabolite, xy = (new_ppms[i], peaks_to_label[i]),
		xytext=(0, 20), textcoords='offset points',
		arrowprops=dict(arrowstyle="->"), rotation = 75,
		va = 'bottom'
		)  
plt.show()

cls
  
#  ax2.annotate(lab, xy=(piccoX,peakY),  xycoords='data',
#                xytext=(0, 20), textcoords='offset points',
#                arrowprops=dict(arrowstyle="->"), rotation = 75,
#		va = 'bottom'
#               )
  