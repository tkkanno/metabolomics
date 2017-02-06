import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats
import os
from numpy import trapz
import statsmodels.sandbox.stats.multicomp as mltc

def scaletointegral(ydata):
        integral = ydata.sum(1) / 100
        newydata = np.transpose(ydata.T / integral)
        return np.array(newydata)
        
def scalepqn(ydata):
      # 1: scale to integral
      newydata = scaletointegral(ydata)
      # 2: reference spectrum is median of all samples
      yref = np.median(ydata, 0)
      # 3: quotient of test spectra with ref spectra
      yquot = newydata / yref
      # 4: median of those quotients
      ymed = np.median(yquot, 0)
      # 5: divide by this median
      newydata = newydata / ymed
      return np.array(newydata)
      
def find_curve(ppm, data, point):
    data_mean  = data.mean(axis = 0)
    idx_og = np.argmin(np.abs(ppm-point))
    print "starting point is %.3f " %ppm[idx_og]
    moment_og = data_mean[idx_og]
    new_moment_right = data_mean[idx_og-1]
    moment = moment_og
    idx = idx_og
    while moment >= new_moment_right:
      idx -= 1
      moment = new_moment_right
      new_moment_right = data_mean[idx]
    print "right side minimum is %.3f " %(ppm[idx])
    idx_right = idx  
    moment = moment_og
    idx = idx_og
    new_moment_left = data_mean[idx+1]
    while moment >= new_moment_left:
      idx+=1
      moment = new_moment_left
      new_moment_left = data_mean[idx]
    print "left side minimum is %.3f " %(ppm[idx])
    idx_left = idx
    
    idx_left -= 1 
    idx_right += 1 
    if np.abs(idx_left-idx_right)<0:
      idx_right = idx_og-1
      idx_left = idx_og+1
    
    return idx_left, idx_right
    
def find_new_peak(ppm, data, point, range_to_search = 0.002):
  idx 	= np.argmin(np.abs(ppm - float(point)))
  old_idx = idx
  idxmin = np.argmin(np.abs(ppm - float(point) - range_to_search)) #signs reversed b/c nmr spectra x-axis is reversed
  idxmax = np.argmin(np.abs(ppm - float(point) + range_to_search))
  sppm = ppm[idxmin:idxmax]
  max_spect = data[:, idxmin:idxmax].mean(axis = 0)
  try:
    b = np.where(float(max_spect.max()) ==max_spect )
    b= b[0][0]
    #take slices of data where peak is highest
    idx = b
    new_ppm = float(sppm[idx])
    new_idx = np.argmin(np.abs(ppm-new_ppm))
    peakY = data[:,new_idx].max()
    return new_ppm
  except ValueError:
    print "something wrong with %.3f" %point
    return point

def get_integrals(ppm, data, ranges):
  return [trapz(ppm[ranges[1]:ranges[0]], i[ranges[1]:ranges[0]]) for i in data ]
      
print """This script will give you boxplots for the metabolites of interest and write the data to a new file \n
	so you can do further analysis of just assigned metabolites. \n
	it will ask you for  a file of the ppm you want to extract and the names of those metabolites"""
directory = raw_input('enter directory ending with /: ')
dat = raw_input('enter pls_input file: ')
mets = raw_input('file name of metabolites and ppm values: ')
print "processing data"
data = np.loadtxt(directory+dat)
clss = data[:,0]
clss = clss[1:]
data = data[:,1:]

try:
  mets = pd.read_csv(directory + mets, sep = ',')
except ValueError:
  mets = pd.read_csv(directory + mets, sep = '\t')
mets.columns = ['metabolites', 'ppm']

#processing of data, separating ppm from dataframe and PQN normalising the data
ppm = data[0]
data = data[1:,:]
data = scalepqn(data)

#finding the proper peaks of the data and adding to mets dataframe
new_ppms = [find_new_peak(ppm, data, i,0.005) for i in mets.ppm]
mets.ppm = new_ppms


#getting the areas under the curves for future box plots
integrals = []
ranges_for_integrals = []
for i in mets.ppm:
  a = find_curve(ppm, data, i)
  ranges_for_integrals.append(a)
  integrals.append(get_integrals(ppm, data, a))
  
ranges_for_integrals = np.array(ranges_for_integrals)
ppm_ranges_for_integrals =ppm[ranges_for_integrals]

integrals = pd.DataFrame(integrals).T
integrals.columns = mets.metabolites
integrals['class'] = clss
integrals['class'][integrals['class']>=1] =1
integrals[integrals < 0 ] = 0
folder_name = '_boxplots_pd_vs_rest/'
try:
  os.mkdir(directory+dat+folder_name)
except OSError:
  pass
print "making and saving the boxplots..."
ticks = ['PD', 'PD spouse', 'Healthy Control', 'Healthy Control spouse', 'IBD']
ticks  = ['PD', 'Rest']

for i in range(len(integrals.columns)):
  if integrals.columns[i] == 'class':
    next
    
  else:
    label = "%s %.3f ppm" %(integrals.columns[i], mets.ppm[i])
    ax = sns.boxplot(x = "class", y = integrals.columns[i], data =integrals, fliersize = 0 )
    ax = sns.stripplot(x = "class", y = integrals.columns[i], data = integrals, jitter = True, color = '.25')
    sns.plt.title(label)
    sns.plt.xticks([0,1], ticks)
    fig = ax.get_figure()
    print "saving boxplot %i" %i
    fig.savefig(directory+dat+folder_name+ integrals.columns[i])
    fig.clear()
print "Done!"
#mean comparison of integrals, ttest and fdr correction
means = integrals.groupby('class', axis = 0).mean()
means = means.T

zeros = integrals[integrals['class']==0]
ones = integrals[integrals['class']==1]
pvalues = np.array([stats.ttest_ind(zeros[zeros.columns[i]], ones[ones.columns[i]])[1] for i in range(len(integrals.columns))])
pvalues = pvalues[:pvalues.shape[0]-1]
fdr_pvalues = np.array(mltc.multipletests(pvalues, alpha  = 0.05, method = 'fdr_bh'))
fdr_pvalues = fdr_pvalues[1]
