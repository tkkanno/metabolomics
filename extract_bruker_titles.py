import os
import pandas as pd
directory  = '/home/louic/Desktop/parkinsons_project/170111_parkserum/cpmg/'

spectra = []
datafiles = []
numspectra = 0

for dir, subdirs, files in os.walk(directory, followlinks = True):
  if '1r' in files:
    spectrumfile = dir + '/' + '1r'
    procsfile    = dir + '/' + 'procs'
    titlefile    = dir + '/' + 'title'
    datafiles.append(titlefile)
datafiles.sort()
titles_list = []
for titlefile in datafiles:
    if os.path.exists(titlefile):
      title = open(titlefile).read().strip()
    else:
      title = "<No title>"
    titles_list.append(title)

titles_list = pd.Series(titles_list)

titles_list.to_csv(directory+'titles_list.csv')
