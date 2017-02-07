import os
import pandas as pd
print "\n \n \n \n \n"
directory  = str(raw_input("Enter directory of Bruker files: \n"))

spectra = []
datafiles = []
numspectra = 0
count = 0
count2 = 0
for dir, subdirs, files in os.walk(directory, followlinks = True):
  if '1r' in files:
    spectrumfile = dir + '/' + '1r'
    procsfile    = dir + '/' + 'procs'
    titlefile    = dir + '/' + 'title'
    datafiles.append(titlefile)
    count+=1
datafiles.sort()
titles_list = []

for titlefile in datafiles:
    if os.path.exists(titlefile):
      count2+=1
      title = open(titlefile).read().strip()
    else:
      title = "<No title>"
    titles_list.append(title)

titles_list = pd.Series(titles_list)
print "Successfully copied %i/%i of the files found" %(count,count2)
titles_list.to_csv(directory+'titles_list.csv')
print "File saved to \n %s \n as 'titles_list.csv'" %directory
