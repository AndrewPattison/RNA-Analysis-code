import pysam
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import re
import pylab
from collections import defaultdict

samfile = pysam.Samfile( "LM2m1.bam", "rb" )
region = samfile.fetch('chrM')

list_of_lists = []
for read in region:
    extracted = str(read)
    AS = extracted.find('AS')
    a_flag_string = extracted[AS:AS+9]
    a_flag_number = re.findall(r'\d+', a_flag_string)
    a_flag_number = int (a_flag_number[0])
    regionpos = read.pos
    list_for_list = []
    list_for_list.append(regionpos)
    list_for_list.append(a_flag_number)
    list_of_lists.append(list_for_list)

    

list_of_lists.sort(key=lambda x:x)

# List I am appending to
newlist=[]
counter = 1
for ls in list_of_lists:
	newlist.append(ls)

sums = defaultdict(int)
counts = defaultdict(int)
    
xaxis = []
yaxis = []
#zaxis = []
for lis in newlist:
    xaxis.append(lis[0])
    yaxis.append(lis[1])
    
xmax = max(xaxis)
ymax = max(yaxis)

print(newlist)

for k in newlist:
    sums[k[0]] += k[1]
    counts[k[0]] += 1
    
    
    
    
positions= []
average_As = []

for x in sums.keys():
    print "average " + str(x) + ": " +  str(float(sums[x])/counts[x])
    positions.append(int(x))
    average_As.append(float(sums[x])/counts[x])
    
print(positions)
print(average_As)

plt.plot(positions, average_As)
plt.show()