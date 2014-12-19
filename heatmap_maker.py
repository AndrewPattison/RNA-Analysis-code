import fnmatch
import os
import pysam
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import re
import pylab
from collections import defaultdict

def get_regions_of_interst():
	
	chromosome = raw_input('What chromosome you you like to explore?')
	start = int(raw_input('What chromosomal position would you like to start at?'))
	end = int(raw_input('What chromosomal position would you like to end at?'))

	return(chromosome,start,end) 


def get_all_bam_vals(coords):
    
    matches = []
    for root, dirnames, filenames in os.walk('.'):
      for filename in fnmatch.filter(filenames, '*bam'):
          matches.append(os.path.join(root, filename))
    
    regions_list = []
    
# For each bam file get the regions corresponding to the give coordinates.  
    for bam_file in matches:
        samfile = pysam.Samfile( bam_file, "rb" )
        region = samfile.fetch(coords[0],coords[1],coords[2])
        regions_list.append(region)
     
    return (regions_list)
        
def get_list_of_lists (regions_list):

    list_of_lists = []
    for region in regions_list:
        for read in region:
            #print(dir(read))
            read_length = read.inferred_length
            direction = read.is_reverse
            extracted = str(read.tags)
            AS = extracted.find('AS')
            a_flag_string = extracted[AS:AS+9]
            a_flag_number = re.findall(r'\d+', a_flag_string)
            a_flag_number = int (a_flag_number[0])
            regionpos = read.pos
            list_for_list = []
            list_for_list.append(regionpos)
            list_for_list.append(a_flag_number)
            list_for_list.append(direction)
            list_for_list.append(read_length)
            list_of_lists.append(list_for_list)
            
    list_of_lists.sort(key=lambda x:x)
    #print(list_of_lists)
    
    return(list_of_lists)

def split_by_strand(list_of_lists):
    plus_reads =[]
    minus_reads = []
    for value_group in list_of_lists:
        if value_group[2] == False: 
            plus_reads.append(value_group)
        if value_group[2] == True: 
                minus_reads.append(value_group)
    return (plus_reads,minus_reads)
    
def plot_matrix_plus(plus_reads):
    

    added= []
    for lis in plus_reads:
        lis[0] += lis [3]
        added.append(lis)
        
    xaxis = []
    yaxis = []
    
    for lis in added:
        yaxis.append(lis[0])
        xaxis.append(lis[1])        

    xmax = max(xaxis)
    ymax = max(yaxis)

    mat = np.zeros((xmax +1 ,ymax + 1))
    print(mat)
    

    for lis in added:
        mat[lis[1], lis[0]] += 1 

    matrix = np.log2(mat +1)
   

    pylab.pcolor(np.array(matrix))
    pylab.colorbar()
    pylab.show()
    print ('done')
    return

def plot_matrix_minus(minus_reads):
    xaxis = []
    yaxis = []
    for lis in minus_reads:
        yaxis.append(lis[1])
        xaxis.append(lis[0])

    xmax = max(xaxis)
    ymax = max(yaxis)
    
    mat = np.zeros((ymax+1 ,xmax + 1))

    for lis in minus_reads:
        mat[lis[1], lis[0]] += 1 
        
    
    matrix = np.log2(mat +1)

    pylab.pcolor(np.array(matrix))
    pylab.colorbar()
    pylab.show()
    print ('done')
    return

start_end = get_regions_of_interst() 
regions_list = get_all_bam_vals(start_end)
list_of_lists = get_list_of_lists(regions_list)
strands = split_by_strand(list_of_lists)
plot_matrix_plus(strands[0])
plot_matrix_minus(strands[1])
