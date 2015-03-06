#!/usr/bin/env python

"""Tool for checking the 3' end of a fastq sequence for 3' uridylation. By Andrew Pattison"""

fastq_file = "231HM-miR242rep1_ACTTGA_L002_R1_001.fastq"

adapter = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

fqfile = open (fastq_file, "rU").readlines()[1::4]
adlen = len(adapter)+6

def get_adapter_reads(adapter_seq,target_file):
	full_adapter_reads = []
	full_length_adapter_reads = []
	for line in target_file:
		# Might be worth putting in a hamming distance calculation here. 
		if adapter_seq[0:20] in line[-22::]:
			full_adapter_reads.append(line)
		if adapter_seq in line:
			full_length_adapter_reads.append(line)

	print '\n','20 bases of the adapter seq are in', len(full_adapter_reads),'reads'  
	print 'The full adapter seq is in', len(full_length_adapter_reads),'reads' ,'\n' 
	return(full_adapter_reads, full_length_adapter_reads)

def check_seq_after_adapter(full_adapter_reads, adapter_seq):
	for line in full_adapter_reads:
		line.split(adapter_seq)
		if len(line[1]) > 1:
			print'More than 1 extra base follows the adapter',(line[1])

def see_pre_seq(adapter_reads, adapter):
	ends_in_pol_A = []
	no_poly_a_end = []
	for line in adapter_reads:
		# clip the adapter and everything after
		split_line = line.rsplit(adapter[0:20],1)[0]
		# grab the last 10 reads
		cut_line = split_line[-10::]
		if 'AAAA' in cut_line:
			ends_in_pol_A.append(cut_line)
		else:
			no_poly_a_end.append(cut_line)
	print 'There are',len(ends_in_pol_A), 'reads with a poly_A tail (4As in the final 10) directly before the adpater'
	print 'There are',len(no_poly_a_end), 'reads without a poly_A tail (4As in the final 10) directly before the adpater', '\n'
	return(ends_in_pol_A, no_poly_a_end)

def check_for_U(poly_a_ends):
	not_A = []
	only_A = 0
	for line in poly_a_ends:
		if line[0:10:] != 'AAAAAAAAAA':
			not_A.append(line)
		else: only_A += 1
	print 'There are', len(not_A), 'reads that end with a sequence other than AAAAAAAAAA directly before the adapter'
	print 'There are', only_A, 'reads that contain only an AAAAAAAAAA directly before the adapter', '\n'
	return(not_A)

def get_end_percents(not_just_A):
	A_reads = []
	C_reads = []
	T_reads = []
	G_reads = []
	for line in not_just_A:
		if line [-1::] == 'A':
			A_reads.append(line)
		elif line [-1::] == 'C':
			C_reads.append(line)
		elif line [-1::] == 'T':
			T_reads.append(line)
		elif line [-1::] == 'G':
			G_reads.append(line)
		else: print 'not a known nucleotide'
		total = len(not_just_A)
	print 'A reads =', float(len (A_reads))/total*100,'%,', 'C reads =', float(len (C_reads))/total*100 ,'%,', 'T reads =', float(len (T_reads))/total*100,'%,', 'G reads =', float(len (G_reads))/total*100, '%'
	return ()



adapter_reads = get_adapter_reads(adapter, fqfile)

post_adapter = check_seq_after_adapter(adapter_reads[1], adapter)

clipped_poly_A = see_pre_seq(adapter_reads[0], adapter)

print '\n', 'of these', len(clipped_poly_A[0]), 'reads with a poly-A tail:', '\n' 

not_just_A = check_for_U(clipped_poly_A[0])

print '\n', 'of these', len(not_just_A), 'reads without a homogenous poly-A sequence:', '\n' 


get_end_percents(not_just_A)