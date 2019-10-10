#!/usr/bin/python

# Copyright 2015 Lina Sieverling, Philip Ginsbach, Lars Feuerbach

# This file is part of TelomereHunter.

# TelomereHunter is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# TelomereHunter is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with TelomereHunter.  If not, see <http://www.gnu.org/licenses/>.



import subprocess
import os
import sys
import re
import csv
import pysam


##################################################################
### get the gc content distribution of the reads in a bam file ###
##################################################################

def get_gc_content_distribution(bam_file, out_dir, pid, sample, remove_duplicates):

	# open input bam_file for reading
	bamfile = pysam.Samfile( bam_file, "rb" )

	# make GC content list
	gc_content_list = { gc_content:0 for gc_content in range(0,100+1) }


	for read in bamfile.fetch(until_eof=True):   

		if read.is_secondary:		#skips all secondary alignments (only needed for RNA analysis!)
			continue
		  
		if remove_duplicates == True and read.is_duplicate:		#if duplicate flag is set: skip all reads that are flagged as optical or PCR duplicate
			continue
		      
		if read.flag >= 2048:						# skip supplementary alignments 
			continue
		  
		sequence = read.seq
		read_length = len(read.seq)
		
		c_count = sequence.count('C')
		g_count = sequence.count('G')
		n_count = sequence.count('N')
		

		if float(n_count)/float(read_length) > 0.2:
			continue
		
		gc_content = int(round(float(c_count + g_count) / float(read_length - n_count) * 100))
		
		gc_content_list[gc_content] += 1
		

	#############################
	### write gc content file ###
	#############################

	# open gc content file for writing
	gc_content_file = open(out_dir + "/" + pid + "_" + sample + "_gc_content.tsv", "w")

	# header
	gc_content_file.write( "gc_content_percent\tread_count\n")

	#write line for each gc content
	for gc_content in range(0,100+1) :
	  
		gc_content_file.write( "%i\t%i\n" % (gc_content, gc_content_list[gc_content]))


	##########################
	### close file handles ###
	##########################

	bamfile.close()
	gc_content_file.close()




####################################################################################################
### estimate the telomere content in telomeric reads per million reads with a similar gc content ###
####################################################################################################

def estimate_telomere_content(input_dir, out_dir, pid, sample, gc_lower, gc_upper):
	
	#gc bins used for correction
	gc_bins_correction = range(gc_lower, gc_upper+1)
	
	# read in total gc content counts
	gc_content_file = input_dir + "/" + pid + "_" + sample + "_gc_content.tsv"

	gc_content_list = {}

	with open(gc_content_file,'rb') as tsvin:
	    tsvin = csv.reader(tsvin, delimiter='\t')
	    
	    next(tsvin, None)  # skip the headers

	    for row in tsvin:
		gc_content_list[int(row[0])]=int(row[1])


	# get total number of reads in these bins
	sum_over_threshold=0

	for gc in gc_bins_correction:
		sum_over_threshold+=gc_content_list[gc]


	# get total read count and number of intratelomeric_reads
	read_count_file = input_dir + "/" + pid + '_readcount.tsv'
	total_read_count = subprocess.Popen('less ' + read_count_file + ' | awk \'{if (NR!=1) {print}}\'  | awk \'{ SUM += $3}  END { print SUM }\'' , shell=True, 
stdout=subprocess.PIPE).stdout.read()
	total_read_count = int(float(total_read_count.strip()))

	filtered_intratelomeric_file = input_dir + "/" + pid + '_filtered_intratelomeric.bam'
	intratel_read_count = subprocess.Popen('samtools view -c ' + filtered_intratelomeric_file, shell=True, stdout=subprocess.PIPE).stdout.read()
	intratel_read_count = int(intratel_read_count.strip())


	#estimate telomere content (intratelomeric reads per million reads with similar gc content)
	tel_content = float(intratel_read_count) / float(sum_over_threshold) *1000000


	##########################
	### write summary file ###
	##########################

	# open summary file for writing
	summary_file = open(out_dir + "/" + pid + "_" + sample + "_summary.tsv", "w")

	# header
	#summary_file.write( "PID\tsample\ttotal_reads\tintratel_reads\tintratel_gc_count_threshold\ttotal_reads_with_tel_gc\ttel_content\n")
	summary_file.write( "PID\tsample\ttotal_reads\tintratel_reads\tgc_bins_for_correction\ttotal_reads_with_tel_gc\ttel_content\n")
	
	# results 
	gc_bins = str(gc_lower) + "-" + str(gc_upper)
	summary_file.write( "%s\t%s\t%i\t%i\t%s\t%i\t%f\n" % (pid, sample, total_read_count, intratel_read_count, gc_bins, sum_over_threshold, tel_content))

	##########################
	### close file handles ###
	##########################

	summary_file.close()
	

