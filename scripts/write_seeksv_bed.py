#!/usr/bin/env python3

import argparse
import sys
import csv
import pdb

def main(args):

	parser = argparse.ArgumentParser(description = "convert seeksv output to bed file")
	parser.add_argument('--seeksv-output', help = 'output from seeksv', required=True)
	parser.add_argument('--chromlist', help = 'list of chromosomes from desired reference, one per line', required=True)
	parser.add_argument('--output', help = 'name of output bed file', default='seeksv.bed')
	args = parser.parse_args()
	
	# read chromsomes into memory
	chroms = set()
	with open(args.chromlist) as handle:
		line = handle.readline().strip()
		chroms = set(line.split(" "))

	# read lines of seeksv output and write bed file
	with open(args.seeksv_output, 'r', newline='') as seeksv, open(args.output, 'w', newline='') as bed:
		reader = csv.DictReader(seeksv, delimiter='\t')
		writer = csv.writer(bed, delimiter='\t')
		
		for row in reader:
			if row['@left_chr'] in chroms:
				# this is a guess - not documented if microhomology comes before or after position
				if row['left_strand'] == "+":
					start = int(row['left_pos'])
					stop = start + int(row['microhomology_length'])
				else:
					start = int(row['left_pos']) -  int(row['microhomology_length'])
					stop = int(row['left_pos'])			
					
				row = (row['@left_chr'], start, stop, row['left_strand'])
				writer.writerow(row)
				
			elif row['right_chr'] in chroms:
				# this is a guess - not documented if microhomology comes before or after position
				if row['right_strand'] == "+":
					start = int(row['right_pos'])
					stop = start + int(row['microhomology_length'])
				else:
					start = int(row['right_pos']) -  int(row['microhomology_length'])
					stop = int(row['right_pos'])		
			
				row = (row['right_chr'], row['right_pos'], row['right_pos'], row['right_strand'])
				writer.writerow(row)
	
	
if __name__ == "__main__":
	main(sys.argv)
