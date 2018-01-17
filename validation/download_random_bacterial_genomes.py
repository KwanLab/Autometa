#!/usr/bin/env python

# Copyright 2018 Ian J. Miller, Evan Rees, Izaak Miller, Jason C. Kwan
#
# This file is part of Autometa.
#
# Autometa is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Autometa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Autometa. If not, see <http://www.gnu.org/licenses/>.

# Program to download a random number of bacterial genomes from NCBI
# Before starting - download all the bacterial genomes from NCBI:
# rsync -avzh --progress ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria ./
# Don't worry - this just downloads a series of symlinks (~300 MB)
# USAGE: download_random_bacterial_genomes.py <number to download>
# By the way - the paths are currently hardcoded below.  Feel free to change.

import sys
import subprocess
import pprint
from random import randint
import os

pp = pprint.PrettyPrinter(indent=4)

n = int(sys.argv[1])

# Store a list of the bacterial species
bacteria_path = '/mnt/not_backed_up/ncbi_bacteria_genomes/bacteria'
ls_output = subprocess.check_output(['ls', '-1', bacteria_path])
bacteria_list = ls_output.split()
highest_index = len(bacteria_list) - 1
print highest_index

# Now make a random dictionary of n bacteria
random_bacteria = dict()
while (len(random_bacteria) < n):
	i = randint(0, highest_index)
	while (i in random_bacteria):
		i = randint(0, highest_index)
	dir_path = bacteria_path + '/' + bacteria_list[i]
	if os.path.isdir(dir_path):
		# Check if representative or latest_assembly_versions dirs exist
		if (os.path.exists(dir_path + '/representative')) or (os.path.exists(dir_path + '/latest_assembly_versions')):
			random_bacteria[i] = 1

# Now download each random bacterial genome (first representative assembly)
destination_path = '/mnt/not_backed_up/ncbi_bacteria_genomes/random_3000_bacteria'
for i in random_bacteria:
	species_directory = bacteria_list[i]
	# Get directory listing for the representative assemblies
	asm_path = bacteria_path + '/' + species_directory + '/representative'
	if not os.path.exists(asm_path):
		asm_path = bacteria_path + '/' + species_directory + '/latest_assembly_versions'
	asm_ls_output = subprocess.check_output(['ls', '-1', asm_path])
	asm_list = asm_ls_output.split()
	
	asm_dir = asm_list[0]
	destination_dir = destination_path + '/' + species_directory
	os.makedirs(destination_dir)

	# Download directory
	print 'Species: ' + species_directory
	print 'Asm: ' + asm_dir
	subprocess.call(['rsync', '-avzh', '--progress', 'ftp.ncbi.nlm.nih.gov::genomes/all/' + asm_dir, destination_dir])