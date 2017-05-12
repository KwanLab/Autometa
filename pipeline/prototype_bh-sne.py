#!/usr/bin/env python

# Prototype code implementation of the Vizbin BH-tSNE algorithm
# See Scientific Reports 2014, 4, 4516

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from scipy import stats
import math
from sklearn import decomposition
from sklearn.manifold import TSNE
from tsne import bh_sne
import numpy as np
import pprint
import pdb

def revcomp( string ):
	trans_dict = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C' }
	complement_list = list()

	for i in range(0, len(string) ):
		if string[i] in trans_dict:
			complement_list.append(trans_dict[string[i]])
		else:
			return -1

	return ''.join(reversed(complement_list))


pp = pprint.PrettyPrinter(indent=4)

input_file = '78.125Mbp.fasta'
length_cutoff = 10000
k_mer_size = 5

print('Filtering sequences')
filtered_sequences = list()
for seq_record in SeqIO.parse(input_file, 'fasta'):
	if len(seq_record.seq) >= length_cutoff:
		filtered_sequences.append(seq_record)

# Now we count k-mers
# First we make a dictionary of all the possible k-mers (discounting revcomps)
# Under each key is an index to be used in the subsequent lists
# The order of the indices depends on the order k-mers were encountered while making the dictionary
print('Counting k-mers')
count = 0
unique_k_mers = dict()

DNA_letters = ['A', 'T', 'C', 'G']
all_k_mers = list(DNA_letters)
for i in range(1, k_mer_size):
	new_list = list()
	for current_seq in all_k_mers:
		for char in DNA_letters:
			new_list.append(current_seq + char)
	all_k_mers = new_list

# Now we trim k-mers and put them in the dictionary
for k_mer in all_k_mers:
	k_mer_reverse = revcomp(k_mer)
	if (k_mer not in unique_k_mers) and (k_mer_reverse not in unique_k_mers):
		unique_k_mers[k_mer] = count
		count += 1

contig_k_mer_counts = list()

for i in range(0, len(filtered_sequences)):
	contig_name = str(filtered_sequences[i].id)
	contig_seq = str(filtered_sequences[i].seq)
	# Initialize the list for the contig to be all 1s - as we can't have zero values in there for CLR later
	list_size = len(unique_k_mers.keys())
	current_contig_k_mer_counts = [1 for k in range(0, list_size)]

	for j in range(0, (len(contig_seq) - k_mer_size)):
		k_mer = contig_seq[j:j + k_mer_size]
		k_mer_reverse = revcomp(k_mer)

		# Find appropriate index
		# Note - this part naturally ignores any k_mers with weird characters
		if k_mer in unique_k_mers:
			index = unique_k_mers[k_mer]
			current_contig_k_mer_counts[index] += 1
		elif k_mer_reverse in unique_k_mers:
			index = unique_k_mers[k_mer_reverse]
			current_contig_k_mer_counts[index] += 1

	contig_k_mer_counts.append(current_contig_k_mer_counts)

# We now remove all the k-mers where all counts are '1'
print ('Trimming k-mers')
columns_to_delete = dict()
for i in range(0, len(unique_k_mers.keys())):
	non_zero_counts = 0
	for j in range(0, len(contig_k_mer_counts)):
		if contig_k_mer_counts[j][i] > 1:
			non_zero_counts += 1

	if non_zero_counts == 0:
		columns_to_delete[i] = 1


filtered_contig_k_mer_counts = list()
for i in range(0, len(contig_k_mer_counts)):
	new_row = list()
	for j in range(0, len(unique_k_mers.keys())):
		if j not in columns_to_delete:
			new_row.append(contig_k_mer_counts[i][j])
	filtered_contig_k_mer_counts.append(new_row)
filtered_contig_k_mer_counts = contig_k_mer_counts

print('Normalizing counts')

k_mer_frequency_matrix = list()

for count_list in filtered_contig_k_mer_counts:
	total_count = 0
	for count in count_list:
		total_count += count

	normalized_list = list()
	for count in count_list:
		normalized_list.append(float(count)/total_count)

	# Now we calculate the Centered log-ratio (CLR) transformation
	# See Aitchison, J. The Statistical Analysis of Compositional Data (1986) and 
	# Pawlowsky-Glahn, Egozcue, Tolosana-Delgado. Lecture Notes on Compositional Data Analysis (2011)
	geometric_mean = stats.mstats.gmean(normalized_list)
	clr_list = list()

	for item in normalized_list:
		intermediate_value = item / geometric_mean
		clr_list.append(math.log(intermediate_value))

	k_mer_frequency_matrix.append(clr_list)


print('Principal component analysis')

pca = decomposition.PCA(n_components=50)
pca_matrix = pca.fit_transform(k_mer_frequency_matrix)

print('BH-tSNE')

# Note TSNE uses 'angle' instead of 'theta'
#model = TSNE(n_components=2, random_state=0, perplexity=30.0, n_iter=1000, angle=0.5, method='barnes_hut')
#bh_tsne_matrix = model.fit_transform(pca_matrix)

X = np.array(pca_matrix)
bh_tsne_matrix = bh_sne(X, d=2, perplexity=30.0, theta=0.5)

print('Outputting file')
output = open('output_table', 'w')
output.write('contig\tbh_tsne_x\tbh_tsne_y\n')

for i in range(0, len(filtered_sequences)):
	contig_name = str(filtered_sequences[i].id)
	bh_tsne_x = str(bh_tsne_matrix[i][0])
	bh_tsne_y = str(bh_tsne_matrix[i][1])

	output.write(contig_name + '\t' + bh_tsne_x + '\t' + bh_tsne_y + '\n')

output.close()
