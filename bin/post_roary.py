import csv
import pandas as pd
import argparse

# Also needs to create the traits.csv file for scoary

#path_to_roary_output = '/Volumes/External/CIDRI_Data/felix_S_pneumoniae/roary/roary/gene_presence_absence.csv'

#path_to_phenotype_data = '/Volumes/External/CIDRI_Data/felix_S_pneumoniae/sample_sheet_FsD_dichotomized.csv'

#path_to_sample_sheet = '/Volumes/External/CIDRI_Data/felix_S_pneumoniae/sample_sheet.csv'

#out_file_name = 'traits.csv'


#print(encoded_sample_table)
# Write to file

'''
input_file = csv.DictReader(open(path_to_roary_output))

for a_gene in input_file:

	gene_results = {}

	for a_sample in a_gene.keys():
		if a_sample.split('_')[0] == 'sample':
			sample_ID = a_sample
			if len(a_gene[a_sample]) > 1:

				# Gene present for this sample

				1 == 1
'''


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='''a description''', epilog="""the epilogue""")

	parser.add_argument('--gene_presence_absence', type=str, help='Path to roary gene presence absence file')

	parser.add_argument('--phenotype_data', type=str, help='CSV file describing the phenotype of the samples')

	parser.add_argument('--sample_sheet', type=str, help='Sample sheet')

	parser.add_argument('--out_file', type=str, default='traits.csv', help='Path for temporary data')

	args = parser.parse_args()

	# ------------------- Making the traits.csv file

	sample_info_obj = pd.read_csv(args.phenotype_data, delimiter=';')

	# Make variables dichotomized
	encoded_sample_table = pd.get_dummies(sample_info_obj)

	# Alter index
	encoded_sample_table['number'] = 'sample_' + encoded_sample_table['number'].astype(str)

	# Remove index column name
	encoded_sample_table = encoded_sample_table.set_index('number')
	encoded_sample_table.index.name = None

	# Write to file
	encoded_sample_table.to_csv(args.out_file)


