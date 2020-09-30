import csv


path_to_roary_output = '/Volumes/External/CIDRI_Data/felix_S_pneumoniae/roary/roary/gene_presence_absence.csv'

path_to_phenotype_data = '/Volumes/External/CIDRI_Data/felix_S_pneumoniae/sample_sheet_FsD.csv'

path_to_sample_sheet = '/Volumes/External/CIDRI_Data/felix_S_pneumoniae/sample_sheet.csv'

input_file = csv.DictReader(open(path_to_roary_output))

for a_gene in input_file:

	gene_results = {}

	for a_sample in a_gene.keys():
		if a_sample.split('_')[0] == 'sample':
			sample_ID = a_sample
			if len(a_gene[a_sample]) > 1:

				# Gene present for this sample

				1 == 1


