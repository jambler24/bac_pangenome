# ![uct-cbio/bacterial_variant_calling](/assets/cbio_logo.png)

#Nextflow pipeline for bacterial genome assembly and pangenomics

## Quickstart 

    nextflow run jambler24/bac_pangenome --reads sample_sheet.csv --genome Streptococcus_pneumoniae.fa 
    
## Basic usage: 
The typical command for running the pipeline is as follows:

    nextflow run jambler24/bac_pangenome --reads sample_sheet.csv --genome refgenome.fa -profile ilifu
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Path to reference genome against which the reads will be aligned (in fasta format) for use in QC steps.
      -profile                      Hardware config to use. Currently profile available for ilifu and UCT's HPC 'uct_hex' - create your own if necessary
      
      
    Other arguments:
      --outdir                      The output directory where the results will be saved
      --SRAdir                      The directory where reads downloaded from the SRA will be stored
      --phenotype_info              Dichotomised CSV file containing traits associated with isolates for bacterial GWAS with scoary
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name      

## Sample file
To allow for both local reads and reads from the [SRA](https://www.ncbi.nlm.nih.gov/sra) to be used, the pipeline has the 
ability to pull reads from the SRA based on the accession number (eg, [SRR5989977](https://www.ncbi.nlm.nih.gov/sra/SRX3145707[accn])). 

The 'number' column must contain a unique value. 

number | origin | replicate | isolate | R1 | R2
------------ | ------------- | ------------- | ------------- | ------------- | -------------
1 | genomic | 1 | wgs_sample_1 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
2 | genomic | 2 | wgs_sample_1 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
3 | genomic | 3 | wgs_sample_1 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
4 | genomic | 1 | wgs_sample_2 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
5 | genomic | 2 | wgs_sample_2 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
6 | genomic | 3 | wgs_sample_2 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
7 | genomic | 1 | wgs_sample_3 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
8 | genomic | 2 | wgs_sample_3 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
9 | genomic | 3 | wgs_sample_3 | path/to/reads/reads_R1.fq | path/to/reads/reads_R2.fq
10 | genomic | 1 | H37Rv | SRR5989977 | 


In the above example, samples 1-9 are locally stored where sample 10 is a control sample from the SRA. 
Including the accession number in the R1 column will result in the reads from the SRA to be downloaded and used in the analysis. 
This must be exported to a csv file, with a comma ',' separating the columns:

    number,origin,replicate,isolate,R1,R2
    1,genomic,1,wgs_sample_1,path/to/reads/reads_R1.fq,path/to/reads/reads_R2.fq
    2,genomic,2,wgs_sample_1,path/to/reads/reads_R1.fq,path/to/reads/reads_R2.fq
    ...
    10,genomic,1,H37Rv,SRR5989977


## Documentation

Read mapping:  [BWA](http://bio-bwa.sourceforge.net)

Genome assembly: [Unicycler](https://github.com/rrwick/Unicycler)

Genome annotation: [Prokka](https://sco.h-its.org/exelixis/web/software/raxml/index.html)

Pan genome analysis: [Roary](https://sanger-pathogens.github.io/Roary/)

Bacterial GWAS: [Scoary](https://github.com/AdmiralenOla/Scoary)

## Other useful info

To create a Singularity image from a Docker image, please make use of 
[Docker to singularity](https://github.com/singularityware/docker2singularity). This is needed to run the pipeline on the
UCT cluster. 

## Built With
* [Nextflow](https://www.nextflow.io/)
* [Docker](https://www.docker.com/what-docker)
* [Singularity](https://singularity.lbl.gov/)

## Credits
This pipeline was developed by members of the Bioinformatics Support Team (BST) at the University of Cape Town. Dr.
Jon Ambler is a member of CIDRI-Africa, and the main developer of this pipeline, using the layout and documentation
 outlined by Dr Katie Lennard and Gerrit Botha. Adapted from the nf-core bacass pipeline. 

Additional thanks to Paolo Di Tommaso, the developer of NextFlow, for their help troubleshooting. 

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details


