# ABO_Allele_Caller
This repo is the (in progress) data for the ABO allele sequence alignment caller.

## About the Repo
Introduction

This R script is designed to analyze genetic data related to ABO alleles. It performs multiple tasks, including data parsing, sequence processing, simulation of patient genotypes, and calling ABO alleles for simulated patients.

The original ABO alleles are extracted from BD MUT Database of allelic variants of Human Blood Group Antigens (doi: 10.1159/000366108). A now defunct website only accessible via FTP (ftp://ftp.ncbi.nlm.nih.gov/pub/mhc/mhc/Final Archive). 


## Organization
Each folder contains a readme document explaining in greater detail the contents. The `code` folder contains scripts to grab the raw data, clean, and output. It contains scripts for the caller to be ran locally and on an HPC. The `raw_data` folder contains the original archived ABO alleles. `Results` folder contains output of the ABO sequence alignment algorithm. `Simulated_genotypes` folder contains the mutated genotypes used to test the performance of the caller on "patients" with known sequences. 

## Contributing

This is an active project, contributions are welcome but may not be appropriate at this time. 

## License

original code and configuration are under the MIT License.