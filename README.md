
---------------------------
# COMPUTE MT DNA ABUNDACE FROM WEIGHTED PROBES 
## also: compute QC measures                      
### Author: Felix Grassmann, UoA, 2020           
### GitHub version 1.0
### this is the companion script to: https://pubmed.ncbi.nlm.nih.gov/33385171/
## Deciphering the genetic and epidemiological landscape of mitochondrial DNA abundance.
### Hägg S, Jylhävä J, Wang Y, Czene K, Grassmann F. Hum Genet. 2020 Dec 31. doi: 10.1007/s00439-020-02249-w. PMID: 33385171 
------------------------------------------


##  1. first steps: 
 a) download the l2r, bim and fam files from UKB with ukbgene (files are QC'ed by UKB i.e. do only contain 265 SNPs and not the around 500 MT SNPs that are on the array)
 b) download and install datamash (version 1.4)
 c) download and install tabix and bcftools
 d) download and install haplogrep
 optional: install mosdepth and samtools if you want to compute the coverage of the MT genome from the exome reads


----------------------------------------

##  2. calculate SDLRR from the lrr values of all autosomes 
 this is used as a quality control measure (excluding samples with SDLRR_auto > 0.36)
 Use MoChA to get this measure: https://github.com/freeseek/mocha
 Just a heads-up, using MoChA might take a while to get the data right, probably best to start from the .cel files and use gtc2vcf

