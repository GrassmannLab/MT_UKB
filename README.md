
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
### a) download the l2r, bim and fam files from UKB with ukbgene (files are QC'ed by UKB i.e. do only contain 265 SNPs and not the around 500 MT SNPs that are on the array)
### b) download and install datamash (version 1.4)
### c) download and install tabix and bcftools
### d) download and install haplogrep
### optional: install mosdepth and samtools if you want to compute the coverage of the MT genome from the exome reads


----------------------------------------

##  2. calculate SDLRR from the lrr values of all autosomes 
### this is used as a quality control measure (excluding samples with SDLRR_auto > 0.36)
### Use MoChA to get this measure: https://github.com/freeseek/mocha
 Just a heads-up, using MoChA might take a while to get the data right, probably best to start from the .cel files and use gtc2vcf


----------------------------------------


##  3. compute mtDNA abundance from all 265 SNPs 
## 3.1 compute the median of the L2R values (mtL2R) on the MT chromosome (no in-house QC, equivalent to MoChA output)
 extract dosages so we do not need to have i/o each step
```shell
zcat _001_ukb_l2r_chrMT_v2.txt.gz > tmp.dos
```

 create dummy file with header
```shell
echo 'medLRRMT' > UKBB_CHRMT_medLRRMT.txt
```


 extract 100k samples at a time, otherwise datamash needs more than 100G of RAM
 also substitute spaces with new line to have the data in columns; there is no missing data but --narm is still used just in case
```shell
cut -f1-100000 -d' ' tmp.dos | /datamash-1.4/datamash --narm -t ' ' median 1-100000 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT.txt
cut -f100001-200000 -d' ' tmp.dos | /datamash-1.4/datamash --narm -t ' ' median 1-100000 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT.txt
cut -f200001-300000 -d' ' tmp.dos | /datamash-1.4/datamash --narm -t ' ' median 1-100000 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT.txt
cut -f300001-400000 -d' ' tmp.dos | /datamash-1.4/datamash --narm -t ' ' median 1-100000 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT.txt
cut -f400001-488377 -d' ' tmp.dos | /datamash-1.4/datamash --narm -t ' ' median 1-88377 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT.txt
```

## 3.2 perform QC and compute mtL2R from QC'ed variants

 EXCLUDE SNPS THAT FAILED QC
```shell
zcat _001_ukb_cal_chrMT_v2.bed.gz > _001_ukb_cal_chrMT_v2.bed
zcat _001_ukb_snp_chrMT_v2.bim.gz > _001_ukb_snp_chrMT_v2.bim
plink2 --bfile _001_ukb_cal_chrMT_v2 --fam ukb2222_cal_v2_s488364.fam --bim _001_ukb_snp_chrMT_v2.bim --missing --hardy --freq  --out UKB_MT_QC
```

 extract lines that will be filtered out from LRR file

```shell
R
```
```R
AF=read.table("UKB_MT_QC.afreq", header=TRUE)
v_missing=read.table("UKB_MT_QC.vmiss", header=TRUE)
QC=cbind(AF, v_missing)

exclude=which(QC$F_MISS>0.1|QC$ALT_FREQS>(1-2.083333e-06)|QC$ALT_FREQS<2.083333e-06)
# construct sed command
paste(paste(exclude, "d", sep=""), collapse=";")
q()
```

```shell
echo 'medLRRMT' > UKBB_CHRMT_medLRRMT_QC.txt
```

 sed will delete rows that contain SNPs that failed QC; then datamash will compute the median
```shell
cut -f1-100000 -d' ' tmp.dos | sed '4d;5d;7d;8d;10d;12d;15d;17d;25d;26d;27d;30d;33d;41d;42d;44d;51d;52d;55d;57d;58d;59d;60d;61d;62d;67d;70d;71d;73d;77d;81d;82d;86d;90d;91d;92d;94d;97d;98d;99d;100d;102d;103d;105d;113d;114d;115d;117d;118d;122d;123d;125d;126d;127d;128d;131d;135d;137d;139d;140d;142d;145d;146d;147d;148d;149d;151d;152d;153d;161d;162d;163d;167d;170d;172d;175d;176d;177d;178d;184d;186d;189d;191d;192d;193d;194d;201d;203d;204d;210d;211d;213d;215d;216d;217d;219d;220d;221d;222d;224d;225d;227d;228d;231d;236d;238d;239d;242d;246d;248d;250d;251d;252d;254d;255d;257d;259d;260d;261d;264d;265d' | /datamash-1.4/datamash --narm -t ' ' median 1-100000 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT_QC.txt
cut -f100001-200000 -d' ' tmp.dos | sed '4d;5d;7d;8d;10d;12d;15d;17d;25d;26d;27d;30d;33d;41d;42d;44d;51d;52d;55d;57d;58d;59d;60d;61d;62d;67d;70d;71d;73d;77d;81d;82d;86d;90d;91d;92d;94d;97d;98d;99d;100d;102d;103d;105d;113d;114d;115d;117d;118d;122d;123d;125d;126d;127d;128d;131d;135d;137d;139d;140d;142d;145d;146d;147d;148d;149d;151d;152d;153d;161d;162d;163d;167d;170d;172d;175d;176d;177d;178d;184d;186d;189d;191d;192d;193d;194d;201d;203d;204d;210d;211d;213d;215d;216d;217d;219d;220d;221d;222d;224d;225d;227d;228d;231d;236d;238d;239d;242d;246d;248d;250d;251d;252d;254d;255d;257d;259d;260d;261d;264d;265d' | /datamash-1.4/datamash --narm -t ' ' median 1-100000 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT_QC.txt
cut -f200001-300000 -d' ' tmp.dos | sed '4d;5d;7d;8d;10d;12d;15d;17d;25d;26d;27d;30d;33d;41d;42d;44d;51d;52d;55d;57d;58d;59d;60d;61d;62d;67d;70d;71d;73d;77d;81d;82d;86d;90d;91d;92d;94d;97d;98d;99d;100d;102d;103d;105d;113d;114d;115d;117d;118d;122d;123d;125d;126d;127d;128d;131d;135d;137d;139d;140d;142d;145d;146d;147d;148d;149d;151d;152d;153d;161d;162d;163d;167d;170d;172d;175d;176d;177d;178d;184d;186d;189d;191d;192d;193d;194d;201d;203d;204d;210d;211d;213d;215d;216d;217d;219d;220d;221d;222d;224d;225d;227d;228d;231d;236d;238d;239d;242d;246d;248d;250d;251d;252d;254d;255d;257d;259d;260d;261d;264d;265d' |  /datamash-1.4/datamash --narm -t ' ' median 1-100000 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT_QC.txt
cut -f300001-400000 -d' ' tmp.dos | sed '4d;5d;7d;8d;10d;12d;15d;17d;25d;26d;27d;30d;33d;41d;42d;44d;51d;52d;55d;57d;58d;59d;60d;61d;62d;67d;70d;71d;73d;77d;81d;82d;86d;90d;91d;92d;94d;97d;98d;99d;100d;102d;103d;105d;113d;114d;115d;117d;118d;122d;123d;125d;126d;127d;128d;131d;135d;137d;139d;140d;142d;145d;146d;147d;148d;149d;151d;152d;153d;161d;162d;163d;167d;170d;172d;175d;176d;177d;178d;184d;186d;189d;191d;192d;193d;194d;201d;203d;204d;210d;211d;213d;215d;216d;217d;219d;220d;221d;222d;224d;225d;227d;228d;231d;236d;238d;239d;242d;246d;248d;250d;251d;252d;254d;255d;257d;259d;260d;261d;264d;265d' | /datamash-1.4/datamash --narm -t ' ' median 1-100000 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT_QC.txt
cut -f400001-488377 -d' ' tmp.dos | sed '4d;5d;7d;8d;10d;12d;15d;17d;25d;26d;27d;30d;33d;41d;42d;44d;51d;52d;55d;57d;58d;59d;60d;61d;62d;67d;70d;71d;73d;77d;81d;82d;86d;90d;91d;92d;94d;97d;98d;99d;100d;102d;103d;105d;113d;114d;115d;117d;118d;122d;123d;125d;126d;127d;128d;131d;135d;137d;139d;140d;142d;145d;146d;147d;148d;149d;151d;152d;153d;161d;162d;163d;167d;170d;172d;175d;176d;177d;178d;184d;186d;189d;191d;192d;193d;194d;201d;203d;204d;210d;211d;213d;215d;216d;217d;219d;220d;221d;222d;224d;225d;227d;228d;231d;236d;238d;239d;242d;246d;248d;250d;251d;252d;254d;255d;257d;259d;260d;261d;264d;265d' | /datamash-1.4/datamash --narm -t ' ' median 1-88377 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT_QC.txt
```

## 3.3 compute mtDNA abundance from high quality SNPs/probes (affy 6.0 snps)
 probes that are useful on Affy 6.0 arrays: http://genvisis.umn.edu/rsrc/MitoCN/AffySnp6/gw6_MT_USE.oneHitWonders_20.06_28_16.txt but only 12 are present on the Axiom arrays

 for us, the rows 169p;31p;33p;49p;150p;157p;38p;202p;236p;183p;201p;166p contain the high quality intensities; those will be extracted by sed


```shell
echo 'medLRRMT' > UKBB_CHRMT_medLRRMT_Genivis.txt
cut -f1-100000 -d' ' tmp.dos | sed -n '169p;31p;33p;49p;150p;157p;38p;202p;236p;183p;201p;166p' | /datamash-1.4/datamash --narm -t ' ' median 1-100000 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT_Genivis.txt
cut -f100001-200000 -d' ' tmp.dos | sed -n  '169p;31p;33p;49p;150p;157p;38p;202p;236p;183p;201p;166p' | /datamash-1.4/datamash --narm -t ' ' median 1-100000 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT_Genivis.txt
cut -f200001-300000 -d' ' tmp.dos | sed -n '169p;31p;33p;49p;150p;157p;38p;202p;236p;183p;201p;166p' | /datamash-1.4/datamash --narm -t ' ' median 1-100000 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT_Genivis.txt
cut -f300001-400000 -d' ' tmp.dos | sed -n '169p;31p;33p;49p;150p;157p;38p;202p;236p;183p;201p;166p' | /datamash-1.4/datamash --narm -t ' ' median 1-100000 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT_Genivis.txt
cut -f400001-488377 -d' ' tmp.dos | sed -n '169p;31p;33p;49p;150p;157p;38p;202p;236p;183p;201p;166p' | /datamash-1.4/datamash --narm -t ' ' median 1-88377 | tr -s ' '  '\n' >> UKBB_CHRMT_medLRRMT_Genivis.txt
```

## 3.4 compute mtDNA abundance from weighted intensities
```shell
echo 'medLRRMT' > UKBB_CHRMT_medLRRMT_QC_Exome.txt
```
 extract the corresponding probes with sed, this may be the same for you but you should check (see supplementary table S1, i.e. here we keep probe 4,19,20,21,24 ...).

 score.txt has the weight for each probe. the intensity of the probe will be multiplied with the weight and then datamash computes the mean out of that

```shell
zcat _001_ukb_l2r_chrMT_v2.txt.gz | cut -f1-100000 -d' '  | sed '1d;2d;3d;5d;6d;7d;8d;9d;10d;11d;12d;13d;14d;15d;16d;17d;18d;22d;24d;25d;26d;27d;29d;30d;33d;34d;35d;36d;37d;39d;40d;41d;42d;44d;45d;46d;47d;50d;51d;52d;53d;54d;55d;57d;58d;59d;60d;62d;64d;66d;67d;69d;70d;71d;72d;73d;74d;75d;76d;77d;78d;80d;81d;82d;84d;86d;87d;88d;89d;90d;91d;92d;94d;95d;96d;97d;98d;99d;100d;101d;102d;103d;104d;105d;106d;107d;108d;110d;112d;113d;114d;115d;117d;118d;119d;120d;121d;122d;123d;124d;125d;126d;127d;128d;129d;133d;134d;135d;136d;137d;138d;139d;140d;141d;142d;143d;145d;146d;147d;148d;149d;150d;151d;152d;153d;155d;157d;159d;160d;161d;162d;163d;167d;168d;169d;170d;172d;173d;174d;175d;176d;177d;178d;179d;180d;181d;182d;183d;184d;185d;186d;187d;188d;189d;190d;191d;192d;193d;194d;195d;196d;197d;198d;199d;200d;201d;202d;203d;204d;207d;208d;209d;210d;211d;212d;213d;215d;216d;218d;219d;220d;221d;222d;223d;224d;225d;226d;227d;228d;229d;230d;231d;232d;233d;235d;236d;238d;239d;241d;242d;245d;246d;247d;248d;249d;250d;251d;252d;253d;254d;255d;256d;257d;258d;259d;260d;261d;262d;264d;265d' | /datamash-1.4/datamash --narm -t ' ' transpose > tmp.dose
awk 'NR==FNR{A[NR]=$1; next} FNR==1{n=NF}{for(i=1; i<=n; i++) $(i)=$i=="NA"?$i:$i*A[i]}1' score.txt tmp.dose > tmp.dose.adj
/datamash-1.4/datamash --narm -t ' ' transpose < tmp.dose.adj | /datamash-1.4/datamash --narm -t ' ' mean 1-100000 | tr -s ' '  '\n'   >> UKBB_CHRMT_medLRRMT_QC_Exome.txt

zcat _001_ukb_l2r_chrMT_v2.txt.gz | cut -f100001-200000 -d' '  | sed '1d;2d;3d;5d;6d;7d;8d;9d;10d;11d;12d;13d;14d;15d;16d;17d;18d;22d;24d;25d;26d;27d;29d;30d;33d;34d;35d;36d;37d;39d;40d;41d;42d;44d;45d;46d;47d;50d;51d;52d;53d;54d;55d;57d;58d;59d;60d;62d;64d;66d;67d;69d;70d;71d;72d;73d;74d;75d;76d;77d;78d;80d;81d;82d;84d;86d;87d;88d;89d;90d;91d;92d;94d;95d;96d;97d;98d;99d;100d;101d;102d;103d;104d;105d;106d;107d;108d;110d;112d;113d;114d;115d;117d;118d;119d;120d;121d;122d;123d;124d;125d;126d;127d;128d;129d;133d;134d;135d;136d;137d;138d;139d;140d;141d;142d;143d;145d;146d;147d;148d;149d;150d;151d;152d;153d;155d;157d;159d;160d;161d;162d;163d;167d;168d;169d;170d;172d;173d;174d;175d;176d;177d;178d;179d;180d;181d;182d;183d;184d;185d;186d;187d;188d;189d;190d;191d;192d;193d;194d;195d;196d;197d;198d;199d;200d;201d;202d;203d;204d;207d;208d;209d;210d;211d;212d;213d;215d;216d;218d;219d;220d;221d;222d;223d;224d;225d;226d;227d;228d;229d;230d;231d;232d;233d;235d;236d;238d;239d;241d;242d;245d;246d;247d;248d;249d;250d;251d;252d;253d;254d;255d;256d;257d;258d;259d;260d;261d;262d;264d;265d' | /datamash-1.4/datamash --narm -t ' ' transpose > tmp.dose
awk 'NR==FNR{A[NR]=$1; next} FNR==1{n=NF}{for(i=1; i<=n; i++) $(i)=$i=="NA"?$i:$i*A[i]}1' score.txt tmp.dose > tmp.dose.adj
/datamash-1.4/datamash --narm -t ' ' transpose < tmp.dose.adj | /datamash-1.4/datamash --narm -t ' ' mean 1-100000 | tr -s ' '  '\n'   >> UKBB_CHRMT_medLRRMT_QC_Exome.txt

zcat _001_ukb_l2r_chrMT_v2.txt.gz | cut -f200001-300000 -d' '  | sed '1d;2d;3d;5d;6d;7d;8d;9d;10d;11d;12d;13d;14d;15d;16d;17d;18d;22d;24d;25d;26d;27d;29d;30d;33d;34d;35d;36d;37d;39d;40d;41d;42d;44d;45d;46d;47d;50d;51d;52d;53d;54d;55d;57d;58d;59d;60d;62d;64d;66d;67d;69d;70d;71d;72d;73d;74d;75d;76d;77d;78d;80d;81d;82d;84d;86d;87d;88d;89d;90d;91d;92d;94d;95d;96d;97d;98d;99d;100d;101d;102d;103d;104d;105d;106d;107d;108d;110d;112d;113d;114d;115d;117d;118d;119d;120d;121d;122d;123d;124d;125d;126d;127d;128d;129d;133d;134d;135d;136d;137d;138d;139d;140d;141d;142d;143d;145d;146d;147d;148d;149d;150d;151d;152d;153d;155d;157d;159d;160d;161d;162d;163d;167d;168d;169d;170d;172d;173d;174d;175d;176d;177d;178d;179d;180d;181d;182d;183d;184d;185d;186d;187d;188d;189d;190d;191d;192d;193d;194d;195d;196d;197d;198d;199d;200d;201d;202d;203d;204d;207d;208d;209d;210d;211d;212d;213d;215d;216d;218d;219d;220d;221d;222d;223d;224d;225d;226d;227d;228d;229d;230d;231d;232d;233d;235d;236d;238d;239d;241d;242d;245d;246d;247d;248d;249d;250d;251d;252d;253d;254d;255d;256d;257d;258d;259d;260d;261d;262d;264d;265d' | /datamash-1.4/datamash --narm -t ' ' transpose > tmp.dose
awk 'NR==FNR{A[NR]=$1; next} FNR==1{n=NF}{for(i=1; i<=n; i++) $(i)=$i=="NA"?$i:$i*A[i]}1' score.txt tmp.dose > tmp.dose.adj
/datamash-1.4/datamash --narm -t ' ' transpose < tmp.dose.adj | /datamash-1.4/datamash --narm -t ' ' mean 1-100000 | tr -s ' '  '\n'   >> UKBB_CHRMT_medLRRMT_QC_Exome.txt

zcat _001_ukb_l2r_chrMT_v2.txt.gz | cut -f300001-400000 -d' '  | sed '1d;2d;3d;5d;6d;7d;8d;9d;10d;11d;12d;13d;14d;15d;16d;17d;18d;22d;24d;25d;26d;27d;29d;30d;33d;34d;35d;36d;37d;39d;40d;41d;42d;44d;45d;46d;47d;50d;51d;52d;53d;54d;55d;57d;58d;59d;60d;62d;64d;66d;67d;69d;70d;71d;72d;73d;74d;75d;76d;77d;78d;80d;81d;82d;84d;86d;87d;88d;89d;90d;91d;92d;94d;95d;96d;97d;98d;99d;100d;101d;102d;103d;104d;105d;106d;107d;108d;110d;112d;113d;114d;115d;117d;118d;119d;120d;121d;122d;123d;124d;125d;126d;127d;128d;129d;133d;134d;135d;136d;137d;138d;139d;140d;141d;142d;143d;145d;146d;147d;148d;149d;150d;151d;152d;153d;155d;157d;159d;160d;161d;162d;163d;167d;168d;169d;170d;172d;173d;174d;175d;176d;177d;178d;179d;180d;181d;182d;183d;184d;185d;186d;187d;188d;189d;190d;191d;192d;193d;194d;195d;196d;197d;198d;199d;200d;201d;202d;203d;204d;207d;208d;209d;210d;211d;212d;213d;215d;216d;218d;219d;220d;221d;222d;223d;224d;225d;226d;227d;228d;229d;230d;231d;232d;233d;235d;236d;238d;239d;241d;242d;245d;246d;247d;248d;249d;250d;251d;252d;253d;254d;255d;256d;257d;258d;259d;260d;261d;262d;264d;265d' | /datamash-1.4/datamash --narm -t ' ' transpose > tmp.dose
awk 'NR==FNR{A[NR]=$1; next} FNR==1{n=NF}{for(i=1; i<=n; i++) $(i)=$i=="NA"?$i:$i*A[i]}1' score.txt tmp.dose > tmp.dose.adj
/datamash-1.4/datamash --narm -t ' ' transpose < tmp.dose.adj | /datamash-1.4/datamash --narm -t ' ' mean 1-100000 | tr -s ' '  '\n'   >> UKBB_CHRMT_medLRRMT_QC_Exome.txt

zcat _001_ukb_l2r_chrMT_v2.txt.gz | cut -f400001-488377 -d' '  | sed '1d;2d;3d;5d;6d;7d;8d;9d;10d;11d;12d;13d;14d;15d;16d;17d;18d;22d;24d;25d;26d;27d;29d;30d;33d;34d;35d;36d;37d;39d;40d;41d;42d;44d;45d;46d;47d;50d;51d;52d;53d;54d;55d;57d;58d;59d;60d;62d;64d;66d;67d;69d;70d;71d;72d;73d;74d;75d;76d;77d;78d;80d;81d;82d;84d;86d;87d;88d;89d;90d;91d;92d;94d;95d;96d;97d;98d;99d;100d;101d;102d;103d;104d;105d;106d;107d;108d;110d;112d;113d;114d;115d;117d;118d;119d;120d;121d;122d;123d;124d;125d;126d;127d;128d;129d;133d;134d;135d;136d;137d;138d;139d;140d;141d;142d;143d;145d;146d;147d;148d;149d;150d;151d;152d;153d;155d;157d;159d;160d;161d;162d;163d;167d;168d;169d;170d;172d;173d;174d;175d;176d;177d;178d;179d;180d;181d;182d;183d;184d;185d;186d;187d;188d;189d;190d;191d;192d;193d;194d;195d;196d;197d;198d;199d;200d;201d;202d;203d;204d;207d;208d;209d;210d;211d;212d;213d;215d;216d;218d;219d;220d;221d;222d;223d;224d;225d;226d;227d;228d;229d;230d;231d;232d;233d;235d;236d;238d;239d;241d;242d;245d;246d;247d;248d;249d;250d;251d;252d;253d;254d;255d;256d;257d;258d;259d;260d;261d;262d;264d;265d' | /datamash-1.4/datamash --narm -t ' ' transpose > tmp.dose
awk 'NR==FNR{A[NR]=$1; next} FNR==1{n=NF}{for(i=1; i<=n; i++) $(i)=$i=="NA"?$i:$i*A[i]}1' score.txt tmp.dose > tmp.dose.adj
/datamash-1.4/datamash --narm -t ' ' transpose < tmp.dose.adj | /datamash-1.4/datamash --narm -t ' ' mean 1-88377 | tr -s ' '  '\n'   >> UKBB_CHRMT_medLRRMT_QC_Exome.txt
```

## 3.5 compute mtDNA abundance from the coverage of the MT genome
 download the CRAM files (around 5 TB in total for 50k participants) of each participant

 index with samtools

```shell
find . -name "*.cram" | xargs -n1 -P20 -I{} samtools index {}
```

 run mosdepth on each file to compute coverage
```shell
find . -name "*.cram" | xargs -n1 -P20 -I{} mosdepth -t 10  -n  --fasta /nfs/home/felgra/bin/mosdepth/hg38.fa.gz --fast-mode {}
```

 aggegrate the coverage of each individual file into one file (I did this with a clunky R script, maybe you can find a better solution)
```shell
R
```
```R
library(data.table)
samples_with_exome=fread("samples_with_exome_data.txt", data.table=FALSE)
samples_with_exome=samples_with_exome[which(samples_with_exome$SPB_CRAM_files!=""),]

samples_with_exome$Total=NA
samples_with_exome$X=NA
samples_with_exome$Y=NA
samples_with_exome$MT=NA
samples_with_exome$CHR1=NA
samples_with_exome$AUTO=NA

file.names=list.files(pattern="summary.txt")
for(file.name in file.names){
    sample_id=gsub("\\..*", "", file.name)
    chr_all=fread(file=paste(sample_id, ".mosdepth.summary.txt", sep=""), data.table=FALSE)
    samples_with_exome[match(sample_id, samples_with_exome$ID), "X"]=chr_all[chr_all$chrom=="X","mean"]
    samples_with_exome[match(sample_id, samples_with_exome$ID), "CHR1"]=chr_all[chr_all$chrom=="1","mean"]
    samples_with_exome[match(sample_id, samples_with_exome$ID), "AUTO"]=mean(chr_all[match(1:22, chr_all$chrom),"mean"])
    samples_with_exome[match(sample_id, samples_with_exome$ID), "Y"]=chr_all[chr_all$chrom=="Y","mean"]
    samples_with_exome[match(sample_id, samples_with_exome$ID), "MT"]=chr_all[chr_all$chrom=="M","mean"]
    samples_with_exome[match(sample_id, samples_with_exome$ID), "Total"]=chr_all[chr_all$chrom=="total","mean"]
    
    }
write.table(samples_with_exome, file="coverage_in_samples_with_exome.txt")
q()
```




----------------------------------------

##   4. Compute haplogroups with haplogrep 

 haplogrep needs an annotated vcf (containing the same rsIDs as haplogrep expects)

 first convert the UKB snp file in plink format to vcf
```shell
plink --bfile _001_ukb_snp_chrMT_v2 --recode vcf --out _001_ukb_snp_chrMT_v2
```

create bed file with the necessary identifiers that are the same as in the vcf supplied by haplogrep
```shell
awk -v OFS='\t' '{print $1, $2, $2, $3}' haplogrep/HG00097.vcf | grep rs > MT_SNPs_new.bed
sed -i 's/MT/26/g' MT_SNPs_new.bed
```
add SNP identifiers to the vcf with bcftools
```shell
bgzip -c MT_SNPs_new.bed > MT_SNPs_new.bed.gz
tabix -p bed MT_SNPs_new.bed.gz
bcftools annotate -c CHROM,FROM,TO,ID -a MT_SNPs_new.bed.gz  -o v --output _001_ukb_snp_chrMT_v2_annotated.vcf  _001_ukb_snp_chrMT_v2.vcf.gz
```
replace CHR 26 with MT
```shell
sed -i 's/^26/MT/g' _001_ukb_snp_chrMT_v2_annotated.vcf 
```
run haplogrep on the annotated vcf
```shell
java -jar /haplogrep/haplogrep-2.1.25.jar --in  _001_ukb_snp_chrMT_v2_annotated.vcf --format vcf --out haplogroups_UKB.txt --chip
```



----------------------------------------

##   5. Quality control: 

### 1. exclude people with sex chromosome anomalies (!is.na(sex_aneu)) and participants whose genetic sex does not match the reported (MOCHA_SEX != UK Biobank phenotype sex)
### 2. exclude individuals that failed genotyping QC: gen_exclude
### 3. exclude individuals with low call rate: gen_missing>0.01
### 4. exclude related individuals (used_in_PCA!=1)
### 5. exclude samples with high SD of autosomal probes SDLLR_auto>0.369986958592406 (i.e. two SD from the mean)
### 6. exclude non-europeans (gen_ethnicity!=1)
### 7. exclude individuals with high BAF_CONCS (>0.52)
### in general, you should  end up with the similar numbers to the one reported in the manuscript. if not, please contact me!


#########################################
##   6. normalize mLRRMT per genotyping plate (variable "22007-0.0" in UKB) 
 this effectively creates the final outcome or exposure for the analyses
 here, we used the weighted mLRRMT value (mLRRMT)

```shell
R
```

```R
phenotypes$mLRRMT_norm=NA
for(plate in unique(phenotypes$genotyping_plate)[order(unique(phenotypes$genotyping_plate))]){
    phenotypes[phenotypes$genotyping_plate==plate, "mLRRMT_norm"]=scale(phenotypes[phenotypes$genotyping_plate==plate, "mLRRMT"])
}


summary(phenotypes$mLRRMT_norm)

     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
 -5.139131 -0.661685 -0.007579  0.000000  0.652375  4.790622 
```


--------------
## -------------- END --------------
 
