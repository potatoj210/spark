# spark
QC and pre-imputation for SPARK 
##keep EUR ancestry and sample QC
```
#BASENAME
inputbed="SPARK.iWES_v2.genotyping_array"
root="wes2_array"
dir="./tempfile"
bed=${dir}/${root}
plink \
--bfile ${inputbed} \
--chr 1-23 \
--keep /home/qm226/spark_qc/script/GWAS_ASD_SPARK/Pre-Imputation/temp/EUR_keep_array.txt \
--make-bed \
--out ${bed}_qc1

plink \
--bfile ${bed}_qc1 \
--geno 0.05 \
--make-bed \
--out ${bed}_qc2
###missingness and frequency by batch 
plink \
--bfile ${bed}_qc1 \
--missing \
--within \
--out ${bed}_qc2

plink \
--bfile ${bed}_qc1 \
--freq \
--within \
--out ${bed}_qc2

plink \
--bfile ${bed}_qc2 \
--mind 0.1 \
--me 0.05 0.1 \
--make-bed \
--out ${bed}_qc3

plink \
--bfile ${bed}_qc2 \
--het \
--out ${bed}_qc3

plink \
--bfile ${bed}_qc2 \
--check-sex \
--out ${bed}_qc3
```
##check het and sex
```
library(data.table)
het = fread("/rds/user/qm226/rds-genetics_hpc-Nl99R8pHODQ/SPARK_data/SPARK_0223_iWES2genotype/Genotyped/tempfile/wes2_array_qc3.het")
het$HET = (het$`N(NM)` - het$`O(HOM)`)/het$`N(NM)` #create heterozygosity stats
mean = mean(het$HET)
sd = sd(het$HET)
het$Z = (het$HET - mean)/sd #create Z scores of heterozygosity
hetoutlier = subset(het, abs(Z) > 3)
failedsample = hetoutlier[,c(1:2)]

sex_check = fread("/rds/user/qm226/rds-genetics_hpc-Nl99R8pHODQ/SPARK_data/SPARK_0223_iWES2genotype/Genotyped/tempfile/wes2_array_qc3.sexcheck")
problem = subset(sex_check, STATUS == "PROBLEM")
problem2 = problem[,c(1:2)]

failedsample2 = rbind(failedsample, problem2)
failedsample2 = unique(failedsample2)

write.table(failedsample2, file = "./temp/wes2_array_failedsample_het_sex.txt", row.names = F, col.names = T, quote = F)
```
##snp QC
```
plink \
--bfile ${bed}_qc4 \
--hwe 0.000001 \
--maf 0.01 \
--make-bed \
--out ${bed}_qc5
```
##exclude twins
```
plink \
--bfile ${bed}_qc5 \
--remove-fam /home/qm226/spark_qc/script/GWAS_ASD_SPARK/Pre-Imputation/temp/twins.txt \
--make-bed \
--out ${bed}_notwin_maf001
```
##batch QC
```
#BASENAME
inputbed="SPARK.iWES_v2.genotyping_array"
root="wes2_array_notwin"
dir="./tempfile"
bed=${dir}/${root}

plink \
--bfile ${bed}_maf001 \
--exclude ./tempfile/snplist.txt \
--make-bed \
--out ${bed}_maf001_batch
```
##generating pseudocontrols
```
plink \
--bfile ${bed}_notwin_maf001_batch \
--keep /home/qm226/spark/script/Pre-Imputation/temp/trio_families2.list \
--make-bed \
--out ${bed}_qcedtrios3

plink \
--bfile ${bed}_qcedtrios3 \
--tucc write-bed \
--out ${bed}_qcedtrios_pseudocontrols3
```
##select cases and controls
```
#FOR duos
plink \
--bfile wes2_array_notwin_maf001_batch \
--keep /home/qm226/spark/script/Pre-Imputation/temp/duolist1.txt \
--make-bed \
--out duo
#FOR singleton
plink \
--bfile wes2_array_notwin_maf001_batch \
--keep /home/qm226/spark/script/Pre-Imputation/temp/singlelist.txt \
--make-bed \
--out sing
#singleton duo and trio merge
plink \
--merge-list /home/qm226/spark/script/Pre-Imputation/temp/merge_list.txt \
--make-bed \
--out case_control
```
##prepare for topmed
```
###update snp names
library(data.table)
bim = fread("/home/qm226/rds/rds-genetics_hpc-Nl99R8pHODQ/SPARK_data/SPARK_0223_iWES2genotype/Genotyped/tempfile/case_control.bim")
bim$update <- gsub("GSA-", "", bim$V2)
write.table(bim[,c("V2", "update")], file = "./temp/updatenames.txt", row.names = F, col.names = F, quote = F)
```
```
plink \
--bfile case_control \
--update-name /home/qm226/spark/script/Pre-Imputation/temp/updatenames.txt \
--make-bed \
--out case_control_update
```
```
perl HRC-1000G-check-bim.pl -b /rds/user/qm226/rds-genetics_hpc-Nl99R8pHODQ/SPARK_data/SPARK_0223_iWES2genotype/Genotyped/tempfile/case_control_update.bim -f /rds/user/qm226/rds-genetics_hpc-Nl99R8pHODQ/SPARK_data/SPARK_0223_iWES2genotype/Genotyped/tempfile/case_control_update.frq -r /rds/user/qm226/rds-genetics_hpc-Nl99R8pHODQ/SPARK_data/TOPMED_referencepanel/PASS.Variantsbravo-dbsnp-all.tab
plink \
--bfile case_control_update \
--exclude /rds/user/qm226/rds-genetics_hpc-Nl99R8pHODQ/SPARK_data/TOPMED_referencepanel/Exclude-case_control_update-HRC.txt \
--make-bed \
--out TEMP1

plink \
--bfile TEMP1 \
--a2-allele /rds/user/qm226/rds-genetics_hpc-Nl99R8pHODQ/SPARK_data/TOPMED_referencepanel/Force-Allele1-case_control_update-HRC.txt \
--make-bed \
--out case_control_update4

for i in {1..23};do
plink \
--bfile case_control_update4 \
--a2-allele /rds/user/qm226/rds-genetics_hpc-Nl99R8pHODQ/SPARK_data/TOPMED_referencepanel/Force-Allele1-case_control_update-HRC.txt \
--chr ${i} \
--recode vcf \
--output-chr chr26 \
--out ./chr_vcf/test_chr${i}
done

for i in {1..23};do
bcftools sort test_chr${i}.vcf -Oz -o test_chr${i}.vcf.gz
done
```
