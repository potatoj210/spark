# spark
QC and pre-imputation for SPARK
1. update snp names
```
library(data.table)
bim = fread("/home/qm226/rds/rds-genetics_hpc-Nl99R8pHODQ/SPARK_data/SPARK_0223_iWES2genotype/Genotyped/tempfile/case_control.bim")
bim$update <- gsub("GSA-", "", bim$V2)
write.table(bim[,c("V2", "update")], file = "./temp/updatenames.txt", row.names = F, col.names = F, quote = F)
```
```

```
