###Code from https://github.com/MorganLevineLab/methylCIPHER
# Paper https://www.biorxiv.org/content/10.1101/2022.07.13.499978v1

#instalar paquetes cargar librerias
# devtools::install_github("MorganLevineLab/prcPhenoAge")
# devtools::install_github("danbelsky/DunedinPoAm38")
# devtools::install_github("MorganLevineLab/methylCIPHER")


library(methylCIPHER)

#Cargar datos
Betas_data=c("met.normal_filter_by_gene_TCGA_BRCA.tsv")
Pheno_data=c("clin.normal_TCGA_BRCA.tsv")
Data.name=c("methylCIPHER_TCGA_BRCA_Normal")


#Read tables
Betas<-read.table(Betas_data)
Pheno<- read.table(Pheno_data,header=T,sep='\t')

# 
# for ( col in 1:ncol(Betas)){
#   colnames(Betas)[col] <-  sub(".01.*", ".01", colnames(Betas)[col])
# }
# Pheno$samples=gsub ("-01.*", "-01", Pheno$samples) 
# Pheno$samples=gsub("-", ".",  Pheno$samples)
# names= colnames(Betas)
# common <- intersect(names, Pheno$samples)  


row.names(Betas)<-Betas$CpGName

#Eliminar la columna CpGname

Betas1 <- Betas[,-1]

#Trasponer como ejemplo
library(data.table)
t_Betas <- transpose(Betas1)
#Devolver los colnames y rownames

rownames(t_Betas) <- colnames(Betas1)
colnames(t_Betas) <- rownames(Betas1)

#PHENODATA#####
#Ordenar columna
library(dplyr)
Pheno_age<-Pheno%>% select (patient, tissue, samples,age_at_index)


###########Saber el porcentaje de CpGs en una muestra para un reloj#########################
CpG_in_Clocks<-getClockProbes(t_Betas)
CpG_in_Clocks
#Reloj
#Pheno_age
Pheno_age<-calcPhenoAge(t_Betas, Pheno_age, imputation = F)
#Running Multiple Epigenetic Clocks Simultaneously
Multiple_clocks<-calcCoreClocks(t_Betas, Pheno_age)
#Multiple_clocks
write.table(Multiple_clocks, paste(Data.name, ".txt",sep=""), row.names = FALSE, dec=".", quote = FALSE,sep = "\t")
write.table(Pheno_age, paste(Data.name, ".txt",sep=""), row.names = FALSE, dec=".", quote = FALSE,sep = "\t")

##### Customize
##### Verificar si hay columnas que contengan NAs para los CpG nombrados como faltantes lo cuales no contarán
#Si obtiene "entero con nombre (0)", entonces no tiene ninguno y ok!!!
#which(apply(t_Betas, 2, function(x)all(is.na(x))))

#In the case that you have CpGs missing from only some samples, 
#Definir si  falten CpG en algunas muestras, Si esta no sale en 0 jecutar la imputación media 
#dentro de sus datos para que los valores de NA para muestras únicas o pocas al menos tengan 
#valores medios en lugar de ser ignorados.

sum(is.na(t_Betas))
##
meanimpute <- function(x){
  apply(x,2,function(z)ifelse(is.na(z),mean(z,na.rm=T),z))
}
mean=meanimpute(t_Betas)
sum(is.na(mean))
Pheno_age2<-calcPhenoAge(mean, Pheno_age, imputation = F)
Multiple_clocks2<-calcCoreClocks(mean, Pheno_age)
#Multiple_clocks
write.table(Multiple_clocks2, paste(Data.name, "_Imputation.txt",sep=""), row.names = FALSE, dec=".", quote = FALSE,sep = "\t")
write.table(Pheno_age2, paste(Data.name, "_Imputation.txt",sep=""), row.names = FALSE, dec=".", quote = FALSE,sep = "\t")
