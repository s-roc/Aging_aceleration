
###https://www.bioconductor.org/packages/release/bioc/vignettes/methylclock/inst/doc/methylclock.html

library(methylclock)
library(ggplot2)

# ELEGIR MATRIZ
#met.normal_filter_by_gene_TCGA_BRCA.tsv
#met.tumor_filter.by_gene_TCGA_BRCA.tsv

name.met= ("met.normal_filter_by_gene_TCGA_BRCA.tsv")
name.clin= ("clin.normal_TCGA_BRCA.tsv")
name.output=("Normal_BRCA_TCGA")
met=read.table(name.met,header=T,sep='\t')
clin=read.table(name.clin,header=T,sep='\t')


colnames(met)<-gsub("\\." ,  "\\-" , colnames(met))

met<-met[order(colnames(met))]
clin<-clin[order(clin$patient),]


##############      RELOJES    #################################
age<-clin$age_at_index 
cpgs.missing <- checkClocks(met)
eClocks<-DNAmAge(met, normalize = TRUE)
plotDNAmAge(eClocks$Levine, age)
#plotCorClocks(eClocks)

#Se puede elegir solo un Reloj
#eClocks<-DNAmAge(met.t, normalize = TRUE, clocks = "Horvath")

#Guardar
write.table(eClocks,file=paste("eclocks",name.output,".tsv",sep=""),sep='\t',quote=TRUE)

######################## AGE ACELERATION #####################################
#   ageAcc: Difference between DNAmAge and chronological age.
#   ageAcc2: Residuals obtained after regressing chronological age and DNAmAge 
#            (similar to IEAA).
#   ageAcc3: Residuals obtained after regressing chronological age and DNAmAge 
#            adjusted for cell counts (similar to EEAA).
###############################################################################


acc<- DNAmAge(met, age=age, normalize = TRUE, clocks = "Horvath")

#acc<- DNAmAge(p, age=age, normalize = TRUE, clocks = "Horvath")



write.table(acc,paste("acc",name.output,".tsv",sep=""),sep='\t',quote=TRUE)
