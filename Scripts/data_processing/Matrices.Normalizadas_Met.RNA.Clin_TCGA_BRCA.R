library(TCGAbiolinks)
library(SummarizedExperiment)
library(gprofiler2)
library(edgeR)
library(tidyverse)
library(data.table)
library(sesameData)
library(sesame)
library(ELMER)
library(tidyverse)
library(tibble)
library(GEOquery)
library(maftools)



#DESCARGAR INFORMACION DE DATOS DE METILACION

mthyltn <-  GDCquery(project = "TCGA-BRCA",
                     data.category = "DNA Methylation",
                     platform="Illumina Human Methylation 450")


mthyltn=getResults(mthyltn)


i=substr(mthyltn$cases,1,19)

#DESCARGAR INFOMRACION DE DATOS DE EXPRESION 

xprssn <- GDCquery(project = "TCGA-BRCA",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = "STAR - Counts")

xprssn=getResults(xprssn)

j=substr(xprssn$cases,1,19)


############## EXTRAER LAS MUESTRAS QUE ESTAN PAREADAS ENTRE METILACION Y EXPRESION ########################
sapply(list(i,j),function(x) length(unique(x)))
#[1]  893 1221 

samples=intersect(i,j)
length(samples)
#[1] 870



samples=data.frame(cbind(samples,
                         sapply(samples,function(x) 
                           unique(as.character(mthyltn$sample_type[i==x])))))            

#Agregar una columna con el tipo de tejido

colnames(samples)[2]="tissue"
samples$patient=substr(samples$samples,1,12)




#Eliminar duplicados 

temp=table(samples[samples$patient%in%samples$patient[duplicated(
  samples$patient)],2:3])
table(apply(temp,2,paste,collapse=""))

#Sacar las muestras metastasicas 

samples<-samples[!(samples$tissue == "Metastatic"),] 

id<-rownames(samples)
samples<-cbind(samples,id)
row.names(samples)<- NULL



###################  INFORMACION CLINICA ##########################

#Descargar la informacion clinica desde TCGA

clin <- GDCquery_clinic("TCGA-BRCA","clinical")

#Extraer las columnas de interes

clin <- clin[,c("submitter_id","ethnicity","age_at_index",
                "ajcc_pathologic_stage","race")]


samples=cbind(samples,t(sapply(samples$patient,function(x) 
  clin[clin$submitter_id==x,2:5])))

#Ordenar las tablas de muestras e informacion clinica 

samples<- samples[order(samples$tissue),]

#Guardar la unformacion de las muetras que se seleccionaron
samples <-data.frame(lapply(samples, as.character), stringsAsFactors=FALSE)





#write.table(samples,"clinical_TCGA_Breast.txt",sep='\t',quote=F)

############# GENERAR MATRIZ DE EXPRESION ######################



# Descargar datos 

exp <- GDCquery(project = "TCGA-BRCA",
                data.category = "Transcriptome Profiling", 
                data.type = "Gene Expression Quantification", 
                workflow.type = "STAR - Counts",
                barcode = samples$id)



GDCdownload(exp)

raw.counts <- GDCprepare(query = exp)


#Extraer la matriz de expresion del objeto Summarized Experiment 

counts <- raw.counts@assays@data@listData$unstranded


# Convertir anotacion de  Ensemble ID a Gene Symbol
ensg<- unlist(strsplit(rownames(raw.counts), split = "[.]"))
ensg<- ensg[c(TRUE, FALSE)]

gconvert <- gconvert(query=ensg, target="HGNC", mthreshold=1, filter_na=FALSE)
class(gconvert[,1])<-"integer"
gconvert<- gconvert[sort.list(gconvert[,1]),]
genesymbol <- gconvert[,5]
genesymbol<-ifelse(is.na(genesymbol), ensg, genesymbol)

rownames(counts) <- genesymbol
colnames(counts) <- colnames(raw.counts)


############### FILTRADO Y NORMALIZACION #################


d<-DGEList(counts)

keep <- filterByExpr(d, min.count = 10, min.total.count = 10, large.n = 10, min.prop = 0.7)

d <- d[keep, , keep.lib.sizes=FALSE]

d$samples$lib.size <- colSums(d$counts)

d <- calcNormFactors(d, method = "TMM")

RNA.norm_cpm<- cpm(d, normalized.sizes = TRUE)



######################## METILACION ##########################

#OBTENER LOS DATOS 
mthyltn_el <-  GDCquery(project = "TCGA-BRCA",
                        data.category = "DNA Methylation",
                        platform="Illumina Human Methylation 450",
                        data.type = "Methylation Beta Value",
                        barcode=samples$id)

GDCdownload(mthyltn_el) 

mthyltn_el=GDCprepare(query= mthyltn_el) 

########## REMOVER NAs ######################################

TCGAprepare_elmer <- function(data,
                              platform,
                              met.na.cut = 0.2,
                              save = FALSE){
  # parameters veryfication
  
  if (missing(data))  stop("Please set the data parameter")
  if (missing(platform))  stop("Please set the platform parameter")
  
  if (grepl("illuminahiseq_rnaseqv2|illuminahiseq_totalrnaseqv2",
            platform, ignore.case = TRUE)) {
    message("============ Pre-pocessing expression data =============")
    message(paste0("1 - expression = log2(expression + 1): ",
                   "To linearize \n    relation between ",
                   "methylation and expression"))
    if(typeof(data) == typeof(SummarizedExperiment())){
      row.names(data) <- paste0("ID",values(data)$entrezgene)
      data <- assay(data)
    }
    
    if(all(grepl("\\|",rownames(data)))){
      message("2 - rownames  (gene|loci) => ('ID'loci) ")
      aux <- strsplit(rownames(data),"\\|")
      GeneID <- unlist(lapply(aux,function(x) x[2]))
      row.names(data) <- paste0("ID",GeneID)
    }
    data <- log2(data+1)
    Exp <- data.matrix(data)
    
    if (save)  save(Exp,file = "Exp_elmer.rda")
    return(Exp)
  }
  
  if (grepl("humanmethylation", platform, ignore.case = TRUE)) {
    message("============ Pre-pocessing methylation data =============")
    if (class(data) == class(data.frame())){
      msg <- paste0("1 - Removing Columns: \n  * Gene_Symbol  \n",
                    "  * Chromosome  \n  * Genomic_Coordinate")
      message(msg)
      data <- subset(data,select = 4:ncol(data))
    }
    if(typeof(data) == typeof(SummarizedExperiment())){
      data <- assay(data)
    }
    msg <- paste0("2 - Removing probes with ",
                  "NA values in more than 20% samples")
    message(msg)
    data <- data[rowMeans(is.na(data)) < met.na.cut,]
    Met <- data.matrix(data)
    if (save)  save(Met,file = "Met_elmer.rda")
    return (Met)
  }
}


met <- TCGAprepare_elmer(mthyltn_el,
                         platform = "HumanMethylation450",
                         save = TRUE,                         met.na.cut = 0.20)
#met2<-met
#met<-met2
#colnames(met)=substr(colnames(met),1,19)


######################### SEPARAR lncRNA"####################################

#Cargar anotacion

aME<-as.data.frame(raw.counts@rowRanges@elementMetadata@listData)

a<-dplyr::filter(aME, gene_type == "lncRNA")


ann=substr(aME$gene_id,1,15)

lncRNA<-RNA.norm_cpm[row.names(RNA.norm_cpm)%in%ann,]


############## SEPARAR GENES   #####################


b<-dplyr::filter(aME, gene_type == "protein_coding")

RNA<-RNA.norm_cpm[row.names(RNA.norm_cpm)%in%b$gene_name,]


####### SEPARAR NORMALES Y TUMORES DE TODAS LAS MATRICES ####################

###### Matriz de RNA
RNA<-as.data.frame(RNA)

RNA.t<-RNA%>% select ( grep ("01A|01B", colnames (RNA)))
RNA.n<-RNA%>% select ( grep ("11A|11B", colnames (RNA)))

RNA.t<-RNA.t%>% select (order (colnames(RNA.t)))
RNA.n<-RNA.n%>% select (order (colnames(RNA.n)))


#############################################################################

lncRNA<-as.data.frame(lncRNA)


lncRNA.t<-lncRNA%>% select ( grep ("01A|01B", colnames (lncRNA)))
lncRNA.n<-lncRNA%>% select ( grep ("11A|11B", colnames (lncRNA)))

lncRNA.t<-lncRNA.t%>% select (order (colnames(lncRNA.t)))
lncRNA.n<-lncRNA.n%>% select (order (colnames(lncRNA.n)))

##########################################################################

met.df<-as.data.frame(met)


met.t<-met.df%>% select ( grep ('01A|01B', colnames (met.df)))
met.n<-met.df%>% select ( grep ('11A|11B', colnames (met.df)))

met.n<-met.n%>% select (order (colnames(met.n)))
met.t<-met.t%>% select (order (colnames(met.t)))

met.t<- tibble::rownames_to_column(met.t, "CpGName") 
met.n<-tibble::rownames_to_column(met.n, "CpGName")



######## ACORTAR EL IDENTIFICADOR 

colnames(RNA.n)=substr(colnames(RNA.n),1,12)
colnames(RNA.t)=substr(colnames(RNA.t),1,12)

colnames(met.n)=substr(colnames(met.n),1,12)
colnames(met.t)=substr(colnames(met.t),1,12)

colnames(lncRNA.n)=substr(colnames(lncRNA.n),1,12)
colnames(lncRNA.t)=substr(colnames(lncRNA.t),1,12)

######################################################################

#REVISAR DUPLICADOS

RNA.n<-RNA.n[,!duplicated(colnames(RNA.n))]
RNA.t<-RNA.t[,!duplicated(colnames(RNA.t)) ]

lncRNA.n<-lncRNA.n[,!duplicated(colnames(lncRNA.n))]
lncRNA.t<-lncRNA.t[,!duplicated(colnames(lncRNA.t)) ]

met.n<-met.n[,!duplicated(colnames(met.n))]
met.t<-met.t[,!duplicated(colnames(met.t)) ]

#Separar clinica

clin.tumor<-filter(samples, tissue == "Primary Tumor")
clin.normal<-filter(samples, tissue == "Solid Tissue Normal")

clin.tumor<-clin.tumor[!duplicated(clin.tumor$patient),]
clin.normal<-clin.normal[!duplicated(clin.normal$patient),]

clin.normal <- clin.normal[with(clin.normal, order(clin.normal$patient)), ] # Orden directo 
clin.tumor <- clin.tumor[with(clin.tumor, order(clin.tumor$patient)), ]

#Verificar orden muestras de la tabla clinica vs df metiliacion
clin.tumor$patient %in% colnames(met.t)
clin.normal$patient %in% colnames(met.n)

######## GUARDAR MATRICES  
#Expresion
write.table(RNA.n,"RNA.normal_cpm_TCGA_breast.tsv",sep='\t',quote=F,row.names=TRUE)
write.table(RNA.t,"RNA.tumor_cpm_TCGA_breast.tsv",sep='\t',quote=F,row.names=TRUE)

#lncRNA

write.table(lncRNA.n,"lncRNA.normal_cpm_TCGA_breast.tsv",sep='\t',quote=F,row.names=TRUE)
write.table(lncRNA.t,"lncRNA.tumor_cpm_TCGA_breast.tsv",sep='\t',quote=F,row.names=TRUE)
write.table(a,"anotacion_lncRNA_TCGA_breast.tsv",sep='\t',quote=F,row.names=TRUE)

#Metilacion

write.table(met.t,"met.tumor_filter.by_gene_TCGA_BRCA.tsv",sep='\t',quote=F)
write.table(met.n,"met.normal_filter_by_gene_TCGA_BRCA.tsv",sep='\t',quote=F)

#Clinica

write.table(clin.normal,"clin.normal_TCGA_BRCA.tsv",sep='\t',quote=F)
write.table(clin.tumor,"clin.tumor_TCGA_BRCA.tsv",sep='\t',quote=F)


