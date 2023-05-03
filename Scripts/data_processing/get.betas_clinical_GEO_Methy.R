
library(sesame)
library(sesameData)
library(parallel)
library(GEOquery)



#Descargar datos 

download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE74214&format=file", "GSE74214_RAW.tar")

#getOption('timeout')
#options(timeout=1500)


#Revisar que la lista de archivos del .tar verificar que haya .idat 
untar("GSE74214_RAW.tar", list = TRUE)


#Descomprimir 
untar("GSE74214_RAW.tar", exdir = "/Users/diana/Desktop/Proyecto Relojes/Sandra/Matrices_Metilacion_Expresion_TCGA.Breast/GEO")





######################################       OBTENER BETAS  ################################
# getBetas= La función toma un SigSet y devuelve un valor beta
#dyeBiasNL = compara los canales verdes y rojos y genera un mapeo para corregir la señal hacia el medio
#noob = normal-exponential out-of-band correccion de fondo con normalizacion para dye-bias
#pOOBAH = Toma un SigDF Y detecta un p-value utilizando la distribución empírica de las sondas fuera de banda 
#         y devuelve un nuevo SigDF con un nuevo slot  enmascardo actualizado.

betas=do.call(cbind,mclapply(searchIDATprefixes("."), function(px) getBetas(dyeBiasNL(noob(pOOBAH(readIDATpair(px))))), mc.cores=2))

#dim(betas)
#class(betas)


write.table(betas,"betas_GEO_GSE74214.tsv",sep='\t',quote=F)

#######################################################################################################

#s= mclapply(searchIDATprefixes("."), readIDATpair, mc.cores=2)
#s2=pOOBAH(s$GSM3931707_200360140015)
#pval = pOOBAH(s2, return.pval=TRUE)
#s3=addMask(s2, pval>0.05)
#s4 = noob(s3)
#s5 = dyeBiasNL(s4)


#sesameQC_plotIntensVsBetas(s$GSM3931707_200360140015) #s -nada-
#sesameQC_plotIntensVsBetas(s2) #+ pOOBAH
#sesameQC_plotIntensVsBetas(s3) # + addMask
#sesameQC_plotIntensVsBetas(s4) # + noob
#sesameQC_plotIntensVsBetas(s5) # + dyeBiasNL



########################## PHENODATA ###########

dd <- getGEO("GSE74214")
geo<- dd[[1]]
#geo

pheno <- pData(geo)


age <- as.numeric(pheno$`subject age:ch1`)
############################################################################
#En el caso de que en un  pheno data tengamos una combinacion de letras y numeros para la edad
pheno <- pData(gse133985)
age<-(pheno$`age:ch1`)


#age<- round(as.numeric(gsub(".*?([0-9]*.[0-9]).*", "\\1", age)),0) 

" .culaquier caracter  *? toma solo los numeros * puede haber . con numeros round para redondear el cero significa que de entero 

############################################################################

id<-pheno$geo_accession

#Hacer data frame datos clinicos 

clin<- data.frame(id,age)

write.table(clin,"clin_GEO_GSE74214.tsv",sep='\t',quote=F)

