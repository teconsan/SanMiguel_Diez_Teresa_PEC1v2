#Script con todo el código utilizado y comentado
#Inicialmente trabajé la actividad en un Rmarkdown que me ha permitido ordenar un poco la secuencia
#PAQUETES INSTALADOS
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(version = "3.20", ask = FALSE) #para Bioconductor

if (!require("metabolomicsWorkbenchR")) {
  BiocManager::install("metabolomicsWorkbenchR") #para metabolómica esencial
}
if (!require("MetaboAnalystR")) {
  BiocManager::install("MetaboAnalystR") #para metabolómica avanzada
}

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse") #para inspección rápida y otras cosas
}
if (!require("VIM")) {
  install.packages("VIM") #para la imputación k-NN
}
if (!requireNamespace("factoextra", quietly = TRUE)) {
  install.packages("factoextra") #para el PCA
}
if (!require("plotly")) {
  install.packages("plotly") #para gráficos interactivos
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

#LIBRERÍAS UTILIZADAS
library(Biobase)
library(readxl) #para abrir xlsx
library(tidyverse) 
library(SummarizedExperiment) #para crear el objeto de tipo SE ideal para omicas
library(VIM) #para la imputación k-NN
library(factoextra) #PCA
library(plotly) #gráficos interactivos
library(pheatmap)#clustering
library(dplyr)
library(ggplot2) #para graficos
library(car)   # Para el test de Levene
library(broom) # Para tidy() si se quiere procesar resultados de test
library(tibble)  # Para trabajar más cómodo con tablas
library(tidyr)   # Para pivotear si hiciera falta

#EXPLORACIÓN DEL FICHERO SELECCIONADO

list.files() #compruebo el fichero del repositorio clonado

files <- list.files(recursive = TRUE, pattern = "GastricCancer.*\\.xlsx$", ignore.case = TRUE)
files <- files[!grepl("~\\$", files)]  # Filtrar archivos temporales de Excel (porque como lo estoy mirando,aparece)
print(files)


excel_sheets("D:/Documentos D/OMICAS/SanMiguel_Diez_Teresa_PEC1v2/GastricCancer_NMR.xlsx") #Inspercciono el excel de interés a ver cuantas hojas tiene

datos <- read_excel("D:/Documentos D/OMICAS/SanMiguel_Diez_Teresa_PEC1v2/GastricCancer_NMR.xlsx", sheet = "Data")
picos <- read_excel("D:/Documentos D/OMICAS/SanMiguel_Diez_Teresa_PEC1v2//GastricCancer_NMR.xlsx", sheet = "Peak")
head(datos) #Tras comprobar que tiene 2, extraigo la información de ambos
head(picos)

glimpse(datos)#inspección rápid para ver NA, si las variables son factor, chr etc
glimpse(picos) #lo mismo con la otra pestaña del excel

matriz_datos <- as.matrix(datos %>% select(starts_with("M"))) #Preparación de la matriz de expresión de datos omicos
matriz_datos <- t(matriz_datos)  # la transpongo 149 metabolitos × 140 muestras

meta_muestras <- datos %>% #Preparación de la información sobre las muestras incluidas
  select(Idx, SampleID, SampleType, Class) %>% 
  column_to_rownames("SampleID")  # Entre IDx y SampleID, elijo SampleID como rownames

colnames(matriz_datos) <- rownames(meta_muestras) #Identificación de muestras en la cabecera

rowData <- picos %>% #*Asignación de la información de los metabolitos como rowData (identificación e indicios de calidad)
  column_to_rownames("Name") %>%  # "Name" es M1, M2...
  .[rownames(matriz_datos), ]     # asegurar orden correcto

se <- SummarizedExperiment( #Creación del objeto SummarizedExperiment
  assays = list(counts = matriz_datos),
  colData = meta_muestras,
  rowData = rowData
)

metadata(se)$fuente <- "https://github.com/nutrimetabolomics/metaboData" #Referencio la fuente
metadata(se)$fecha <- Sys.Date() #Incluyo la fecha

save(se, file = "SE_GastCancer.Rda")#guardo el SummarizedExperiment

#Filtrado por calidad
# Extraer los valores de QC_RSD y Perc_missing desde rowData
rsd <- rowData(se)$QC_RSD
perc_missing <- rowData(se)$Perc_missing

# Crear un vector lógico con los metabolitos que <10% de datos faltantes y <20% de variabilidad entre réplicas
filtro_buena_calidad <- (rsd < 20) & (perc_missing < 10)

# Filtrar el objeto para conservar solo esos metabolitos
se_filtrado <- se[filtro_buena_calidad, ]

# Ver cuántos quedaron
dim(se_filtrado)

#guardo también el SE ya filtrado por calidad
save(se, file = "SE_GastCancer_filtrado.Rda")

#Guardado de datos en formato txt
#Guardado de elementos 1:Guardar la matriz de datos
write.table(assay(se_filtrado), file = "datos_matriz.txt", sep = "\t", quote = FALSE, row.names = TRUE)

#Guardado de elementos 2: Guardar los metadatos de las muestras (colData)
write.table(as.data.frame(colData(se_filtrado)), file = "metadatos_muestras.txt", sep = "\t", quote = FALSE, row.names = TRUE)

#Guardado de elementos 3: Guardar metadatos de variables(rowData)
write.table(as.data.frame(rowData(se_filtrado)), file = "metadatos_variables.txt", sep = "\t", quote = FALSE, row.names = TRUE)

