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

if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")#Machinelearning



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
library(caret)  #Para machinelearning-partición estratificada
library(mixOmics) #Para crear el modelo de machinelearning

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

#ANÁLISIS DE LOS DATOS
#1. ANÁLISIS PCA
#1.1 Transformación log10 para reducir efecto de valores extremos
#Filtramos/extraemos las muestras de la matriz del se_filtado, que están en assay
datos_raw <- assay(se_filtrado)  # 52 metabs × 140 muestras
#Transponemos muestras como filas
datos_raw <- t(datos_raw)  # 140 filas (muestras), 52 columnas (metabolitos)
#Hacemos el log10
# Log-transformación (añadiendo +1 para evitar log(0))
datos_log <- log10(datos_raw + 1)

#1.2 Imputación k-NN mejor que imputación medianas k=3
datos_knn <- kNN(as.data.frame(datos_log), k = 3, imp_var = FALSE)
datos_knn <- as.matrix(datos_knn)  # Para asegurarnos de que sea matriz

#1.3 Análisis PCA
res.pca <- prcomp(datos_knn, scale. = TRUE)

#1.4 Visualización PCA
library(factoextra) 
fviz_pca_ind(res.pca,
             geom.ind = "point",
             habillage = colData(se_filtrado)$SampleType,
             addEllipses = TRUE,
             title = "PCA con log10 + imputación k-NN (QC vs Sample)")
#1.5 Tabla para ver la carga de cada metabolito
# Paso 1: Extraer las cargas del PCA
loadings <- res.pca$rotation
df_loadings <- as.data.frame(loadings)
df_loadings$M_ID <- rownames(df_loadings)  # Aquí están los M1, M2, ...

# Paso 2: Calcular contribuciones
df_loadings$contrib_PC1 <- df_loadings$PC1^2
df_loadings$contrib_PC2 <- df_loadings$PC2^2
df_loadings$contrib_total <- df_loadings$contrib_PC1 + df_loadings$contrib_PC2

# Paso 3: Crear tabla con nombres de metabolitos
df_nombres <- as.data.frame(rowData(se_filtrado)) %>%
  mutate(M_ID = rownames(.)) %>%
  select(M_ID, Nombre = Label)

# Paso 4: Unir ambas tablas por M_ID
df_completo <- left_join(df_loadings, df_nombres, by = "M_ID")

# Paso 5: Ordenar por contribución total
df_ordenado <- df_completo %>%
  arrange(desc(contrib_total)) %>%
  select(Nombre, M_ID, PC1, PC2, contrib_PC1, contrib_PC2, contrib_total)

# Mostrar los 10 principales
head(df_ordenado, 10)

#1.6 Gráficos
#1.6.1 Gráfico estático de cargas
df_cargas <- df_ordenado

ggplot(df_cargas, aes(x = PC1, y = PC2, label = Nombre)) +
  geom_point(color = "green", size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5, color = "black", size = 3) +
  labs(title = "PCA Loadings",
       x = "PC1",
       y = "PC2") +
  theme_minimal()

#1.6.2 Gráfico "interactivo" para verlo online
plot_ly(data = df_cargas,
        x = ~PC1,
        y = ~PC2,
        type = 'scatter',
        mode = 'markers',
        text = ~Nombre,
        hoverinfo = 'text',
        marker = list(size = 8, color = 'green')) %>%
  layout(title = "PCA Loadings (interactivo)",
         xaxis = list(title = "PC1"),
         yaxis = list(title = "PC2"))


#2. ANÁLISIS ESTADÍSTICO
#2.1 Filtrado de datos de pacientes con cáncer y sanos
# Extraer datos del SE filtrado
datos <- assay(se_filtrado) #ojo que con esto sobreescribo "datos" quer era como había llamado inicialmente las cosas
metainfo <- colData(se_filtrado)

# Filtrar solo GC y HE (excluye QC y BN)
idx_gc_he <- metainfo$Class %in% c("GC", "HE")
datos_gc_he <- datos[, idx_gc_he]
metainfo_gc_he <- metainfo[idx_gc_he, ]

#2.2 Test de Wilcoxon
# Función para comparar cada metabolito
comparar_metabolito <- function(x, grupo) {
  grupo1 <- x[grupo == "GC"]
  grupo2 <- x[grupo == "HE"]
  
  # Test de Wilcoxon
  test <- wilcox.test(grupo1, grupo2)
  
  # Fold change (GC / HE)
  fc <- median(grupo1, na.rm = TRUE) / median(grupo2, na.rm = TRUE)
  
  data.frame(p_value = test$p.value, 
             fold_change = fc)
}
# Aplicar a cada metabolito
grupo <- metainfo_gc_he$Class
resultados <- apply(datos_gc_he, 1, comparar_metabolito, grupo = grupo)
# Convertir a tabla
tabla_resultados <- do.call(rbind, resultados)
tabla_resultados <- as.data.frame(tabla_resultados)
tabla_resultados$metabolito <- rownames(tabla_resultados)
# Ajustar p-valores por FDR
tabla_resultados$adj_p <- p.adjust(tabla_resultados$p_value, method = "fdr")
# Ordenar por significancia
tabla_ordenada <- tabla_resultados %>% arrange(p_value)
head(tabla_ordenada)
# Añadir el nombre del metabolito para que sea más legible
# Convertir rowData a data.frame que llamo metanombres
metanombres <- as.data.frame(rowData(se_filtrado))
# Añadir identificador como fila para hacer merge
metanombres$metabolito <- rownames(metanombres)
# Juntar con la tabla ordenada
tabla_con_nombres <- tabla_ordenada %>%
  left_join(metanombres %>% select(metabolito, Label), by = "metabolito") %>%
  relocate(Label, .after = metabolito)  # mover el nombre al lado del código
head(tabla_con_nombres)

# 2.3 VOLCANO PLOT
# Añadir columnas a la tabla
tabla_con_nombres <- tabla_con_nombres %>%
  mutate(log2FC = log2(fold_change),
         negLogP = -log10(p_value),
         Significativo = case_when(
           adj_p < 0.05 & log2FC > 1 ~ "GC ↑",
           adj_p < 0.05 & log2FC < -1 ~ "HE ↑",
           TRUE ~ "No significativo"
         ))
# señalo en rojo los aumentados en cáncer y en azul los aumentados en sanos 
# (en realidad son los disminuidos en cáncer, que es lo mismo que decir los aumentados 
# en sanos)
ggplot(tabla_con_nombres, aes(x = log2FC, y = negLogP, color = Significativo)) +
  geom_point() +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("GC ↑" = "red", "HE ↑" = "blue", "No significativo" = "grey")) +
  labs(title = "Volcano plot: GC vs HE",
       x = "log2(Fold Change)",
       y = "-log10(p-valor)") +
  theme_minimal()
# Guardo la tabla
write.csv(tabla_con_nombres, "Resultados_GC_vs_HE.csv", row.names = FALSE)

#Confirmación de que el test no-paramétrico era la opción adecuada
#Compruebo normalidad y homogeneidad de varianzas
# Crear función que evalúe por metabolito
evaluar_metabolito <- function(x, grupo) {
  # Separar los valores por grupo
  valores_gc <- x[grupo == "GC"]
  valores_he <- x[grupo == "HE"]
  
  # Calcular medias
  media_gc <- mean(valores_gc, na.rm = TRUE)
  media_he <- mean(valores_he, na.rm = TRUE)
  
  # IC 95% (asumiendo t-distribución)
  n_gc <- sum(!is.na(valores_gc))
  n_he <- sum(!is.na(valores_he))
  se_gc <- sd(valores_gc, na.rm = TRUE) / sqrt(n_gc)
  se_he <- sd(valores_he, na.rm = TRUE) / sqrt(n_he)
  ic95_gc <- c(media_gc - 1.96 * se_gc, media_gc + 1.96 * se_gc)
  ic95_he <- c(media_he - 1.96 * se_he, media_he + 1.96 * se_he)
  
  # Shapiro para normalidad
  shapiro_gc <- if (length(valores_gc) >= 3) shapiro.test(valores_gc)$p.value else NA
  shapiro_he <- if (length(valores_he) >= 3) shapiro.test(valores_he)$p.value else NA
  
  # Levene test
  df_tmp <- data.frame(valores = c(valores_gc, valores_he),
                       grupo = rep(c("GC", "HE"), times = c(length(valores_gc), length(valores_he))))
  p_levene <- if (nrow(df_tmp) >= 4) leveneTest(valores ~ grupo, data = df_tmp)$`Pr(>F)`[1] else NA
  
  # Indicadores lógicos
  normalidad_gc <- !is.na(shapiro_gc) && shapiro_gc >= 0.05
  normalidad_he <- !is.na(shapiro_he) && shapiro_he >= 0.05
  homogeneidad <- !is.na(p_levene) && p_levene >= 0.05
  
  # Resultado
  data.frame(
    media_gc = media_gc,
    ic95_gc_inf = ic95_gc[1],
    ic95_gc_sup = ic95_gc[2],
    media_he = media_he,
    ic95_he_inf = ic95_he[1],
    ic95_he_sup = ic95_he[2],
    p_shapiro_gc = shapiro_gc,
    p_shapiro_he = shapiro_he,
    p_levene = p_levene,
    cumple_shapiro_gc = normalidad_gc,
    cumple_shapiro_he = normalidad_he,
    cumple_levene = homogeneidad
  )
}

# Aplicar a todos los metabolitos
resultados_supuestos <- apply(datos_gc_he, 1, evaluar_metabolito, grupo = grupo)

# Convertir lista en tabla
tabla_supuestos <- bind_rows(resultados_supuestos, .id = "metabolito")

# Añadir nombres de metabolitos
tabla_supuestos <- tabla_supuestos %>%
  left_join(metanombres %>% select(metabolito, Label), by = "metabolito") %>%
  relocate(Label, .after = metabolito)

# Ver los primeros resultados
head(tabla_supuestos)

# MACHINE LEARNING
# 1.Preparación
#Extraer los datos de la matriz de expresión y de la metainformación
datos <- assay(se_filtrado)
metainfo <- colData(se_filtrado)

# Filtrar GC y HE
idx_gc_he <- metainfo$Class %in% c("GC", "HE")
datos_gc_he <- datos[, idx_gc_he]
metainfo_gc_he <- metainfo[idx_gc_he, ]

# 1.1 Transponer y transformar los datos
#A seguramos que estén en formato muestras × metabolitos, y luego imputamos
# Transponer: filas = muestras, columnas = metabolitos
datos_gc_he_t <- t(datos_gc_he)
# Log10 (+1 para evitar log(0))
datos_log <- log10(datos_gc_he_t + 1)
# Imputar con kNN
datos_knn <- kNN(as.data.frame(datos_log), k = 3, imp_var = FALSE)
datos_knn <- as.matrix(datos_knn)

# 1.2 Dividir entre train y test 75/25
set.seed(123)  # reproducibilidad
# Dividir
particion <- createDataPartition(metainfo_gc_he$Class, p = 0.75, list = FALSE)
# Crear subconjuntos
datos_train <- datos_knn[particion, ]
datos_test  <- datos_knn[-particion, ]

metainfo_train <- metainfo_gc_he[particion, ]
metainfo_test  <- metainfo_gc_he[-particion, ]

# 1.3 Escalado Z-score. Esto se hace con medias y SD del train.  
# Escalar train
datostrain_escalado <- scale(datos_train)
# Guardar medias y SD
media_train <- attr(datostrain_escalado, "scaled:center")
sd_train <- attr(datostrain_escalado, "scaled:scale")
# Escalar test usando media y SD del train
datostest_escalado <- scale(datos_test, center = media_train, scale = sd_train)

# 2. Crear modelo PLS-DA y entrenarlo
# Me baso en el tutorial para elegir dos componentes. Todavía no conozco estas herramientas
# Convertir clases a factor
Y_train <- factor(metainfo_train$Class, levels = c("HE", "GC"))
# Entrenar modelo PLS-DA
modelo_plsda_final <- plsda(datostrain_escalado, Y_train, ncomp = 2)


# 3. Evaluar la validez del modelo creado
set.seed(123)
perf_plsda <- perf(modelo_plsda_final, validation = "Mfold", folds = 5, nrepeat = 10, progressBar = TRUE)
plot(perf_plsda)  # ver tasa de error por número de componentes

# Extraer métricas de error
# Error 1: global usando 2 componentes y distancia por centroides
perf_plsda$error.rate$overall["comp2", "centroids.dist"]

# Error 2: Error por clase usando 2 componentes y distancia por centroides
perf_plsda$error.rate.class$centroids.dist[, "comp2"]

# 3. Número óptimo de componentes según error global (distancia centroides)
perf_plsda$choice.ncomp["overall", "centroids.dist"]