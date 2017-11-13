rm(list=ls())
# Display the current working directory
getwd()
# If necessary, change the working directory to the one containing the data
workingDir<-"M:/Code for Projects/Project 1"
setwd(workingDir)
options(stringsAsFactors=FALSE)
### Source file for functions used
source("Functions.R")
source("https://bioconductor.org/biocLite.R")
### Set the CRAN mirror for installing packages
local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cloud.r-project.org/"
  options(repos = r)
})
### Packages that will be used for this analysis
packages<-c("igraph", "Matrix", "WGCNA", "flashClust", "gdata", "MASS", "corrplot", "huge")
### Function "ipak" is defined in the end of the file. Check to see if multiple packages are installed.
### If they are not, it installs and loads them. It takes a while...
# biocLite(c("GO.db", "preprocessCore", "impute"),suppressUpdates=TRUE) 
ipak(packages)
# Load the DILGOM 2007 metabolite data and give rownames
metabolomics<-as.data.frame(read.table("u7671950_export.v6uUbn.txt",sep=",",header=T))
rownames(metabolomics)=metabolomics$ID
# Load the DILGOM 2007 phenotype data and select Sex, Age, BMI, 
# Cholesterol, Fasting Glucose, Mean Diastolic Blood Pressure, Blood Pressure 
# Medication, Diabetes Diagnosis, Cholesterol Medication. Give rownames
phenotypes<-as.data.frame(read.csv("exls.csv",header=T,sep=","))
rownames(phenotypes)=phenotypes$ID
phenotypes<-phenotypes[,c(2,3,5,9,12,14,17,18,21)]
# Select only one measurement per particle size in metabolites.
metabolomics<-metabolomics[,c(2,4,8,15,22,28,32,36,42,47,50,57,63,69,72,77:84,88:90,92:96,99,100,102:115,118,119,121,122,124:130) ]
# Screen samples 
# NAs in FR07_38 are coded as not diagnosed with diabetes
phenotypes$FR07_38[is.na(phenotypes$FR07_38)] <- 0
# NAs in Q29 are coded as not receiving medication for blood pressure
phenotypes$Q29[is.na(phenotypes$Q29)] <- 1
# Merge Metabolite and Phenotype data 
data <- merge(metabolomics,phenotypes, by=0)
rm(phenotypes,metabolomics)
#  Remove individuals diagnosed with diabetes (FR07_38>1)
data<-data[ which(data$FR07_38==1|data$FR07_38==0),]
# Remove individuals with outlying fasting glucose (SL_GLUK_0H>10)
data<-data[ which(data$SL_GLUK_0H<=10), ]
# Remove individuals with cholesterol medication (K34=1)
data<-data[ which(data$K34==0), ]
# Remove samples with any NAs
data<-na.omit(data)
# Give rownames
rownames(data)<-data$Row.names
# Remove variables used in screening (Cholesterol, Glucose, diabetes diagnosis, cholesterol medication)
data<-data[,-c(1,63,64,67,68)]
# Rename Phenotypes
colnames(data)[59]<-"SEX"
colnames(data)[60]<-"AGE"
colnames(data)[61]<-"BMI "
colnames(data)[62]<-"DIASTOLIC"
colnames(data)[63]<-"BLOOD_PRESSURE_MEDICATION"
data$SEX<-factor(data$SEX,levels=c(1,2),labels = c("male","female"))
data$BLOOD_PRESSURE_MEDICATION<-factor(data$BLOOD_PRESSURE_MEDICATION,levels=c(1,2),labels = c("no","yes"))
# Adjust metabolites for Diastolic blood pressure and Blood pressure medication
for (i in 1:58) {
  Z<-lm((data[,i])~data$DIASTOLIC+data$BLOOD_PRESSURE_MEDICATION)
  data[,i]<-Z$residuals
}
# Transform BMI to categorical variables (3 levels)
data$BMI_class<-cut(data$BMI, c(min(data$BMI),quantile(data$BMI,1/3),quantile(data$BMI,2/3),max(data$BMI)),labels=c("low","medium","high"),include.lowest = TRUE)
data <- within(data, BMI_class <- relevel(BMI_class, ref = 2))
# Data containing the metabolite part related to BMI
data.bmi.info<-c()
for (i in 1:58) {
  Z<-lm((data[,i])~BMI_class+AGE+SEX+AGE:SEX+BMI_class:AGE+BMI_class:SEX+BMI_class:SEX:AGE,data=data)
  # Keeping the predicted values from the BMI related part
  val<-rowSums((predict(Z,type="terms"))[,c(1,5,6,7)])
  data.bmi.info<-cbind(data.bmi.info,val)
}
colnames(data.bmi.info)<-colnames(data)[1:58]
rownames(data.bmi.info)<-rownames(data)
# Data containing all information
data.all.info<-data[,1:58]

#########################################################################
################ ANALYSIS FOR DATA CONTAINING ALL INFO ##################
#########################################################################
#########################################################################
#########################################################################
################################ WGCNA ##################################
#########################################################################
#########################################################################
### WGCNA network estimation
heatmap.data<- data.all.info
### Tuning parameters for correlation plot
softPower = soft.thres(heatmap.data);
### What is the minimum allowed number of metabolites in a module 
minModuleSize = 8;
### Define the similarity and dissimilarity matrix 
sim.mat<-(abs(cor(heatmap.data, method="pearson",use="pairwise.complete.obs")))^softPower;diag(sim.mat)<-0
diss.mat = 1-sim.mat
### Metabolites dendrogram
metabolTree = flashClust(as.dist(diss.mat ), method = "average");
### Identification of modules
dynamicMods = cutreeDynamic(dendro = metabolTree, distM = diss.mat,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)
### Visualization of correlation plot
TOMplot(diss.mat, metabolTree, dynamicColors, main = "Network correlation plot")
### How many metabolites does every cluster contain
table(dynamicColors)
### Compute PC1 of each of the clusters
MEList = moduleEigengenes(heatmap.data, colors = dynamicColors)
### The values of the PC1
MEs = MEList$eigengenes
### For keeping the names of the samples
rownames(MEs)<-rownames(heatmap.data)
### Dissimilarity of the PC1s
MEDiss = 1-cor(MEs,use="pairwise.complete.obs");
### Dendrogram of clusters
METree = hclust(as.dist(MEDiss), method = "average");
sizeGrWindow(7, 6)
### Plotting of the dendrogram
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.2
abline(h=MEDissThres, col = "red")
### Merge 2 clusters if their correlation is over 0.8
merge = mergeCloseModules(heatmap.data, dynamicColors, cutHeight = MEDissThres, verbose = 2)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
dynamicColors<-mergedColors
### Visualize the final correlation matrix
TOMplot(diss.mat, metabolTree, dynamicColors, main = "Network correlation plot")

### Save the plot as image
png(file = "Heatmap.png")
TOMplot(diss.mat, metabolTree, dynamicColors, main = "Network correlation plot")
dev.off()

### For network visualization keep the top 10% of the edges
q<-quantile(upperTriangle(sim.mat, diag=FALSE),0.9);q
### Compute adjacency matrix
adj<-topCpercent(x=sim.mat,c=0.1)
### Conformity based quantities
length(conformityBasedNetworkConcepts(sim.mat)$conformityBasedNCs$Conformity)
conformityBasedNetworkConcepts(sim.mat)$fundamentalNCs$Density
conformityBasedNetworkConcepts(sim.mat)$fundamentalNCs$Centralization
conformityBasedNetworkConcepts(sim.mat)$fundamentalNCs$Heterogeneity
color="black"
length(conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$conformityBasedNCs$Conformity)
conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$fundamentalNCs$Density
conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$fundamentalNCs$Centralization
conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$fundamentalNCs$Heterogeneity
sort(colnames(sim.mat)[which(dynamicColors==color)])

### Alternatively, compute Density, Centralization, and Heterogeneity from here
quantities<-DensCentrHet1(sim.mat,dynamicColors);quantities

### Top 20 connected metabolites
mostConnected(topCmatrix,nr=20);

### Metabolites within each module
for (color in rownames(quantities)[-dim(quantities)[1]]){
  one<-sort(colnames(heatmap.data)[which(dynamicColors==color)])
  print(color,quote=FALSE)
  cat(one, sep =", ")
  print("\n")
}


### For saving the edges to vizualize them in Gephi
write.edges(adj)

#########################################################################
#########################################################################
#########################################################################
############################### GLASSO ##################################
#########################################################################
#########################################################################

data<- data.all.info
### You have to omit the NA data here
data<-as.matrix(na.omit(data))

### Grid for selecting the regularization parameter
sequ<-seq(from=1,to=0.001, by=-0.001)
### Fitting GL for every lambda from the grid
GL<-huge(data,method="glasso",lambda=sequ)
### Selecting the optimal lambda
GL.net<-huge.select(est=GL,criterion="stars",rep.num=100,stars.thresh=0.02)
### Visual inspection of the network
huge.plot(GL.net$refit)
### Optimal regularization parameter
GL.net$opt.lambda
### Intensity matrix of network
GL.int<-as.matrix((GL.net$refit)*GL.net$merge[[GL.net$opt.index]]);diag(GL.int)<-0;colnames(GL.int)<-rownames(GL.int)<-colnames(data)

### What is the minimum allowed number of metabolites in a module 
minModuleSize = 5;
### Define the similarity and dissimilarity matrix 
sim.mat<-GL.int
diss.mat = 1-sim.mat
### Metabolites dendrogram
metabolTree = flashClust(as.dist(diss.mat ), method = "average");
### Identification of modules
dynamicMods = cutreeDynamic(dendro = metabolTree, distM = diss.mat,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
### Visualization of correlation plot
TOMplot(diss.mat, metabolTree, dynamicColors, main = "Network correlation plot")
### How many metabolites does every cluster contain
table(dynamicColors)
### Compute PC1 of each of the clusters
MEList = moduleEigengenes(data, colors = dynamicColors)
### The values of the PC1
MEs = MEList$eigengenes
### For keeping the names of the samples
rownames(MEs)<-rownames(data)
### Dissimilarity of the PC1s
MEDiss = 1-cor(MEs,use="pairwise.complete.obs");
### Dendrogram of clusters
METree = hclust(as.dist(MEDiss), method = "average");
sizeGrWindow(7, 6)
### Plotting of the dendrogram
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.1
abline(h=MEDissThres, col = "red")
### Merge 2 clusters if their correlation is over 0.8
merge = mergeCloseModules(data, dynamicColors, cutHeight = MEDissThres, verbose = 2)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
dynamicColors<-mergedColors
### Visualize the final correlation matrix
TOMplot(diss.mat, metabolTree, dynamicColors, main = "Network correlation plot")

### Save the plot as image
png(file = "Heatmap.png")
TOMplot(diss.mat, metabolTree, dynamicColors, main = "Network correlation plot")
dev.off()

### For network visualization keep the top 10% of the edges
q<-quantile(upperTriangle(sim.mat, diag=FALSE),0.9);q
### Compute adjacency matrix
adj<-GL.int

### Conformity based quantities
length(conformityBasedNetworkConcepts(sim.mat)$conformityBasedNCs$Conformity)
conformityBasedNetworkConcepts(sim.mat)$fundamentalNCs$Density
conformityBasedNetworkConcepts(sim.mat)$fundamentalNCs$Centralization
conformityBasedNetworkConcepts(sim.mat)$fundamentalNCs$Heterogeneity
color="black"
length(conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$conformityBasedNCs$Conformity)
conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$fundamentalNCs$Density
conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$fundamentalNCs$Centralization
conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$fundamentalNCs$Heterogeneity
sort(colnames(sim.mat)[which(dynamicColors==color)])

### Alternatively, compute Density, Centralization, and Heterogeneity from here
quantities<-DensCentrHet1(sim.mat,dynamicColors);quantities

### Top 20 connected metabolites
mostConnected(adj>0,nr=20);

### Metabolites within each module
for (color in rownames(quantities)[-dim(quantities)[1]]){
  one<-sort(colnames(data)[which(dynamicColors==color)])
  print(color,quote=FALSE)
  cat(one, sep =", ")
  print("\n")
}


### For saving the edges to vizualize them in Gephi
write.edges(adj)


#########################################################################
################ ANALYSIS FOR DATA CONTAINING BMI INFO ##################
#########################################################################
#########################################################################
#########################################################################
################################ WGCNA ##################################
#########################################################################
#########################################################################
### WGCNA network estimation
heatmap.data<- data.bmi.info
### Tuning parameters for correlation plot
softPower = soft.thres(heatmap.data);
### What is the minimum allowed number of metabolites in a module 
minModuleSize = 5;
### Define the similarity and dissimilarity matrix 
sim.mat<-(abs(cor(heatmap.data, method="pearson",use="pairwise.complete.obs")))^softPower;diag(sim.mat)<-0
diss.mat = 1-sim.mat
### Metabolites dendrogram
metabolTree = flashClust(as.dist(diss.mat ), method = "average");
### Identification of modules
dynamicMods = cutreeDynamic(dendro = metabolTree, distM = diss.mat,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
### Visualization of correlation plot
TOMplot(diss.mat, metabolTree, dynamicColors, main = "Network correlation plot")
### How many metabolites does every cluster contain
table(dynamicColors)
### Compute PC1 of each of the clusters
MEList = moduleEigengenes(heatmap.data, colors = dynamicColors)
### The values of the PC1
MEs = MEList$eigengenes
### For keeping the names of the samples
rownames(MEs)<-rownames(heatmap.data)
### Dissimilarity of the PC1s
MEDiss = 1-cor(MEs,use="pairwise.complete.obs");
### Dendrogram of clusters
METree = hclust(as.dist(MEDiss), method = "average");
sizeGrWindow(7, 6)
### Plotting of the dendrogram
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.1
abline(h=MEDissThres, col = "red")
### Merge 2 clusters if their correlation is over 0.8
merge = mergeCloseModules(heatmap.data, dynamicColors, cutHeight = MEDissThres, verbose = 2)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
dynamicColors<-mergedColors
### Visualize the final correlation matrix
TOMplot(diss.mat, metabolTree, dynamicColors, main = "Network correlation plot")

### Save the plot as image
png(file = "Heatmap.png")
TOMplot(diss.mat, metabolTree, dynamicColors, main = "Network correlation plot")
dev.off()

### For network visualization keep the top 10% of the edges
q<-quantile(upperTriangle(sim.mat, diag=FALSE),0.9);q
### Compute adjacency matrix
adj<-topCpercent(x=sim.mat,c=0.1)

### Conformity based quantities
length(conformityBasedNetworkConcepts(sim.mat)$conformityBasedNCs$Conformity)
conformityBasedNetworkConcepts(sim.mat)$fundamentalNCs$Density
conformityBasedNetworkConcepts(sim.mat)$fundamentalNCs$Centralization
conformityBasedNetworkConcepts(sim.mat)$fundamentalNCs$Heterogeneity
color="black"
length(conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$conformityBasedNCs$Conformity)
conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$fundamentalNCs$Density
conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$fundamentalNCs$Centralization
conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$fundamentalNCs$Heterogeneity
sort(colnames(sim.mat)[which(dynamicColors==color)])

### Alternatively, compute Density, Centralization, and Heterogeneity from here
quantities<-DensCentrHet1(sim.mat,dynamicColors);quantities

### Top 20 connected metabolites
mostConnected(topCmatrix,nr=20);

### Metabolites within each module
for (color in rownames(quantities)[-dim(quantities)[1]]){
  one<-sort(colnames(heatmap.data)[which(dynamicColors==color)])
  print(color,quote=FALSE)
  cat(one, sep =", ")
  print("\n")
}


### For saving the edges to vizualize them in Gephi
write.edges(adj)

#########################################################################
#########################################################################
#########################################################################
############################### GLASSO ##################################
#########################################################################
#########################################################################

data<- data.bmi.info
### You have to omit the NA data here
data<-as.matrix(na.omit(data))

### Grid for selecting the regularization parameter
sequ<-seq(from=1,to=0.001, by=-0.001)
### Fitting GL for every lambda from the grid
GL<-huge(data,method="glasso",lambda=sequ)
### Selecting the optimal lambda
GL.net<-huge.select(est=GL,criterion="stars",rep.num=100,stars.thresh=0.02)
### Visual inspection of the network
huge.plot(GL.net$refit)
### Optimal regularization parameter
GL.net$opt.lambda
### Intensity matrix of network
GL.int<-as.matrix((GL.net$refit)*GL.net$merge[[GL.net$opt.index]]);diag(GL.int)<-0;colnames(GL.int)<-rownames(GL.int)<-colnames(data)

### What is the minimum allowed number of metabolites in a module 
minModuleSize = 5;
### Define the similarity and dissimilarity matrix 
sim.mat<-GL.int
diss.mat = 1-sim.mat
### Metabolites dendrogram
metabolTree = flashClust(as.dist(diss.mat ), method = "average");
### Identification of modules
dynamicMods = cutreeDynamic(dendro = metabolTree, distM = diss.mat,
                            pamStage = TRUE,deepSplit=4, pamRespectsDendro = TRUE,method="hybrid",
                            minClusterSize = minModuleSize);plot(metabolTree)
dynamicColors = labels2colors(dynamicMods)
### Visualization of correlation plot
TOMplot(diss.mat, metabolTree, dynamicColors, main = "Network correlation plot")
### How many metabolites does every cluster contain
table(dynamicColors)
### Compute PC1 of each of the clusters
MEList = moduleEigengenes(data, colors = dynamicColors)
### The values of the PC1
MEs = MEList$eigengenes
### For keeping the names of the samples
rownames(MEs)<-rownames(data)
### Dissimilarity of the PC1s
MEDiss = 1-cor(MEs,use="pairwise.complete.obs");
### Dendrogram of clusters
METree = hclust(as.dist(MEDiss), method = "average");
sizeGrWindow(7, 6)
### Plotting of the dendrogram
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.1
abline(h=MEDissThres, col = "red")
### Merge 2 clusters if their correlation is over 0.8
merge = mergeCloseModules(data, dynamicColors, cutHeight = MEDissThres, verbose = 2)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
dynamicColors<-mergedColors
### Visualize the final correlation matrix
TOMplot(diss.mat, metabolTree, dynamicColors, main = "Network correlation plot")

### Save the plot as image
png(file = "Heatmap.png")
TOMplot(diss.mat, metabolTree, dynamicColors, main = "Network correlation plot")
dev.off()

### For network visualization keep the top 10% of the edges
q<-quantile(upperTriangle(sim.mat, diag=FALSE),0.9);q
### Compute adjacency matrix
adj<-GL.int

### Conformity based quantities
length(conformityBasedNetworkConcepts(sim.mat)$conformityBasedNCs$Conformity)
conformityBasedNetworkConcepts(sim.mat)$fundamentalNCs$Density
conformityBasedNetworkConcepts(sim.mat)$fundamentalNCs$Centralization
conformityBasedNetworkConcepts(sim.mat)$fundamentalNCs$Heterogeneity
color="blue"
length(conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$conformityBasedNCs$Conformity)
conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$fundamentalNCs$Density
conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$fundamentalNCs$Centralization
conformityBasedNetworkConcepts(sim.mat[which(dynamicColors==color),which(dynamicColors==color)])$fundamentalNCs$Heterogeneity
sort(colnames(sim.mat)[which(dynamicColors==color)])

### Alternatively, compute Density, Centralization, and Heterogeneity from here
quantities<-DensCentrHet1(sim.mat,dynamicColors);quantities

### Top 20 connected metabolites
mostConnected(adj>0,nr=20);

### Metabolites within each module
for (color in rownames(quantities)[-dim(quantities)[1]]){
  one<-sort(colnames(data)[which(dynamicColors==color)])
  print(color,quote=FALSE)
  cat(one, sep =", ")
  print("\n")
}


### For saving the edges to vizualize them in Gephi
write.edges(adj)