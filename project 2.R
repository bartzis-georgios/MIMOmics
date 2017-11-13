rm(list=ls())
# Display the current working directory
getwd()
# If necessary, change the working directory to the one containing the data
workingDir<-"M:/Code for Projects/Project 2"
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
packages<-c("psych","rgl","lme4","flashClust","WGCNA","gdata","igraph","huge")
### Function "ipak" is defined in the end of the file. Check to see if multiple packages are installed.
### If they are not, it installs and loads them. It takes a while...
# biocLite(c("GO.db", "preprocessCore", "impute"),suppressUpdates=TRUE) 
ipak(packages)

### Load metabolite data 2007 
data07<-read.table("DILGOM 2007.txt",sep=",",header=T)
rownames(data07)=data07[,1]
data07<-data07[,2:151]
### Metabolites to be kept
a<-c("XXL.VLDL.L","XL.VLDL.L","L.VLDL.L","M.VLDL.L",
     "S.VLDL.L","XS.VLDL.L","IDL.L","IDL.C","L.LDL.L","M.LDL.L",
     "S.LDL.L","XL.HDL.L","L.HDL.L","M.HDL.L","S.HDL.L",
     "VLDL.D","LDL.D","HDL.D","Serum.C","VLDL.C",
     "Remnant.C","LDL.C","HDL.C","HDL2.C","HDL3.C",
     "Serum.TG","TotPG","PC","SM","TotCho",
     "ApoA1","ApoB","TotFA","DHA","LA",
     "FAw3","FAw6","PUFA","MUFA","SFA","FAw3.FA",
     "FAw6.FA","Glc","Lac","Pyr","Cit",
     "Glol","Ala","Gln","Gly","His",
     "Ile","Leu","Val","Phe","Tyr",
     "Ace","AcAce","bOHBut","Crea","Alb",
     "Gp")
### Keeping the metabolites we need
metabolite_list<-c()
for(i in 1:length(a)){
  metabolite_list<-c(metabolite_list,which(colnames(data07)==a[i]))
}
data07<-data07[,metabolite_list]
### Remove some metabolites since we do not have any information on the mGWAS
met.rem07<-c()
met.rem07<-c(met.rem07,which(colnames(data07)=="FAw6"))
met.rem07<-c(met.rem07,which(colnames(data07)=="FAw6.FA"))
met.rem07<-c(met.rem07,which(colnames(data07)=="Gly"))
met.rem07<-c(met.rem07,which(colnames(data07)=="MUFA"))
met.rem07<-c(met.rem07,which(colnames(data07)=="Remnant.C"))
met.rem07<-c(met.rem07,which(colnames(data07)=="SFA"))
met.rem07<-c(met.rem07,which(colnames(data07)=="VLDL.C"))
data07<-data07[,-c(met.rem07)]

### Load metabolite data 2014
data14<-read.table("DILGOM 2014.txt",sep=",",header=T)
rownames(data14)=data14[,1]
data14<-data14[,2:151]
### Keeping the metabolites we need
metabolite_list<-c()
for(i in 1:length(a)){
  metabolite_list<-c(metabolite_list,which(colnames(data14)==a[i]))
}
data14<-data14[,metabolite_list]
### Remove some metabolites since we do not have any information on the mGWAS
met.rem14<-c()
met.rem14<-c(met.rem14,which(colnames(data14)=="FAw6"))
met.rem14<-c(met.rem14,which(colnames(data14)=="FAw6.FA"))
met.rem14<-c(met.rem14,which(colnames(data14)=="Gly"))
met.rem14<-c(met.rem14,which(colnames(data14)=="MUFA"))
met.rem14<-c(met.rem14,which(colnames(data14)=="Remnant.C"))
met.rem14<-c(met.rem14,which(colnames(data14)=="SFA"))
met.rem14<-c(met.rem14,which(colnames(data14)=="VLDL.C"))
data14<-data14[,-c(met.rem14)]
### Number of metabolites
p=dim(data14)[2]
### Names of metabolites
mcn<-colnames(data07)

### Load phenotypes of 2007
pheno07<-read.csv("Pheno 2007.csv",header=T,sep=",")
### Keep the phenotypes we need
pheno07<-pheno07[,c(636,2,3,5,12,18,21)]
rownames(pheno07)=pheno07[,1]
pheno07<-pheno07[,2:7]
### Rename the Phenotypes
colnames(pheno07)[1]<-"Sex"
colnames(pheno07)[2]<-"Age"
colnames(pheno07)[3]<-"BMI"
colnames(pheno07)[4]<-"SL_GLUK"
colnames(pheno07)[5]<-"FR07_38"
colnames(pheno07)[6]<-"K34"
### Filter the individuals
pheno07$FR07_38[is.na(pheno07$FR07_38)] <- 0
pheno07<-pheno07[which(pheno07$FR07_38<=1),]
pheno07<-pheno07[which(pheno07$SL_GLUK!="NA"),]
pheno07<-pheno07[which(pheno07$SL_GLUK<=10),]
### Keep only Age and Sex
pheno07<-pheno07[,c(1:2)]
pheno07<-na.omit(pheno07)

### Load phenotypes of 2014
pheno14<-read.table("Pheno 2014.txt",sep=",",header=T)
rownames(pheno14)<-pheno14[,1]
### Keep the phenotypes we need
pheno14<-pheno14[,c(67,66,68,71,11,8)]
### Rename the Phenotypes
colnames(pheno14)[1]<-"Sex"
colnames(pheno14)[2]<-"Age"
colnames(pheno14)[3]<-"BMI"
colnames(pheno14)[4]<-"SL_GLUK"
colnames(pheno14)[5]<-"FR07_38"
colnames(pheno14)[6]<-"K34"
### Filter the individuals
pheno14<-pheno14[which(pheno14$FR07_38==1),]
pheno14<-pheno14[which(pheno14$SL_GLUK!="NA"),]
pheno14<-pheno14[which(pheno14$SL_GLUK<=10),]
pheno14<-pheno14[,c(1:2)]
pheno14<-na.omit(pheno14)

### Names of Phenotypes
pcn<-colnames(pheno14)
### For finding which people were in the 1st or both study waves
### ID numbers of people in wave 2007
x<-rownames(pheno07)
### ID numbers of people in wave 2014
y<-rownames(pheno14)
### Pool of all IDs
common.samples<-intersect(x, y)
### People in Wave 2007
unique07samples<-x[which(is.element(x,y)=="FALSE")]
### People in wave 2014
unique14samples<-y[which(is.element(y,x)=="FALSE")]
### One Age and sex per individual of the people being in both waves
ph.common<-(pheno07[common.samples,]+pheno14[common.samples,])/2
### Putting for the unique people of 2007 one age which is from the midst of the period 2007-2014
ph07unique<-pheno07[unique07samples,];ph07unique$Age<-ph07unique$Age+3.5
### Putting for the unique people of 2014 one age which is from the midst of the period 2007-2014
ph14unique<-pheno14[unique14samples,];ph14unique$Age<-ph14unique$Age-3.5
### Median age
median(rbind(ph.common,ph07unique,ph14unique)[,2])
### Phenotypes of people at 2007
pheno07<-rbind(ph.common,ph07unique)
### Phenotypes of people at 2014
pheno14<-rbind(ph.common,ph14unique)
### Sex of people of 2007
pheno07$Sex<-as.factor(pheno07$Sex)
### Sex of people of 2014
pheno14$Sex<-as.factor(pheno14$Sex)
### Make Age binary at 50
pheno07$Age<-as.numeric(pheno07$Age>=50);pheno07$Age<-as.factor(pheno07$Age)
pheno14$Age<-as.numeric(pheno14$Age>=50);pheno14$Age<-as.factor(pheno14$Age)

### Load SNP data
genotypes<-read.csv("Genotypes.csv",header=T,sep=";",row.names=1)
### Load the matrix with the effect of each SNP to the metabolites
eff.matrix<-read.csv("SNPs-metabolites.csv",header=T,sep=";")

### Load FFQ
food<-read.csv("Pheno 2007.csv",header=T,sep=",");rownames(food)<-food$ID
### Keep the following variables
food<-(food[,c(282:321,335,337,339,341,343,345,347,349,351,353,355,357,359,361,363)])
### Name the variables
colnames(food)<-c("Sweet coffeebread or pies", "Sweet biscuits","Other sweet pastry","Salty pies and pastry","Pizza",
                  "Burgers","Pasta or rice","Porridge","Cereals or muesli","Non-flavoured yoghurt","Flavoured yoghurt","Low fat cheese",
                  "Other cheese", "Ice cream, puddings","Cooked or smashed potatoes","Roasted potatoes or french fries","Vegetarian food",
                  "Cooked vegetables", "Fresh salad, fresh vegetables","Salad dressing or oil","Fruits","Fresh or frosen berries",
                  "Fruit and berry juices","Fish and other fishfood combined", "Salmon, rainbow trout","Herring", "Other fish",
                  "Meat","Poultry meat","Sausages","Cutlets","Cold cuts","Eggs","Chocolate","Other candies","Salty snacks","Store boughtready meal",
                  "Rye bread, rye crisp","Yeast bread,graham bread", "White bread","Coffee/day","Tea/day","Chocolate milk/day","Milk/day",
                  "Sour milk/day","Tap water/day","Well water/day","Bottled water/day","Juice/day","Energy drink/day",
                  "Low alcohol /day","Cola/day","Cola light/day","Soft drink/day","Low calory soft drink/day")
### Remove cases with NAs
food<-na.omit(food)

### Factor analysis in food questionnaire
### Treat food frequencies as numeric
food[] <- lapply(food, as.numeric)
### Scale the FFQ
food<-scale(food)
### Store the correlation matrix of the FFQ
corMat<-cor(food)
### How many factors to be asked
factors=6
### Apply Factor Analysis
solution <- fa(r = corMat, nfactors = factors, rotate = "oblimin", fm = "pa",scores="tenBerge")
### Make the Scree Plot
scree.p<-cbind(1:55,solution$e.values)
### Visualize the Scree Plot
plot(scree.p,pch = 19,type = "b",main="Scree plot",xlab="Number of factors",ylab="Eigenvalues of factors")
### Store the Loadings
f.loadings<-solution$loadings
### Print the Loadings sorted
ICLUST.sort (f.loadings)
### Store the ten Berge Factor Scores
F1<-factor.scores(food, solution, method = "tenBerge")$score[,1];
F2<-factor.scores(food, solution, method = "tenBerge")$score[,2];
F3<-factor.scores(food, solution, method = "tenBerge")$score[,3];
F4<-factor.scores(food, solution, method = "tenBerge")$score[,4];
F5<-factor.scores(food, solution, method = "tenBerge")$score[,5];
F6<-factor.scores(food, solution, method = "tenBerge")$score[,6];
### 3D plot of the first 3 Factors
plot3d(F1, F2, F3, col="red", size=3)

### Merge the data of 2007 so they are in the same order
merged <- merge(data07,pheno07,by=0);rownames(merged) <- merged[,1];merged <- merged[,-1];
### Merge all diets to the merged file
merged <-merge(merged ,F1,by=0);rownames(merged ) <- merged [,1];merged <- merged [,-1];colnames(merged )[p+2+1]<-"F1"
merged <-merge(merged ,F2,by=0);rownames(merged ) <- merged [,1];merged <- merged [,-1];colnames(merged )[p+2+2]<-"F2"
merged <-merge(merged ,F3,by=0);rownames(merged ) <- merged [,1];merged <- merged [,-1];colnames(merged )[p+2+3]<-"F3"
merged <-merge(merged ,F4,by=0);rownames(merged ) <- merged [,1];merged <- merged [,-1];colnames(merged )[p+2+4]<-"F4"
merged <-merge(merged ,F5,by=0);rownames(merged ) <- merged [,1];merged <- merged [,-1];colnames(merged )[p+2+5]<-"F5"
merged <-merge(merged ,F6,by=0);rownames(merged ) <- merged [,1];merged <- merged [,-1];colnames(merged )[p+2+6]<-"F6"
### Separate the metabolite data 
data07<-merged[,1:p];colnames(data07)<-mcn
### Separate the Phenotype data
pheno07<-merged[,(p+1):(p+2)];colnames(pheno07)<-pcn;pheno07[,1]<-as.factor(pheno07[,1])
### Separate the Diet data
F07<-merged[,(p+3):(p+3+5)];


### Merge the data of 2014 so they are in the same order
merged <- merge(data14,pheno14,by=0);rownames(merged) <- merged[,1];merged <- merged[,-1];
### Merge all diets to the merged file
merged <-merge(merged ,F1,by=0);rownames(merged ) <- merged [,1];merged <- merged [,-1];colnames(merged )[p+2+1]<-"F1"
merged <-merge(merged ,F2,by=0);rownames(merged ) <- merged [,1];merged <- merged [,-1];colnames(merged )[p+2+2]<-"F2"
merged <-merge(merged ,F3,by=0);rownames(merged ) <- merged [,1];merged <- merged [,-1];colnames(merged )[p+2+3]<-"F3"
merged <-merge(merged ,F4,by=0);rownames(merged ) <- merged [,1];merged <- merged [,-1];colnames(merged )[p+2+4]<-"F4"
merged <-merge(merged ,F5,by=0);rownames(merged ) <- merged [,1];merged <- merged [,-1];colnames(merged )[p+2+5]<-"F5"
merged <-merge(merged ,F6,by=0);rownames(merged ) <- merged [,1];merged <- merged [,-1];colnames(merged )[p+2+6]<-"F6"
### Separate the metabolite data 
data14<-merged[,1:p];colnames(data14)<-mcn
### Separate the Phenotype data
pheno14<-merged[,(p+1):(p+2)];colnames(pheno14)<-pcn;pheno14[,1]<-as.factor(pheno14[,1])
### Separate the Diet data
F14<-merged[,(p+3):(p+3+5)];

### Store the IDs of samples used
samples<-union(rownames(data07),rownames(data14))
### Prepare data for Linear Mixed Model (all data)
### Make the subject's ID factor for 2007
subject<-factor(rownames(data07))
### Put all the data together for 2007
data07<-cbind(data07,subject,pheno07,F07)
### Put time=0 for everyone in 2007
data07$time<-factor(c(rep(0,dim(data07)[1])))
### Number of variables in 2007
p07<-dim(data07)[2]
### Number of samples in 2007
n07<-dim(data07)[1]
### Make the subject's ID factor for 2014
subject<-factor(rownames(data14))
### Put all the data together for 2014
data14<-cbind(data14,subject,pheno14,F14)
### Put time=1 for everyone in 2014
data14$time<-factor(c(rep(1,dim(data14)[1])))
### Number of variables in 2014
p14<-dim(data14)[2]
### Number of samples in 2014
n14<-dim(data14)[1]

### Put both waves together
metabol<-rbind(data07,data14)
### Calculate how many people are in total
n<-length(unique(metabol$subject))

### Give the same score in the diet to each person for both waves (People keeping the same dietary patterns)
diets<-matrix(0,1,6);colnames(diets)<-c("F1","F2","F3","F4","F5","F6");
for (r in 1:n){
  id<-which(rownames(metabol)==unique(metabol$subject)[r])
  diets<-rbind(diets,metabol[id,(p+4):(p+4+5)])
}
diets<-diets[-1,]

### Make matrix for computing the Polygenic effects
data<-matrix(0,n,p)
rownames(data)<-unique(metabol$subject)

### Compute polygenic scores
Pol.sc<-matrix(0,dim(metabol)[1],p);dim(Pol.sc)

### For every metabolite do the following
for (i in 1:p){
### Name of the metabolite  
  mcn[i]
### Which are the SNPs associated  
  rsid<-eff.matrix$rsid[which(eff.matrix$trait==mcn[i])];rsid
### How many SNPs are associated  
  nsnps<-length(rsid);nsnps
### Take the metabolite concentrations
  y<-metabol[,i];y
### If there is at least 1 SNP associated  
  if(nsnps!=0){
### Recover the position of the SNP the effect, the allele frequency and check if it needs recoding    
    snp.info<-matrix(0,nsnps,4);snp.info
    colnames(snp.info)<-c("position","beta","effect allele frequency","recode")
    rownames(snp.info)<-rsid;snp.info
### Find these info
      for (j in 1:nsnps){
      snp.info[j,1]<-grep(rsid[j],colnames(genotypes),ignore.case = TRUE)
      snp.info[j,2]<-eff.matrix$beta[which((eff.matrix$trait==mcn[i])&(eff.matrix$rsid==rsid[j]))]
      snp.info[j,3]<-eff.matrix$coded_af[which((eff.matrix$trait==mcn[i])&(eff.matrix$rsid==rsid[j]))]
      snp.info[j,4]<-2*(snp.info[j,3]>0.5)
      }
### Order all samples    
    ordering<-c()
    for (k in 1:length(metabol$subject)){
      ordering<-c(ordering,which(rownames(genotypes)==metabol$subject[k]))
    }
### Get the SNP data for it 
    x<-genotypes[ordering,snp.info[,1]];x
### Compute the PRS    
    G<-rowSums(t(snp.info[,2]*t(x)));G
### Scale the PRS    
    G<-scale(G)
  }
### Store the PRS  
  Pol.sc[,i]<-G
### If no SNPs are associated  
  if(nsnps==0){
    G<-rep(0,times=dim(metabol)[1]);G
  }
}
### Take the matrix of the PRSs
G<-Pol.sc
### Put names for the variables
colnames(G)<-mcn
### Put names for the samples
rownames(G)<-metabol$subject


### Fitting the Linear mixed effects model
### For storing all coefficients of all terms of the model
coefficients<-matrix(0,p,40);colnames(coefficients)<-c("Age","G","S","F1","F2","F3","F4","F5","F6","time","F1G","F2G","F3G","F4G","F5G","F6G","F1A","F2A","F3A","F4A","F5A","F6A","F1S","F2S","F3S","F4S","F5S","F6S","AS","GS","F1T","F2T","F3T","F4T","F5T","F6T","ST","AT","GT","AG")
### For storing the coefficient of the Genetic part
betaG<-c()
### For storing all random intercepts 
random.effects<-matrix(0,n,1)
### Names of the samples
rownames(random.effects)<-unique(metabol$subject)
### Storing the variance of the model
vASGF<-c()
### Storing the residual variance of the model
res.vASGF<-c()

### For every metabolite fit a linear mixed effects model
for (i in 1:p){
### Which is the metabolite  
  y=metabol[,i]
### What are the SNPs associated
  rsid<-eff.matrix$rsid[which(eff.matrix$trait==mcn[i])];rsid
### How many SNPs are these  
  nsnps<-length(rsid);nsnps
### If there are associated SNPs you need to put the PRS  
  if(nsnps!=0){
### MIxed Model    
    modelASGF<-lmer(y ~ 1 + Age + G[,i] +Sex+F1 + F2 + F3 + F4 + F5 + F6 + F1:G[,i] + F2:G[,i] + F3:G[,i] + F4:G[,i] +
                      F5:G[,i] + F6:G[,i] + F1:Age + F2:Age + F3:Age + F4:Age + F5:Age + F6:Age + 
                      F1:Sex + F2:Sex + F3:Sex + F4:Sex + F5:Sex + F6:Sex + 
                      Sex:Age + Sex:G[,i] +F1:time + F2:time + F3:time + F4:time + F5:time + F6:time + 
                      Sex:time + Age:time + time:G[,i]+Age:G[,i]+
                      time + (1 | subject), data=metabol)
### Store variance components    
    vASGF<-c(vASGF,as.data.frame(VarCorr(modelASGF))[1,4])
### Store residual variance    
    res.vASGF<-c(res.vASGF,as.data.frame(VarCorr(modelASGF))[2,4])
### Store the coefficients    
    coefficients[i,]<-summary(modelASGF)$coef[2:41,1];
### Store the random effects    
    random.effects<-merge(random.effects,as.matrix((coef(modelASGF))$subject)[,1],by=0)
    rownames(random.effects) <- random.effects[,1];random.effects<- random.effects[,-1];
    colnames(random.effects)[dim(random.effects)[2]]<-mcn[i]
  }
### If there are no associated SNPs you should omit the G term
  if(nsnps==0){
### Mixed model in this case    
    modelASGF<-lmer(y ~ 1 + Age + Sex+F1 + F2 + F3 + F4 + F5 + F6 + 
                      F1:Age + F2:Age + F3:Age + F4:Age + F5:Age + F6:Age + 
                      F1:Sex + F2:Sex + F3:Sex + F4:Sex + F5:Sex + F6:Sex + 
                      Sex:Age + F1:time + F2:time + F3:time + F4:time + F5:time + F6:time + 
                      Sex:time + Age:time + 
                      time + (1 | subject), data=metabol)
### Variance components    
    vASGF<-c(vASGF,as.data.frame(VarCorr(modelASGF))[1,4])
### Residual variance    
    res.vASGF<-c(res.vASGF,as.data.frame(VarCorr(modelASGF))[2,4])
### Store coefficients    
    coefficients[i,1]<-summary(modelASGF)$coef[2,1]
    coefficients[i,3]<-summary(modelASGF)$coef[3,1]
    coefficients[i,4:10]<-summary(modelASGF)$coef[c(4:10),1]
    coefficients[i,17:28]<-summary(modelASGF)$coef[c(11:22),1]
    coefficients[i,29]<-summary(modelASGF)$coef[c(23),1]
    coefficients[i,31:38]<-summary(modelASGF)$coef[c(24:31),1]
### Store random effects    
    random.effects<-merge(random.effects,as.matrix((coef(modelASGF))$subject)[,1],by=0)
    rownames(random.effects) <- random.effects[,1];random.effects<- random.effects[,-1];
    colnames(random.effects)[dim(random.effects)[2]]<-mcn[i]
  }
}
### Put as rownames the names of the metabolites
rownames(coefficients)<-mcn
### Remove the empty 1st column of the matrix
random.effects<-random.effects[,-1]
### Store the coefficients of the PRS
betaG<-coefficients[,2]


### Visualize the intraclass correlation
plot(vASGF/(res.vASGF+vASGF),type="b",xaxt = "n",xlab="",ylab="Intraclass Correlation",ylim=c(0,1))
axis(1, at=1:p, labels=mcn)


### Visualize the genetic influence
plot(betaG^2*apply(G, 2, var),type="b",ylab="% of explained variance",xlab="metabolites",xaxt = "n")
axis(1, at=1:p, labels=colnames(G))

### Keeping diet related part
### The final samples to be used
samples<-rownames(metabol)
### Make a data matrix
data<-matrix(0,length(samples),1)
### Sample names of the matrix
rownames(data)<-samples
### The random effects 
random.effects<-random.effects[metabol$subject,]
### Sample names of the random effects
rownames(random.effects)<-samples

### For every metabolite find the dietary and lifestyle part
for (i in 1:p){
### Covariates for all samples  
  Fs<-metabol[,(p+2):(p+10)]
### Age for all samples  
  As<-matrix(0,length(samples),1);rownames(As)<-samples;As[,1]<-Fs[samples,2];colnames(As)<-"Age";As<-As-1
### Sex for all samples  
  Ss<-matrix(0,length(samples),1);rownames(Ss)<-samples;Ss[,1]<-Fs[samples,1];colnames(Ss)<-"Sex";Ss<-Ss-1
### Time for all samples  
  Ts<-matrix(0,length(samples),1);rownames(Ts)<-rownames(metabol);Ts[,1]<-metabol$time;colnames(Ts)<-"Time";Ts<-Ts-1
###   Diets for all samples
  Fs<-Fs[samples,3:8]
### PRS for all samples  
  Gs<-G[metabol$subject,i];names(Gs)<-samples
### Age by  Sex as covariate 
  ASs<-merge(As,Ss,by=0);rownames(ASs) <- ASs[,1];ASs <- ASs [,-1];Ns<-rownames(ASs)
  ASs<-(data.matrix(ASs[,1:length(ASs)-1], rownames.force = NA))*ASs[,length(ASs)];rownames(ASs)<-Ns
### Diet by time as covariate  
  FTs<-merge(Fs,Ts,by=0);rownames(FTs) <- rownames(metabol);FTs <- FTs [,-1];colnames(FTs)[dim(FTs)[2]]<-colnames(metabol)[i]
  FTs<-(data.matrix(FTs[,1:length(FTs)-1], rownames.force = NA))*FTs[,length(FTs)]
### Diet by Genetics as covariate  
  FGs<-merge(Fs,Gs,by=0);rownames(FGs) <- FGs[,1];FGs <- FGs [,-1];colnames(FGs)[dim(FGs)[2]]<-colnames(metabol)[i]
  FGs<-(data.matrix(FGs[,1:length(FGs)-1], rownames.force = NA))*FGs[,length(FGs)]
### Diet by Age as covariate  
  FAs<-merge(Fs,As,by=0);rownames(FAs) <- FAs[,1];FAs <- FAs [,-1];
  FAs<-(data.matrix(FAs[,1:length(FAs)-1], rownames.force = NA))*FAs[,length(FAs)]
### Diet by Sex as covariate  
  FSs<-merge(Fs,Ss,by=0);rownames(FSs) <- FSs[,1];FSs <- FSs [,-1];
  FSs<-(data.matrix(FSs[,1:length(FSs)-1], rownames.force = NA))*FSs[,length(FSs)]
### Time by Genetics as covariate  
  TGs<-merge(Gs,Ts,by=0);rownames(TGs) <- TGs[,1];TGs <- TGs [,-1];Ns<-rownames(TGs)
  TGs<-(data.matrix(TGs[,1:length(TGs)-1], rownames.force = NA))*TGs[,length(TGs)];rownames(TGs)<-Ns
### Time by Age as covariate  
  TAs<-merge(Ts,As,by=0);rownames(TAs) <- TAs[,1];TAs <- TAs [,-1];Ns<-rownames(TAs)
  TAs<-(data.matrix(TAs[,1:length(TAs)-1], rownames.force = NA))*TAs[,length(TAs)];rownames(TAs)<-Ns
### Time by Sex as covariate
  TSs<-merge(Ts,Ss,by=0);rownames(TSs) <- TSs[,1];TSs <- TSs [,-1];Ns<-rownames(TSs)
  TSs<-(data.matrix(TSs[,1:length(TSs)-1], rownames.force = NA))*TSs[,length(TSs)];rownames(TSs)<-Ns
### Genetics by Age as covariate
  GAs<-merge(Gs,As,by=0);rownames(GAs) <- GAs[,1];GAs <- GAs [,-1];colnames(GAs)[dim(GAs)[2]]<-colnames(metabol)[i];Ns<-rownames(GAs)
  GAs<-(data.matrix(GAs[,1:length(GAs)-1], rownames.force = NA))*GAs[,length(GAs)];rownames(GAs)<-Ns
### Genetics by Sex as covariate  
  GSs<-merge(Gs,Ss,by=0);rownames(GSs) <- GSs[,1];GSs <- GSs [,-1];colnames(GSs)[dim(GSs)[2]]<-colnames(metabol)[i];Ns<-rownames(GSs)
  GSs<-(data.matrix(GSs[,1:length(GSs)-1], rownames.force = NA))*GSs[,length(GSs)];rownames(GSs)<-Ns
### Merge all covariates into one data matrix  
  X<-merge(As,Gs,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-merge(X,Ss,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-merge(X,Fs,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-merge(X,Ts,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-merge(X,FGs,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-merge(X,FAs,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-merge(X,FSs,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-merge(X,ASs,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-merge(X,GSs,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-merge(X,FTs,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-merge(X,TSs,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-merge(X,TAs,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-merge(X,TGs,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-merge(X,GAs,by=0);rownames(X) <- X[,1];X <- X [,-1];colnames(X)<-colnames(coefficients)[1:length(X)]
  X<-data.matrix(X)
### Puting one and zero if we want or not to take into account that term  in the fitted values
  coefficients<-as.data.frame(coefficients)
  coefficients$Age <- coefficients$Age  *0
  coefficients$G   <- coefficients$G    *0
  coefficients$S   <- coefficients$S    *0
  coefficients$F1  <- coefficients$F1   *01
  coefficients$F2  <- coefficients$F2   *01
  coefficients$F3  <- coefficients$F3   *01
  coefficients$F4  <- coefficients$F4   *01
  coefficients$F5  <- coefficients$F5   *01
  coefficients$F6  <- coefficients$F6   *01
  coefficients$time<- coefficients$time *0
  coefficients$F1G <- coefficients$F1G  *01
  coefficients$F2G <- coefficients$F2G  *01
  coefficients$F3G <- coefficients$F3G  *01
  coefficients$F4G <- coefficients$F4G  *01
  coefficients$F5G <- coefficients$F5G  *01
  coefficients$F6G <- coefficients$F6G  *01
  coefficients$F1A <- coefficients$F1A  *01
  coefficients$F2A <- coefficients$F2A  *01
  coefficients$F3A <- coefficients$F3A  *01
  coefficients$F4A <- coefficients$F4A  *01
  coefficients$F5A <- coefficients$F5A  *01
  coefficients$F6A <- coefficients$F6A  *01
  coefficients$F1S <- coefficients$F1S  *01
  coefficients$F2S <- coefficients$F2S  *01
  coefficients$F3S <- coefficients$F3S  *01
  coefficients$F4S <- coefficients$F4S  *01
  coefficients$F5S <- coefficients$F5S  *01
  coefficients$F6S <- coefficients$F6S  *01
  coefficients$AS  <- coefficients$AS   *0
  coefficients$GS  <- coefficients$GS   *0
  coefficients$F1T <- coefficients$F1T  *1
  coefficients$F2T <- coefficients$F2T  *1
  coefficients$F3T <- coefficients$F3T  *1
  coefficients$F4T <- coefficients$F4T  *1
  coefficients$F5T <- coefficients$F5T  *1
  coefficients$F6T <- coefficients$F6T  *1
  coefficients$ST  <- coefficients$ST   *0
  coefficients$AT  <- coefficients$AT   *0
  coefficients$GT  <- coefficients$GT   *0
  coefficients$AG  <- coefficients$AG   *0
### Fitted values
  data<-merge(data,rowSums(X%*%diag(coefficients[i,])),by=0);rownames(data) <- data[,1];data <- data [,-1]
  colnames(data)[length(data)]<-colnames(metabol)[i]
}
data<-data[,-1]

### merge diets, metabolites so they are in same order
merged<-cbind(data[samples,],Fs[samples,])
data<-merged[,1:p];diets<-merged[,(p+1):length(merged)]

### Take the diet part
dietary<-data
### Take liestyle part
lifestyle<-random.effects
### Take the diet and lifestyle part
dietLife<-dietary+lifestyle
### Remove NAs if any
dietLife<-na.omit(dietLife)
dietary<-na.omit(dietary)
lifestyle<-na.omit(lifestyle)


#########################################################################
########## ANALYSIS FOR DATA CONTAINING DIET AND LIFESTYLE ##############
#########################################################################
#########################################################################
#########################################################################
################################ WGCNA ##################################
#########################################################################
#########################################################################
### WGCNA network estimation
heatmap.data<- dietLife
### Tuning parameters for correlation plot
softPower = soft.thres(heatmap.data);
### What is the minimum allowed number of metabolites in a module 
minModuleSize = 3;
### Define the similarity and dissimilarity matrix 
sim.mat<-(abs(cor(heatmap.data, method="pearson",use="pairwise.complete.obs")))^softPower;diag(sim.mat)<-0
diss.mat = 1-sim.mat
### Metabolites dendrogram
metabolTree = flashClust(as.dist(diss.mat ), method = "average");
### Identification of modules
dynamicMods = cutreeDynamic(dendro = metabolTree, distM = diss.mat,
                            pamStage = TRUE,deepSplit=4, pamRespectsDendro = TRUE,method="hybrid",
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

### Compute Density, Centralization, and Heterogeneity from here
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
########################################################
### Profile plots of modules by diet
par(mfrow=c(7,6))
diets

for (i in rownames(quantities)[1:7]){
  for (diet in 1:6){
    y=(scale(diets[,diet]))
    plot(c(-1.6:1.6),c(-1.6:1.6),type="n",xlab="",ylab="",ylim=range(-1.1:1.1))
    for (j in which(dynamicColors==i)){
      x=(scale(dietary[,j]))
      a<-cbind(x,y)
      a<-a[which(abs(a[,2])<1.5),]
      smoothingSpline = smooth.spline(a[,2], a[,1], spar=0.99)
      lines(smoothingSpline,col=i)
      #  lab <- colnames(data)[j]
      #  text(locator(1), lab, adj=0) 
    }
  }
}


#########################################################################
#########################################################################
#########################################################################
############################### GLASSO ##################################
#########################################################################
#########################################################################

data<- dietLife
### You have to omit the NA data here
data<-as.matrix(na.omit(data))

### Grid for selecting the regularization parameter
sequ<-seq(from=1,to=0.001, by=-0.001)
### Fitting GL for every lambda from the grid
GL<-huge(data,method="glasso",lambda=sequ)
### Selecting the optimal lambda
GL.net<-huge.select(est=GL,criterion="stars",rep.num=100,stars.thresh=0.05)
### Visual inspection of the network
huge.plot(GL.net$refit)
### Optimal regularization parameter
GL.net$opt.lambda
### Intensity matrix of network
GL.int<-as.matrix((GL.net$refit)*GL.net$merge[[GL.net$opt.index]]);diag(GL.int)<-0;colnames(GL.int)<-rownames(GL.int)<-colnames(data)

### What is the minimum allowed number of metabolites in a module 
minModuleSize = 3;
### Define the similarity and dissimilarity matrix 
sim.mat<-GL.int
diss.mat = 1-sim.mat
### Metabolites dendrogram
metabolTree = flashClust(as.dist(diss.mat ), method = "average");
### Identification of modules
dynamicMods <- cutreeDynamic(dendro = metabolTree,method="hybrid", distM = diss.mat,pamRespectsDendro = FALSE, 
                             pamStage = TRUE,deepSplit=4, minClusterSize = minModuleSize);
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
MEDissThres = 0.2
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


### Profile plots of modules by diet
par(mfrow=c(8,6))

for (i in rownames(quantities)[1:8]){
  for (diet in 1:6){
    y=(scale(diets[,diet]))
    plot(c(-1.6:1.6),c(-1.6:1.6),type="n",xlab="",ylab="",ylim=range(-1.1:1.1))
    for (j in which(dynamicColors==i)){
      x=(scale(dietary[,j]))
      a<-cbind(x,y)
      a<-a[which(abs(a[,2])<1.5),]
      smoothingSpline = smooth.spline(a[,2], a[,1], spar=0.99)
      lines(smoothingSpline,col=i)
      #  lab <- colnames(data)[j]
      #  text(locator(1), lab, adj=0) 
    }
  }
}
