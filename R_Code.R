#Steps to download the dataset :

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install()

#BiocManager::install("antiProfilesData")

dataset <- antiProfilesData::apColonData

#Question 1:
#Explore the data:
#Show the type of the data and assign each of feature, phenotype and expression data to different variables
print(dataset)
print(class(dataset))
print(typeof(dataset))
print(str(dataset))
print(summary(dataset))

pdata=pData(dataset)#phenotype data (samples): sample X information
edata=exprs(dataset)#expression data (genomics data): features(genes) X sample, [i,j] = count gene/sample
fdata = fData(dataset)# features data: genomic features: features(genes) X features , feature1: gene name

pdataDim = dim(pdata)
edataDim = dim(edata)
fdataDim = dim(fdata)

pdataCols = colnames(pdata)
edataCols = colnames(edata)
fdataCols = colnames(fdata)

pdataRows = rownames(pdata)
edataRows = rownames(edata)
fdataRows = rownames(fdata)

#1.a Show the type of each column
for (i in pdataCols){
  temp = pdata[,i]
  cat(i)
  cat(" : ")
  cat(class(temp))
  cat("\n")
}
for (i in edataCols){
  temp = edata[,i]
  cat(i)
  cat(" : ")
  cat(class(temp))
  cat("\n")
}
for (i in fdataCols){
  temp = fdata[,i]
  cat(i)
  cat(" : ")
  cat(class(temp))
  cat("\n")
}

#1.b Show column names and rows name
print(pdataCols)
print(pdataRows)
print(edataCols)
print(edataRows)
print(fdataCols)
print(fdataRows)

#1.c Calculate summary of each column
for (i in pdataCols){
  temp = pdata[,i]
  print(i)
  print(summary(temp))
}
for (i in edataCols){
  temp = edata[,i]
  print(i)
  print(summary(temp))
}
for (i in fdataCols){
  temp = fdata[,i]
  print(i)
  print(summary(temp))
}

#1.d Show frequency of categorical data, taking into the consideration, NA values frequency if any.
for (i in pdataCols){
  temp = pdata[,i]
  if(class(temp) == "character"){
    print(i)
    print(table(factor(temp),useNA="ifany"))
  }
}
for (i in edataCols){
  temp = edata[,i]
  if(class(temp) == "character"){
    print(i)
    print(table(factor(temp),useNA="ifany"))
  }
}
for (i in fdataCols){
  temp = fdata[,i]
  if(class(temp) == "character"){
    print(i)
    print(table(factor(temp),useNA="ifany"))
  }
}

#1.e Calculate the correlation and covariance between the first 10 columns only of our data set and draw full correlation matrix.
print(cor(edata[,c(1:10)]))
print(cov(edata[,c(1:10)]))
res = cor(edata)
library(corrplot)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res, col = col, symm = T)

#1.f For both samples: GSM95478,GSM95473 show the plot with a line of their relation.
par(mfrow =c(1,2))

lm = lm(edata[,'GSM95473']~edata[,'GSM95478'] )
plot(edata[,'GSM95478'],edata[,'GSM95473'],col=2,main = "Relation between two samples and best fit line")
abline(lm,col=5,lwd=3)

lm = lm(edata[,'GSM95478']~edata[,'GSM95473'] )
plot(edata[,'GSM95473'],edata[,'GSM95478'],col=2,main = "Relation between two samples and best fit line")
abline(lm,col=5,lwd=3)

#Question 2:
#Using PCA and SVD, Prove by plotting and values that both can return the same result by
#suitable normalization.

#Calculate the singular vectors,Here we calculate the singular vectors: 
edata_log = log2(edata + 5)
edata_centered = edata_log - rowMeans(edata_log)
svd1 = svd(edata_centered)
names(svd1)

#Plot top two principal components
par(mfrow=c(1,2))
plot(svd1$v[,1],col=5,ylab="1st PC")
plot(svd1$v[,2],col=6,ylab="2nd PC")

#Plot PC1 vs. PC2,A very common plot is to plot PC1 versus PC2 to see if you can see any "clusters" or "groups".
#One thing you can do is color them by different variables to see if clusters stand out. 
par(mfrow=c(1,1))
plot(svd1$v[,1],svd1$v[,2],ylab="2nd PC", xlab="1st PC",col=c(5,6))

#PCs versus SVs,What we have been plotting is not exactly the principal components. 
#Here before normalization, each gives different outputs
par(mfrow=c(1,2))
pc1 = prcomp(edata_log)
plot(pc1$rotation[,1],svd1$v[,1],col=c(5,6),main= "PC1 VS SVD1 before Normalization")
plot(pc1$rotation[,2],svd1$v[,2],col=c(5,6),main= "PC2 VS SVD2 before Normalization")

#To get the actual PCs you have to subtract the column means rather than the row means when normalizing. 
#Here after normalization, each give same results
par(mfrow=c(1,2))
edata_centered2 = t(t(edata_log) - colMeans(edata_log))
svd2 = svd(edata_centered2)
plot(pc1$rotation[,1],svd2$v[,1],col=c(5,6),main= "PC1 VS SVD1 after Normalization")
plot(pc1$rotation[,2],svd2$v[,2],col=c(5,6),main= "PC2 VS SVD2 after Normalization")

#Question 3:
#256 visual artists were surveyed to find out their zodiac sign. The results were: Aries (29),
#Taurus (24), Gemini (22), Cancer (19), Leo (21), Virgo (18), Libra (19), Scorpio (20), Sagittarius
#(23), Capricorn (18), Aquarius (20), Pisces (23).
#3.1) Test the hypothesis that zodiac signs are evenly distributed across visual artists.
#3.2) Explicitly mention your H1 and Ho assumption

#It needs to approve hypothesis H1 that a certain assumption is founded in the population by disapprove hypothesis
#Ho (is called null Hypothesis) which is the opposite of H1. we uses significance level to accept or reject our assumption.

#What does H0 mean in statistics?
#The null hypothesis (H0) is a statement of "no difference," "no association," or "no treatment effect." 
#. The alternative hypothesis, Ha is a statement of "difference," "association," or "treatment effect." 
#Ho is assumed to be true until proven otherwise. However, Ha is the hypothesis the researcher hopes to bolster.
#Ho :Prototypes occur at a ratio of 3:2:2:1:2:1:1:1:2:1:1:2

#What does H1 represent in statistics?
#The alternative hypothesis, H1 or Ha, is a statistical proposition stating that there is a significant difference 
#between a hypothesized value of a population parameter and its estimated value.
#H1 :Prototypes occur at a ratio of 29:24:22:19:21:18:19:20:23:18:20:23


# Null Hypothesis(Ho): Zodiac signs are evenly distributed across visual artist.

# Alternative Hypothesis(H1): Zodiac signs are not evenly distributed  across visual artist.

Phenotypes <- as.factor(c(rep("Aries" ,29),rep("Taurus",24) , rep("Gemini",22),rep("Cancer",19),rep("Leo",21),
rep("Virgo",18),rep("Libra",19),rep("Scorpio",20),rep("Sagittarius",23),rep("Capricorn",18),rep("Aquarius",20),
rep("Pisces",23)))
p <- c(1,1,1,1,1,1,1,1,1,1,1,1)
p <-  p/sum(p)
table(Phenotypes)
chisq.test(table(Phenotypes), p=p)

#Question 4:
#Plot hierarchical clusters on our first 10 columns of edata and apply the kmeans to all the edata
#columns and show the centroid of the result.

# By default calculates the distance between rows so we transpose to take the samples distances
par(mfrow =c(1,2))
dist1 = dist(t(edata[,c(1:10)]))
#Hierarchical clustering
hclust1 = hclust(dist1)
plot(hclust1,hang = -1) # hang = - 1 to make labels written on the same level
plot(hclust1)

#K-Means Clustering:
kmeans1 = kmeans(edata,centers=7)
#Comparing number of data points assigned to each cluster
print(table(kmeans1$cluster))
