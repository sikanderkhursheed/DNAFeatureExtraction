#Setup and library imports 
library(dplyr)
library(stringr)
library(ngram)
library(stringr)
library(e1071)
library(klaR)
library(caret)
library(mlbench)
library(Boruta)
library(randomForest)
library(qualV)
library(tcR)
library(reshape2)

#Loading the dataset

setwd("C:\\Users\\HP\\Documents\\")
Data<-read.csv("AM_repeats_clean.csv", header=TRUE,stringsAsFactors=FALSE)

#Structure of the dataset
str(Data)


#Feature Extraction
Following are the feature extraction approaches.

####1. Bag-of approach
u <- unique(unlist(strsplit(Data$Sequence,split="")))
bag <- Data
count_Char <- function(string, char){
  count = 0
  for (i in string){
    if (i == char){
      count = count + 1
    }
  }
  return(count)
}

bagof <- function(dataset, colnum){
  x = 1
  while (x <= nrow(dataset)){
    for (alpha in u){
      dataset[x,alpha] <- count_Char(unlist(strsplit(dataset[x,colnum],"")), alpha)
    }
    x = x + 1
  }
  return(dataset)
}
bag<-bagof(bag, 2)

#For plotting bag of approach
d <- melt(bag)
ggplot(d, aes(variable,value, col=variable)) +
  geom_point()  +
  stat_smooth()



####2. N-Gram approach

#####2.1 Bi-Gram approach
bi <- get.phrasetable(ngram(Data$Sequence, n=2, sep = ""))
bigramdf<-Data

#For removing the space between the ngrams
for (i in 1:nrow(bi)) {
  bi[i,1] <- gsub(' ', '',bi[i,1])
}

#Create a dataframe with initially all values = 0
for (n in bi[,1]){
  bigramdf[,n]<-0
}

x = 1
while (x <= nrow(bigramdf)){
  for (n in 1:nrow(bi)){
    
    bigramdf[x,n+2] <- str_count(bigramdf[x,2], bi[n,1])
    
  }
  x = x + 1
}
dim(bigramdf)


#####2.2 Tri-Gram approach

tri <- get.phrasetable(ngram(Data$Sequence, n=3, sep = ""))
trigramdf <- Data

#For removing the space between the ngrams
for (i in 1:nrow(tri)) {
  tri[i,1] <- gsub(' ', '',tri[i,1])
}

#Create a dataframe with initially all values = 0
for (n in tri[,1]){
  trigramdf[,n]<-0
}

x = 1
while (x <= nrow(trigramdf)){
  for (n in 1:nrow(tri)){
    
    trigramdf[x,n+2] <- str_count(trigramdf[x,2], tri[n,1])
    
  }
  x = x + 1
}
dim(trigramdf)

#####2.3 Four-Gram approach
four <- get.phrasetable(ngram(Data$Sequence, n=4, sep = ""))
gramdf4 <- Data

#For removing the space between the ngrams
for (i in 1:nrow(four)) {
  four[i,1] <- gsub(' ', '',four[i,1])
}

#Create a dataframe with initially all values = 0
for (n in four[,1]){
  gramdf4[,n]<-0
}

x = 1
while (x <= nrow(gramdf4)){
  for (n in 1:nrow(four)){
    
    gramdf4[x,n+2] <- str_count(gramdf4[x,2], four[n,1])
    
  }
  x = x + 1
}

dim(gramdf4)

####3. Finding the first character of each sequence 
ch <- unlist(substr(Data$Sequence, 1 , 1))

i=1
for(i in ch){
  ch <- str_replace(ch,'a', '1')
  ch <- str_replace(ch,'d', '2')
  ch <- str_replace(ch,'g', '3')
  ch <- str_replace(ch,'l', '4')
  ch <- str_replace(ch,'n', '5')
  ch <- str_replace(ch,'s', '6')
  ch <- str_replace(ch,'t', '7')
}
ch<-data.frame(ch)
ch <- as.numeric(unlist(ch))

####4. Counting length of each sequence
len_seq<-data.frame(names=Data$Sequence,chr=apply(Data,2,nchar)[,2])

rep_char<-sapply(letters, function(x) x<-sum(x==unlist(strsplit(Data$Sequence,""))))


####6. kmers approach
#####6.1 2-mers approach
#General function for counting kmers
kmers.count <- function(dataset, colnum, kdata){
  x = 1
  while (x <= nrow(dataset)){
    for (k in kdata[,1]){
      dataset[x,k] <- str_count(dataset[x,colnum],k)
    }
    x = x + 1
  }
  return(dataset)
}

kmers2<-get.kmers(Data[,2],.k=2)
kmers.df<-Data
for (k in kmers2[,1]){
  kmers.df[,k]<-0
}
kmers.df<-kmers.count(kmers.df, 2, kmers2)
dim(kmers.df)


#####6.2 3-Mers approach
kmers3<-get.kmers(Data[,2],.k=3)
kmers.df3<-Data
for (k in kmers3[,1]){
  kmers.df3[,k]<-0
}
kmers.df3<-kmers.count(kmers.df3, 2, kmers3)
dim(kmers.df3)


#####6.3 4-Mers approach
kmers4<-get.kmers(Data[,2],.k=4)
kmers.df4<-Data
for (k in kmers4[,1]){
  kmers.df4[,k]<-0
}
kmers.df4<-kmers.count(kmers.df4, 2, kmers4)
dim(kmers.df4)


####7. Longest Common Subsequence
r <- nrow(Data)
names <- Data[,1]

LCSMatrix <- matrix(0,nrow = r, ncol= r, dimnames = list( c(names), c(names)))

compare <- function(){
  
  current_row = 1
  Comparing_row = 1
  LCS<-c()
  LengthLCS<-c()
  S1<-c()
  S2<-c()
  QualitySimilarityIndex<-c()
  
  while(current_row <= nrow(Data)){
    string1 = Data[current_row,2]
    
    while(Comparing_row <= nrow(Data)){
      string2 = Data[Comparing_row,2]
      
      l=LCS(unlist(strsplit(string1,"")),unlist(strsplit(string2,"")))
      S1<-append(S1,string1)
      S2<-append(S2,string2)
      
      LCS<-append(LCS,paste(l$LCS,collapse=""))
      
      LengthLCS<-append(LengthLCS,l$LLCS)
      
      QualitySimilarityIndex<-append(QualitySimilarityIndex,l$QSI)
      
      Comparing_row = Comparing_row+1
    }
    
    current_row = current_row+1
    Comparing_row = current_row
  }
  LCS.data<<-data.frame(S1,S2,LCS,LengthLCS,QualitySimilarityIndex)
}
compare()
LCS.data[,1]<-sapply(LCS.data[,1],as.character)
LCS.data[,2]<-sapply(LCS.data[,2],as.character)
LCS.data[,3]<-sapply(LCS.data[,3],as.character)


#Now merging all the features that we extracted
total<-cbind(bag,rep_char,len_seq$chr,bigramdf[c(-1,-2)],trigramdf[c(-1,-2)],gramdf4[c(-1,-2)], ch)


#The dimension of the features that we extracted are as follow
dim(total)


#Feature Selection
####1.	Zero-and Near Zero-Variance Predictors
nzv <- nearZeroVar(total)
filtered <- total[, -nzv]


#The dimensions of the dataset after removing near zero-variance predictors are
dim(filtered)


#####2. Removing Highly Correlated Predictors
descrCor <- cor(filtered[c(-1,-2)])
highlyCorDescr <- findCorrelation(descrCor, cutoff = .75)
filtered <- filtered[,-highlyCorDescr]


#In order to have class labels we areusing Levenshtein distance and generating hierarchical clustering model of the sequences.

#Finding distance (Levenshtein)
levDist  <- adist(as.list(filtered[,2]))
colnames(levDist) <- filtered[,1]
rownames(levDist) <- filtered[,1]

#hierarchical clustering method
hc <- hclust(as.dist(levDist))
plot(hc, main="Levenshtein on SSRs", cex = 0.6, col = "blue")
rect.hclust(hc,k=6)
filtered <-data.frame(filtered, as.factor(cutree(hc,k=6)))
names(filtered)[63] <- "Class"

filtered$Class <- as.factor(filtered$Class)


#####3. Boruta 
set.seed(111)
Boruta<-Boruta(Class~.,data = filtered, doTrace=2,maxRuns=500)
plot(Boruta,las=2,cex.axis=0.7)
plotImpHistory(Boruta)

Finding the non rejected formula
getNonRejectedFormula(Boruta)


#Finding confirmedFormula
getConfirmedFormula(Boruta)

#Tentative fix
TentativeRoughFix(Boruta)

#finding more mean importance
head(attStats(Boruta))


######selecting training and testing set

tsize <- floor(0.75 * nrow(filtered))
dt <- sample(x = nrow(filtered), size = tsize)

dt.train = filtered[dt,]
dt.test = filtered[-dt,]

train<-filtered[dt,c(-1,-2)]
test<-filtered[-dt,c(-1,-2)]

data.train = as.factor(dt.train$Class)
data.test = as.factor(dt.test$Class)


#Classification Models
###1. Naive Bayes Classification Model
m <- naiveBayes(as.matrix(dt.train), data.train)
p = predict(m, as.matrix(dt.test))
 
confusionMatrix(p, data.test)


###2. Random Forest Classification Model
set.seed(333)
rf60 <- randomForest(Class~.,data =train)

p<-predict(rf60,test)
confusionMatrix(p,test$Class)

####2.1 Random Forest on reduced data set
rf4<-randomForest(Class ~ + d + s + a + q + t + l + p + len_seq.chr + 
                    ss + qq + sq + as + sg + qe + sv + ls + dq + gd  + da + 
                    lp + ssq + sgq + sag + vlp + lpq + tssq + sgqq + sqsg + adss + 
                    svss + qsgq + gdqq + agdq + aggq + sagg + sgvs + gvss + svlp, data=train)

p1<-predict(rf4,test)
confusionMatrix(p1,test$Class)

####2.2 Random Forest on confirmed data set
rfc<-randomForest(Class ~ + d + s + q + t + l + p + len_seq.chr + ss + 
                    qq + sq + as + sg + qe + sv + ls + dq + gd  + da + lp + 
                    ssq + sgq + sag + vlp + tssq + sgqq + sqsg + adss + svss + 
                    qsgq + agdq + aggq + sagg + sgvs + gvss,data = train)


p2<-predict(rfc,test)
confusionMatrix(p2,test$Class)
