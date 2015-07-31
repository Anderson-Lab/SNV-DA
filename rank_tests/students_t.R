

fishers = read.csv("~/repos//SNV-DA/fishers/SNVM.fishers_rank.csv", header=F, na.string="Na",stringsAsFactors=FALSE)



nonsyn <- read.csv("SNVM.unsup_total.csv", header=F, na.string="Na",stringsAsFactors=FALSE)
topNon <- read.delim("~/repos//SNV-DA/spls2/intronic_parallel_top_SNVs.txt",stringsAsFactors=FALSE,header=F)
topNon <- read.delim("~/repos//SNV-DA/spls2/nonsyn_parallel_top_SNVs.txt",stringsAsFactors=FALSE,header=F)
topNon <-  unlist(topNon[,1])
tmp <- topNon[1]
head(tmp)
topNon <- topNon[2:length(topNon)]
topNon[length(topNon)+1] = tmp


topData <- cbind.data.frame(topNon, c(1:length(topNon)))
colnames(topData) = c("SNV", "Rank")


SNVs <- nonsyn[,1:2]
nonsyn <- nonsyn[,c(-1,-2)]

ps = c()
for(i in 1:nrow(nonsyn)){
  print(i)
  test <- tryCatch(t.test(na.omit(as.numeric(nonsyn[i,1:11])), na.omit(as.numeric(nonsyn[i,12:20]))), error = function(e) {"NaN"})
  if (test != "NaN"){
    ps = c(ps, test$p.value)
  }
  else {
    ps = c(ps, test)
  }
}
test = cbind.data.frame(SNVs, ps)
head(test)
test2 = test[with(test, order(ps)),]
colnames(test2) = c("SNV", "type", "ps")
#write.table(test2, file="SNVM.students_t_rank.txt", sep="\t")
test2 = test2[which(grepl("nonsyn", test2$type)),]

colnames(fishers) = c("SNV", "p-value", "tab", "type")
test2 = fishers[which(grepl("nonsyn", fishers$type)),]

saved = topData

Rank1 = saved$Rank



topData = saved
loci = matrix(,nrow=length(topData$SNV), ncol=1)
for (i in 1:length(topData$SNV)){
 loci[i] = substr(toString(topData$SNV[i]), 1, nchar(toString(topData$SNV[i])) - 5) 
  
}

loci = topData$SNV

Rank1 = saved$Rank

Rank2 = c()
for (i in 1:length(test2$SNV)){
  found = F
  #print(iter)
  #iter = iter + 1
  if(test2$SNV[i] %in% loci){
    print("!")
  for(x in 1:length(topData$SNV)){
   if(unlist(test2$SNV)[i] == loci[x]){
    Rank2 = c(Rank2, topData$Rank[x])
    found = T
    topData = topData[c(-x),] 
    break
    }
   }
  }
   if(found == F){
    Rank2 = c(Rank2, 489)
    index = index + 1
     
   }
   
}
  
Rank1 = c(Rank1, rep(489, length(Rank2) - length(Rank1)))


intronsPLSRRank = Rank1
intronFishersRank = Rank2
intronStudentsRank = Rank2

newRank1 =  c(intronsPLSRRank, rep(489, length(intronFishersRank) - length(intronsPLSRRank)))
corTest <- cor.test(newRank1, intronFishersRank, method="spearman")
corTest
corTest <- cor.test(newRank1, intronFishersRank, method="pearson")
corTest
corTest <- cor.test(newRank1, intronFishersRank, method="kendall")
corTest

friedman.test(cbind(newRank1, intronFishersRank))

newRank1 =  c(intronsPLSRRank, rep(489, length(intronStudentsRank) - length(intronsPLSRRank)))
corTest <- cor.test(newRank1, intronStudentsRank, method="spearman")
corTest
corTest <- cor.test(newRank1, intronStudentsRank, method="pearson")
corTest
corTest <- cor.test(newRank1, intronStudentsRank, method="kendall")
corTest

friedman.test(cbind(newRank1, intronStudentsRank))


nonsynsPLSRRank = Rank1
nonsynFishersRank = Rank2
nonsynStudentsRank = Rank2

newRank1 =  c(nonsynsPLSRRank, rep(489, length(nonsynFishersRank) - length(nonsynsPLSRRank)))
corTest <- cor.test(newRank1, nonsynFishersRank, method="spearman")
corTest
corTest <- cor.test(newRank1, nonsynFishersRank, method="pearson")
corTest
corTest <- cor.test(newRank1, nonsynFishersRank, method="kendall")
corTest

friedman.test(cbind(newRank1, nonsynFishersRank))

newRank1 =  c(nonsynsPLSRRank, rep(489, length(nonsynStudentsRank) - length(nonsynsPLSRRank)))
corTest <- cor.test(newRank1, nonsynStudentsRank, method="spearman")
corTest
corTest <- cor.test(newRank1, nonsynStudentsRank, method="pearson")
corTest
corTest <- cor.test(newRank1, nonsynStudentsRank, method="kendall")
corTest

friedman.test(cbind(newRank1, nonsynStudentsRank))



