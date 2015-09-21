################################################
################################################
# Part 1: Set up the environment
library(qtl)
cross<-read.cross(format="csv",
                 file="LerxCvix.csv",
                 na.strings="NA",
                 genotypes=c("A","B"))
class(cross)[1]<-"riself"
cross<-calc.genoprob(cross, step=2, error.prob=0.0001, map.function="kosambi")
save(cross, file="lerxcviCross.RData")
########################






################################################
################################################
# Part 2: Initial QTL analysis of all traits
rm(list=ls())
load("lerxcviCross.RData")
phenames(cross)

########################
# Part 2.1: make a vector of phenotypes you are interested in
phes<-phenames(cross)[-1]

########################
# Part 2.2: plot all traits - using a for loop
plot(scanone(cross, pheno.col=1), type="n", ylim=c(0,30))
title("scanone of all phenotypes")
count<-1
for(i in phes){
  #scanone
  s1<-scanone(cross, pheno.col=i,  method="hk", model="normal")
  plot(s1,col=count, add=T)
  #permutations
  perms<-scanone(cross, pheno.col=i,  method="hk", model="normal", n.perm=100)
  add.threshold(s1,perms=perms, alpha=0.05, col=count, lty=2)
  count<-count+1
#   summary - if you want
#   cat("result for",i,":\n")
#   print(summary(s1, perms=perms, alpha=.05))
}
legend("topright", phes, col=1:6, lty=1)

########################
# Part 2.3: alternatively, all traits using internal R/qtl
s1s<-scanone(cross, pheno.col=phes,  method="hk", model="normal")
# R/qtl can only plot 3 traits at once
plot(s1s, lodcolumn=1:3, ylim=c(0,30)) # 1st 3
plot(s1s, lodcolumn=4:6, add=T, col=c(4:6)) # 2nd 3
legend("topright", phes, col=1:6, lty=1)
########################





################################################
################################################
# part 3: multiple QTL Modeling - single trait
phe<-"SD.FT"
form1<- "y ~ Q1"
scan1<-scanone(cross, pheno.col=phe,  method="hk", model="normal")
print(maxS1<-max(scan1))
scan1chr<-as.numeric(as.character(maxS1$chr))

########################
# part 3.1: make a QTL model that fits the marker under the highest QTL peak
?makeqtl
qtl1<-makeqtl(cross, chr=scan1chr, pos=maxS1$pos, what="prob")

#refine the position of this QTL
?refineqtl
summary(ref1 <- refineqtl(cross, pheno.col = phe, formula = form1, qtl=qtl1,
                          method="hk", model="normal"))
#compare to original - 
print(maxS1) #no difference

########################
# part 3.2 search for another QTL controlling for variation at this position
?addqtl
form2 <- "y ~ Q1 + Q2"
scan2<-addqtl(cross, pheno.col = phe, qtl = ref1, formula= form2, 
              method="hk", model="normal")

# plot result
plot(scan2,col="red", main="multi QTL scan comparisons")
plot(scan1, add=T, col="black")

maxS2<-max(scan2)
scan2chr<-as.numeric(as.character(maxS2$chr))

# make new qtl model
qtl2<-makeqtl(cross, 
              chr=c(scan1chr,scan2chr),
              pos=c(maxS1$pos,maxS2$pos), what="prob")

# refine again
summary(ref2 <- refineqtl(cross, pheno.col = phe, formula = form2, qtl=qtl2,
                          method="hk", model="normal"))

########################
# part 3.3: repeat for 3rd and 4th QTL
form3 <- "y ~ Q1 + Q2 + Q3"
scan3<-addqtl(cross, pheno.col = phe, qtl = ref2, formula= form3, 
              method="hk", model="normal")
maxS3<-max(scan3)
scan3chr<-as.numeric(as.character(maxS3$chr))
qtl3<-makeqtl(cross, 
              chr=c(ref2$chr,scan3chr),
              pos=c(ref2$pos,maxS3$pos), what="prob")
summary(ref3 <- refineqtl(cross, pheno.col = phe, formula = form3, qtl=qtl3,
                          method="hk", model="normal"))

form4 <- "y ~ Q1 + Q2 + Q3 + Q4"
scan4<-addqtl(cross, pheno.col = phe, qtl = ref3, formula= form4, 
              method="hk", model="normal")
maxS4<-max(scan4)
scan4chr<-as.numeric(as.character(maxS4$chr))
qtl4<-makeqtl(cross, 
              chr=c(ref3$chr,scan4chr),
              pos=c(ref3$pos,maxS4$pos), what="prob")
summary(ref4 <- refineqtl(cross, pheno.col = phe, formula = form4, qtl=qtl4,
                          method="hk", model="normal"))

########################
# part 3.4: plot result
plot(scan3, add=T, col="blue")
plot(scan4, add=T, col="green")
perms<-scanone(cross, pheno.col=phe,  method="hk", model="normal", n.perm=100)
add.threshold(scan1,perms=perms, alpha=0.05, lty=2)
legend(100,c("initial scanone", "search for 2nd QTL", "3rd", "4th"), 
       col=c("black","red","blue","green"), lty=1)

########################
# part 3.5: Fit a model that looks appropriate - drop the 4th QTL
?fitqtl
summary(fit3<-fitqtl(cross, pheno.col=phe, qtl=qtl3, formula=form3, 
                     method="hk", model="normal",
                     dropone=TRUE, get.ests=TRUE))

########################
# part 3.6: Look for epistasis
?addint
?scantwo
plot(twoway<-scantwo(cross, pheno.col=phe, 
                     method="hk", model="normal"))
epi3<-addint(cross, pheno.col=phe,
             qtl=qtl3, formula=form3, 
             method="hk", model="normal")
plot(twoway<-scantwo(cross, pheno.col=phe, 
                     method="hk", model="normal"))

########################
# part 3.7: fit final model, report statistics, make interaction plots
finalFormula<- "y ~ Q1 + Q2 + Q3 + Q2:Q3"
summary(finalRef <- refineqtl(cross, pheno.col = phe, formula = finalFormula, qtl=qtl3,
                          method="hk", model="normal"))
finalModel<-finalRef
summary(fitqtl(cross, pheno.col=phe,
               qtl=finalModel, formula=finalFormula, 
               method="hk", model="normal",
               dropone=TRUE, get.ests=TRUE))
?effectplot
?find.marker
mars<-find.marker(cross, chr=finalModel$chr, pos=finalModel$pos)
effectplot(cross, pheno.col=phe,
           mname1=mars[2], mname2=mars[3], add.legend=F)
legend("topleft", title=mars[2], col=c("red","blue"), c("AA","BB"), lty=1, pch=1)

?plotLodProfile
plot(finalModel)
plotLodProfile(finalRef, ylim=c(0,60), showallchr=T, 
               main=paste("LOD Profile of 3-QTL model for",phe))
########################






################################################
################################################
# Part 4: Automated stepwise model selection 

########################
# part 4.1: run scantwo permutations (these would normally be done on TACC)
s2perms<-scantwo(cross, pheno.col=phe,  method="hk", model="normal", n.perm=20)

########################
# part 4.2: calculate penalties
?calc.penalties
pens<-calc.penalties(s2perms)
print(pens)

########################
# part 4.3: run stepwise model selection
?stepwiseqtl
stepout<-stepwiseqtl(cross, pheno.col=phe,  
                     method="hk", model="normal", penalties=pens, max.qtl=6,
                     additive.only=FALSE, refine.locations=TRUE, 
                     keeptrace=TRUE, keeplodprofile=TRUE)

########################
# part 4.4: plot results
plot(stepout)
plotModel(stepout)
thetrace <- attr(stepout, "trace")
par(mfrow=c(4,3), mar=c(2, 4, 2, 2))
for(i in seq(along=thetrace))
  plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

########################
# part 4.5: fit results via anova
summary(fitqtl(cross, qtl=stepout, formula=formula(stepout), pheno.col=phe,
               method="hk", model="normal", get.ests=T, dropone=T))


