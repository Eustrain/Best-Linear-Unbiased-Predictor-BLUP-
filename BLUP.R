
########################################
#####Variance components and BLUP#############

library(lme4)
dat <- read.csv("data.csv",header = T)
names(dat) <- c("geno","rep","ears" ,"len"  ,   "weight"  ,"yield" )
head(dat)

trait <- c("ears" ,"len"  ,   "weight"  ,"yield" )


names(fw)

fw <- read.csv("dat.csv",header = T)
head(fw)
fw <- dplyr::filter(fw,env==1)
fw$PUDMZ <- as.numeric(fw$PUDMZ)
fw$ALTMZ <- as.numeric(fw$ALTMZ)
fw <- as.data.frame(fw)
names(fw)
trait <- c("ALTMZ")

y <- fw$ALTMZ

dd <- Blups(fw$geno,fw$rep,fw$BLOCK,trait,model = "Lattice",fw)
dd$VarianceComponets

Blups <- function(Entry,Rep,Block=NULL,Traits,model=c("RCB","Lattice"),data){


  
if (!is.null(Block)){
    Block <- as.factor(fw$BLOCK)
}  
  
geno<-as.factor(fw$entry)
rep <-as.factor(fw$rep)
nrep <- length(levels(rep))

###########################################  
nrep <- length(levels(rep))
ng <-length(levels(geno))
leng_traits<- length(trait)

###########################################



#########*Traits
VarComp <- data.frame()
H2 <- data.frame()
Blups<- data.frame(matrix(vector(),ng,1, dimnames=list(c(), c("Entry"))))
Blups$Entry <- levels(geno)

##########################################

for (i in 1: leng_traits) {
  
  y <- fw[,trait[i]]
  
  #if(model=="RCB"){
  #fm <- lmer(y ~ (1|geno)+(1|rep))
  #}
  
  
#else if (model=="Lattice"){
  fm <- lmer(y~(1|geno)+(1|rep)+(1|Block:rep))
#}
  
    #*BLUPS
  blup <- coef(fm)$geno
  colnames(blup) <- trait[1]
  Blups<-cbind(Blups,blup)
  #*H2
  vc <- VarCorr(fm,comp="Variance")
  VG <- vc$geno[1]
  VP<- vc$geno[1] + attr(vc, "sc")^2/nrep 
  h2 <-data.frame( h2=VG/VP)
  h2$trait <- trait[1]
  H2<- rbind(H2,h2)
  #* Variance Components
  varComp<-as.data.frame(vc)
  drops_cc<-  c("Block:rep")
  drops <- c("var1","var2","sdcor") 
  varComp<-varComp[ , !(names(varComp) %in% drops )]
  varComp<-varComp[ !(varComp$grp %in% drops_cc), ]
  varComp$Trait<-trait[1]
 VarComp <- rbind(VarComp,varComp)
  
}


VarComp_Output <-reshape(VarComp, idvar = "Trait", timevar = "grp", direction = "wide")
VarComp_Output<- round(VarComp_Output[2:4],4)
row.names(VarComp_Output) <-trait[1]
VarComp_Output <- cbind(VarComp_Output ,H2$h2)
names(VarComp_Output) <- c("Genotype","Rep","Error","H2")

result <- list(VarianceComponets=VarComp_Output,Blups=Blups)

return(result)

}




