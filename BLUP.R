############******BLUP******************#############################
#### all effects are  taken as random
#### Blup(Entry,Rep,Block,Env,Traits,model=c("RCB","Lattice"),data)
#### if there not BLOCK then Block=NULL equals ENV
####  Traits = names of traits (one o more), it has be  numeric 

####################################################################
####################**Across Enviroment RCB **######################
####################      Example 1     ###########################

dat <- read.csv("dat_RCBD Data.csv",head=T)
names(dat)
dat$YLD <- as.numeric(dat$YLD)
dat$AD <- as.numeric(dat$AD)
dat$SD <- as.numeric(dat$SD)
dat$PH <- as.numeric(dat$PH)
dat$EH <- as.numeric(dat$EH)
dat$rEPH <- as.numeric(dat$rEPH)
trait <- c("ears" ,"len"  ,   "weight"  ,"yield")

Ex_1<- Blup_2(dat$Entry,dat$Replication,dat$BLK,dat$Location,trait,model = "RCB",dat)
Ex_1$VarianceComponets# <- variance componets and H2
Ex_1$Blups ## <- BLUPS

####################################################################
####################**One Enviroment  RCB **#######################
####################      Example 2    ###########################

dat <- read.csv("data_2.csv",header = T)
names(dat) <- c("geno","rep","ears" ,"len"  ,   "weight"  ,"yield" )
head(dat)
trait <- c("ears" ,"len"  ,   "weight"  ,"yield")
Ex_2<- Blup_2(dat$geno,dat$rep,Block = NULL,Env=NULL,trait,model = "RCB",dat)
Ex_2$VarianceComponets# <- variance componets and H2
Ex_2$Blups ## <- BLUPS


#################Scrip


Blup<-function(Entry,Rep,Block=NULL,Env= NULL,Traits,model=c("RCB","Lattice"),data){
  
  if (!is.null(Block)){
    Block <- as.factor(Block)
                            }  
  if (!is.null(Env)){
    Block <- as.factor(Env)
    Loc <- as.factor(Env) 
    nLoc <- length(levels(Loc))
                                 }  
  
  
  
  Entry<-as.factor(Entry)
  Rep <-as.factor(Rep)
  nrep <- length(levels(Rep))  
  nE <-length(levels(Entry))
  leng_traits<- length(trait)
  
  ############Data frame out
  VarComp <- data.frame()
  H2 <- data.frame()
  Blups<- data.frame(matrix(vector(),nE,1, dimnames=list(c(), c("Entry"))))
  Blups$Entry <- levels(Entry)
  ############Loop
  
  if (!is.null(Env)) {
    for (i in 1: leng_traits) {
      
      y <- data[,trait[i]]
      
      if (model=="Lattice"){
        fm<- lmer(y~(1|Entry)+(1|Rep:Loc)+(1|Block:Rep:Loc)+(1|Entry:Loc)+(1|Loc))
      }
      
      else if (model=="RCB"){
        fm<- lmer(y~(1|Entry)+(1|Rep:Loc)+(1|Entry:Loc)+(1|Loc))
      }
      
      blup <- coef(fm)$Entry
      colnames(blup) <- trait[i]
      Blups<-cbind(Blups,blup)
      
      
      ### h2
      
      vc <- VarCorr(fm,comp="Variance")
      VG <- vc$Entry[1]
      VGE <- vc$`Entry:Loc`[1]
      VP<- VG + VGE/nLoc + attr(vc, "sc")^2/(nrep*nLoc)
      h2 <-data.frame( h2=VG/VP)
      h2$trait <- trait[i]
      H2<- rbind(H2,h2)
      
      
      #### Variance Componets
      
      varComp<-as.data.frame(vc)
      drops_cc<-  c("Rep:Loc","Block:Rep:Loc")
      drops_c <- c("var1","var2","sdcor")
      varComp<-varComp[ , !(names(varComp) %in% drops_c )]
      varComp<-varComp[ !(varComp$grp %in% drops_cc), ]
      varComp$Trait<-trait[i]
      VarComp <- rbind(VarComp,varComp)
      
    } #end loop
    
    
    
    VarComp_Output <-reshape(VarComp, idvar = "Trait", timevar = "grp", direction = "wide")
    names(VarComp_Output) <- c("Trait","GxE","Env","Entry","Error")                             
    VarComp_Output$H2 <- H2$h2
    
  }##End across enviroment
  
#####################################################################
  if (is.null(Env)) {
    
    for (i in 1: leng_traits) {
      
      y <- data[,trait[i]]
      
      if(model=="RCB"){
        fm <- lmer(y ~ (1|Entry)+(1|Rep))
      }
      
      else if (model=="Lattice"){
        fm <- lmer(y~(1|Entry)+(1|Rep)+(1|Block:Rep))
      }
      
      blup <- coef(fm)$Entry
      colnames(blup) <- trait[i]
      Blups<-cbind(Blups,blup)
      
      ######h2
      vc <- VarCorr(fm,comp="Variance")
      VG <- vc$Entry[1]
      VP<- vc$Entry[1] + attr(vc, "sc")^2/nrep 
      h2 <-data.frame( h2=VG/VP)
      h2$trait <- trait[i]
      H2<- rbind(H2,h2)  
      
      #* Variance Components
      varComp<-as.data.frame(vc)
      drops_cc<-  c("Block:Rep")
      drops <- c("var1","var2","sdcor") 
      varComp<-varComp[ , !(names(varComp) %in% drops )]
      varComp<-varComp[ !(varComp$grp %in% drops_cc), ]
      varComp$Trait<-trait[i]
      VarComp <- rbind(VarComp,varComp)              
    }##end for loop for one enviroment
    
    VarComp_Output <-reshape(VarComp, idvar = "Trait", timevar = "grp", direction = "wide")
    VarComp_Output<- round(VarComp_Output[2:4],4)
    row.names(VarComp_Output) <-trait
    VarComp_Output <- cbind(VarComp_Output ,H2$h2)
    names(VarComp_Output) <- c("Genotype","Rep","Error","H2")
    
  } # end for one enviroment
  
  result <- list(VarianceComponets=VarComp_Output,Blups=Blups)
  
} ##end function




















