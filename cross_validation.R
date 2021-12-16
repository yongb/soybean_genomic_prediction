rm(list=ls())
##################### Prep ##################################################################
setwd("C:/Dropbox/Ph.D/GS/Data/Pheno")
library(rrBLUP)
library(kernlab)
library(BLR)
library(brnn)
library(randomForest)
cores = 4
install.packages("brnn")
load("SNP_TP.RData")
yield <- read.delim("SCN.txt", header=T)
lines <- as.character(yield$name)  # has to be character to extract the columns of dataframes!!!!!!!!!!!!!!!!!!!!!!!!
yieldlines <- YY[,lines]
row.names(yieldlines) <-YY$name
N <- yieldlines[M,]
rrSNP <- rr(N)  
Markers <- data.matrix(rrSNP[1:(dim(rrSNP)[1]),1:(dim(rrSNP)[2])],rownames.force = NA)
class(Markers) <- "numeric"
impute = A.mat(Markers, min.MAF=0.01, max.missing=0.5, impute.method="mean", n.core=4,return.imputed=T)
A = impute$A
Markers_impute=impute$imputed
tM <- t(Markers_impute)
X <- t(na.omit(tM))
n <- nrow(X); p <- ncol(X)
y <- yield[,2]

##################### Cross Validation Split ##################################################################
set.seed(1234)
folds<-5
sets<-rep(1:5,57)
rep = 10
COR.CV<-matrix(NA,rep,folds)
modellist=c("RR","BLR","BCP","RF","NN")
model=modellist[1]

for(r in 1:rep){
  sets<-sets[order(runif(n))]  
   for(fold in 1:folds) {
    tst<-which(sets==fold)
    XTRN<-X[-tst,] ; yTRN<-y[-tst]
    XTST<-X[tst,] ;  yTST<-y[tst]
    
     if (model == "RF") {
        cat("=================== Random Forest is running now! ===============") 
        penalty = 0
        if(penalty){
        pValues<-numeric()
        for(i in 1:p){
        fm<-lm(yTRN~XTRN[,i])
        pValues[i]<-summary(fm)$coef[2,4]
        print(paste('Fitting Marker ',i,'!',sep=''))
        }
        nMarkers<-100
       selSNPs<-order(pValues)[1:nMarkers]
       XTRN<-XTRN[,selSNPs]
       XTST<-XTST[,selSNPs]
        }
       NF<-randomForest(y=yTRN,x=XTRN,ntree=500,mtry=4)
        yHatNF<-predict(NF,XTST)
        COR.CV[r,fold] <- cor(yHatNF,yTST)
        cat("Accuracy is ", COR.CV[r,fold], " in ", r," Rep, and in ", fold, " fold! \n\n")
     }
        
  
  if (model == "RR") {
    cat("=================== Ridge regression is running now! =============== \n\n")
      rrBLUP <- mixed.solve(yTRN, Z=XTRN, K=NULL, SE = FALSE, return.Hinv=FALSE)
      SNP_Effect = as.matrix(rrBLUP$u)     
      yHatrr=((XTST %*% SNP_Effect)[,1] + rrBLUP$beta)    
      COR.CV[r,fold] <- cor(yHatrr,yTST)
     cat("Accuracy is ", COR.CV[r,fold], " in ", r," Rep, and in ", fold, " fold! \n\n")
    }
  
  if (model == "BLR") {
      cat("=================== Bayesian lasso is running now! ===============\n\n")
      nIter<-5000    
      burnIn<-1000     
      thin<-10
      w<-rep(1/nrow(A),folds) 
      yHatCV<-numeric()
      priorBL<-list(
        varE=list(df=3,S=2.5),
        varU=list(df=3,S=0.63),
        lambda = list(shape=0.52,rate=1e-5,value=20,type='random')
      )
        prefix<-paste('PM_BL','_fold_',fold,'_',sep='')
        fm<-BLR(y=yTRN,XL=XTRN,prior=priorBL,
                nIter=nIter,burnIn=burnIn,thin=thin)
        yHatCV[tst]<-fm$yHat[fm$tst]
        w[fold]<-w[fold]*length(fm$tst)
        COR.CV[r,fold]<-cor(yHatCV[tst],yTST)
        cat("Accuracy is ", COR.CV[r,fold], " in ", r," Rep, and in ", fold, " fold! \n\n")
      }
    }
    
  if (model == "NN") {
      cat("=================== Nueral Networks is running now! ===============\n\n")
      penalty = 0
      if(penalty){
        pValues<-numeric()
        for(i in 1:p){
          fm<-lm(yTRN~XTRN[,i])
          pValues[i]<-summary(fm)$coef[2,4]
          print(paste('Fitting Marker ',i,'!',sep=''))
        }
        nMarkers<-100
        selSNPs<-order(pValues)[1:nMarkers]
        XTRN<-XTRN[,selSNPs]
        XTST<-XTST[,selSNPs]
      }
        NN<-brnn(y=yTRN,X=XTRN,neurons=4, epochs=100, cores=cores)
        yHatNN<-predictions.nn(X=XTST,theta=NN$theta, neurons=4)
        COR.CV[r,fold] <- cor(yHatNN,yTST)
       
      cat("Accuracy is ", COR.CV[r,fold], " in ", r," Rep, and in ", fold, " fold! \n\n")
    }
  
    if (model == "BCP") {  
      cat("=================== Bayesian CPi is running now! ===============\n\n")
      numiter <- 5000    
      burnin <- 1000
      thin <- 10
    
      markerMeans = colMeans(X)
      X = t(t(X) - markerMeans)
      check <- X[,colMeans(X) == 0]
      todelete <- colnames(X[,colMeans(X) == 0])
      X<-X[,!colnames(X) %in% todelete]
      for(fold in 1:folds) {
        pi = 0.95; logPi = log(pi); logPiComp = log(1-pi); nuVar = 4; scaleRes = 1; nuRes = 4
        tst<-which(sets==fold)
        XTRN<-X[-tst,] ; yTRN<-y[-tst]
        XTST<-X[tst,] ;  yTST<-y[tst]  
        nmarkers = ncol(XTRN);nrecords = nrow(XTRN);size = ncol(XTRN)
        vara = 2.26; storePi = array(0.0,numiter); p = markerMeans/2.0;mean2pq = mean(2*p*(1-p))
        varEffects  = vara/(nmarkers*(1-pi)*mean2pq); scaleVar = varEffects*(nuVar-2)/nuVar
        b = array(0.0,size); meanb = b; mu = mean(y); var = array(0.0,nmarkers); ppa = array(0.0,nmarkers)
        meanMu = 0; piMean = 0.0; meanVar = 0.0  
        # adjust y
        ycorr = yTRN - mu 
        # mcmc sampling
        for (iter in 1:numiter){ 
          # sample vare
          vare = ( t(ycorr)%*%ycorr + nuRes*scaleRes )/rchisq(1,nrecords + nuRes)   
          # sample intercept
          ycorr = ycorr + mu
          rhs    = sum(ycorr)
          invLhs = 1.0/(nrecords)
          mean = rhs*invLhs                            
          mu = rnorm(1,mean,sqrt(abs(invLhs*vare)))                                                  
          ycorr = ycorr - mu
          meanMu = meanMu + mu  
          # sample delta and effect for each locus
          nLoci = 0
          for (locus in 1:nmarkers) {
            ycorr = ycorr + XTRN[,locus]*b[locus];
            rhs = t(XTRN[,locus])%*%ycorr;
            xpx = t(XTRN[,locus])%*%XTRN[,locus];
            v0  =  xpx*vare;
            v1  =  (xpx^2*varEffects + xpx*vare);
            logDelta0 = -0.5*(log(v0) + rhs^2/v0) + logPi;     
            logDelta1 = -0.5*(log(v1) + rhs^2/v1) + logPiComp;
            probDelta1 = 1.0/(1.0 + exp(logDelta0-logDelta1));
            u = runif(1);
            if (u < probDelta1) {
              nLoci = nLoci + 1
              lhs = xpx + vare/varEffects
              invLhs = 1.0/lhs
              mean = invLhs*rhs
              b[locus]= rnorm(1,mean,sqrt(invLhs*vare))
              ycorr = ycorr - XTRN[,locus]*b[locus]
              meanb[locus] = meanb[locus] + b[locus]
              ppa[locus] = ppa[locus] + 1
              var[locus] = varEffects
            } else {
              b[locus] = 0.0
              var[locus] = 0.0
            } #probDelta
          }   #locus
          # sample common variance 
          countLoci = sum(var!=0)
          varEffects = (scaleVar*nuVar + sum(b^2))/rchisq(1,nuVar+countLoci); #
          meanVar = meanVar + varEffects;    
          # sample Pi
          aa = nmarkers-countLoci + 1;
          bb = countLoci + 1;
          pi = rbeta(1, aa, bb);  # 
          storePi[iter] = pi;
          piMean = piMean + pi;
          logPi     = log(pi);
          logPiComp = log(1-pi);
          # tracking and monitoring  &(iter %% 100 == 0) )
          if ((iter > burnin)&(iter %% thin == 0) ) {   
            cat ("iteration ",iter," number of loci in model = ", nLoci," in fold ", fold, "!\n");
            meanbb  = meanb/(iter-burnin)
            yHat   = XTST %*% meanbb
            cat ("corr = ",cor(yTST,yHat), "\n\n")
          } 
        } # iter
        meanbb  = meanb/(numiter-burnin)
        yHat   = XTST %*% meanbb
        COR.CV[fold] <- cor(yTST,yHat)
        cat("Accuracy is ", COR.CV[fold], " in ", fold, " fold! \n\n")
      } # fold
      hist(storePi,prob=T,ylim=c(0,50),xlab="Pi",main="Posterior Distribution of Pi")
      COR.CV[fold+1]<-mean(COR.CV[1:fold]); COR.CV[fold+2]<-sd(COR.CV[1:fold])/sqrt(length(5))
      cat("The average accuracy is ", COR.CV[fold+1], "\n\n")
      write.csv(COR.CV, file="Cross_Validation_BayesC_seedqt.csv")
    }  
}

    
    if (model == "RKHS") {
      cat("=================== Reproducing Kernel is running now! ===============\n\n")
      load("RKHSw.rda")
      D<-as.matrix(dist(X,method="euclidean"))^2
      D<-D/mean(D)
      #h<-c(1e-2,.1,.4,.8,1.5,3,5)
      h<-c(3)
      for(fold in 1:folds){
        yNA<-y
        tst<-which(sets==fold)
        yNA[tst]<-NA
        ### FITS MODELS #################################
        PMSE<-numeric() ; VARE<-numeric(); VARU<-numeric() 
        pD<-numeric(); DIC<-numeric()
        fmList<-list()
        for(i in 1:length(h)){
          print(paste('Working with h=',h[i],sep=''))
          # COMPUTES THE KERNEL
          K<-exp(-h(i)*D)
          # FITS THE MODEL
          prefix<- paste(h[i], "_",sep="")
          fm<-RKHS(y=yNA,K=list(list(K=K,df0=5,S0=2)),
                   nIter=15000,burnIn=5000,df0=5,S0=2,saveAt=prefix)
          fmList[[i]]<-fm
          PMSE[i]<-mean((y[tst]-fm$yHat[tst])^2)
          VARE[i]<-fm$varE
          VARU[i]<-fm$K[[1]]$varU
          DIC[i]<-fm$fit$DIC
          pD[i]<-fm$fit$pD
        }
        #R2<-1-PMSE/mean((y[tst]-mean(y[-tst]))^2)
        COR.CV[fold] <- cor(y[tst],fm$yHat[tst])
        cat("Accuracy is ", COR.CV[fold], " in ", fold, " fold! \n\n")
      }
      ## PLOTS ############################### 
      plot(VARE~h,xlab="Bandwidth", ylab="Residual Variance",type="o",col=4)
      plot(I(VARE/VARU)~h,xlab="Bandwidth",
           ylab="variance ratio (noise/signal)",type="o",col=4)
      plot(pD~h,xlab="Bandwidth", ylab="pD",type="o",col=2)
      plot(DIC~h,xlab="Bandwidth", ylab="DIC",type="o",col=2)
      plot(R2~h,xlab="Bandwidth", ylab="R-squared",type="o",col=2)
    }
  }
}
COR.CV[fold+1]<-mean(COR.CV[1:fold]); COR.CV[fold+2]<-sd(COR.CV[1:fold])/sqrt(length(5))
cat("The average accuracy is ", COR.CV[fold+1], "\n\n")
write.csv(COR.CV, file="Cross_Validation_RR_yield100.csv")
 
############ SVM #########################################################
for(fold in 1:folds) {
  tst<-which(sets==fold)
  XTRN<-X[-tst,] ; yTRN<-y[-tst]
  XTST<-X[tst,] ;  yTST<-y[tst]
  svp <- ksvm(y=yTRN,x=XTRN,type="eps-svr",kernel='rbfdot',C=100,scaled=c())
  yHatSV<-predict(svp,XTST)
  COR.CV[fold] <- cor(yHatSV,yTST)
  cat("Accuracy is ", COR.CV[fold], " in ", fold, " fold! \n\n")
}
COR.CV[fold+1]<-mean(COR.CV[1:fold]); COR.CV[fold+2]<-sd(COR.CV[1:fold])/sqrt(length(5))
cat("The average accuracy is ", COR.CV[fold+1], "\n\n")
write.csv(COR.CV, file="Cross_Validation_SVM_oil.csv")

##################### SNP Effect Summary ####################################################

library(ggplot2)
library(plyr)
d <- density(GH$V8)
plot(d, main="SNP Effect Distribution",type="o",col="blue")
lines(d,type="o",pch=22, lty=2,col="red")
legend("topright", c("FD","TN"), cex=0.8, 
       col=c("blue","red"), pch=21:22, lty=1:2);
legend("topright", inset=.05, col=c("blue","red"), pch=21:22, lty=1:2)

polygon(d, col="red", border="blue") 
abs=abs(SNP_Effect)
CDF <-data.frame(cdf=(ecdf(abs))(abs)*length(abs), blup=abs)
qplot(blup, length(abs)-cdf, data=CDF, geom="point",color="Red", xlab="Marker Effects", 
      ylab="Marker Counts", main="Empirical distribution of marker effects for yiled", 
      cex.lab=2, cex.axis=20)
CDF_sorted <- CDF[with(CDF, order(-blup)), ]
SNP200 <- row.names(CDF_sorted[1:200,])
CDF$cdf=NULL
as.numeric(CDF$blup)
normalize(Protein)
names(CDF)<-"ST"
write.csv(SNP_Effect,file="SNP_Effect_Yield_RR2.csv")
rm(list=ls())
Yield <- read.csv("SNP_Effect_Yield_RR2.csv", header=T)
Y <- Yield[order(-Yield$V1),] 
Y200 <- Y[1:500,]
Protein <- read.csv("SNP_Effect_Protein.csv", header=T)
P <- Protein[order(-Protein$Protein),] 
P200 <- P[1:500,]
Oil <- read.csv("SNP_Effect_Oil.csv", header=T)
O <- Oil[order(-Oil$Oil),] 
O200 <- O[1:500,]
SeedWeight <- read.csv("SNP_Effect_SeedWeight.csv", header=T)
SW <- SeedWeight[order(-SeedWeight$SW),] 
SW200 <- SW[1:500,]
SeedQuality <- read.csv("SNP_Effect_SeedQuality.csv", header=T)
SQ<- SeedQuality[order(-SeedQuality$ST),] 
SQ200 <- SQ[1:500,]
YieldBL <- read.csv("SNP_Effect_Yield_BLR.csv", header=T)
PI <- merge(PI,SQ200, all=F, by="X")
write.csv(PI,file="SNP_Effect_Summary.csv")
PI$X=NULL
PI1 <- as.matrix(PI)
as.numeric(PI1)
P <- (na.omit(PI1))
colnames(P)[5]<- "SeedQuality"
normalize(P)
install.packages("som")
library(som)
Pt <- t(P)
Q <- normalize(Pt, byrow=TRUE)
Qt <- t(Q)
heatmap_yield <- heatmap(P, col = cm.colors(256), scale="column", margins=c(5,10), Rowv=NA, Colv=NA)
names(PI)
names(PI)[3] <- "BLR"
attach(PI)
qplot(SNP_Effect,P_Value, color = "red", cex.lab=200, cex.axis=20)


##################### Cross Validation with population subsets ##################################################################
accuracy <- NULL
cycles<-500
subset_n <- round(c(1:9)*0.1*n)
for (c in 1:9) { 
   for (r in 1:cycles) {    
    train=as.matrix(sample(1:n, subset_n[c]))
    test<-setdiff(1:n, train)
    pheno=yield[train,2]
    g=Marker[train,]
    pheno_valid=yield[test,2]
    gt=Marker[test,]
    pheno_answer<-mixed.solve(pheno, Z=g, K=NULL, SE = FALSE, return.Hinv=FALSE)
    PHE = pheno_answer$u
    e4 = as.matrix(PHE)
    pred_pheno_valid = gt %*% e4
    pred_pheno=(pred_pheno_valid[,1] + pheno_answer$beta)
    #pred_yield
    accuracy[r,c] <- cor(pred_pheno_valid, pheno_valid, use="complete")
    cat("This is cycle ", r, " with ", subset_n[c], " lines! \n\n")
  }
}
write.csv(e4, file="Marker_effect_seedqt.csv")
#write.csv(accuracy, file="accuracy_cross_validation_popn_seedquality2.csv")

##################### Cross Validation with marker subsets ##################################################################
subset_p <- round(c(1:9)*0.1*p)
for (c in 1:9) { 
  for (r in 1:cycles) {
    train=as.matrix(sample(1:n, 0.9*n))
    test<-setdiff(1:n, train)
    pheno=yield[train,2]
    submarker<-as.vector(sample(3:p,subset_p[c]))
    g=Marker[train,submarker]
    pheno_valid=yield[test,2]
    gt=Marker[test,submarker]
    pheno_answer<-mixed.solve(pheno, Z=g, K=NULL, SE = FALSE, return.Hinv=FALSE)
    PHE = pheno_answer$u
    e = as.matrix(PHE)
    pred_pheno_valid = gt %*% e
    pred_pheno=(pred_pheno_valid[,1] + pheno_answer$beta)
    #pred_yield
    accuracy[r,c] <- cor(pred_pheno_valid, pheno_valid, use="complete")
    cat("This is cycle ", r, " with ", subset_p[c], " SNPs! \n\n")
  }
}
write.csv(accuracy, file="accuracy_cross_validation_marker_seedquality.csv")


# Cumulative distribution of SNP blups
library(ggplot2)
library(plyr)
blup_cdf<-data.frame(cdf=(ecdf(e))(e)*1129, blup=e, abs_blup=abs(e))
qplot(blup, cdf, data=blup_cdf, geom="step", xlab="Marker Effect", ylab="Marker Counts",
      main="Empirical distribution of marker effects for ")
write.csv(blup_cdf,file="cdf output of snp blups abs.csv")

qplot(blup, 1129-cdf, data=blup_cdf, geom="point", color="Red", xlab="Marker Effect", 
      ylab="Marker Counts", main="Empirical distribution of marker effects for oil",
      cex.lab=20, cex.axis=20)

abs_blup=abs(e)
blup_cdf_abs<-data.frame(cdf=(ecdf(abs_blup))(abs_blup)*1129, blup=abs_blup)
qplot(blup, 1129-cdf, data=blup_cdf_abs, geom="point",color="Red", xlab="Marker Effects", 
      ylab="Marker Counts", main="Empirical distribution of marker effects for oil", 
      cex.lab=2, cex.axis=20)

blup_cdf_abs_sorted <- arrange(blup_cdf_abs,desc(blup))
blup20_abs <- blup_cdf_abs_sorted[1:50,]
qplot(blup, 1129-cdf, data=blup20_abs, geom="point", size=1, color="Red", xlab="", 
      ylab="", main="", cex.lab=2, cex.axis=20)
         

dat <- data.frame(dens = c(rnorm(100), rnorm(100, 10, 5))
                  , lines = rep(c("a", "b"), each = 100))


RR <- read.csv("Summary_RR.csv")
boxplot(formula= Accuracy ~ Trait, data=RR, notch=TRUE,
        col=(c("gold","darkgreen","darkblue","pink","black")),
        main="Accuracy among traits", xlab="Traits") 
library(doBy)
points(x=1:5, y=summaryBy(formula= Accuracy ~ Trait, data=RR, FUN=mean, na.rm=T)$dat.mean, 
       pch=8)
}


me <- maker_effect_summary
densityplot(~ME, me, groups = Trait) 
                        auto.key = list(columns + xlab = "Marker effect")
            library(sm)
            attach(me)
            
            # create value labels
            trait.f <- factor(Trait, levels= c("Yield", "Protein", "Oil"),
                            labels = c("Yield", "Protein", "Oil"))
            
            # plot densities
            sm.density.compare(ME, Trait, xlab="Marker effect")
            title(main="Marker effect among traits")
            # add legend via mouse click
            colfill<-c(2:(2+length(levels(trait.f))))
            legend(locator(1), levels(trait.f), fill=colfill) 
            library(sm)
            sm.density.compare(bayes, newdata)
            # Add a legend (the color numbers start from 2 and go up)
            legend("topright", levels(data$cond), fill=2+(0:nlevels(data$cond)))