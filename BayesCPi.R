# Parameters
nmarkers    = 2000;    #number of markers
numiter     = 2000;    #number of iterations
pi          = 0.5;
vara        = 1.0;
logPi     = log(pi);
logPiComp = log(1-pi);
mean2pq     = 0.5;
nua       = 4;

cat("Rohan Fernando's implementation of Bayes CPi\n")
# input data
data     = matrix(scan("trainData.out0"),ncol=nmarkers+2,byrow=TRUE);
nrecords = dim(data)[1];
startMarker = 1901;
x = cbind(1,data[,startMarker:nmarkers]);
y = data[,nmarkers+1];
a =  data[,nmarkers+2];
storePi = array(0.0,numiter);
# inital values
nmarkers = nmarkers - startMarker + 1;
varEffects  = vara/(nmarkers*(1-pi)*mean2pq);
scalec      = varEffects*(nua-2)/nua; 

b = array(0.0,nmarkers+1);
meanb   = b;
b[1]   = mean(y);
var    = array(0.0,nmarkers);
ppa    = array(0.0,nmarkers); 
piMean = 0.0;
meanVar = 0.0;

# adjust y
ycorr = y - x%*%b;

# mcmc sampling
for (iter in 1:numiter){
	
# sample vare
	vare = ( t(ycorr)%*%ycorr )/rchisq(1,nrecords + 3);
	
# sample intercept
	ycorr = ycorr + x[,1]*b[1];
	rhs    = sum(ycorr)/vare;
	invLhs = 1.0/(nrecords/vare);
	mean = rhs*invLhs;
	b[1] = rnorm(1,mean,sqrt(invLhs));
	ycorr = ycorr - x[,1]*b[1];
	meanb[1] = meanb[1] + b[1];
	
# sample delta (slide 48)  and effect for each locus
    nLoci = 0;
	for (locus in 1:nmarkers){
		ycorr = ycorr + x[,1+locus]*b[1+locus];
		rhs = t(x[,1+locus])%*%ycorr;
		xpx = t(x[,1+locus])%*%x[,1+locus];
		v0  =  xpx*vare;
		v1  =  (xpx^2*varEffects + xpx*vare);
		logDelta0 = -0.5*(log(v0) + rhs^2/v0) + logPi;        
		logDelta1 = -0.5*(log(v1) + rhs^2/v1) + logPiComp;
		probDelta1 = 1.0/(1.0 + exp(logDelta0-logDelta1));
		u = runif(1);
		if(u < probDelta1) {
			nLoci = nLoci + 1;
			lhs = xpx/vare + 1.0/varEffects;
			invLhs = 1.0/lhs;
			mean = invLhs*rhs/vare;
			b[1+locus]= rnorm(1,mean,sqrt(invLhs));
			ycorr = ycorr - x[,1+locus]*b[1+locus];
			meanb[1+locus] = meanb[1+locus] + b[1+locus];
			ppa[locus] = ppa[locus] + 1;
			var[locus] = varEffects;
		}
		else {
			b[1+locus] = 0.0;
			var[locus] = 0.0;
		}
		
	}
	if(iter %% 100 == 0) cat ("iteration ",iter," number of loci in model = ", nLoci,"\n");
# sample common variance 
	countLoci = 0;
	sum = 0.0;
	for (locus in 1:nmarkers){
		if(var[locus]>0.0){
			countLoci = countLoci + 1;
			sum = sum + b[1+locus]^2;
		}
	}
	varEffects = (scalec*nua + sum)/rchisq(1,nua+countLoci); # slide 50
	meanVar = meanVar + varEffects;
# sample Pi
	aa = nmarkers-countLoci + 1;
	bb = countLoci + 1;
	pi = rbeta(1, aa, bb);  # slide 52
	storePi[iter] = pi;
	piMean = piMean + pi;
	scalec =  (nua-2)/nua*vara/((1-pi)*nmarkers*mean2pq)
	logPi     = log(pi);
	logPiComp = log(1-pi);
}
meanb  = meanb/numiter;
ppa    = ppa/numiter;
piMean = piMean/numiter;
meanVar = meanVar/numiter;
aHat   = x %*% meanb;
corr = cor(a,aHat);
cat ("corr = ",corr, "\n");

nmarkers = 2000;

data     = matrix(scan("testData.out0"),ncol=nmarkers+2,byrow=TRUE);
nrecords = dim(data)[1];
x = cbind(1,data[,startMarker:nmarkers]);
y = data[,nmarkers+1];
a =  data[,nmarkers+2];

aHat   = x %*% meanb;
corr = cor(a,aHat);
cat ("corr = ",corr, "\n");


