# Parameters
nmarkers    = 2000;   #number of markers
numiter     = 2000;   #number of iterations
numMHIter   = 100;    #use 1 for Bayes A 
pi          = 0.95;   #Change this to run Bayes B rather than Bayes A
vara        = 1.0; 

cat("Rohan Fernando's implementation of Bayes B (or A if pi=0)\n")
# input data
data     = matrix(scan("trainData.out0"),ncol=nmarkers+2,byrow=TRUE);
nrecords = dim(data)[1];
startMarker = 1901;
x = cbind(1,data[,startMarker:nmarkers]);
y = data[,nmarkers+1];
a =  data[,nmarkers+2];

# inital values
nmarkers = nmarkers - startMarker + 1;
mean2pq     = 0.5;
scaleb  = 0.5*vara/(nmarkers*(1-pi)*mean2pq);  

b = array(0.0,nmarkers+1);
meanb = b;
b[1] = mean(y);
var  = array(0.0,nmarkers);
ppa  = array(0.0,nmarkers);

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
	
# sample variance and effect for each locus
    nLoci = 0;
	for (locus in 1:nmarkers){
		ycorr = ycorr + x[,1+locus]*b[1+locus];
		rhs = t(x[,1+locus])%*%ycorr;
		totalSS = sum(ycorr^2)/vare;
		xpx = t(x[,1+locus])%*%x[,1+locus];
		v1  =  (xpx^2*var[locus] + xpx*vare);           # slide 47
		v2  =  xpx*vare;
		logDataNullModel = -0.5*(log(v2) + rhs^2/v2);   # slide 47
		if (var[locus] > 0.0){
			logDataOld   = -0.5*(log(v1) + rhs^2/v1);
		}
		else {
			logDataOld = logDataNullModel;            
		}
		for (mhiter in 1:numMHIter){
			u = runif(1);
			varCandidate = 0;
			if (u > pi){
				varCandidate = scaleb*4/rchisq(1,4);
			}
			if (varCandidate > 0.0){
				v1  =  (xpx^2*varCandidate + xpx*vare);
				logDataNew =  -0.5*(log(v1) + rhs^2/v1);  
			}
			else{
				logDataNew = logDataNullModel;
			}
			acceptProb = exp(logDataNew-logDataOld); # slide 45
			u = runif(1);
			if(u <acceptProb) {                      
				var[locus] = varCandidate;
				logDataOld = logDataNew;
			}
		}
		if(var[locus]) {
			nLoci = nLoci + 1;
			lhs = xpx/vare + 1.0/var[locus];
			invLhs = 1.0/lhs;
			mean = invLhs*rhs/vare;
			b[1+locus]= rnorm(1,mean,sqrt(invLhs));          
			ycorr = ycorr - x[,1+locus]*b[1+locus];
			meanb[1+locus] = meanb[1+locus] + b[1+locus];
			ppa[locus] = ppa[locus] + 1;
		}
		else b[1+locus] = 0.0;
	}
	if ((iter %% 100)==0) cat ("iteration ",iter," number of loci in model = ", nLoci,"\n");
}
meanb = meanb/numiter;
ppa   = ppa/numiter;
aHat  = x %*% meanb;
corr = cor(a,aHat);
cat ("corr = ",corr, "\n");

nmarkers=2000
data     = matrix(scan("testData.out0"),ncol=nmarkers+2,byrow=TRUE);
nrecords = dim(data)[1];
x = cbind(1,data[,startMarker:nmarkers]);
y = data[,nmarkers+1];
a =  data[,nmarkers+2];

aHat   = x %*% meanb;
corr = cor(a,aHat);
cat ("corr = ",corr, "\n");



