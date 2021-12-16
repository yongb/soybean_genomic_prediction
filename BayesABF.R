# slides 51

# Parameters
nmarkers    = 2000;   #number of markers
numiter     = 2000;    #number of iterations
numMHIter   = 10;
pi          = 0.9;
vara        = 1.0; 

logPi     = log(pi);
logPiComp = log(1-pi);

cat("Rohan Fernando's implementation of Bayes B Fast uses new MH proposal(or A if pi=0)\n")
# input data
data     = matrix(scan("trainData.out0"),ncol=nmarkers+2,byrow=TRUE);
nrecords = dim(data)[1];
startMarker = 1901;
x = cbind(1,data[,startMarker:nmarkers]);
y = data[,nmarkers+1];
a =  data[,nmarkers+2];

#x = matrix(scan("Rgenos.012"),ncol=453,byrow=TRUE);
#nrecords = dim(x)[1];
#x = cbind(1,x);
#y = matrix(scan("Rphenos"),ncol=1,byrow=TRUE);


# inital values
nmarkers = nmarkers - startMarker + 1;
#nmarkers    = 453;
mean2pq     = 0.2072;
scaleb   =  0.5*vara/(nmarkers*(1-pi)*mean2pq);
consta   = 2.0*log(scaleb*2.0);

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
		xpx = t(x[,1+locus])%*%x[,1+locus];
		v1  =  (xpx^2*var[locus] + xpx*vare);
		v2  =  xpx*vare;
		logDataNullModel = -0.5*(log(v2) + rhs^2/v2);
		if (var[locus] > 0.0){
			logDataOld  = -0.5*(log(v1) + rhs^2/v1);
			sv = scaleb;
			logPriorOld = consta -3.0*log(var[locus]) - sv*4/(2.0*var[locus]) + logPiComp;
		}
		else {
			logDataOld  = logDataNullModel;
			logPriorOld = logPi;
		}
		logPosteriorOld = logDataOld + logPriorOld;
		for (mhiter in 1:numMHIter){
			u = runif(1);
			varCandidate = 0;
			if (u > 0.5){
				if (var[locus] > 0){
					varCandidate = var[locus]*2 /rchisq(1,4);
				}
				else {
					varCandidate = scaleb*4/rchisq(1,4);
				}
			}
			if (varCandidate > 0.0){
				v1  =  (xpx^2*varCandidate + xpx*vare);
				logDataNew =  -0.5*(log(v1) + rhs^2/v1);
				sv = scaleb;
				logPriorNew = consta -3.0*log(varCandidate) - sv*4/(2.0*varCandidate) + logPiComp;
				logPosteriorNew = logDataNew + logPriorNew;
				if(var[locus]>0){
					sv = varCandidate*0.5;
					logProposalOld = 2.0*log(sv*2.0) -3.0*log(var[locus])   - sv*4/(2.0*var[locus]);
					sv = var[locus]*0.5;
					logProposalNew = 2.0*log(sv*2.0) -3.0*log(varCandidate) - sv*4/(2.0*varCandidate);				
				}
				else {
					logProposalOld = 0.0;
					sv = scaleb;
					logProposalNew = consta -3.0*log(varCandidate) - sv*4/(2.0*varCandidate);	
				}
			}
			else {
				logDataNew = logDataNullModel;
				logPriorNew = logPi;
				logPosteriorNew = logDataNew + logPriorNew;
				if (var[locus]>0){
					sv = scaleb;
					logProposalOld = consta -3.0*log(var[locus]) - sv*4/(2.0*var[locus]);
					logProposalNew = 0.0;
				}
				else {
					logProposalOld = 0.0;
					logProposalNew = 0.0;
				}
			}
			acceptProb = (exp(logPosteriorNew+logProposalOld-logPosteriorOld-logProposalNew))[1,1];
			u = runif(1);
#cat(locus,var[locus],varCandidate,u,acceptProb,"\n");
			if(u < acceptProb) {
				var[locus]      = varCandidate;
				logPosteriorOld = logPosteriorNew;			
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


