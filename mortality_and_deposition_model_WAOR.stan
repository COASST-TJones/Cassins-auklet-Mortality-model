// Carcass deposition model
// Estimates daily deposition rate via a spatio-temporal
// function with spatial and temporal double-logistic functions
// Carcass counts are related to deposition rates via
// survey timing, and carcass persistence and detection
// Detection: p_val
// Persistence: psi 

// Split into 3 events
// Event 1: October central
// Event 2: December central
// Event 3: December northern

functions {
  // Function for creating a vector of persistence rates
	// according to 1:ndays gap
	// Modelled as psi^gap
	vector Survgen(real psival, int ndays);
	vector Survgen(real psival, int ndays){
		vector[ndays] psivec;
		for(i in 1:ndays){
			psivec[i] = psival^(ndays-i);
		}
		return(psivec);
	}
	
	// Function for extracting a specific row vector
	// From the deposition matrix
	vector AGrab(matrix A, int strt, int ndays, int rownumb);
	vector AGrab(matrix A, int strt, int ndays, int rownumb){
		vector[ndays] Grab;
		for(i in 1:ndays){
			Grab[i] = A[rownumb,strt+i];
		}
		return(Grab);
	}
	
	// Function for evaluating the dot product of two vectors
	// mostly a convenience function (still need this?)
	real dopro(vector A, vector B);
	real dopro(vector A, vector B){
		real outdp;
		outdp = dot_product(A,B);
		return(outdp);
	}
	
	// Function for evaluating the double logistic function
	// For a given time/location (D) with parameters
	// St - start location/time, En - end location/time
	// B_St - Rate at start, B_En - Rate at end
	real DBLlog(real D, real St, real En, real B_St, real B_En);
	real DBLlog(real D, real St, real En, real B_St, real B_En){
		real D_out;
		D_out = 1/((1 + exp(-(D-St)/B_St)) * (1 + exp(-(En-D)/B_En)));
		return(D_out);
	}
	
	// Function for evaluating the single logistic function
	// For a given time/location (D) with parameters
	// St - start location/time, B_St - Rate at start
	real SINlog(real D, real St, real B_St);
	real SINlog(real D, real St, real B_St){
		real D_out;
		D_out = 1/(1 + exp(-(D-St)/B_St));
		return(D_out);
	}
}

data {
	int<lower=1> Nbeach;          	// number of beaches
	int<lower=1> Nday;				      // Number of days
	int<lower=1> NAll;              // number of segments (for scaling)
	real Allspace[NAll];            // all beach segments (for scaling)
	real day[Nday];				          // list of days (for scaling)
	int nsurv[Nbeach];				      // Number of surveys on each beach
	vector[Nbeach] BeachMid;        // surveyed segment midpoints
	vector[Nbeach] BeachLength;     // surveyed segment lengths
	int<lower=1> nsurvmax;          // maximum number of surveys
	int daynum[Nbeach,nsurvmax];    // matrix of survey dates (beaches by row)
	real<lower=0> p_val;            // detection rate
	vector[Nbeach] psi;             // persistence rate
	int C[Nbeach, nsurvmax];        // observed counts
	matrix[Nday,213] BigDMat;       // deposition matrix
}

parameters {
	real<lower=0> D0;				               // background deposition
	real<upper=80> T1s;	                   // Start time of Event 1
	real<lower=T1s> T1e;                   // End time of Event 1
	real<lower=0> B_T1e;                   // Rate at start of Event 1
	real<lower=0> B_T1s;                   // Rate at end of Ev 1
	real<lower=0> M1;                      // Magnitude of Ev 1
	real<lower=100> T2s;                   // Start time of Ev 2	
	real<lower=T2s, upper=180> T2e;        // End time of Ev 2
	real<lower=0> B_T2e;                   // Rate at end of Ev 2
	real<lower=0> B_T2s;                   // Rate at start of Ev 2
	real<lower=0> M2;                      // Magnitude of Ev 2
	real<lower=300, upper=900> S1s;	       // Southern extent Ev 1
	real<lower=S1s> S1e;                   // Northern extent Ev 1
	real<lower=0> B_S1e;                   // Rate at northern limit Ev 1
	real<lower=0> B_S1s;                   // Rate at southern limit Ev 1
	real<lower=1050, upper=1150> S3s;      // Southern Extent Ev 3
	real<lower=0> B_S3s;                   // Rate at southern limit Ev 3
	real<lower=0, upper=2> M3;             // Magnitude of Ev 3
	real<lower=300, upper=900> S2s;	       // Southern limit of Ev 2
	real<lower=S2s, upper=S3s> S2e;        // Northern limit of Ev 2
	real<lower=0> B_S2e;                   // Rate at northern limit Ev 2
	real<lower=0> B_S2s;                   // Rate at southern limit Ev 2
	real<lower=0> nu;	                     // Dispersion param			
}

model{
	// Model components
	vector[Nday] Timeprof1;		             // Temporal profile Ev1 (scaled max=1)	
	vector[Nday] Timeprof2;	               // Temporal profile Ev2 (scaled max=1)
  vector[Nday] TPmag1;	                 // DBL (for scaling) Ev 1	time	
	vector[Nday] TPmag2;                   // DBL (for scaling) Ev 2 time
	vector[NAll] SPmag1;                   // DBL (for scaling) Ev 1 space
	vector[NAll] SPmag2;                   // DBL (for scaling) Ev 2 space
	vector[NAll] SPmag3;                   // DBL (for scaling) Ev 3 space
	vector[Nbeach] Spaceprof1;				     // Spatial profile Ev1 (scaled max=1)
	vector[Nbeach] Spaceprof2;				     // Spatial profile Ev2 (scaled max=1)
	vector[Nbeach] Spaceprof3;             // Spatial profile Ev3 (scaled max=1)
	matrix[Nbeach,Nday] ABD;			         // Beach and day-specific deposition
	matrix[Nbeach,nsurvmax] Sum_AS;		     // Cumulative deposition (intermediate step)
	matrix[Nbeach,nsurvmax] C_hat;         // Expected counts
	real TPscale1;                         // The following are all the max of the unscaled profiles                
	real TPscale2;
	real SPscale1;
	real SPscale2;
	real SPscale3;
	matrix[Nday,213] MortMat1;             // Mortality matrix
	matrix[Nday,213] MortMat2;
	vector[Nday] Temp1;
	vector[Nday] Temp2;
	
	// Priors
	D0 ~ normal(0,10);
	T1s ~ normal(75,30);	
	T1e ~ normal(100,30);
	B_T1e ~ normal(0,100);
	B_T1s ~ normal(0,100);
	
	T2s ~ normal(110,30);	
	T2e ~ normal(150,50);
	B_T2e ~ normal(0,100);
	B_T2s ~ normal(0,100);
	
	S1s ~ normal(500,300);	
	S1e ~ normal(900,300);
	B_S1e ~ normal(0,100);
	B_S1s ~ normal(0,100);
	
	S2s ~ normal(650,300);	
	S2e ~ normal(950,300);
	B_S2e ~ normal(0,100);
	B_S2s ~ normal(0,100);
		
	M1 ~ normal(0,100);
	M2 ~ normal(0,100);
	
	S3s ~ normal(1100,200);
	B_S3s ~ normal(0,100);
	M3 ~ normal(0,100);
	
	nu ~ normal(0,100);					
	
	// Scaling
	// Evaluate the double-logistic functions across all time to find max for scaling
	// Scaling is required so that event magnitudes (M1, M2, M3) are properly defined
	for(i in 1:Nday){
		TPmag1[i] = DBLlog(day[i], T1s, T1e, B_T1s, B_T1e);
		TPmag2[i] = DBLlog(day[i], T2s, T2e, B_T2s, B_T2e);
	}
	// Similar scaling eval for spatial functions
	for(i in 1:NAll){
		SPmag1[i] = DBLlog(Allspace[i], S1s, S1e, B_S1s, B_S1e);
		SPmag2[i] = DBLlog(Allspace[i], S2s, S2e, B_S2s, B_S2e);
		SPmag3[i] = SINlog(Allspace[i], S3s, B_S3s);
	}	
	// Find maximum values for scaling	
	TPscale1 = max(TPmag1);
	TPscale2 = max(TPmag2);
	SPscale1 = max(SPmag1);
	SPscale2 = max(SPmag2);
	SPscale3 = max(SPmag3);
	
	// Time profile of mortality = the relative number of birds dying on day i
	// These are scaled
	for(i in 1:Nday){
		Timeprof1[i] = M1*DBLlog(day[i], T1s, T1e, B_T1s, B_T1e)/TPscale1;
		Timeprof2[i] = M2*DBLlog(day[i], T2s, T2e, B_T2s, B_T2e)/TPscale2;
	}
	// Spatial profiles (scaled) same as above
	for(i in 1:Nbeach){
		Spaceprof1[i] = DBLlog(BeachMid[i], S1s, S1e, B_S1s, B_S1e)/SPscale1;
		Spaceprof2[i] = DBLlog(BeachMid[i], S2s, S2e, B_S2s, B_S2e)/SPscale2;
		Spaceprof3[i] = M3*SINlog(BeachMid[i], S3s, B_S3s)/SPscale3;
	}	
	
	// Calculate deposition time series
	// The time-profile is converted to deposition by being passed through
	// a deposition matrix
	// If N birds die on day d, then these are commuted through the deposition
	// matrix onto birds deposited on days d to d+14 (length of particle sims)
	// So MortMat is the weighted deposition
	for(i in 1:213){
	  MortMat1[,i] = Timeprof1[i]*BigDMat[,i];
		MortMat2[,i] = Timeprof2[i]*BigDMat[,i];
	}
	
	// Now take the sum of those to get expected temporal deposition
	// on each day as the sum across previous days mortality that result in
	// deposition on day i
	for(i in 1:Nday){
	  Temp1[i] = sum(MortMat1[i,]);
		Temp2[i] = sum(MortMat2[i,]);
	}
		
	// Beach and day-specific deposition matrix
	// Combine the daily deposition (time profile * dep matrix) with
	// space to account for spatial variation
	// ABD is then the expected deposition each day on a specific beach
	for(i in 1:Nbeach){
		for(j in 1:Nday){
			ABD[i,j] = (D0 + (Temp1[j]*Spaceprof1[i]) +	(Temp2[j]*(Spaceprof2[i] + Spaceprof3[i])))*BeachLength[i];
	  }
	}
		
	// Summation over previous days which feed into Cbd
	// Effectively, if there was x carcasses deposited d days ago
  // then we expect x*psi^d of those to remain
  // So if we have a vector of x and multiply by a vector of lag-specific 
  // persistence rates and take the sum we get the number of carcasses
  // expected to remain for a particular survey
	for(i in 1:Nbeach){
		Sum_AS[i,1] = 0;
		for(j in 2:nsurvmax){
			Sum_AS[i,j] = dopro(AGrab(ABD, daynum[i,j-1], daynum[i,j]-daynum[i,j-1], i), Survgen(psi[i],daynum[i,j]-daynum[i,j-1])); 
		}
	}
	
	// Calculation of expected count for a given survey
	// Consists of birds that were present, but uncounted on previous survey
	// + new deposition (minus persistence loss) in intervening time
	// and including missed detections
	for(i in 1:Nbeach){
	  C_hat[i,1] = 0;
		for(j in 2:nsurvmax){
			C_hat[i,j] = ((C_hat[i,j-1]*(1-p_val)*(psi[i]^(daynum[i,j]-daynum[i,j-1]))) + p_val*Sum_AS[i,j]);
		}
	}
	
	// Likelihood
	for(i in 1:Nbeach){											// loops over beaches
		for(j in 2:nsurv[i]){									// loops over surveys within beaches
			C[i,j] ~ neg_binomial_2(C_hat[i,j],nu);
		}
	}
}