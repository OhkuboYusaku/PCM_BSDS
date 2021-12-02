#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include "/Users/SuperOyasumi/Library/Mobile Documents/com~apple~CloudDocs/simulate/ohkubo_simu/MT.h"
#include <Rcpp.h>
using namespace Rcpp;

#define N 6	// # branch
//#define PI 3.1415926535897932384626433832795
long  seed = 5;    // [time?]



double Uniform( void ){
  return genrand_real3();
}

double rcauchy( double mu, double gamma){
  return mu + gamma*tan(M_PI*( Uniform()-0.5 ));
}


double rnorm2( double mu, double sigma ){
  double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
  return mu + sigma*z;
}

float expdev(long idum){
  return -log(Uniform())/ idum;
  
}

int rbinom(double p){
  double j;
  j=genrand_real1();
  if (j<p) {
    return 1;
  }else{
    return 0  ;  }
  
}

//double get_LH(double, double, double, double);
// get_LH(MRCA, ev_rate_1, alpha, beta);


// phylogeny & phenotype stable //

double phenotype[N];
int branch_a[N];
int branch_d[N];
double er[N];	//evolutionary rate
double bl[N];	//blanch length

// data //

double human = 1321.27;		//phenotype[2]
double chimp = 348.08;		//phenotype[3]
double gorilla = 467.12;	//phenotype[4]
double orang = 334.76;		//phenotype[5]

double sd_human = 123.84;
double sd_chimp = 46.26;
double sd_gorilla = 73.84;
double sd_orang = 52.04;
//end

//transformed parameters
//log((1/(sqrt(2*PI)*sd_human)) * (exp( -(phenotype[2]-human)*(phenotype[2]-human) / (2*sd_human*sd_human) )) )





// parameter set //
double MRCA;			//333.0-1330.0, step_min = 25.0

//double step_MRCA = 200.;


double ev_rate_1;			// 1.-10000., step_min = 25.0
//double step_ev_1 = 1000.;

double k;


//double step_beta = 25.0;

double mutation_step = 1.;




// end //


// [[Rcpp::export]]
NumericVector rcpp_LH(double MRCA, double ev_rate_1, double k){
  
  bl[0]=6.48; // (MRCA, A1)
  bl[1]=2.48; // (A1, A2)
  bl[2]=6.18; // (A2, human)
  bl[3]=6.18; // (A2, chimp)
  bl[4]=8.65; // (A1, gorilla)
  bl[5]=15.13; // (MRCA, orang)
  
  branch_a[0]=0; branch_d[0]=1; // (MRCA, A1)
  branch_a[1]=1; branch_d[1]=2; // (A1, A2)
  branch_a[2]=2; branch_d[2]=3; // (A2, human)
  branch_a[3]=2; branch_d[3]=4; // (A2, chimp)
  branch_a[4]=1; branch_d[4]=5; // (A1, gorilla)
  branch_a[5]=0; branch_d[5]=6; // (MRCA, orang)
  
  
  er[0]=ev_rate_1;
  er[1]=ev_rate_1;
  er[2]=ev_rate_1;
  er[3]=ev_rate_1;
  er[4]=ev_rate_1;
  er[5]=ev_rate_1;
  
  int mu_plus, mu_minus;
  int i, j;
  
  //double ln_prob_norm;
  
  
  phenotype[0] = 0;
  phenotype[1] = 0;
  phenotype[2] = 0;
  phenotype[3] = 0;
  phenotype[4] = 0;
  phenotype[5] = 0;
  
  for (i=0; i<N; i++) {
    
    
    
    //evolution of sp[2]
    if (i==2) {
      
      if (k>0){
        //mu_plus = poisson_rnd( er[i]*bl[i]*k );
        //mu_minus = poisson_rnd( er[i]*bl[i]/k );
        
        mu_plus  = rnorm2( er[i]*bl[i]*k , sqrt(er[i]*bl[i]*k));
        mu_minus = rnorm2( er[i]*bl[i]/k , sqrt(er[i]*bl[i]/k));
        
      }else if (k<0){
        
        mu_plus  = rnorm2( er[i]*bl[i]/(-k) , sqrt(er[i]*bl[i]/(-k)));
        mu_minus = rnorm2( er[i]*bl[i]*(-k) , sqrt(er[i]*bl[i]*(-k)));
        
      }else{
        mu_plus  = rnorm2( er[i]*bl[i] , sqrt(er[i]*bl[i]));
        mu_minus = rnorm2( er[i]*bl[i] , sqrt(er[i]*bl[i]));
        
      }
      
      phenotype[i]= rnorm2(mu_plus/seed-mu_minus/seed,sqrt((mu_plus+mu_minus)/pow(seed,2)));
      
      
    }else {
      mu_plus  = rnorm2( er[i]*bl[i] , sqrt(er[i]*bl[i]));
      mu_minus = rnorm2( er[i]*bl[i] , sqrt(er[i]*bl[i]));
      
      phenotype[i]= rnorm2(mu_plus/seed-mu_minus/seed,sqrt((mu_plus+mu_minus)/pow(seed,2)));
    }
    
    for (j=0; j<N; j++) {
      if (branch_d[i]==branch_a[j]) {
        phenotype[j] += phenotype[i];
      }
    }
    
  }
  
  for (i=0; i<N; i++){ phenotype[i] += MRCA;}
  //printf("%lf\n",phenotype[0]);
  
  
  NumericVector v = NumericVector::create(phenotype[2],phenotype[3],phenotype[4],phenotype[5]);
  
    return (v);
    
}
