#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include "/Users/SuperOyasumi/Library/Mobile Documents/com~apple~CloudDocs/simulate/ohkubo_simu/MT.h"
#include <Rcpp.h>
using namespace Rcpp;
float inf = std::numeric_limits<float>::infinity(); // C++

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

double get_LH(double, double, double, double);
// get_LH(MRCA, ev_rate_1, alpha, beta);


// phylogeny & phenotype stable //

double phenotype[N];
int branch_a[N];
int branch_d[N];
double er[N];	//evolutionary rate
double bl[N];	//blanch length

// data //


//end

//transformed parameters
//log((1/(sqrt(2*PI)*sd_human)) * (exp( -(phenotype[2]-human)*(phenotype[2]-human) / (2*sd_human*sd_human) )) )





// parameter set //
double MRCA;			//333.0-1330.0, step_min = 25.0
double MIN_MRCA = 0.0;
double MAX_MRCA = 10000.0;
//double step_MRCA = 200.;


double ev_rate_1;			// 1.-10000., step_min = 25.0
double MIN_ev_rate_1 = 1.;
double MAX_ev_rate_1 = 10000.;
//double step_ev_1 = 1000.;

double k;
double MIN_k = 0.0000; // should be larger than zero
double MAX_k = 30.;

//double step_beta = 25.0;

double mutation_step = 1.;




// end //



double get_LH(double MRCA, double ev_rate_1, double k, double human, double chimp, double gorilla, double orang){
  
  
  
  double sd_human = 123.84;
  double sd_chimp = 46.26;
  double sd_gorilla = 73.84;
  double sd_orang = 52.04;
  
  er[0]=ev_rate_1;
  er[1]=ev_rate_1;
  er[2]=ev_rate_1;
  er[3]=ev_rate_1;
  er[4]=ev_rate_1;
  er[5]=ev_rate_1;
  
  int mu_plus, mu_minus;
  int i, j;
  
  double ln_prob_norm;
  
  
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
  
  for (i=0; i<N; i++) phenotype[i] += MRCA;
  //printf("%lf\n",phenotype[3]);
  
  ln_prob_norm =
    log((1/(sqrt(2*PI)*sd_human)) * (exp( -(phenotype[2]-human)*(phenotype[2]-human) / (2*sd_human*sd_human) )) )
    + log((1/(sqrt(2*PI)*sd_chimp)) * (exp( -(phenotype[3]-chimp)*(phenotype[3]-chimp) / (2*sd_chimp*sd_chimp) )) )
    + log((1/(sqrt(2*PI)*sd_gorilla)) * (exp( -(phenotype[4]-gorilla)*(phenotype[4]-gorilla) / (2*sd_gorilla*sd_gorilla) )))
    + log((1/(sqrt(2*PI)*sd_orang)) * (exp( -(phenotype[5]-orang)*(phenotype[5]-orang) / (2*sd_orang*sd_orang) )));
    
    
    return (ln_prob_norm);
    
}

//return ln_prob_norm;





// [[Rcpp::export]]
DataFrame rcpp_ABC_Rhat(int n, double phenotype2, double phenotype3, double phenotype4, double phenotype5)
{    // time count
  // clock_t start, end;
  double human = phenotype2;		//phenotype[2]
  double chimp = phenotype3;		//phenotype[3]
  double gorilla = phenotype4;	//phenotype[4]
  double orang = phenotype5;		//phenotype[5]
  
  
  int j;
  int chain =4;
  int perc = n/10;  
  double rhat=0;
  
  //double abc_result[n][7];
  NumericMatrix result_MRCA (n, chain);
  NumericMatrix result_k (n, chain);
  NumericMatrix result_ev_rate (n, chain);
  NumericMatrix result_loglik (n, chain);
  // time count
  // start = clock();
  // timer;
  // time(&timer);
  
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
  
  float inf = INFINITY;  
  
  //fprintf(fp, "\n*** %s*** Phylogenetic model of phenotypic evolution version. 1.0 by Nobuyuki Kutsukake ***\n\n", ctime(&timer));
  //end
  
  // parameter set //
  double MRCA;            //333.0-1330.0, step_min = 25.0
  double MIN_MRCA = 0.0;
  double MAX_MRCA = 10000.0;
  //double step_MRCA = 200.;
  
  
  double ev_rate_1;            // 1.-10000., step_min = 25.0
  double MIN_ev_rate_1 = 1.;
  double MAX_ev_rate_1 = 10000.;
  //double step_ev_1 = 1000.;
  
  double k;
  double MIN_k = 0; // should be larger than zero
  double MAX_k = 30.;
  
  
  int i;
  double temp_LH;
  
  do{
    //printf("Sampling Started...\n");    
    
    for (j=0; j<chain; j++ ){
      //Generating Initial Value
      do {
        
        ev_rate_1 = (MAX_ev_rate_1 - MIN_ev_rate_1) * genrand_real2() + MIN_ev_rate_1;
        MRCA = (MAX_MRCA-MIN_MRCA)*genrand_real2()+MIN_MRCA;
        k = (MAX_k - MIN_k)*genrand_real2()+MIN_k;
        
        temp_LH =  get_LH(MRCA, ev_rate_1, k, human, chimp, gorilla, orang);
      }while (temp_LH < -200);
      
      
      result_MRCA(0,j)=MRCA;
      result_k(0,j)=k;
      result_ev_rate(0,j)=ev_rate_1;
      result_loglik(0,j)=temp_LH;
      // Rcout << "The value of k : " << result_k(0,j) << "\n";
      //printf("%1.2lf,%1.2lf,%1.2lf,%1.2lf\n", abc_result[0][0],  abc_result[0][1],  abc_result[0][2],  abc_result[0][3]);
      
      
      //Sampling from the conditional simulation
      for (i=1; i<n; ) {
        
        //提案??????
        MRCA = rnorm2(result_MRCA(i-1, j), 20);
        k = result_k(i-1, j);
        ev_rate_1 = result_ev_rate(i-1, j);
        
        
        //事後???度=尤度*区間付き一様事前??????
        if ((MIN_MRCA < MRCA && MRCA < MAX_MRCA)) {
          temp_LH =  get_LH(MRCA, ev_rate_1, k, human, chimp, gorilla, orang);;
        }else{
          temp_LH = -inf;
          
        }
        
        if (temp_LH > result_loglik(i-1, j)) {
          result_MRCA(i, j)=MRCA;
          result_loglik(i, j)=temp_LH;
        }else if (rbinom(exp(temp_LH - result_loglik(i-1, j)))==1){
          result_MRCA(i, j)=MRCA;
          result_loglik(i ,j)=temp_LH;
        }else{
          result_MRCA(i, j) = result_MRCA(i-1, j);
          result_loglik(i, j)= result_loglik(i-1, j);
        }
        
        
        //提案??????
        MRCA = result_MRCA(i, j);
        k =  rnorm2(result_k(i-1, j),1);
        ev_rate_1 = rcauchy(result_ev_rate(i-1, j), 50);
        
        //事後???度=尤度*区間付き一様事前??????
        if (
            (MIN_k< k && k <MAX_k)&&
              (MIN_ev_rate_1 < ev_rate_1 && ev_rate_1 < MAX_ev_rate_1)) {
          temp_LH =  get_LH(MRCA, ev_rate_1, k, human, chimp, gorilla, orang);;
        }else{
          temp_LH = -inf;
          
        }
        
        
        if (temp_LH > result_loglik(i,j)) {
          result_k(i, j)=k;
          result_ev_rate(i, j)=ev_rate_1;
          result_loglik(i, j)=temp_LH;
        }else if (rbinom(exp(temp_LH - result_loglik(i, j)))==1){
          result_k(i, j)=k;
          result_ev_rate(i, j)=ev_rate_1;
          result_loglik(i, j)=temp_LH;
        }else{
          result_k(i, j) = result_k(i-1, j);
          result_ev_rate(i, j) = result_ev_rate(i-1, j);
          
        }
        if (i%perc==0){
          //printf("%8.3lf loglik\n",result_loglik(i, j));
          // printf("%d samples collected\n",i);
        }
        i++;
      }
      //printf("%d chain collected\n",j);
    }
    
    //calculation of Rhat, convergence monitor.
    double b=0;
    double w=0;
    
    double tmp;
    double var_hat =0;
    
    double within_chain_mean[chain];
    double power_s[chain];
    
    double between_chain_mean=0;
    
    for (j=0; j<chain; j++){
      tmp=0;
      for(i=10001; i< n; i++ ){
        tmp+= result_loglik(i, j);
      }
      within_chain_mean[j] = tmp/(n-10000);
      Rcout << "The value of within chain mean: " << within_chain_mean [j]<< "\n";
      between_chain_mean +=  within_chain_mean[j];
    }
    between_chain_mean= between_chain_mean/chain;
    Rcout << "The value of between chain mean: " << between_chain_mean<< "\n";
    
    // printf("%lf,%lf,%lf,%lf, %lf\n", within_chain_mean[1],within_chain_mean[2],within_chain_mean[3],within_chain_mean[4],between_chain_mean);
    
    
    b=0;
    for (j=0; j< chain; j++ )
    {
      b+= pow(within_chain_mean[j]-between_chain_mean,2);
    }
    b = (n-10000)*b/(chain-1);
    Rcout << "The value of b: " << b<< "\n";
    
    for (j=0; j<chain; j++){
      for (i=10001; i<n; i++){
        power_s[j]+= pow(result_loglik(i,j)-within_chain_mean[j],2);
        //Rcout << "The value of power s: " << power_s[j]<< "\n";
        
      }
      power_s[j]=power_s[j]/(n-1);
      //Rcout << "The value of power s: " << power_s[j]<< "\n";
    }

    w=0;
    for (j=0; j<chain;j++){w+=power_s[j];}
    w=w/chain;
    Rcout << "The value of power w: " << w<< "\n";
    double w_prep =(n-10000-1);
    Rcout << "The value of power w_prep: " << w_prep << "\n";
    w_prep=w_prep/(n-10000);
    Rcout << "The value of power w_prep: " << w_prep << "\n";
    var_hat = w_prep*w + b/(n-10000);
    
    Rcout << "The value of power b/n: " << b/(n-10000)<< "\n";
    
    Rcout << "The value of power var_hat: " << var_hat<< "\n";
    
    rhat = sqrt (var_hat/w);
    Rcout << "The value of Rhat: " << rhat<< "\n";
    
  }while (!( 1 < rhat &&  rhat < 1.01));
  
    DataFrame df = DataFrame::create( Named("MRCA") = result_MRCA,
                                    Named("k") = result_k,
                                    Named("ev_rate") = result_ev_rate,
                                    Named("loglik") = result_loglik,
                                    Named("Rhat") = rhat) ;
  
  //printf("Sampling finished...\n");
  return df;
  //return rhat;
}
