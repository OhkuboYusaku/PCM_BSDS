data {
  // import data sizes //
    int<lower=0> N; // total sample size
    int<lower=0> N_sp; // number of species
    int<lower=0> D; // dimension of explanatory variables
    int D_edge; // location of the phyloenetic edge, where directional selection occured
  
  // import phylogenetic data //
    int len_phylo; // length of phylogenetic tree
    vector[len_phylo] branch_len; // vector of branch length of all the edges 
    int tree_obj[len_phylo, 2]; // tree structure
    int MRCA_ij[N_sp, N_sp]; // i,j elements correspond to the location of their MRCA in the tree
  
  // import regressin model data //
    vector[N] y;        // objective variable
    matrix[N, D] X;   // explanatory variable
    matrix[N, N_sp] Z;  // spices id for random effect
}

 parameters {
  real MRCA; // most common recent ancestor
  real<lower=0> ev; // evolution rate
  //real<lower=1> k[len_phylo];
  real sel; // strength of directional selection
  vector[D] beta; // regression coefficients
  real<lower=0> sigma_y; // standard deviation of objective variable
  vector[N_sp] b; // species specific random effect
  
}

transformed parameters{
  vector[len_phylo] k;
  vector[len_phylo+1] sim_mean;
  vector[len_phylo+1] sim_var;
  
  cov_matrix[N_sp] vcv_BSDE;

  sim_mean[tree_obj[1,1]]= MRCA;
  sim_var [tree_obj[1,1]]= 0;
  
  for(i in 1:len_phylo) k[i] = 1;
  k[D_edge] = exp(sel);
  
  for(i in 1:len_phylo){
    sim_mean[tree_obj[i,2]] =  sim_mean[tree_obj[i,1]] + (branch_len[i]*ev*(k[i]^2-1))/k[i];
    sim_var [tree_obj[i,2]] =  sim_var[tree_obj[i,1]] + (2*branch_len[i]*ev*(k[i]^2+1))/k[i];
  }
  
  for(i in 1: N_sp){
    for(j in i: N_sp){
      if(i != j){
      vcv_BSDE[i, j] = sim_var[MRCA_ij[i, j]];
      vcv_BSDE[j, i] = sim_var[MRCA_ij[i, j]];
      }else{
        vcv_BSDE[i, j] = sim_var[i];
        vcv_BSDE[j, i] = sim_var[i];
      }
    }
  }
 
}

model {
  //target+= uniform_lpdf(ev| 1, 10000);
  //target+= uniform_lpdf(MRCA| 0, 1000);
  
  // add prob. density for realized random-effects
  target+= multi_normal_lpdf(b| sim_mean[1:N_sp], vcv_BSDE);
  
  // define likelihood
  target+= normal_lpdf(y| X * beta + Z * b, sigma_y);
  //target+= poisson_log(y| mu);
  //target+= bernoulli_logit(y| mu);
}


generated quantities {
  vector[N] log_likelihood;
  // obtain log likelihood
  for(n in 1:N){log_likelihood[n] = normal_lpdf(y[n]| X[n,] * beta + Z[n,] * b, sigma_y);}
}
