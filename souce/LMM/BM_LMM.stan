data {
  // import data sizes //
    int<lower=0> N;
    int<lower=0> N_sp;
    int<lower=0> D;
  
  // import phylogenetic data //
    int len_phylo;
    vector[len_phylo] branch_len;
    int tree_obj[len_phylo, 2];
    int MRCA_ij[N_sp, N_sp]; // i,j elements correspond to the location of their MRCA in the tree
  
  // import regressin model data //
    vector[N] y;        // objective variable
    matrix[N, D] X;   // explanatory variable
    matrix[N, N_sp] Z;  // spices id for random effect
}

 parameters {
  real MRCA;
    real<lower=0> ev_base;
  //real<lower=1> k[len_phylo];
  vector[N_sp] b;
  vector[D] beta;
  real<lower=0> sigma_y;

}

transformed parameters{
  //real<lower=0> sigma[N];

  vector[len_phylo+1] sim_mean;
  vector[len_phylo+1] sim_var;
  
  cov_matrix[N_sp] vcv_BSDE;

  sim_mean[tree_obj[1,1]]= MRCA;
  sim_var [tree_obj[1,1]]= 0;
  
  
  for(i in 1:len_phylo){
    sim_mean[tree_obj[i,2]] = sim_mean[tree_obj[i,1]];
    sim_var [tree_obj[i,2]] =  sim_var[tree_obj[i,1]] + (2*branch_len[i]*ev_base);
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
  
  //target+= multi_normal_cholesky_lpdf(b| sim_mean[1:N_sp], vcv_BSDE);
  target+= multi_normal_lpdf(b| sim_mean[1:N_sp], vcv_BSDE);
  
  target+= normal_lpdf(y| X * beta + Z * b, sigma_y);

  //target+= poisson_log(y| mu);
  //target+= bernoulli_logit(y| mu);
}

generated quantities {
vector[N] log_likelihood;

  for(n in 1:N){log_likelihood[n] = normal_lpdf(y[n]| X[n,] * beta + Z[n,] * b, sigma_y);}
}

