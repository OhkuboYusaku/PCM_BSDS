 data {
  // import data sizes //
    int<lower=0> N;
    int<lower=0> N_sp;
    
  // import phylogenetic data //
    int DE_edge;
    int len_phylo;
    vector[len_phylo] branch_len;
    int tree_obj[len_phylo, 2];
    int MRCA_ij[N_sp, N_sp]; // i,j elements correspond to the location of their MRCA in the tree
  
  // import regressin model data //
    vector[N] y;        // objective variable
    matrix[N, N_sp] Z;  // spices id for random effect
}

 parameters {
  real MRCA;
  real<lower=0> ev_base;
  //real<lower=1> k[len_phylo];
  vector[N_sp] b;
  real<lower=0> sigma_y;
  real acc_DE_edge;
}

transformed parameters{
  //real<lower=0> sigma[N];

  vector[len_phylo+1] sim_mean;
  vector[len_phylo+1] sim_var;
  vector[len_phylo] acc;

  cov_matrix[N_sp] vcv_BSDE;
  
  for(i in 1:N_sp){ sim_mean[i] = MRCA;}
  
  sim_var [tree_obj[1,1]]= 0;

  for(i in 1:len_phylo) acc[i] = 1;
  acc[DE_edge] = exp(acc_DE_edge); // set the edge number where directional evolution occured.
  
  for(i in 1:len_phylo){
    sim_var [tree_obj[i,2]] =  sim_var[tree_obj[i,1]] + (2*branch_len[i]*ev_base*acc[i]);
  }
  
  for(i in 1: N_sp){
    for(j in i: N_sp){
      if(i != j){
      vcv_BSDE[i, j] = sim_var[MRCA_ij[i, j]];
      vcv_BSDE[j, i] = sim_var[MRCA_ij[i, j]];
      }else{
        vcv_BSDE[i, j] = sim_var[i];
      }
    }
  }
  
}

model {
  target+= uniform_lpdf(acc_DE_edge| -10E06, 10E06);
  //target+= uniform_lpdf(MRCA| 0, 1000);
  
  //target+= multi_normal_cholesky_lpdf(b| sim_mean[1:N_sp], vcv_BSDE);
  target+= multi_normal_lpdf(b| sim_mean[1:N_sp], vcv_BSDE);
  

  target+= normal_lpdf(y| Z * b, sigma_y);
  //target+= poisson_log(y| mu);
  //target+= bernoulli_logit(y| mu);
}


generated quantities {
vector[N] log_likelihood;

  for(n in 1:N){log_likelihood[n] = normal_lpdf(y[n]| Z[n,] * b, sigma_y);}
} 

