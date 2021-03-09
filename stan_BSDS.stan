data {
    // import data sizes //
      int<lower=0> N; // total sample size
      int<lower=0> N_sp; // number of species
      int D_edge; // location of the phyloenetic edge, where directional selection occured

    // import phylogenetic data //
      int len_phylo; // length of phylogenetic tree
      vector[len_phylo] branch_len; // vector of branch length of all the edges 
      int tree_obj[len_phylo, 2];
      int MRCA_ij[N_sp, N_sp]; // i,j elements correspond to the location of their MRCA in the tree
      
    // import regressin model data //
      vector[N] y;        // objective variable
      int Z[N];  // spices id for random effect

      }

  parameters {
    real MRCA;
    real<lower=0> ev;
    real<lower=1> sel;
  }

  transformed parameters{
    real<lower=0> sigma[N];
    vector<lower=1>[len_phylo] k;

    vector[len_phylo+1] sim_mean;
    vector[len_phylo+1] sim_var;
  
    cov_matrix[N_sp] vcv_BSDE;
    vector[N] mu;
  
    sim_mean[tree_obj[1,1]]= MRCA;
    sim_var [tree_obj[1,1]]= 0;
  
    for(i in 1:len_phylo) k[i] = 1;
    k[DS_edge] = sel; // set the edge number where directional evolution occured.
  
    for(i in 1:len_phylo){
      sim_mean[tree_obj[i,2]] = sim_mean[tree_obj[i,1]] + (branch_len[i]*ev*(k[i]^2-1))/k[i];
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
  
    for(i in 1:N){
      mu[i] = sim_mean[Z[i]];
      sigma[i] = sqrt(fabs(vcv_BSDE[Z[i],Z[i]]));
    }
  
  }

  model {
    // set the prior
      //target+= uniform_lpdf(ev| 0, 10000); // prior on ev
      //target+= uniform_lpdf(sel| 1, 30);   // prior on sel
      //target+= uniform_lpdf(MRCA| 1, 30);  // prior on MRCA
    
    // set the likelihood
     target+= normal_lpdf(y| mu, sigma);

  }

  generated quantities {
  vector[N] log_likelihood;
  // obtain log likelihood
    for(i in 1:N){log_likelihood[i] = normal_lpdf(y[i]| mu[i], sigma[i]);}
  }


