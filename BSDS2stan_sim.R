########################################## BSDS_sim #####################################################
#R-code to conduct simulation experiments for BSDS model
#Article Ref: 
#  Authors: Ohkubo et al. (unpub)
#Title: "Evaluating a strength of directional selection 
#         using a novel branch-specific directional selection (BSDS) model of phylogenetic comparative method", 
#Volume: XXX.

#Softwere Information:
#  Import: "phylo" object of sps., and evolutionary parameters.
#Outputs: list of mean vector and variancecovariance matrix of tip sp. 
#Dependencies: "ape"

#Source Author: Ph.D. Ohkubo Yusaku123
#1.	Center for Data Assimilation Research and Applications, 
#  Joint Support Center for Data Science Research, Research Organization of Information and Systems, Tokyo, Japan
#2.	Institute of Statistical Mathematics, Tokyo, Japan
#3.	Graduate School of Environmental Science, Hokkaido University, Sapporo, Japan.

# Email: ohkubo.yusaku1989[at]gmail.com
############################################################################################################

BSDE2_stan_sim<-function(N_sample, N_sp, MRCA, ev_rate, sel, n_sim){
  library(ape)
  library(rstan)
  library(loo)
  # define stan model
  # Model 1: BSDE
  stancode_BSDE <- '
   data {
    // import data sizes //
      int<lower=0> N;
      int<lower=0> N_sp;
      int DE_edge;

    // import phylogenetic data //
      int len_phylo;
      vector[len_phylo] branch_len;
      int tree_obj[len_phylo, 2];
      int MRCA_ij[N_sp, N_sp]; // i,j elements correspond to the location of their MRCA in the tree
  
    // import regressin model data //
      vector[N] y;        // objective variable
      int<lower=1> Z[N];  // spices id for random effect
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
    k[DE_edge] = sel; // set the edge number where directional evolution occured.
  
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
      //target+= uniform_lpdf(ev| 0, 10000);
      //target+= uniform_lpdf(sel| 1, 30);
    
    // set the likelihood
     target+= normal_lpdf(y| mu, sigma);

  }

  generated quantities {
  vector[N] log_likelihood;

    for(i in 1:N){log_likelihood[i] = normal_lpdf(y[i]| mu[i], sigma[i]);}
  }
  '
  
  stan_BSDE<- stan_model(model_code = stancode_BSDE, verbose = F)
  
  ## Model2: variable rate Brownian Motion model
  
  stancode_vBM <- '
  data {
  // import data sizes //
    int<lower=0> N;
    int<lower=0> N_sp;
    int DE_edge;
  // import phylogenetic data //
    int len_phylo;
    vector[len_phylo] branch_len;
    int tree_obj[len_phylo, 2];
    int MRCA_ij[N_sp, N_sp]; // i,j elements correspond to the location of their MRCA in the tree
  
  // import regressin model data //
    vector[N] y;        // objective variable
    int<lower=1> Z[N];  // spices id
}

 parameters {
  real MRCA;
  real<lower=0> ev_base;
  real<lower=1> acc_DE_edge;

}

transformed parameters{
  real<lower=0> sigma[N];
  vector[len_phylo] ev;

  vector[len_phylo+1] sim_mean;
  vector<lower=0>[len_phylo+1] sim_var;
  
  cov_matrix[N_sp] vcv_BSDE;
  vector[N] mu;
  
  sim_mean[tree_obj[1,1]]= MRCA;
  sim_var [tree_obj[1,1]]= 0;
  
  for(i in 1:len_phylo) ev[i] = ev_base;
  ev[DE_edge] = ev_base * acc_DE_edge; // set the edge number where directional evolution occured.
  
  for(i in 1:len_phylo){
    sim_mean[tree_obj[i,2]] = sim_mean[tree_obj[i,1]];
    sim_var [tree_obj[i,2]] =  sim_var[tree_obj[i,1]] + (2*branch_len[i]*ev[i]);
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
   // mu[i] = b[Z[i]];
  }
  
}

model {
  // set the prior
    //target+= uniform_lpdf(ev_base| 0, 10000);
    //target+= uniform_lpdf(acc_DE_edge| 1, 1000);

  // set the likelihood
    target+= normal_lpdf(y| mu, sigma);

}

  generated quantities {
  vector[N] log_likelihood;
    for(i in 1:N){log_likelihood[i] = normal_lpdf(y[i]| mu[i], sigma[i]);}
  }
  '

stan_vBM<- stan_model(model_code = stancode_vBM, verbose = F)

## Model3: Brownian Motion model

stancode_BM <- '
 data {
  // import data sizes //
    int<lower=0> N;
    int<lower=0> N_sp;
    int DE_edge;
    //int<lower=0> D;
  
  // import phylogenetic data //
    int len_phylo;
    vector[len_phylo] branch_len;
    int tree_obj[len_phylo, 2];
    int MRCA_ij[N_sp, N_sp]; // i,j elements correspond to the location of their MRCA in the tree
  
  // import regressin model data //
    vector[N] y;        // objective variable
    int<lower=1> Z[N];  // spices id for random effect
}

 parameters {
  real MRCA;
  real<lower=0> ev_base;

}

transformed parameters{
  real<lower=0> sigma[N];

  vector[len_phylo+1] sim_mean;
  vector[len_phylo+1] sim_var;
  
  cov_matrix[N_sp] vcv_BSDE;
  vector[N] mu;
  
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
  
  for(i in 1:N){
    mu[i] = sim_mean[Z[i]];
    sigma[i] = sqrt(fabs(vcv_BSDE[Z[i],Z[i]]));
   // mu[i] = b[Z[i]];
  }
  
}

model {
  // set the prior
    //target+= uniform_lpdf(ev_base| 0, 10000);

  // set the likelihood
    target+= normal_lpdf(y| mu, sigma);

}

generated quantities {
vector[N] log_likelihood;

  for(i in 1:N){log_likelihood[i] = normal_lpdf(y[i]| mu[i], sigma[i]);}
}

  '

stan_BM<- stan_model(model_code = stancode_BM, verbose = F)
  
  # start definition of R-functions
  phylo2ABC<- function(phylo, MRCA, ev, k){
    #define data
    tree_obj<- as.matrix(phylo[["edge"]])
    simulated_trait<- numeric(length(phylo$edge.length))
    len_phylo<- length(phylo$edge.length)
    
    ## conduct evolution simulation
    simulated_trait[tree_obj[1,1]]<- MRCA
    for(i in 1:len_phylo){
      branch_len=phylo$edge.length[i]
      mut_plus<-  max(0, rnorm(1, branch_len*ev*k[i], sqrt(branch_len*ev*k[i]))) # the number of positive mutation
      mut_minus<- max(0, rnorm(1, branch_len*ev/k[i], sqrt(branch_len*ev/k[i]))) # the number of negative mutation
      simulated_trait[tree_obj[i,2]]<- rnorm(1, (mut_plus - mut_minus), sqrt(mut_plus + mut_minus)) + simulated_trait[tree_obj[i,1]]
    }
    
    # output only tips
    N_tip<- len_phylo - phylo$Nnode +1
    tip<- numeric(N_tip)
    simulated_trait<- simulated_trait[1:N_tip] # delete inner nodes
    
    return(simulated_trait)
  }
  
  
  phylo2stan_data<- function(phylo, y, Z, DE_edge){
    len_phylo<- length(phylo$edge.length)
    N_tip<- len_phylo - phylo$Nnode +1
    branch_len<- phylo$edge.length
    tree_obj<- as.matrix(phylo$edge)
    MRCA_ij<- matrix(0, N_tip, N_tip) ## i,j elements correspond to the location of their MRCA in the tree
    
    for(i in 1:N_tip){
      for(j in i:N_tip){
        MRCA_ij[i,j]<- MRCA_ij[j,i]<- getMRCA(phylo, tip=c(i,j))
      }
    }
    
    dat<- list(N=length(y), N_sp=N_tip, 
               len_phylo=len_phylo, branch_len=branch_len, tree_obj=tree_obj, MRCA_ij=MRCA_ij,
               y=y, Z=Z, DE_edge=DE_edge)
    
    return (dat)
  }
  
  map_mcmc <- function(z){ density(z)$x[which.max(density(z)$y)]} 
  
  # end definition of R-functions
  
  sim_date<- as.Date(Sys.time()) #extract the experiment date
  result_BSDE<- matrix(0, n_sim, 18) # results container
  result_vBM<-  matrix(0, n_sim, 6)  # results container
  result_BM<-  matrix(0, n_sim, 6)

  
 
  for(i_sim in 1:n_sim){
    ##### generate data ####
    phylo<- rtree(N_sp) # random sampling of phylogenetic tree
    len_phylo<- length(phylo$edge.length)
    phylo$edge.length<- phylo$edge.length*100
    N_tip<- len_phylo - phylo$Nnode + 1
    
    k<- numeric(length(phylo$edge.length))
    DE_edge<- floor(runif(1, min=1, max=len_phylo+1))
    
    k[1:length(phylo$edge.length)]<- 1
    k[DE_edge]<- sel
    
    y<- matrix(0, N_sample, N_tip)
    for(j in 1: N_sample){
      y[j,]<-phylo2ABC(phylo, MRCA, ev_rate, k=k)
    }
    y<- as.numeric(y)
    Z<- numeric(N_tip*N_sample)
    for(j in 1:N_tip){
      for(k in 1:N_sample){
        Z[(j-1)*N_sample + k]<- j
      }
    }
    
    #### conduct Bayesian inference ###
    dat<- phylo2stan_data(phylo = phylo, y=y, Z=Z, DE_edge=DE_edge)
    
    war<- 5000
    ite<- 10000
    cha<- 2
    
    # BSDE model
    par = c("MRCA",  "ev", "sel", "log_likelihood")
    Rhat<- FALSE
    fit_BSDE<- sampling(stan_BSDE, data = dat, chains = cha, pars = par,
                        warmup = war, iter=ite,thin=1, verbose=F, open_progress=F,
                        refresh = 0, control = list(adapt_delta=0.9))  
    Rhat_below_1.05<- all(summary(fit_BSDE)$summary[,"Rhat"] <= 1.01, na.rm = T)
    
    # if MCMC not converged, then multiply the iterations and retry
      while(Rhat_below_1.05==FALSE) {
        war<- 2 * war
        ite<- 2 * ite
        fit_BSDE<- sampling(stan_BSDE, data = dat, chains = cha, pars = par,
                            warmup = war, iter=ite,thin=1, verbose=F, open_progress=F,
                            refresh = 0, control = list(adapt_delta=0.9))  
        Rhat_below_1.05<- all(summary(fit_BSDE)$summary[,"Rhat"] <= 1.01, na.rm = T)
      }
    
    a<- extract(fit_BSDE)
    result_BSDE[i_sim,1]<- mean(a$MRCA)
    result_BSDE[i_sim,2]<- mean(a$ev)
    result_BSDE[i_sim,3]<- mean(a$sel)

    result_BSDE[i_sim,4]<- median(a$MRCA)
    result_BSDE[i_sim,5]<- median(a$ev)
    result_BSDE[i_sim,6]<- median(a$sel)

    result_BSDE[i_sim,7]<- sd(a$MRCA)
    result_BSDE[i_sim,8]<- sd(a$ev)
    result_BSDE[i_sim,9]<- sd(a$sel)

    full_CI<- sort(a$sel) # sort it.
    
    result_BSDE[i_sim,10]<- full_CI[0.025 * length(full_CI)] ##CI of k
    result_BSDE[i_sim,11]<- full_CI[0.975 * length(full_CI)]
    if ((result_BSDE[i_sim,10] < sel) && (sel < result_BSDE[i_sim,11])) {
      result_BSDE[i_sim,12]<- 1
    }else{
      result_BSDE[i_sim,12]<- 0  
    }
    
    result_BSDE[i_sim,13]<- waic(a$log_likelihood)[["waic"]]
    result_BSDE[i_sim,14]<- mean(a$lp__)
    result_BSDE[i_sim,15]<- ite
    
    result_BSDE[i_sim,16]<- map_mcmc(a$sel)
    
    # variable-rate Brownian Motion model
    war<- 5000
    ite<- 10000

    par = c("MRCA", "log_likelihood")
    Rhat<- FALSE
    fit_vBM<- sampling(stan_vBM, data = dat, chains = cha, pars = par,
                       warmup = war, iter=ite,thin=10, verbose=F, open_progress=F,
                       refresh = 0, control = list(adapt_delta=0.9))
    Rhat_below_1.05<- all(summary(fit_vBM)$summary[,"Rhat"] <= 1.01, na.rm = T)
    
    # if MCMC not converged, then multiply the iterations and retry
    while(Rhat_below_1.05==FALSE) {
      war<- 2 * war
      ite<- 2 * ite
      fit_vBM<- sampling(stan_vBM, data = dat, chains = cha, 
                         warmup = war, iter=ite,thin=10, verbose=F, open_progress=F,
                         refresh = 0, control = list(adapt_delta=0.9))
      Rhat_below_1.05<- all(summary(fit_vBM)$summary[,"Rhat"] <= 1.01, na.rm = T)
    }

    
    a<- extract(fit_vBM)
    result_vBM[i_sim,1]<- mean(a$MRCA)
    result_vBM[i_sim,2]<- median(a$MRCA)
    result_vBM[i_sim,3]<- sd(a$MRCA)

    result_vBM[i_sim,4]<- waic(a$log_likelihood)[["waic"]]
    result_vBM[i_sim,5]<- mean(a$lp__)
    result_vBM[i_sim,6]<- ite
    
    # ordinal Brownian Motion model
    war<- 5000
    ite<- 10000
    
    par = c("MRCA", "log_likelihood")
    Rhat<- FALSE
    fit_BM<- sampling(stan_BM, data = dat, chains = cha, pars = par,
                       warmup = war, iter=ite,thin=10, verbose=F, open_progress=F,
                       refresh = 0, control = list(adapt_delta=0.9))
    Rhat_below_1.05<- all(summary(fit_BM)$summary[,"Rhat"] <= 1.01, na.rm = T)
    
    # if MCMC not converged, then multiply the iterations and retry
    while(Rhat_below_1.05==FALSE) {
      war<- 2 * war
      ite<- 2 * ite
      fit_BM<- sampling(stan_BM, data = dat, chains = cha, 
                         warmup = war, iter=ite,thin=10, verbose=F, open_progress=F,
                         refresh = 0, control = list(adapt_delta=0.9))
      Rhat_below_1.05<- all(summary(fit_BM)$summary[,"Rhat"] <= 1.01, na.rm = T)
    }
    
    
    a<- extract(fit_BM)
    result_BM[i_sim,1]<- mean(a$MRCA)
    result_BM[i_sim,2]<- median(a$MRCA)
    result_BM[i_sim,3]<- sd(a$MRCA)
    
    result_BM[i_sim,4]<- waic(a$log_likelihood)[["waic"]]
    result_BM[i_sim,5]<- mean(a$lp__)
    result_BM[i_sim,6]<- ite
    
    if((result_BSDE[i_sim,13] < result_vBM[i_sim,4])){
      result_BSDE[i_sim, 17]<- 0
    }else{
      result_BSDE[i_sim, 17]<- 1
    }
    
    if((result_BSDE[i_sim,14] > result_vBM[i_sim,5])){
      result_BSDE[i_sim, 18]<- 0
    }else{
      result_BSDE[i_sim, 18]<- 1
    }
  
    print(i_sim)
  }
  
  dat1 = data.frame(
    mean_MRCA = result_BSDE[,1], vBM_mean_MRCA = result_vBM[,1], BM_mean_MRCA = result_BM[,1], 
    med_MRCA = result_BSDE[,4],  vBM_med_MRCA = result_vBM[,2], BM_med_MRCA = result_BM[,2], 
    sd_MRCA = result_BSDE[,7],   vBM_sd_MRCA = result_vBM[,3],  as_BM_MRCA = result_BM[,3], 
    
    mean_ev = result_BSDE[,2], #vBM_mean_ev_common = result_vBM[,2], 
    med_ev = result_BSDE[,5],  #vBM_med_ev_common = result_vBM[,5], 
    sd_ev = result_BSDE[,8],   #vBM_sd_ev_common = result_vBM[,8],  
    
    mean_sel = result_BSDE[,3], med_sel= result_BSDE[,6], MAP_sel= result_BSDE[,16], sd_sel= result_BSDE[,9], 
    lo_sel = result_BSDE[,10],  up_sel = result_BSDE[,11], inc_CI = result_BSDE[,12],
    
    waic = result_BSDE[,13], vBM_waic = result_vBM[,4], BM_waic = result_BM[,4],
    mar_lik = result_BSDE[,14], vBM_mar_lik = result_vBM[,5], BM_mar_lik = result_BM[,5], 
    max_ite_BSDE = result_BSDE[,15], max_ite_vBM = result_vBM[,6], max_ite_BM = result_BM[,6],
    
    waic_error = result_BSDE[,17], mar_lik_error = result_BSDE[,18] 
  )
  
  return(dat1)
}

# example 
N_sample<- 10 # number of sample per species
N_sp<- 10 # number of species
sel<- 2 # strength of the directional selection

exp1<- BSDE2_stan_sim(N_sample=N_sample, N_sp=N_sp, MRCA=100, ev_rate=10, sel=sel, n_sim=10)
