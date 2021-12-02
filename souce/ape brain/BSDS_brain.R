###########################################################################################################
#R-codes to analysis brain volume size by Stan, a HMC sampler.
#Article Ref: 
#  Authors: Ohkubo et al.
#  Title: (Unpublished)
#  Volume: 

#Source Author: Ph.D. Ohkubo Yusaku (CHAIN, Hokkaido Univ, Japan. ohkubo.yusaku1989[at]gmail.com)
#  Date: 2020/Mar.8
############################################################################################################
library(ape)
library(rstan)
library(loo)

data<- read.csv("brain.csv")

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

######### brain size analysis #########  
# data prep
s<- "ape((1:8.65,(2:6.18,3:6.18):2.48):6.48,4:15.13);" # import phylogenetic tree in newick format
cat(s, file = "ex.tre", sep = "\n")
# convert into phylo object
tree.ape <- read.tree("ex.tre")
tree.ape$tip.label<- c("Gorilla gorilla", "Homo sapiens", "Pan troglodytes", " Pongo pygmaeus")

N_tip<- nlevels(data$spp)
sample<- nrow(data)

y<- data$brain
Z<- numeric(sample)
for(i in 1:sample){
  if(data$spp[i]=="gorilla")Z[i]<- 1
  if(data$spp[i]=="human")Z[i]<- 2
  if(data$spp[i]=="chimp")Z[i]<- 3
  if(data$spp[i]=="orang")Z[i]<- 4 
}

DE_edge = 4 # the edge number of directional evolution occured
dat<- phylo2stan_data(phylo = tree.ape, y=y, Z=dummies::dummy(Z), DE_edge = 4)
cha<- 2
#options(mc.cores = parallel::detectCores())
war<- 25000
ite<- 50000

# conduct HMC sampling by Stan
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
    matrix[N, N_sp] Z;  // spices id for random effect
}

 parameters {
  real MRCA;
  real<lower=0> ev;
  //real<lower=1> k[len_phylo];
  real sel;
  real<lower=0> sigma_y;
  vector[N_sp] b;
  
}

transformed parameters{
  //real<lower=0> sigma[N];
  vector[len_phylo] k;
  vector[len_phylo+1] sim_mean;
  vector[len_phylo+1] sim_var;
  
  cov_matrix[N_sp] vcv_BSDE;

  sim_mean[tree_obj[1,1]]= MRCA;
  sim_var [tree_obj[1,1]]= 0;
  
  for(i in 1:len_phylo) k[i] = 1;
  k[DE_edge] = exp(sel);
  
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
  
}

model {
  //target+= uniform_lpdf(ev| 1, 10000);
  //target+= uniform_lpdf(MRCA| 330, 1330);
  //target+= uniform_lpdf(sel| 0, log(20));
  
  //target+= multi_normal_cholesky_lpdf(b| sim_mean[1:N_sp], vcv_BSDE);
  target+= multi_normal_lpdf(b| sim_mean[1:N_sp], vcv_BSDE);
  

  target+= normal_lpdf(y| Z * b, sigma_y);
  //target+= poisson_log(y| mu);
  //target+= bernoulli_logit(y| mu);
}


generated quantities {
vector[N] log_likelihood;
vector[N_sp] b_pred;

real k_;
    
    k_ = exp(sel);
    for(n in 1:N){log_likelihood[n] = normal_lpdf(y[n]| Z[n,] * b, sigma_y);}
  }
  '

postpredict_BSDE<- stan_model(model_code = stancode_BSDE, verbose = F)


#par = c("Gorilla","Human", "Chimp", "Orang")
init1<- list(MRCA=500, ev=500, sigma_y=1, sel=2, b=rep(500, 4))
init2<- list(MRCA=500, ev=500, sigma_y=1, sel=2, b=rep(500, 4))

init<- list(init1, init2)
fit_BSDE<- sampling(postpredict_BSDE, data = dat, chains = cha, init=init,
                    warmup = war, iter=ite, thin=10, verbose=F, open_progress=F,
                    control = list(adapt_delta=0.99)) 

plot(fit_BSDE,plotfun="hist", par=c("MRCA","ev", "k_", "sigma_y"))
plot(fit_BSDE,plotfun="trace", par=c("MRCA","ev", "sel", "sigma_y"))

print(fit_BSDE)
# Model2: variable rate Brownian Motion model

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
    matrix[N, N_sp] Z;  // spices id for random effect
}

 parameters {
  real MRCA;
  real<lower=0> ev_base;
  real<lower=0> acc_DE_edge;
  real<lower=0> sigma_y;
  vector[N_sp] b;

}

transformed parameters{
  vector[len_phylo] ev;

  vector[len_phylo+1] sim_mean;
  vector<lower=0>[len_phylo+1] sim_var;
  
  cov_matrix[N_sp] vcv_vBM;

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
      vcv_vBM[i, j] = sim_var[MRCA_ij[i, j]];
      vcv_vBM[j, i] = sim_var[MRCA_ij[i, j]];
      }else{
        vcv_vBM[i, j] = sim_var[i];
        vcv_vBM[j, i] = sim_var[i];
      }
    }
  }
  
}

model {
  // set the prior
    //target+= uniform_lpdf(ev_base| 0, 10000);
    target+= uniform_lpdf(acc_DE_edge| 1, 10E10);

  // set the likelihood
  target+= multi_normal_lpdf(b| sim_mean[1:N_sp], vcv_vBM);
  target+= normal_lpdf(y| Z * b, sigma_y);
}

 generated quantities {
  vector[N] log_likelihood;

  for(i in 1:N){log_likelihood[i] = normal_lpdf(y[i]| Z[i,] * b, sigma_y);}

  }

  '

postpredict_vBM<- stan_model(model_code = stancode_vBM, verbose = F)
fit_vBM<- sampling(postpredict_vBM, data = dat, chains = cha,
                   warmup = war, iter=ite, init=0,
                   control=list(adapt_delta=0.99))
print(fit_vBM)

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
    matrix[N, N_sp] Z;  // spices id for random effect
}

 parameters {
  real MRCA;
  real<lower=0> ev_base;
  real<lower=0> sigma_y;
  vector[N_sp] b;

}

transformed parameters{
  vector[len_phylo+1] sim_mean;
  vector[len_phylo+1] sim_var;
  
  cov_matrix[N_sp] vcv_BM;

  sim_mean[tree_obj[1,1]]= MRCA;
  sim_var [tree_obj[1,1]]= 0;
  
  
  for(i in 1:len_phylo){
    sim_mean[tree_obj[i,2]] = sim_mean[tree_obj[i,1]];
    sim_var [tree_obj[i,2]] =  sim_var[tree_obj[i,1]] + (2*branch_len[i]*ev_base);
  }
  
  for(i in 1: N_sp){
    for(j in i: N_sp){
      if(i != j){
      vcv_BM[i, j] = sim_var[MRCA_ij[i, j]];
      vcv_BM[j, i] = sim_var[MRCA_ij[i, j]];
      }else{
        vcv_BM[i, j] = sim_var[i];
        vcv_BM[j, i] = sim_var[i];
      }
    }
  }
  
}

model {
  // set the prior
    //target+= uniform_lpdf(ev_base| 0, 10000);

  // set the likelihood
  target+= multi_normal_lpdf(b| sim_mean[1:N_sp], vcv_BM);
  target+= normal_lpdf(y| Z * b, sigma_y);
}

generated quantities {
  vector[N] log_likelihood;

  for(i in 1:N){log_likelihood[i] = normal_lpdf(y[i]| Z[i,] * b, sigma_y);}

}
  '

stan_BM<- stan_model(model_code = stancode_BM, verbose = F)
#options(mc.cores = parallel::detectCores())
fit_BM<- sampling(stan_BM, data = dat, chains = cha,
                  warmup = war, iter=ite, init=0,
                  control=list(adapt_delta=0.95))
print(fit_BM)
