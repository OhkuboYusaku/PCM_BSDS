install.packages(c("ape", "rstan", "dummies"))

# load packages
library(ape)
library(rstan)
library(dummies)

# define data-arrangement function
BSDSLMM_data<- function(phylo, y, X, Z, D_edge){
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
             y=y, X=X, Z=dummy(Z), D_edge=D_edge, D=ncol(X))
  
  return (dat)
}

# load data
data<- read.csv("BSDS_LMM_sample.csv")

Y<- dat$y
X<- as.matrix(data$X)
sp_ID<- (data$sp_ID)

phylo<- read.tree("BSDS_LMM_tree")
plot(phylo)
axisPhylo()
phylo$edge # tree構造:[,1]の親種から[,2]の子孫種へエッジが伸びている

D_edge<- 18

dat<- BSDSLMM_data(phylo, y, X, Z, D_edge = D_edge)

# config. Stan setting
par<- c("MRCA", "beta", "log_likelihood")
scr<-"BSDS_LMM.stan"

war<- 5000
ite<- 25000
cha<- 2
options(mc.cores = parallel::detectCores())

# sampling
fit_BSDS<-stan(file = scr, model_name = scr, data = dat, pars = par, chains = cha, 
          warmup = war, iter = ite, thin = 10, control = list(adapt_delta=0.95))

print(fit_BSDS)

# calculate Widely-Applicable Information Criterion (WAIC)
library(loo)
waic(extract(fit_BSDS)$log_likelihood)
