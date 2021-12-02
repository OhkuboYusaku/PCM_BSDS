setwd("~/Google ドライブ/multi_noX")

BSDE2_stan_sim<-function(N_sample, N_sp, MRCA, ev, sel, n_sim){
  library(ape)
  library(rstan)
  library(loo)
  # define stan model
  stan_BSDE_unif<- stan_model(file = "BSDS_LMM_noX.stan")
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
  
  map_mcmc <- function(z){ density(z)$x[which.max(density(z)$y)]} 
  
  # end definition of R-functions
  
  sim_date<- as.Date(Sys.time()) #extract the experiment date
  result_BSDE<- matrix(NA, n_sim, 30) # results container
 
  for(i_sim in 1:n_sim){
    ##### generate data ####
    phylo<- rtree(N_sp) # random sampling of phylogenetic tree
    len_phylo<- length(phylo$edge.length)
    phylo$edge.length<- phylo$edge.length
    N_tip<- len_phylo - phylo$Nnode + 1
    
    k<- numeric(length(phylo$edge.length))
    DE_edge<- sample(1:len_phylo, 2)
    
    k[1:length(phylo$edge.length)]<- 1
    k[DE_edge]<- sel
    
    branch_len<- phylo$edge.length
    tree_obj<- as.matrix(phylo$edge)
    MRCA_ij<- matrix(0, N_tip, N_tip) ## i,j elements correspond to the location of their MRCA in the tree
    
    for(i in 1:N_tip){
      for(j in i:N_tip){
        MRCA_ij[i,j]<- MRCA_ij[j,i]<- getMRCA(phylo, tip=c(i,j))
      }
    }
    
    Z<- matrix(0, N_tip*N_sample, N_tip)
    for(j in 1:N_tip){
      for(i in 1:N_sample){
        Z[(j-1)*N_sample + i, j]<- 1 
      }
    }
    
    latent<- numeric(N_tip)
    latent<- phylo2ABC(phylo, MRCA, ev, k=(k))
    
    y<- as.numeric(Z%*%latent + rnorm(N_tip*N_sample, 0, 1))
    
    dat<- list(N=length(y), N_sp=N_tip,
               len_phylo=len_phylo, branch_len=branch_len, tree_obj=tree_obj, MRCA_ij=MRCA_ij,
               y=y, Z=(Z), DE_edge=DE_edge)
    
    #### conduct Bayesian inference ###

    war<- 5000
    ite<- 10000
    cha<- 2
    
    # BSDE model
    par = c("MRCA",  "ev", "sel1", "sel2", "sigma_y")
    Rhat<- FALSE
    fit_BSDE<- sampling(stan_BSDE_unif, data = dat, chains = cha, pars = par,
                        warmup = war, iter=ite,thin=1, verbose=F, open_progress=F,
                        refresh = 0, control = list(adapt_delta=0.9))  
    Rhat_below_1.05<- all(summary(fit_BSDE)$summary[,"Rhat"] <= 1.01, na.rm = T)
    
    # if MCMC not converged, then multiply the iterations and retry
      while(Rhat_below_1.05==FALSE) {
        war<- 2 * war
        ite<- 2 * ite
        fit_BSDE<- sampling(stan_BSDE_unif, data = dat, chains = cha, pars = par,
                            warmup = war, iter=ite,thin=1, verbose=F, open_progress=F,
                            refresh = 0, control = list(adapt_delta=0.9))  
        Rhat_below_1.05<- all(summary(fit_BSDE)$summary[,"Rhat"] <= 1.01, na.rm = T)
      }
    
    a<- extract(fit_BSDE)
    
    result_BSDE[i_sim,1]<- summary(fit_BSDE)$summary["MRCA",1]
    result_BSDE[i_sim,2]<- summary(fit_BSDE)$summary["ev",1]
    result_BSDE[i_sim,3]<- summary(fit_BSDE)$summary["sel1",1]
    result_BSDE[i_sim,4]<- summary(fit_BSDE)$summary["sel2",1]
    
    result_BSDE[i_sim,5]<- summary(fit_BSDE)$summary["MRCA",6]
    result_BSDE[i_sim,6]<- summary(fit_BSDE)$summary["ev",6]
    result_BSDE[i_sim,7]<- summary(fit_BSDE)$summary["sel1",6]
    result_BSDE[i_sim,8]<- summary(fit_BSDE)$summary["sel2",6]
    
    result_BSDE[i_sim,9]<- summary(fit_BSDE)$summary["MRCA",4]
    result_BSDE[i_sim,10]<- summary(fit_BSDE)$summary["ev",4]
    result_BSDE[i_sim,11]<- summary(fit_BSDE)$summary["sel1",4]
    result_BSDE[i_sim,12]<- summary(fit_BSDE)$summary["sel2",4]
    
    result_BSDE[i_sim,13]<- summary(fit_BSDE)$summary["MRCA",8]
    result_BSDE[i_sim,14]<- summary(fit_BSDE)$summary["ev",8]
    result_BSDE[i_sim,15]<- summary(fit_BSDE)$summary["sel1",8]
    result_BSDE[i_sim,16]<- summary(fit_BSDE)$summary["sel2",8]
    
    if ((result_BSDE[i_sim,9] < MRCA) && (MRCA < result_BSDE[i_sim,13])) {
      result_BSDE[i_sim,17]<- 1
    }else{
      result_BSDE[i_sim,17]<- 0  
    }
    
    if ((result_BSDE[i_sim,10] < ev) && (ev < result_BSDE[i_sim,14])) {
      result_BSDE[i_sim,18]<- 1
    }else{
      result_BSDE[i_sim,18]<- 0  
    }    
    
    if ((result_BSDE[i_sim,11] < log(sel[1])) && (log(sel[1]) < result_BSDE[i_sim,15])) {
      result_BSDE[i_sim,19]<- 1
    }else{
      result_BSDE[i_sim,19]<- 0  
    }
    
    if ((result_BSDE[i_sim,12] < log(sel[2])) && (log(sel[2]) < result_BSDE[i_sim,16])) {
      result_BSDE[i_sim,20]<- 1
    }else{
      result_BSDE[i_sim,20]<- 0  
    }
    
    result_BSDE[i_sim,21]<- map_mcmc(a$sel1)
    result_BSDE[i_sim,22]<- map_mcmc(a$sel2)
    
    result_BSDE[i_sim,23]<- summary(fit_BSDE)$summary["sigma_y",1]
    
    result_BSDE[i_sim,24]<- summary(fit_BSDE)$summary["sigma_y",4]
    
    result_BSDE[i_sim,25]<- summary(fit_BSDE)$summary["sigma_y",8]
    
    if ((result_BSDE[i_sim,24] < 1) && (1 < result_BSDE[i_sim,25])) {
      result_BSDE[i_sim,26]<- 1
    }else{
      result_BSDE[i_sim,26]<- 0  
    }    
    
    
    print(i_sim)
  }
  
  dat1 = data.frame(
    mean_MRCA = result_BSDE[,1], mean_ev = result_BSDE[,2], 
    mean_k1 = result_BSDE[,3], mean_k2 = result_BSDE[,4], mean_sigma_y = result_BSDE[,23], 
    
    med_MRCA = result_BSDE[,5],  med_ev = result_BSDE[,6], 
    med_k1 = result_BSDE[,7], med_k2 = result_BSDE[,8], 

    MAP_k1 = result_BSDE[,21], lo_sel = result_BSDE[,11], up_sel = result_BSDE[,15],
    MAP_k2 = result_BSDE[,22], lo_sel = result_BSDE[,12], up_sel = result_BSDE[,16],
    
    len_MRCA = result_BSDE[,13] -result_BSDE[,9],
    len_ev = result_BSDE[,14] -result_BSDE[,10],
    len_sel1 = result_BSDE[,15]- result_BSDE[,11],
    len_sel2 = result_BSDE[,16]- result_BSDE[,12],
    
    len_sigma_y = result_BSDE[,25]- result_BSDE[,24],
    
    CI_MRCA = result_BSDE[,17],  CI_ev = result_BSDE[,18], 
    CI_k1 = result_BSDE[,19], CI_k1 = result_BSDE[,20],
    CI_sigma_y = result_BSDE[,26] 

  )
  
  file.name1 <- sprintf("%.1f_%.1f_sel1&sel2_%.1f_%.2f_%.2f_%.1f.csv",N_sample, N_sp, MRCA, sel[1], sel[2], ev)
  
  write.csv(dat1, file.name1)
  #print(file.name1)
  
  return(dat1)
  
}

# example 
N_sample<- 10
N_sp<- 50
sel<- c(1.5, 2)

#exp1<- BSDE2_stan_sim(N_sample=N_sample, N_sp=N_sp, MRCA=100, ev_rate=10, sel=sel, n_sim=1000)
#exp2<- BSDE2_stan_sim(N_sample=N_sample, N_sp=N_sp, MRCA=100, ev_rate=50, sel=sel , n_sim=1000)
exp3<- BSDE2_stan_sim(N_sample=N_sample, N_sp=N_sp, MRCA=100, ev=10, sel=sel, n_sim=1000)
#exp4<- BSDE2_stan_sim(N_sample=N_sample, N_sp=N_sp, MRCA=100, ev_rate=500, sel=sel, n_sim=1000)
exp5<- BSDE2_stan_sim(N_sample=N_sample, N_sp=N_sp, MRCA=100, ev=50, sel=sel, n_sim=1000)
