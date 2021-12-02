BSDE_sim<- function(N_sample, N_tip, MRCA, ev, sel, n_sim=100){
  library(ape)
  library(rstan)
  library(loo)
  library(dummies)

  stan_BSDE_unif<- stan_model(file = "BSDS_noX2.stan")
  stan_vBM_unif<- stan_model(file = "vBM_noX2.stan")
  stan_BM_unif<- stan_model(file = "BM_noX2.stan")

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
  
  result_BSDE<- matrix(0, n_sim, 26)
  result_vBM<- matrix(0, n_sim, 26)
  result_BM<- matrix(0, n_sim, 26)
  
  for(i_sim in 1:n_sim){
    cha<- 1
    
    phylo<- rtree(N_tip)
    len_phylo<- length(phylo$edge.length)
    phylo$edge.length<- phylo$edge.length

    ## import phylogenetic data //
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
    
    k<- numeric(length(phylo$edge.length))
    DE_edge<- floor(runif(1, min=1, max=len_phylo+1))
    
    k[1:length(phylo$edge.length)]<- 1
    k[DE_edge]<- sel
    
    latent<- numeric(N_tip)
    latent<- phylo2ABC(phylo, MRCA, ev, k=(k))
    
    y<- as.numeric(Z%*%latent + rnorm(N_tip*N_sample, 0, 1))
    
    dat<- list(N=length(y), N_sp=N_tip,
               len_phylo=len_phylo, branch_len=branch_len, tree_obj=tree_obj, MRCA_ij=MRCA_ij,
               y=y, Z=(Z), DE_edge=DE_edge)
    
    war<-5000
    ite<-10000
    cha<-1
    
    fit_BSDE<- sampling(stan_BSDE_unif, data = dat, chains = cha, refresh=ite/2,
                        par=c("MRCA", "ev", "sel","sigma_y","b", "log_likelihood"),
                        warmup = war, iter=ite,thin=1, verbose=F, open_progress=F,
                        control = list(adapt_delta=0.9))  
    Rhat_below_1.01<- all(summary(fit_BSDE)$summary[,"Rhat"] <= 1.01, na.rm = T)
    
    # if MCMC not converged, then multiply the iterations and retry
    while(Rhat_below_1.01==FALSE) {
      war<- 2 * war
      ite<- 2 * ite
      fit_BSDE<- sampling(stan_BSDE_unif, data = dat, chains = cha, refresh=ite/2,
                          par=c("MRCA", "ev", "sel","sigma_y","b", "log_likelihood"),
                          warmup = war, iter=ite,thin=1, verbose=F, open_progress=F,
                          control = list(adapt_delta=0.9))  
      Rhat_below_1.01<- all(summary(fit_BSDE)$summary[,"Rhat"] <= 1.01, na.rm = T)
    }
    
    a<- extract(fit_BSDE)
    result_BSDE[i_sim,1]<- summary(fit_BSDE)$summary["MRCA",1]
    result_BSDE[i_sim,2]<- summary(fit_BSDE)$summary["ev",1]
    result_BSDE[i_sim,3]<- summary(fit_BSDE)$summary["sel",1]
    
    result_BSDE[i_sim,4]<- summary(fit_BSDE)$summary["MRCA",6]
    result_BSDE[i_sim,5]<- summary(fit_BSDE)$summary["ev",6]
    result_BSDE[i_sim,6]<- summary(fit_BSDE)$summary["sel",6]
    
    result_BSDE[i_sim,7]<- summary(fit_BSDE)$summary["MRCA",4]
    result_BSDE[i_sim,8]<- summary(fit_BSDE)$summary["ev",4]
    result_BSDE[i_sim,9]<- summary(fit_BSDE)$summary["sel",4]
    
    result_BSDE[i_sim,10]<- summary(fit_BSDE)$summary["MRCA",8]
    result_BSDE[i_sim,11]<- summary(fit_BSDE)$summary["ev",8]
    result_BSDE[i_sim,12]<- summary(fit_BSDE)$summary["sel",8]
    
    if ((result_BSDE[i_sim,7] < MRCA) && (MRCA < result_BSDE[i_sim,10])) {
      result_BSDE[i_sim,13]<- 1
    }else{
      result_BSDE[i_sim,13]<- 0  
    }
    
    if ((result_BSDE[i_sim,8] < ev) && (ev < result_BSDE[i_sim,11])) {
      result_BSDE[i_sim,14]<- 1
    }else{
      result_BSDE[i_sim,14]<- 0  
    }    
    
    if ((result_BSDE[i_sim,9] < log(sel)) && (log(sel) < result_BSDE[i_sim,12])) {
      result_BSDE[i_sim,15]<- 1
    }else{
      result_BSDE[i_sim,15]<- 0  
    }
    
    result_BSDE[i_sim,16]<- map_mcmc(a$sel)
    
    result_BSDE[i_sim,19]<- summary(fit_BSDE)$summary["sigma_y",1]

    result_BSDE[i_sim,21]<- summary(fit_BSDE)$summary["sigma_y",4]

    result_BSDE[i_sim,23]<- summary(fit_BSDE)$summary["sigma_y",8]

    if ((result_BSDE[i_sim,21] < 1) && (1 < result_BSDE[i_sim,23])) {
      result_BSDE[i_sim,25]<- 1
    }else{
      result_BSDE[i_sim,25]<- 0  
    }    
    
    war<-10000
    ite<-20000
    cha<-1
    
    fit_vBM<- sampling(stan_vBM_unif, data = dat, chains = 1, init=0, refresh=ite/2,
                       par=c("MRCA","acc_DE_edge", "ev_base","sigma_y","b", "log_likelihood"),
                       warmup = war, iter=ite,thin=10, verbose=F, open_progress=F,
                       control = list(adapt_delta=0.9))
    Rhat_below_1.01<- all(summary(fit_vBM)$summary[,"Rhat"] <= 1.01, na.rm = T)
    
    # if MCMC not converged, then multiply the iterations and retry
    while(Rhat_below_1.01==FALSE) {
      war<- 2 * war
      ite<- 2 * ite
      fit_vBM<- sampling(stan_vBM_unif, data = dat, chains = cha, refresh=ite/2,
                         par=c("MRCA","acc_DE_edge", "ev_base","sigma_y","b", "log_likelihood"),
                         warmup = war, iter=ite,thin=10, verbose=F, open_progress=F, refresh=0, 
                         control = list(adapt_delta=0.9))
      Rhat_below_1.01<- all(summary(fit_vBM)$summary[,"Rhat"] <= 1.01, na.rm = T)
    }    
    
    a<- extract(fit_vBM)
    
    result_vBM[i_sim,1]<- summary(fit_vBM)$summary["MRCA",1]
    result_vBM[i_sim,2]<- summary(fit_vBM)$summary["ev_base",1]
    result_vBM[i_sim,3]<- summary(fit_vBM)$summary["acc_DE_edge",1]
    
    result_vBM[i_sim,4]<- summary(fit_vBM)$summary["MRCA",6]
    result_vBM[i_sim,5]<- summary(fit_vBM)$summary["ev_base",6]
    result_vBM[i_sim,6]<- summary(fit_vBM)$summary["acc_DE_edge",6]
    
    result_vBM[i_sim,7]<- summary(fit_vBM)$summary["MRCA",4]
    result_vBM[i_sim,8]<- summary(fit_vBM)$summary["ev_base",4]
    result_vBM[i_sim,9]<- summary(fit_vBM)$summary["acc_DE_edge",4]
    
    result_vBM[i_sim,10]<- summary(fit_vBM)$summary["MRCA",8]
    result_vBM[i_sim,11]<- summary(fit_vBM)$summary["ev_base",8]
    result_vBM[i_sim,12]<- summary(fit_vBM)$summary["acc_DE_edge",8]
    
    if ((result_vBM[i_sim,7] < MRCA) && (MRCA < result_vBM[i_sim,10])) {
      result_vBM[i_sim,13]<- 1
    }else{
      result_vBM[i_sim,13]<- 0  
    }
    
    if ((result_vBM[i_sim,8] < ev) && (ev < result_vBM[i_sim,11])) {
      result_vBM[i_sim,14]<- 1
    }else{
      result_vBM[i_sim,14]<- 0  
    }    
    
    result_vBM[i_sim,19]<- summary(fit_vBM)$summary["sigma_y",1]
    
    result_vBM[i_sim,21]<- summary(fit_vBM)$summary["sigma_y",4]
    
    result_vBM[i_sim,23]<- summary(fit_vBM)$summary["sigma_y",8]
    
    if ((result_vBM[i_sim,21] < 1) && (1 < result_vBM[i_sim,23])) {
      result_vBM[i_sim,25]<- 1
    }else{
      result_vBM[i_sim,25]<- 0  
    }    
    
    war<-10000
    ite<-20000
    
    fit_BM<- sampling(stan_BM_unif, data = dat, chains = cha, refresh=ite/2,
                      par=c("MRCA", "ev_base","sigma_y","b", "log_likelihood"),
                      warmup = war, iter=ite,thin=10, verbose=F, open_progress=F, 
                      control = list(adapt_delta=0.9))
    Rhat_below_1.01<- all(summary(fit_BM)$summary[,"Rhat"] <= 1.01, na.rm = T)
    
    # if MCMC not converged, then multiply the iterations and retry
    while(Rhat_below_1.01==FALSE) {
      war<- 2 * war
      ite<- 2 * ite
      fit_BM<- sampling(stan_BM_unif, data = dat, chains = cha, refresh=ite/2,
                        par=c("MRCA", "ev_base","sigma_y","b", "log_likelihood"),
                        warmup = war, iter=ite,thin=10, verbose=F, open_progress=F, 
                        control = list(adapt_delta=0.9))
      Rhat_below_1.01<- all(summary(fit_BM)$summary[,"Rhat"] <= 1.01, na.rm = T)
    }
    
    a<- extract(fit_BM)
    
    result_BM[i_sim,1]<- summary(fit_BM)$summary["MRCA",1]
    result_BM[i_sim,2]<- summary(fit_BM)$summary["ev_base",1]

    result_BM[i_sim,4]<- summary(fit_BM)$summary["MRCA",6]
    result_BM[i_sim,5]<- summary(fit_BM)$summary["ev_base",6]

    result_BM[i_sim,7]<- summary(fit_BM)$summary["MRCA",4]
    result_BM[i_sim,8]<- summary(fit_BM)$summary["ev_base",4]

    result_BM[i_sim,10]<- summary(fit_BM)$summary["MRCA",8]
    result_BM[i_sim,11]<- summary(fit_BM)$summary["ev_base",8]

    if ((result_BM[i_sim,7] < MRCA) && (MRCA < result_BM[i_sim,10])) {
      result_BM[i_sim,13]<- 1
    }else{
      result_BM[i_sim,13]<- 0  
    }
    
    if ((result_BM[i_sim,8] < ev) && (ev < result_BM[i_sim,11])) {
      result_BM[i_sim,14]<- 1
    }else{
      result_BM[i_sim,14]<- 0  
    }    
    

    result_BM[i_sim,19]<- summary(fit_BM)$summary["sigma_y",1]

    result_BM[i_sim,21]<- summary(fit_BM)$summary["sigma_y",4]

    result_BM[i_sim,23]<- summary(fit_BM)$summary["sigma_y",8]

    if ((result_BM[i_sim,21] < 1) && (1 < result_BM[i_sim,23])) {
      result_BM[i_sim,25]<- 1
    }else{
      result_BM[i_sim,25]<- 0  
    }    
    
    print(i_sim)
  }
  
  dat1 = data.frame(
    mean_MRCA = result_BSDE[,1], mean_ev = result_BSDE[,2], mean_k = result_BSDE[,3],
    mean_sigma_y = result_BSDE[,19], 
    med_MRCA = result_BSDE[,4],  med_ev = result_BSDE[,5], med_k = result_BSDE[,6], 
    MAP_k = result_BSDE[,16], lo_sel = result_BSDE[,9], up_sel = result_BSDE[,12],

    len_MRCA = result_BSDE[,10] -result_BSDE[,7],
    len_ev = result_BSDE[,11] -result_BSDE[,8],
    len_sel = result_BSDE[,12]- result_BSDE[,9],
    len_sigma_y = result_BSDE[,23]- result_BSDE[,21],

    CI_MRCA = result_BSDE[,13],  CI_ev = result_BSDE[,14], CI_k = result_BSDE[,15],
    CI_sigma_y = result_BSDE[,25], 
    #########################################    #########################################
    vBM_mean_MRCA = result_vBM[,1], vBM_mean_ev = result_vBM[,2], vBM_mean_acc = result_vBM[,3],
    vBM_mean_sigma_y = result_vBM[,19], 
    vBM_med_MRCA = result_vBM[,4],  vBM_med_ev = result_vBM[,5], vBM_med_acc = result_vBM[,6], 
    vBM_MAP_acc = result_vBM[,16], vBM_lo_acc = result_vBM[,9], vBM_up_acc = result_vBM[,12],
    
    vBM_len_MRCA = result_vBM[,10] -result_vBM[,7],
    vBM_len_ev = result_vBM[,11] -result_vBM[,8],
    vBM_len_acc = result_vBM[,12]- result_vBM[,9],
    vBM_len_sigma_y = result_vBM[,23]- result_vBM[,21],
    
    vBM_CI_MRCA = result_vBM[,13],  vBM_CI_ev = result_vBM[,14],
    vBM_CI_sigma_y = result_vBM[,25], 
    #########################################    #########################################
    BM_mean_MRCA = result_BM[,1], BM_mean_ev = result_BM[,2], 
    BM_mean_sigma_y = result_BM[,19], 
    BM_med_MRCA = result_BM[,4],  BM_med_ev = result_BM[,5],

    BM_len_MRCA = result_BM[,10] -result_BM[,7],
    BM_len_ev = result_BM[,11] -result_BM[,8],
    BM_len_sigma_y = result_BM[,23]- result_BM[,21],

    BM_CI_MRCA = result_BM[,13],  BM_CI_ev = result_BM[,14], 
    BM_CI_sigma_y = result_BM[,25]
    
    #########################################    #########################################
  )
  
  file.name1 <- sprintf("2021onepath.comp%.1f_%.1f_BSDE_noX2_%.1f_%.1f_%.1f.csv", N_sample, N_tip, MRCA, sel, ev)

  write.csv(dat1, file.name1)
  #print(file.name1)

  return(dat1)
  
}

    
exp1<- BSDE_sim(N_sample=10, N_tip=50, MRCA=100, ev=100, sel=3, n_sim=1000)