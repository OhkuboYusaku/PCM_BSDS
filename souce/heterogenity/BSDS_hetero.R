BSDE_onepath<- function(N_sample, N_tip, MRCA, ev, sel, n_sim=100){
  library(ape)
  library(rstan)
  library(loo)
  library(dummies)
  
  stan_BSDE_unif2<- stan_model(file = "BSDS_noX_hetero_cauchy.stan")
  
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
  
  result_BSDE<- matrix(0, n_sim, 24)

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
    
    #Z<- matrix(0, N_tip*sample, N_tip)
    #for(j in 1:N_tip){
    #  for(i in 1:sample){
    #    Z[(j-1)*sample + i, j]<- 1 
    #  }
    #}
    
    k<- numeric(length(phylo$edge.length))
    DE_edge<- floor(runif(1, min=1, max=len_phylo+1))
    
    k[1:length(phylo$edge.length)]<- 1
    k[DE_edge]<- sel
    
    latent<- numeric(N_tip)
    latent<- phylo2ABC(phylo, MRCA, ev, k=exp(k))
    
    #sigma<- rgamma(N_sp, 3, 1) 
    sigma<- abs(rcauchy(N_sp, 0, 2))
    
    y<- as.numeric(Z%*%latent + rnorm(N_sample*N_sp, 0, Z%*%sqrt(sigma)))
    
    dat<- list(N=length(y), N_sp=N_tip,
               len_phylo=len_phylo, branch_len=branch_len, tree_obj=tree_obj, MRCA_ij=MRCA_ij,
               y=y, Z=(Z), DE_edge=DE_edge)
    
    war<-5000
    ite<-25000
    cha<-1
    
    #fit_BSDE<- sampling(stan_BSDE_unif1, data = dat, chains = cha, 
    #                    par=c("MRCA", "ev", "sel","sigma_y", "shape", "rate", "log_likelihood"),
    #                    warmup = war, iter=ite,thin=1, verbose=F, open_progress=F,
    #                    control = list(adapt_delta=0.9))  
    fit_BSDE<- sampling(stan_BSDE_unif2, data = dat, chains = cha, 
                        par=c("MRCA", "ev", "sel","sigma_y", "rate", "log_likelihood"),
                        warmup = war, iter=ite,thin=1, verbose=F, open_progress=F,
                        control = list(adapt_delta=0.9))  
    
        Rhat_below_1.01<- all(summary(fit_BSDE)$summary[,"Rhat"] <= 1.01, na.rm = T)
    
    # if MCMC not converged, then multiply the iterations and retry
    while(Rhat_below_1.01==FALSE) {
      war<- 2 * war
      ite<- 2 * ite
      fit_BSDE<- sampling(stan_BSDE_unif2, data = dat, chains = cha, 
                          par=c("MRCA", "ev", "sel","", "rate", "log_likelihood"),
                          warmup = war, iter=ite,thin=1, verbose=F, open_progress=F,refresh=0,
                          control = list(adapt_delta=0.9))  
      Rhat_below_1.01<- all(summary(fit_BSDE)$summary[,"Rhat"] <= 1.01, na.rm = T)
    }
    
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
    
    #result_BSDE[i_sim,16]<- map_mcmc(a$sel)

    #result_BSDE[i_sim,17]<- summary(fit_BSDE)$summary["shape",1]
    result_BSDE[i_sim,18]<- summary(fit_BSDE)$summary["rate",1]

    #result_BSDE[i_sim,19]<- summary(fit_BSDE)$summary["shape",4]
    result_BSDE[i_sim,20]<- summary(fit_BSDE)$summary["rate",4]
    
    #result_BSDE[i_sim,21]<- summary(fit_BSDE)$summary["shape",8]
    result_BSDE[i_sim,22]<- summary(fit_BSDE)$summary["rate",8]
    
    #if ((result_BSDE[i_sim,19] < 3) && (3 < result_BSDE[i_sim,21])) {
    #  result_BSDE[i_sim,23]<- 1
    #}else{
    #  result_BSDE[i_sim,23]<- 0  
    #}
    
    if ((result_BSDE[i_sim,20] < 1) && (1 < result_BSDE[i_sim,22])) {
      result_BSDE[i_sim,24]<- 1
    }else{
      result_BSDE[i_sim,24]<- 0  
    }
    print(i_sim)
  }
  
  dat1 = data.frame(
    mean_MRCA = result_BSDE[,1], mean_ev = result_BSDE[,2], mean_k = result_BSDE[,3], 
    med_MRCA = result_BSDE[,4],  med_ev = result_BSDE[,5], med_k = result_BSDE[,6], 
    #MAP_k = result_BSDE[,16], 
    lo_sel = result_BSDE[,9], up_sel = result_BSDE[,12],
    #mean_shape = result_BSDE[,17], 
    mean_rate = result_BSDE[,18],
    
    len_MRCA = result_BSDE[,10] -result_BSDE[,7],
    len_ev = result_BSDE[,11] -result_BSDE[,8],
    len_sel = result_BSDE[,12]- result_BSDE[,9],
    
    CI_MRCA = result_BSDE[,13], CI_ev = result_BSDE[,14], CI_k = result_BSDE[,15],
    #CI_shape = result_BSDE[,23], 
    CI_rate = result_BSDE[,24]
  )
  
  file.name1 <- sprintf("Onepath%.1f_%.1f_BSDE_hetero_%.1f__%.1f_%.1f.csv", N_sample, N_tip, MRCA, sel, ev)

  write.csv(dat1, file.name1)
  #print(file.name1)

  return(dat1)
  
}


# example 
N_sample<- 10
N_tip<- 10
sel<- 2

#exp1<- BSDE_onepath(N_sample=N_sample, N_tip=N_tip, MRCA=100, ev=10, sel=sel, n_sim=2)

exp1<- BSDE_onepath(N_sample=N_sample, N_tip=N_tip, MRCA=100, ev=10, sel=sel, n_sim=100)
exp2<- BSDE_onepath(N_sample=N_sample, N_tip=N_tip, MRCA=100, ev=50, sel=sel, n_sim=100)
#exp3<- BSDE_MLE_sim(N_sample=N_sample, N_sp=N_sp, MRCA=100, ev_rate=100, sel=sel, n_sim=1000)
#exp4<- BSDE2_stan_sim(N_sample=N_sample, N_sp=N_sp, MRCA=100, ev_rate=500, sel=sel, n_sim=1000)
#exp5<- BSDE2_stan_sim(N_sample=N_sample, N_sp=N_sp, MRCA=100, ev_rate=1000, sel=sel, n_sim=1000)
