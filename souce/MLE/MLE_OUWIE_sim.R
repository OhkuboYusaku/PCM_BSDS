OU_sim<- function(N_sample, N_sp, MRCA, ev, sel, n_sim=100){
  
library(ape)
library(OUwie)
library(mvnfast)
  
  loglik<- function(phylo, theta, y, DS_edge){
    MRCA<- theta[1]
    ev<- (theta[2])
    k<- numeric(length(phylo$edge.length)) + 1
    for(i_DS in 1:length(DS_edge)){
      k[DS_edge[i_DS]]<- exp(theta[3])
    }
    
    #define data
    tree_obj<- as.matrix(phylo[["edge"]])
    simulated_mean_trait<- numeric(length(phylo$edge.length))
    simulated_var_trait<- numeric(length(phylo$edge.length))
    len_phylo<- length(phylo$edge.length)
    ###### NOTE ######
    # simulated_mean_trait[i] and simulated_var_trait[i] contains simulated trait of sp i.
    
    ## conduct evolution simulation
    simulated_mean_trait[tree_obj[1,1]]<- MRCA
    simulated_var_trait[tree_obj[1,1]]<- 0
    
    for(i in 1:len_phylo){
      branch_len=phylo$edge.length[i]
      simulated_mean_trait[tree_obj[i,2]]<- simulated_mean_trait[tree_obj[i,1]] + (branch_len*ev*(k[i]^2-1))/k[i]
      simulated_var_trait[tree_obj[i,2]]<- simulated_var_trait[tree_obj[i,1]] + (2*branch_len*ev*(k[i]^2+1))/k[i]
    }
    
    # output vcov matrix ot tip sps.
    N_tip<- len_phylo - phylo$Nnode +1
    vcv_BSDE<- matrix(0, N_tip, N_tip)
    
    for(i in 1: N_tip){
      for(j in i: N_tip){
        if(i != j){
          MRCA_ij<- getMRCA(phylo, tip=c(i,j))
          vcv_BSDE[i, j]<- vcv_BSDE[j, i]<- simulated_var_trait[MRCA_ij]
        }else{
          vcv_BSDE[i, j]<- vcv_BSDE[j, i]<- simulated_var_trait[i]
        }
      }
    }
    
    if((sum(is.na(vcv_BSDE))>0)||(min(eigen(vcv_BSDE)$values)<=0)){return(-Inf)}
    loglik<- try(dmvn(y, simulated_mean_trait[1:N_sp],  Matrix::nearPD(vcv_BSDE)$mat, log=T))
    if(is.numeric(loglik)==F){return(-10000)}
    
    return(loglik)
  }
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


sim_date<- as.Date(Sys.time()) #extract the experiment date
result_OUMA<- matrix(NA, n_sim, 30) # results container
result_OU1<- matrix(NA, n_sim, 30) # results container
result_BMS<- matrix(NA, n_sim, 30) # results container
result_BM1<- matrix(NA, n_sim, 30) # results container
result_BSDS<- matrix(NA, n_sim, 30) # results container

for(i_sim in 1:n_sim){
  ##### generate data ####
  phylo<- rtree(N_sp) # random sampling of phylogenetic tree
  len_phylo<- length(phylo$edge.length)
  phylo$edge.length<- phylo$edge.length
  N_tip<- len_phylo - phylo$Nnode + 1
  
  Z1<- numeric(N_tip*N_sample)
  for(j in 1:N_tip){
    for(k in 1:N_sample){
      Z1[(j-1)*N_sample + k]<- j
    }
  }
  
  Z<- matrix(0, N_tip*N_sample, N_tip)
  for(j in 1:N_tip){
    for(i in 1:N_sample){
      Z[(j-1)*N_sample + i, j]<- 1 
    }
  }
  
  k<- numeric(length(phylo$edge.length))
  DE_edge<- sample(1:len_phylo, 1)
  
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
  
  hist_phylo<- numeric(len_phylo)
  for(i in 1:len_phylo){
    hist_phylo[tree_obj[i,2]] = hist_phylo[tree_obj[i,1]] +branch_len[i]
  }
  
  latent<- numeric(N_tip)
  latent<- phylo2ABC(phylo, MRCA, ev, k=(k))
  
  y<- as.numeric(Z%*%latent + rnorm(N_tip*N_sample, 0, 1))
  model<- lm(y~as.factor(Z1)+0)
  
  sp<- rep(1:N_sp)
  ans.DE<- phylo$edge[DE_edge,1]
  dec.DE<- phylo$edge[DE_edge,2]
  
  Reg<- rep(1, N_sp)
  
  if(dec.DE<N_sp+1){
    Reg[dec.DE]<- 2
    phylo$node.label<- rep(1, phylo$Nnode)
    phylo$node.label[ans.DE-N_sp]<- 2
  }else{
    phylo$node.label<- rep(1, phylo$Nnode)
    phylo$node.label[ans.DE-N_sp]<- 2
    phylo$node.label[dec.DE-N_sp]<- 2
  }
  
  dat1<- data.frame(Genus_species=phylo$tip.label)
  dat1[,2]<- data.frame(Reg=Reg)
  dat1[,3]<- data.frame(X=model$coefficients)
  
  # BSDE model
  ########  ########  ########  ########  ########  ########  ########
  BSDS_MLE<- try(optim(par=c(MRCA, ev, log(sel)), fn=loglik, gr = NULL,
                       phylo=phylo, y=model$coefficients, DS_edge=DE_edge,
                       method = "BFGS",
                       control=list(fnscale=-1, trace=0, maxit=1000000), hessian=T))
  
  if (class(BSDS_MLE) == "try-error") {
    result_BSDS[i_sim,26]<-　1
    
    BSDS_MLE<- (optim(par=c(MRCA, ev, log(sel)), fn=loglik, gr = NULL,
                         phylo=phylo, y=model$coefficients, DS_edge=DE_edge,
                         method = "SANN",
                         control=list(fnscale=-1, trace=1000, maxit=100000), hessian=T))
    
    BSDS_MLE<- (optim(par=BSDS_MLE$par, fn=loglik, gr = NULL,
                      phylo=phylo, y=model$coefficients, DS_edge=DE_edge,
                      method = "BFGS",
                      control=list(fnscale=-1, trace=1000, maxit=1000000), hessian=T))
      }
  
  result_BSDS[i_sim,25]<- BSDS_MLE$convergence
  
  
  #print(BSDE_MLE$value)

  result_BSDS[i_sim,1:3]<- BSDS_MLE$par[1:3] #MLE
  result_BSDS[i_sim,4:6]<- sqrt(diag(solve(-BSDS_MLE$hessian))) # SE
  result_BSDS[i_sim,7:9]<- BSDS_MLE$par[1:3]/result_BSDS[i_sim,4:6] # Z-value

  if(result_BSDS[i_sim,4]!="NaN"){
    result_BSDS[i_sim,11]<- result_BSDS[i_sim,1] - 1.96*result_BSDS[i_sim,4]#CI of k
    result_BSDS[i_sim,12]<- result_BSDS[i_sim,1] + 1.96*result_BSDS[i_sim,4]
    if ((result_BSDS[i_sim,11] < MRCA) && (MRCA < result_BSDS[i_sim,12])) {
      result_BSDS[i_sim,13]<- 1
    }else{
      result_BSDS[i_sim,13]<- 0  
    }
  }
  
  
  if(result_BSDS[i_sim,5]!="NaN"){
    result_BSDS[i_sim,14]<- result_BSDS[i_sim,2] - 1.96*result_BSDS[i_sim,5]#CI of k
    result_BSDS[i_sim,15]<- result_BSDS[i_sim,2] + 1.96*result_BSDS[i_sim,5]
    if ((result_BSDS[i_sim,14] < ev) && (ev < result_BSDS[i_sim,15])) {
      result_BSDS[i_sim,16]<- 1
    }else{
      result_BSDS[i_sim,16]<- 0  
    }
  }
  
  if(result_BSDS[i_sim,6]!="NaN"){
    result_BSDS[i_sim,17]<- result_BSDS[i_sim,3] - 1.96*result_BSDS[i_sim,6]#CI of k
    result_BSDS[i_sim,18]<- result_BSDS[i_sim,3] + 1.96*result_BSDS[i_sim,6]
    if ((result_BSDS[i_sim,17] < log(sel)) && (log(sel) < result_BSDS[i_sim,18])) {
      result_BSDS[i_sim,19]<- 1
    }else{
      result_BSDS[i_sim,19]<- 0  
    }
  }
  
  result_BSDS[i_sim,20]<- -BSDS_MLE$value 
  #result_BSDE[i_sim,21:23]<- eigen(-BSDE_MLE$hessian)$values[1:3]
  result_BSDS[i_sim,24]<-　det(-BSDS_MLE$hessian)
  
  ########  ########  ########  ########  ########  ########  ########
  OUMA<- try(OUwie(phylo, dat1, model=c("OUMA"), algorithm="invert", get.root.theta=T,
                   starting.vals=c(0.5,20),scaleHeight=FALSE, 
                   check.identify=T,warn=F,quiet=T, root.age=max(hist_phylo)))
  
  #mean
  result_OUMA[i_sim,1]<- OUMA$solution[1,1] #alpha1
  result_OUMA[i_sim,2]<- OUMA$solution[1,2] #alpha2
  result_OUMA[i_sim,3]<- OUMA$solution[2,1] #sig1
  result_OUMA[i_sim,4]<- OUMA$solution[2,2] #sig2
  result_OUMA[i_sim,5]<- OUMA$theta[1,1] # theta_MRCA
  result_OUMA[i_sim,6]<- OUMA$theta[2,1] # theta_sel
  result_OUMA[i_sim,7]<- OUMA$theta[3,1] # theta_BM
  
  # SE
  result_OUMA[i_sim,8]<- OUMA$theta[1,2] # theta_MRCA
  result_OUMA[i_sim,9]<- OUMA$theta[2,2] # theta_sel
  result_OUMA[i_sim,10]<- OUMA$theta[3,2] # theta_BM
  
  #CI 
  if(is.nan(result_OUMA[i_sim,8])==F)result_OUMA[i_sim,11]<- try(result_OUMA[i_sim,5] - 1.96*result_OUMA[i_sim,8])
  if(is.nan(result_OUMA[i_sim,9])==F) result_OUMA[i_sim,12]<- try(result_OUMA[i_sim,6] - 1.96*result_OUMA[i_sim,9])
  if(is.nan(result_OUMA[i_sim,10])==F)result_OUMA[i_sim,13]<- try(result_OUMA[i_sim,7] - 1.96*result_OUMA[i_sim,10])
  
  result_OUMA[i_sim,14]<- result_OUMA[i_sim,5] + 1.96*result_OUMA[i_sim,8]
  result_OUMA[i_sim,15]<- result_OUMA[i_sim,6] + 1.96*result_OUMA[i_sim,9]
  result_OUMA[i_sim,16]<- result_OUMA[i_sim,7] + 1.96*result_OUMA[i_sim,10]
  
  if(is.nan(result_OUMA[i_sim,8])==F){
    if ((result_OUMA[i_sim,11] < MRCA) && (MRCA < result_OUMA[i_sim,14])) {
      result_OUMA[i_sim,17]<- 1
    }else{
      result_OUMA[i_sim,17]<- 0  
    }
  }

  
  if(is.nan(result_OUMA[i_sim,10])==F){
    if ((result_OUMA[i_sim,13] < MRCA) && (MRCA < result_OUMA[i_sim,16])) {
      result_OUMA[i_sim,18]<- 1
    }else{
      result_OUMA[i_sim,18]<- 0  
    }
  }

  result_OUMA[i_sim,19]<- OUMA$loglik
  
  ########  ########  ########  ########  ########  ########  ########
  OU1<- try(OUwie(phylo, dat1, model=c("OU1"), algorithm="invert", get.root.theta=T,
                  starting.vals=c(0.5,20), scaleHeight=FALSE, 
                  check.identify=T,warn=F,quiet=T, root.age=max(hist_phylo)))
  
  #mean
  result_OU1[i_sim,1]<- OU1$solution[1,1] #alpha1
  result_OU1[i_sim,3]<- OU1$solution[2,1] #sig1
  result_OU1[i_sim,5]<- OU1$theta[1,1] # theta_MRCA
  result_OU1[i_sim,6]<- OU1$theta[2,1] # theta_other

  # SE
  result_OU1[i_sim,8]<- OU1$theta[1,2] # theta_MRCA
  result_OU1[i_sim,9]<- OU1$theta[2,2] # theta_other

  #CI 
  if(is.nan(result_OU1[i_sim,8])==F)result_OU1[i_sim,11]<- result_OU1[i_sim,5] - 1.96*result_OU1[i_sim,8]
  if(is.nan(result_OU1[i_sim,9])==F)result_OU1[i_sim,12]<- result_OU1[i_sim,6] - 1.96*result_OU1[i_sim,9]

  result_OU1[i_sim,14]<- result_OU1[i_sim,5] + 1.96*result_OU1[i_sim,8]
  result_OU1[i_sim,15]<- result_OU1[i_sim,6] + 1.96*result_OU1[i_sim,9]

  if(is.nan(result_OU1[i_sim,8])==F){
    if ((result_OU1[i_sim,11] < MRCA) && (MRCA < result_OU1[i_sim,14])) {
      result_OU1[i_sim,17]<- 1
    }else{
      result_OU1[i_sim,17]<- 0  
    }
  }

  
  if(is.nan(result_OU1[i_sim,9])==F){
    if ((result_OU1[i_sim,12] < MRCA) && (MRCA < result_OU1[i_sim,15])) {
      result_OU1[i_sim,18]<- 1
    }else{
      result_OU1[i_sim,18]<- 0  
    }
  }

  result_OU1[i_sim,19]<- OU1$loglik
  
  ########  ########  ########  ########  ########  ########  ########
  BMS<- try(OUwie(phylo, dat1, model=c("BMS"), algorithm="invert", get.root.theta=T,
                  starting.vals=c(0.5), scaleHeight=FALSE, 
                  check.identify=T,warn=F,quiet=T, root.age=max(hist_phylo)))
  
  #mean
  result_BMS[i_sim,1]<- BMS$solution[2,1] #sig1
  result_BMS[i_sim,3]<- BMS$solution[2,2] #sig2
  result_BMS[i_sim,5]<- BMS$theta[1,1] # theta_MRCA

  # SE
  result_BMS[i_sim,8]<- BMS$theta[1,2] # theta_MRCA

  #CI 
  result_BMS[i_sim,11]<- result_BMS[i_sim,5] - 1.96*result_BMS[i_sim,8]
  result_BMS[i_sim,14]<- result_BMS[i_sim,5] + 1.96*result_BMS[i_sim,8]

  if ((result_BMS[i_sim,11] < MRCA) && (MRCA < result_BMS[i_sim,14])) {
    result_BMS[i_sim,17]<- 1
  }else{
    result_BMS[i_sim,17]<- 0  
  }
  
  result_BMS[i_sim,19]<- BMS$loglik
  
  ########  ########  ########  ########  ########  ########  ########
  BM1<- try(OUwie(phylo, dat1, model=c("BM1"), algorithm="invert", get.root.theta=T,
                  starting.vals=c(0.5), scaleHeight=FALSE, 
                  check.identify=T,warn=F,quiet=T, root.age=max(hist_phylo)))
  
  #mean
  result_BM1[i_sim,1]<- BM1$solution[2] #sig1
  result_BM1[i_sim,5]<- BM1$theta[1,1] # theta_MRCA
  
  # SE
  result_BM1[i_sim,8]<- BM1$theta[1,2] # theta_MRCA
  
  #CI 
  result_BM1[i_sim,11]<- result_BM1[i_sim,5] - 1.96*result_BM1[i_sim,8]
  result_BM1[i_sim,14]<- result_BM1[i_sim,5] + 1.96*result_BM1[i_sim,8]
  
  if ((result_BM1[i_sim,11] < MRCA) && (MRCA < result_BM1[i_sim,14])) {
    result_BM1[i_sim,17]<- 1
  }else{
    result_BM1[i_sim,17]<- 0  
  }
  
  result_BM1[i_sim,19]<- BM1$loglik
  
  print(i_sim)
}
  
  dat1 = data.frame(
    MLE_MRCA = result_BSDS[,1], MLE_ev = result_BSDS[,2], MLE_k = exp(result_BSDS[,3]), 
    SE_MRCA = result_BSDS[,4],  SE_ev = result_BSDS[,5], SE_k = result_BSDS[,6], 
    z_MRCA = result_BSDS[,7],   z_ev = result_BSDS[,8],  z_k = result_BSDS[,9], 

    lo_MRCA = result_BSDS[,11], up_MRCA = result_BSDS[,12],CI_MRCA = result_BSDS[,13], 
    lo_ev = result_BSDS[,14], up_ev = result_BSDS[,15],CI_ev = result_BSDS[,16], 
    lo_k = result_BSDS[,17], up_k = result_BSDS[,18], CI_k = result_BSDS[,19], 
    
    len_MRCA = result_BSDS[,12] -result_BSDS[,11],
    len_ev = result_BSDS[,15] -result_BSDS[,14],
    len_k= result_BSDS[,18] -result_BSDS[,17],
    
   loglik_BSDS = result_BSDS[,20], det_BSDS_Hess =  result_BSDS[,24],
   vcv_not_positive = result_BSDS[,26],   conv = result_BSDS[,25],
   #mean
    OUMA_alpha1 = result_OUMA[,1], #alpha1
    OUMA_alpha2 = result_OUMA[,2], #alpha2
    OUMA_sig1 = result_OUMA[,3],  #sig1
    OUMA_sig2 = result_OUMA[,4],  #sig2
    OUMA_theta1 = result_OUMA[,5], # theta_MRCA
    OUMA_theta_sel = result_OUMA[,6], # theta_sel
    OUMA_theta_BM = result_OUMA[,7],   # theta_BM
    
    # SE
    OUMA_SE_MRCA = result_OUMA[,8],# theta_MRCA
    OUMA_SE_sel = result_OUMA[,9], # theta_sel
    OUMA_SE_BM = result_OUMA[,10], # theta_BM
    
    #CI 
    OUMA_lo_MRCA = result_OUMA[,11],
    OUMA_lo_sel = result_OUMA[,12],
    OUMA_lo_BM = result_OUMA[,13],
    
    OUMA_up_MRCA = result_OUMA[,14],
    OUMA_up_sel = result_OUMA[,15],
    OUMA_up_BM = result_OUMA[,16],
    
    OUMA_CI_MRCA = result_OUMA[,17],
    
    OUMA_CI_BM = result_OUMA[,18],
   
   OUMA_len_MRCA = result_OUMA[,14] -result_OUMA[,11],
   OUMA_len_BM = result_OUMA[,16] -result_OUMA[,13],

    OUMA_loglik = result_OUMA[,19],
    
    ##########################################    ##########################################
    OU_alpha1 = result_OU1[,1], #alpha1
    OU_sig1 = result_OU1[,3],  #sig1
    OU_theta1 = result_OU1[,5], # theta_MRCA
    OU_theta_sel = result_OU1[,6], # theta_sel

    # SE
    OU_SE_MRCA = result_OU1[,8],# theta_MRCA
    OU_SE_sel = result_OU1[,9], # theta_sel

    #CI 
    OU_lo_MRCA = result_OU1[,11],
    OU_lo_sel = result_OU1[,12],

    OU_up_MRCA = result_OU1[,14],
    OU_up_MRCA = result_OU1[,15],

    OU_CI_MRCA = result_OU1[,17],
    OU_CI_BM = result_OU1[,18],
   
   OU_len_MRCA = result_OU1[,14] -result_OU1[,11],
   OU_len_BM = result_OU1[,15] -result_OU1[,12],
    
    OU_loglik = result_OU1[,19],
    
    ##########################################    ##########################################
    BMS_sig1 = result_BMS[,1], #alpha1
    BMS_sig2 = result_BMS[,3],  #sig1
    BMS_MRCA = result_OU1[,5], # theta_MRCA

    # SE
    BMS_SE_MRCA = result_BMS[,8],# theta_MRCA
    #CI 
    BMS_lo_MRCA = result_BMS[,11],
    BMS_up_MRCA = result_BMS[,14],
    BMS_CI_MRCA = result_BMS[,17],
   
    BMS_len_MRCA = result_OU1[,14] -result_OU1[,11],

    BMS_loglik = result_BMS[,19],
    
    ##########################################    ##########################################
    BM_sig1 = result_BM1[,1], #alpha1
    BM_MRCA = result_BM1[,5], # theta_MRCA
    
    # SE
    BM_SE_MRCA = result_BM1[,8],# theta_MRCA
    #CI 
    BM_lo_MRCA = result_BM1[,11],
    BM_up_MRCA = result_BM1[,14],
    BM_CI_MRCA = result_BM1[,17],
    BM_len_MRCA = result_OU1[,14] -result_OU1[,11],
   
    BM_loglik = result_BM1[,19]
  )
  
  file.name1 <- sprintf("%.1f_%.1f_allstars_%.1f_%.2f_%.1f.csv",N_sample, N_sp, MRCA, sel, ev)
  
  write.csv(dat1, file.name1)
  #print(file.name1)
  return(dat1)
  
}

MRCA<- 100
N_sp<- 10
N_sample<- 10 
ev<- 100 

exp1<- OU_sim(N_sample=N_sample, N_sp=N_sp, MRCA=MRCA, ev=ev, sel=0.8, n_sim=1000)
exp2<- OU_sim(N_sample=N_sample, N_sp=N_sp, MRCA=MRCA, ev=ev, sel=1, n_sim=1000)
exp3<- OU_sim(N_sample=N_sample, N_sp=N_sp, MRCA=MRCA, ev=ev, sel=1.5, n_sim=1000)
exp4<- OU_sim(N_sample=N_sample, N_sp=N_sp, MRCA=MRCA, ev=ev, sel=2, n_sim=1000)
exp5<- OU_sim(N_sample=N_sample, N_sp=N_sp, MRCA=MRCA, ev=ev, sel=3, n_sim=1000)

N_sp<- 50
exp6<- OU_sim(N_sample=N_sample, N_sp=N_sp, MRCA=MRCA, ev=ev, sel=0.8, n_sim=1000)
exp7<- OU_sim(N_sample=N_sample, N_sp=N_sp, MRCA=MRCA, ev=ev, sel=1, n_sim=1000)
exp8<- OU_sim(N_sample=N_sample, N_sp=N_sp, MRCA=MRCA, ev=ev, sel=1.5, n_sim=1000)
exp9<- OU_sim(N_sample=N_sample, N_sp=N_sp, MRCA=MRCA, ev=ev, sel=2, n_sim=1000)
exp10<- OU_sim(N_sample=N_sample, N_sp=N_sp, MRCA=MRCA, ev=ev, sel=3, n_sim=1000)
