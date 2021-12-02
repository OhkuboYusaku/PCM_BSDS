library("Rcpp")
sourceCpp('get_LH.cpp')
sourceCpp('ABC_gen_Rhat.cpp')

ABC_MAP_Rhat<-function(MRCA, ev_rate, k, n){
  tmp<-matrix(0,n,12)
  for(i in 1:n){
    sim.tip<-rcpp_LH(MRCA, ev_rate, k) ##Generate simulated sp.
    MCMC<-rcpp_ABC_Rhat(110000, sim.tip[1], sim.tip[2], sim.tip[3] , sim.tip[4]) ##Conduct ABC-based estimation
    MCMC<-MCMC[10000:110000,] ##Extract MCMC samples after warm-up
    MCMC<-MCMC[order(MCMC[,13]), ] # ordering by loglik
    
    tmp[i,1]<-mean(MCMC[,1]) ## mean MRCA
    tmp[i,2]<-mean(MCMC[,5]) ## mean k
    tmp[i,3]<-mean(MCMC[,9]) ## mean ev_rate
    tmp[i,4]<-MCMC[1,5] ## MAP of k
    tmp[i,5]<-MCMC[1,9] ## MAP of ev_rate
    
    tmp[i,6]<-(MCMC[2,5]) ## median k
    tmp[i,7]<-(MCMC[2,9]) ## median ev_rate
    
    tmp[i,11]<-(MCMC[1,17])

    print(i)
  }
  return(tmp)
}

result<-ABC_MAP_Rhat(100,50,1.5,10000)
write.csv(result,"ABC_MAP_Rhat(100,50,1.5,10000)).csv")

#......