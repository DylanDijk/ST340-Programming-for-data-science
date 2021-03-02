library(extraDistr)



#### This function is going to calculate the probability of generating a data point when given certain mus ####

prob_bernoulli = function(x, mu){
  
  p = length(x)
  
  prob_of_elements = vector(length = p)
  
  for(i in 1:p){
  prob_of_elements[i] = (mu[i]^x[i])*((1-mu[i])^(1-x[i]))
  }
  
  prod = prod(prob_of_elements)
  
  return(prod)
}








#### EM bernoulli algorithm ####

my_em_mix_bernoulli = function(xs, K, numit = 50){
  
  n = nrow(xs)
  p = ncol(xs)
  
  
  weights = rep(1/K, K)
  mus <- .2 + .6*xs[sample(n,K),]
  #mus = matrix(c(0.4,0.6), nrow = K, ncol = p, byrow = TRUE)
  

iterations = 1 

while(iterations < (numit + 1)){

denom_to_sum = matrix(0, nrow = n, ncol = K)

    
  for(i in 1:n){
    for(k in 1:K){
      
      denom_to_sum[i,k] = (weights[k])*(prob_bernoulli(xs[i,], mus[k,]))
      
    }
  }
  
  denom = rowSums(denom_to_sum)
  

# gamma
  
  gamma = matrix(0, nrow = n, ncol = K)
  
  for(i in 1:n){
    for(k in 1:K){
      
      gamma[i,k] = ((weights[k])*(prob_bernoulli(xs[i,], mus[k,])))
      gamma[i,k] = (gamma[i,k])/(denom[i])
      
    }
  }
  

# mus and weights  
  
  i = 1:n
for(k in 1:K){
  mus[k,] = colSums((gamma[i,k])*(xs[i,]))
  mus[k,] = (mus[k,])/(sum(gamma[i,k]))
  sum(mus[k,])

  weights[k] = (sum(gamma[i,k]))/n
}  

  
iterations = iterations + 1     
  
}  
  
return(list(W = weights,Mu =  mus, Gammas = gamma))
  
}  
  






# simulated_matrices = vector(mode = 'list', length = 100)
# 
# for(j in 1:100){
# 
# simul = matrix(nrow = 500, ncol = 2)
# 
# mu_1 = c(0.1,0.9)
# mu_2 = c(0.8,0.2)
# mus = rbind(mu_1, mu_2)
# 
# for(i in 1:500){
#   
#   mu = sample(2,1,prob=c(0.5,0.5))
#   mu = mus[mu,]
#   
#   simul[i,] = rbern(2, mu)
#   
# }
# 
# simulated_matrices[[j]] = simul
# 
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# results = vector(mode = 'list', length = 100)
# 
# 
# for(i in 1:100){
# 
# results[[i]] = my_em_mix_bernoulli(xs = simulated_matrices[[i]], K = 2, numit = 50)$Mu
# 
# }
# 
# 
# 
# apply(simplify2array(results), 1:2, mean)
