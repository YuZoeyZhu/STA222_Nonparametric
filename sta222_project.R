# generate the data

output3 = list()
output4 = list()

n0 = n1 = 500
x0 = rep(NA, n0)
x1 = rep(NA, n1)
set.seed(1234)
for(j in 1:n0){
  r = runif(1)
  if(r<0.5){
    x0[j] = rnorm(1, 2, 5)
  }else if(r>0.9){
    x0[j] = rnorm(1, -3, 5)
  }else{
    x0[j] = rnorm(1, 5, 1)
  }
  
  if(r<0.2){
    x1[j] = rnorm(1, 4, 2)
  }else if(r>0.7){
    x1[j] = rnorm(1, 1, 5)
  }else{
    x1[j] = rnorm(1, -5, 3)
  }
}

hist(log(x0), probability = TRUE)
lines(density(log(x0)))
hist(x1, probability = TRUE)
plot(density(log(x0)))
lines(density(log(x1)))

M = 10000
B = 5000
sample1 = list()
sample1$theta = matrix(NA, nrow = B+M, ncol = n0+n1)
sample1$theta.0 = matrix(NA, nrow = B+M, ncol = n0)
sample1$theta.1 = matrix(NA, nrow = B+M, ncol = n1)
sample1$phi = matrix(NA, nrow = B+M, ncol = n1)
sample1$sigma2 = rep(NA, B+M)
sample1$mu.H = rep(NA, B+M)
sample1$mu.G = rep(NA, B+M)
sample1$tao2.H = rep(NA, B+M)
sample1$tao2.G = rep(NA, B+M)
sample1$alpha.H = rep(NA, B+M)
sample1$alpha.G = rep(NA, B+M)
# initialization
case = 2
a_alpha.H = a_alpha.G = 5
b_alpha.H = b_alpha.G = 1
a_sigma2 = 3
b_sigma2 = 30
a_mu.H = a_mu.G = 0
b_mu.H = b_mu.G = 10
a_tao2.H = a_tao2.G = 3
b_tao2.H = b_tao2.G = 10
## alpha
alpha.H = rgamma(1, shape = a_alpha.H, rate = b_alpha.H)
alpha.G = rgamma(1, shape = a_alpha.G, rate = b_alpha.G)
## sigma2
sigma2 = 1/rgamma(1, shape = a_sigma2, rate = b_sigma2)
## mu, tao2
mu.H = rnorm(1, a_mu.H, b_mu.H)
mu.G = rnorm(1, a_mu.G, b_mu.G)
tao2.H = 1/rgamma(1, shape = a_tao2.H, rate = b_tao2.H)
tao2.G = 1/rgamma(1, shape = a_tao2.G, rate = b_tao2.G)
## theta
theta = rep(NA, n0+n1)
theta.0  = rep(NA, n0)
theta.1  = rep(NA, n1)
theta[1] = rnorm(1, mean = mu.H, sd = sqrt(tao2.H))
for (i in 2:(n0+n1)){
  # unique_theta = unique(na.omit(theta))
  table_theta = table(na.omit(theta))
  table_theta_values = as.numeric(table_theta)
  table_theta_names = as.numeric(names(table_theta))
  prob = c(alpha.H/(alpha.H + i - 1), rep(1/(alpha.H + i - 1), length(table_theta_values))*table_theta_values)
  new_theta = rnorm(1, mean = mu.H, sd = sqrt(tao2.H))
  sample_pool = c(new_theta, table_theta_names)
  theta[i] = sample(x = sample_pool, size = 1, prob = prob)
}

table(theta)
prior_theta = theta 
theta.0 = theta[1:n0]
theta.1 = theta[(n0+1):(n0+n1)]
prior_theta.0 = theta.0
prior_theta.1 = theta.1

## phi
phi = rep(NA, n1)
phi[1] = rnorm(1, mean = mu.G, sd = sqrt(tao2.G))
for (i in 2:n1){
  # unique_theta = unique(na.omit(theta))
  table_phi = table(na.omit(phi))
  table_phi_values = as.numeric(table_phi)
  table_phi_names = as.numeric(names(table_phi))
  prob = c(alpha.G/(alpha.G + i - 1), rep(1/(alpha.G + i - 1), length(table_phi_values))*table_phi_values)
  new_phi = rnorm(1, mean = mu.G, sd = sqrt(tao2.G))
  sample_pool = c(new_phi, table_phi_names)
  phi[i] = sample(x = sample_pool, size = 1, prob = prob)
}

table(phi)
prior_phi = phi 

alpha.H_cur = alpha.H
alpha.G_cur = alpha.G
sigma2_cur = sigma2
mu.H_cur = mu.H
mu.G_cur = mu.G
tao2.H_cur = tao2.H
tao2.G_cur = tao2.G

acc_num = 0

for (iter in 1:(B+M)) {
  if (iter %% 1000 == 0) {
    print(paste("at iteration", iter))
    print(paste("The acceptance rate is: ", round(acc_num/iter, 2)))
  }
  # sample theta
  for(i in 1:n0){
    ######### sample theta.0 ###################
    # remove theta i
    removed_i_theta.0 = theta.0[-i]
    
    # get n_star, theta_star
    table_removed_i_theta.0 = table(removed_i_theta.0)
    nj_removed_i.0 = as.numeric(table_removed_i_theta.0)
    theta_star_removed_i.0 = as.numeric(names(table_removed_i_theta.0))
    n_star_removed_i.0 = length(theta_star_removed_i.0)
    
    # get qj and q0.theta.0
    qj.theta.0 = dnorm(x0[i], mean = theta_star_removed_i.0, sd = sqrt(sigma2_cur))
    q0.theta.0 = (1/sqrt(2*pi)) * (1/sqrt(sigma2_cur + tao2.H_cur)) * exp(-0.5*(x0[i]^2/sigma2_cur + mu.H_cur^2/tao2.H_cur - (sigma2_cur*mu.H_cur + tao2.H_cur*x0[i])^2/(sigma2_cur*tao2.H_cur*(sigma2_cur + tao2.H_cur))))
    
    # get vector of probs
    prob_h.theta.0 = (alpha.H_cur*q0.theta.0)/(alpha.H_cur*q0.theta.0 + sum(nj_removed_i.0*qj.theta.0))
    prob_existed.theta.0 = rep(NA, n_star_removed_i.0)
    for (j in 1:n_star_removed_i.0){
      prob_existed.theta.0[j] = (nj_removed_i.0[j]*qj.theta.0[j])/(alpha.H_cur*q0.theta.0 + sum(nj_removed_i.0*qj.theta.0))
    }
    probs.theta.0 = c(prob_h.theta.0, prob_existed.theta.0)
    sum(probs.theta.0)
    # new theta from h
    var_h.theta.0 = 1/(1/tao2.H_cur + 1/sigma2_cur)
    mean_h.theta.0 = var_h.theta.0*(mu.H_cur/tao2.H_cur + x0[i]/sigma2_cur)
    new_theta_h.theta.0 = rnorm(1, mean = mean_h.theta.0, sd = sqrt(var_h.theta.0))
    
    # sample pool
    sample_pool_sim.theta.0 = c(new_theta_h.theta.0, theta_star_removed_i.0)
    
    # update theta i
    theta.0[i] = sample(x = sample_pool_sim.theta.0, size = 1, prob = probs.theta.0)
  }
  
  for(i in 1:n1){
    #############################################
    ######### sample theta.1 ###################
    # remove theta i
    removed_i_theta.1 = theta.1[-i]
    
    # get n_star, theta_star
    table_removed_i_theta.1 = table(removed_i_theta.1)
    nj_removed_i.1 = as.numeric(table_removed_i_theta.1)
    theta_star_removed_i.1 = as.numeric(names(table_removed_i_theta.1))
    n_star_removed_i.1 = length(theta_star_removed_i.1)
    
    # get qj and q0.theta.1
    qj.theta.1 = dnorm(x1[i], mean = theta_star_removed_i.1, sd = sqrt(sigma2_cur))
    q0.theta.1 = (1/sqrt(2*pi)) * (1/sqrt(sigma2_cur + tao2.G_cur)) * exp(-0.5*(x1[i]^2/sigma2_cur + mu.G_cur^2/tao2.G_cur - (sigma2_cur*mu.G_cur + tao2.G_cur*x1[i])^2/(sigma2_cur*tao2.G_cur*(sigma2_cur + tao2.G_cur))))
    
    # get vector of probs
    prob_h.theta.1 = (alpha.G_cur*q0.theta.1)/(alpha.G_cur*q0.theta.1 + sum(nj_removed_i.1*qj.theta.1))
    prob_existed.theta.1 = rep(NA, n_star_removed_i.1)
    for (j in 1:n_star_removed_i.1){
      prob_existed.theta.1[j] = (nj_removed_i.1[j]*qj.theta.1[j])/(alpha.G_cur*q0.theta.1 + sum(nj_removed_i.1*qj.theta.1))
    }
    probs.theta.1 = c(prob_h.theta.1, prob_existed.theta.1)
    
    # new theta from h
    var_h.theta.1 = 1/(1/tao2.G_cur + 1/sigma2_cur)
    mean_h.theta.1 = var_h.theta.1*(mu.G_cur/tao2.G_cur + x1[i]/sigma2_cur)
    new_theta_h.theta.1 = rnorm(1, mean = mean_h.theta.1, sd = sqrt(var_h.theta.1))
    
    # sample pool
    sample_pool_sim.theta.1 = c(new_theta_h.theta.1, theta_star_removed_i.1)
    
    # update theta i
    theta.1_update = sample(x = sample_pool_sim.theta.1, size = 1, prob = probs.theta.1)
    
    #############################################
    ######### sample phi ###################
    # remove theta i
    removed_i_phi = phi[-i]
    
    # get n_star, theta_star
    table_removed_i_phi = table(removed_i_phi)
    nj_removed_i = as.numeric(table_removed_i_phi)
    phi_star_removed_i = as.numeric(names(table_removed_i_phi))
    n_star_removed_i = length(phi_star_removed_i)
    
    # get qj and q0.phi
    qj.phi = dnorm(x1[i], mean = phi_star_removed_i, sd = sqrt(sigma2_cur))
    q0.phi = (1/sqrt(2*pi)) * (1/sqrt(sigma2_cur + tao2.G_cur)) * exp(-0.5*(x1[i]^2/sigma2_cur + mu.G_cur^2/tao2.G_cur - (sigma2_cur*mu.G_cur + tao2.G_cur*x1[i])^2/(sigma2_cur*tao2.G_cur*(sigma2_cur + tao2.G_cur))))
    
    # get vector of probs
    prob_h.phi = (alpha.G_cur*q0.phi)/(alpha.G_cur*q0.phi + sum(nj_removed_i*qj.phi))
    prob_existed.phi = rep(NA, n_star_removed_i)
    for (j in 1:n_star_removed_i){
      prob_existed.phi[j] = (nj_removed_i[j]*qj.phi[j])/(alpha.G_cur*q0.phi + sum(nj_removed_i*qj.phi))
    }
    probs.phi = c(prob_h.phi, prob_existed.phi)
    
    # new theta from h
    var_h.phi = 1/(1/tao2.G_cur + 1/sigma2_cur)
    mean_h.phi = var_h.phi*(mu.G_cur/tao2.G_cur + x1[i]/sigma2_cur)
    new_phi_h = rnorm(1, mean = mean_h.phi, sd = sqrt(var_h.phi))
    
    # sample pool
    sample_pool_sim.phi = c(new_phi_h, phi_star_removed_i)
    
    # update theta i
    phi_update = sample(x = sample_pool_sim.phi, size = 1, prob = probs.phi)
    
    ############## Metropolis Step #################
    ### evaluate with the proposed value
    # match_index_theta.1 = match(theta.1_update, sample_pool_sim.theta.1)
    # p_theta.1 = probs.theta.1[match_index_theta.1]
    
    # match_index_phi = match(phi_update, sample_pool_sim.phi)
    # p_phi = probs.phi[match_index_phi]
    
    # rho_pro <- log(dnorm(x1[i], max(theta.1_update, phi_update), sigma2_cur)) + log(p_theta.1) + log(p_phi)
    rho_pro <- log(dnorm(x1[i], max(theta.1_update, phi_update), sigma2_cur)) 
    
    ### evaluate with the current value
    theta.1_cur = theta.1[i]
    phi_cur = phi[i]
    
    # match_index_theta.1_cur = match(theta.1_cur, theta_star_removed_i.1)
    # p_theta.1_cur = ifelse(is.na(match_index_theta.1_cur), prob_h.theta.1, prob_existed.theta.1[match_index_theta.1])

    # match_index_phi_cur = match(phi_cur, phi_star_removed_i)
    # p_phi_cur = ifelse(is.na(match_index_phi_cur), prob_h.phi, prob_existed.phi[match_index_phi_cur])
    
    # rho_cur <- log(dnorm(x1[i], max(theta.1_cur, phi_cur), sigma2_cur)) + log(p_theta.1_cur) + log(p_phi_cur)
    rho_cur <- log(dnorm(x1[i], max(theta.1_cur, phi_cur), sigma2_cur)) 
    
    ## update the parameter
    if(log(runif(1)) < (rho_pro - rho_cur))  ## accept w/p rho
    {
      theta.1[i] <- theta.1_update
      phi[i] <- phi_update
      acc_num <- acc_num + 1
    }
  }
  
  
  # save theta, phi
  sample1$theta.0[iter, ] = theta.0
  sample1$theta.1[iter, ] = theta.1
  sample1$theta[iter, ] = c(theta.0, theta.1)
  sample1$phi[iter, ] = phi
  
  # update n_star and theta_star
  theta = c(theta.0, theta.1)
  table_theta = table(theta)
  nj.theta = as.numeric(table_theta)
  theta_star = as.numeric(names(table_theta))
  n_star.theta = length(theta_star)
  
  table_phi = table(phi)
  nj.phi = as.numeric(table_phi)
  phi_star = as.numeric(names(table_phi))
  n_star.phi = length(phi_star)
  
  # sample sigma2
  sigma2_cur = 1/rgamma(1, shape = a_sigma2 + (n0+n1)/2, rate = b_sigma2 + sum((x0 - theta.0)^2)/2 + sum((x1 - max(theta.1, phi))^2)/2)
  # save sigma2
  sample1$sigma2[iter] = sigma2_cur
  
  # sample mu
  var_mu.H = 1/(1/b_mu.H + n_star.theta/tao2.H_cur)
  mean_mu.H = var_mu.H*(a_mu.H/b_mu.H + sum(theta_star)/tao2.H_cur)
  mu.H_cur = rnorm(1, mean = mean_mu.H, sd = sqrt(var_mu.H))
  
  var_mu.G = 1/(1/b_mu.G + n_star.phi/tao2.G_cur)
  mean_mu.G = var_mu.G*(a_mu.G/b_mu.G + sum(phi_star)/tao2.G_cur)
  mu.G_cur = rnorm(1, mean = mean_mu.G, sd = sqrt(var_mu.G))
  
  # save mu
  sample1$mu.H[iter] = mu.H_cur
  sample1$mu.G[iter] = mu.G_cur
  
  # sample tao2
  tao2.H_cur = 1/rgamma(1, shape = a_tao2.H + n_star.theta/2, rate = b_tao2.H + sum((theta_star - mu.H_cur)^2)/2)
  tao2.G_cur = 1/rgamma(1, shape = a_tao2.G + n_star.phi/2, rate = b_tao2.G + sum((phi_star - mu.G_cur)^2)/2)
  # save tap2
  sample1$tao2.H[iter] = tao2.H_cur
  sample1$tao2.G[iter] = tao2.G_cur
  
  # sample alpha
  eta.H = rbeta(1, alpha.H_cur+1, n0+n1)
  gamma_mix1.H = rgamma(1, a_alpha.H+n_star.theta, b_alpha.H-log(eta.H))
  gamma_mix2.H = rgamma(1, a_alpha.H+n_star.theta-1, b_alpha.H-log(eta.H))
  epsilon.H = (a_alpha.H+n_star.theta-1) / ((n0+n1)*(b_alpha.H-log(eta.H)) + a_alpha.H + n_star.theta - 1)
  alpha.H_cur = sample(c(gamma_mix1.H, gamma_mix2.H), size = 1, prob = c(epsilon.H, 1-epsilon.H))
  
  eta.G = rbeta(1, alpha.G_cur+1, n1)
  gamma_mix1.G = rgamma(1, a_alpha.G+n_star.phi, b_alpha.G-log(eta.G))
  gamma_mix2.G = rgamma(1, a_alpha.G+n_star.phi-1, b_alpha.G-log(eta.G))
  epsilon.G = (a_alpha.G+n_star.phi-1) / (n1*(b_alpha.G-log(eta.G)) + a_alpha.G + n_star.phi - 1)
  alpha.G_cur = sample(c(gamma_mix1.G, gamma_mix2.G), size = 1, prob = c(epsilon.G, 1-epsilon.G))
  # save alpha
  sample1$alpha.H[iter] = alpha.H_cur
  sample1$alpha.G[iter] = alpha.G_cur
}

save(sample1, file=paste("/Users/zoeyyuzhu/Desktop/ucsc/courses/STA 222/project/polya_urn_case", case, "_a.alpha.H", a_alpha.H, "_b.alpha.H", b_alpha.H, "_a.mu.H", a_mu.H , "_b.mu.H", b_mu.H, "_a.tao2.H", a_tao2.H, "_b.tao2.H", b_tao2.H, ".RData", sep=""))


# Posterior of alpha, mu, tao2, and phi
par(mfrow = c(2, 4))
plot(sample1$alpha.H[(B+1):(B+M)], type = "l", ylab = "", main = expression(alpha[H]))
plot(sample1$alpha.G[(B+1):(B+M)], type = "l", ylab = "", main = expression(alpha[G]))
plot(sample1$mu.H[(B+1):(B+M)], type = "l",  ylab = "", main = expression(mu[H]))
plot(sample1$mu.G[(B+1):(B+M)], type = "l",  ylab = "", main = expression(mu[G]))
plot(sample1$tao2.H[(B+1):(B+M)], type = "l", ylab = "", main = expression(tau[H]^2))
plot(sample1$tao2.G[(B+1):(B+M)], type = "l", ylab = "", main = expression(tau[G]^2))
plot(sample1$sigma2[(B+1):(B+M)], type = "l",ylab = "", main = expression(sigma^2))

par(mfrow = c(2, 2))
hist(sample1$alpha[(B+1):(B+M)],  ylab="",xlab = "", main = expression(alpha))
abline(v = mean(sample1$alpha[(B+1):(B+M)]), col = "red")
abline(v = quantile(sample1$alpha[(B+1):(B+M)], 0.025), col = "blue")
abline(v = quantile(sample1$alpha[(B+1):(B+M)], 0.975), col = "blue")

hist(sample1$mu[(B+1):(B+M)],  ylab="",xlab = "", main = expression(mu))
abline(v = mean(sample1$mu[(B+1):(B+M)]), col = "red")
abline(v = quantile(sample1$mu[(B+1):(B+M)], 0.025), col = "blue")
abline(v = quantile(sample1$mu[(B+1):(B+M)], 0.975), col = "blue")


hist(sample1$tao2[(B+1):(B+M)],  ylab="",xlab = "", main = expression(tau^2))
abline(v = mean(sample1$tao2[(B+1):(B+M)]), col = "red")
abline(v = quantile(sample1$tao2[(B+1):(B+M)], 0.025), col = "blue")



abline(v = quantile(sample1$tao2[(B+1):(B+M)], 0.975), col = "blue")
hist(sample1$phi[(B+1):(B+M)],  ylab="",xlab = "", main = expression(phi))
abline(v = mean(sample1$phi[(B+1):(B+M)]), col = "red")
abline(v = quantile(sample1$phi[(B+1):(B+M)], 0.025), col = "blue")
abline(v = quantile(sample1$phi[(B+1):(B+M)], 0.975), col = "blue")

# Posterior of Theta
par(mfrow = c(2,2))
theta.0_poster_median = apply(sample1$theta.0[(B+1):(B+M), ], 2, median)
theta.0_poster_quantile1 = apply(sample1$theta.0[(B+1):(B+M), ], 2, function(x) quantile(x, 0.025))
theta.0_poster_quantile2 = apply(sample1$theta.0[(B+1):(B+M), ], 2, function(x) quantile(x, 0.975))

hist(theta.0_poster_median, main = "Posterior Median of Theta", xlab = "")
abline(v = 2, col = "blue", lwd = 2)
abline(v = -3, col = "blue", lwd = 2)
abline(v = 5, col = "blue", lwd = 2)

hist(theta.0_poster_quantile1, main = "Posterior Lower Quantile of Theta", xlab = "")
hist(theta.0_poster_quantile2, main = "Posterior Upper Quantile of Theta", xlab = "")


# Posterior Predictive of new Theta0
theta.new =  rep(NA, M)
for (i in 1:M){
  # unique_theta = unique(na.omit(theta))
  table_theta_post = table(sample1$theta[B+i, ])
  table_theta_post_values = as.numeric(table_theta_post)
  table_theta_post_names = as.numeric(names(table_theta_post))
  prob_new = c(sample1$alpha.H[B+i]/(sample1$alpha.H[B+i] + n0 + n1), rep(1/(sample1$alpha.H[B+i] + n0 + n1), length(table_theta_post_values))*table_theta_post_values)
  new_theta = rnorm(1, mean = sample1$mu.H[B+i], sd = sqrt(sample1$tao2.H[B+i]))
  sample_pool_new = c(new_theta, table_theta_post_names)
  theta.new[i] = sample(x = sample_pool_new, size = 1, prob = prob_new)
}

phi.new =  rep(NA, M)
for (i in 1:M){
  # unique_phi = unique(na.omit(phi))
  table_phi_post = table(sample1$phi[B+i, ])
  table_phi_post_values = as.numeric(table_phi_post)
  table_phi_post_names = as.numeric(names(table_phi_post))
  prob_new = c(sample1$alpha.G[B+i]/(sample1$alpha.G[B+i] + n1), rep(1/(sample1$alpha.G[B+i] + n1), length(table_phi_post_values))*table_phi_post_values)
  new_phi = rnorm(1, mean = sample1$mu.G[B+i], sd = sqrt(sample1$tao2.G[B+i]))
  sample_pool_new = c(new_phi, table_phi_post_names)
  phi.new[i] = sample(x = sample_pool_new, size = 1, prob = prob_new)
}
# par(mfrow = c(1, 1))
density_theta.new <- density(theta.new) # returns the density data
plot(density_theta.new, main = "Density Plot of Posterior Predictive Theta_new") # plots the results

# par(mfrow = c(1, 1))
density_phi.new <- density(phi.new) # returns the density data
plot(density_phi.new, main = "Density Plot of Posterior Predictive Phi_new") # plots the results


# Posterior Predictive of new y0
x0.new= rep(NA, M)
for (i in 1:M){
  x0.new[i]  = rnorm(1, mean = theta.new[i], sd = sqrt(sample1$sigma2[B+i]))
}

# Prior Predictive of y
x0_p = rep(NA, n0)
for (i in 1:n0){
  x0_p[i]  = rnorm(1, mean = prior_theta.0[i], sd = sqrt(sigma2))
}

par(mfrow = c(2,2))
hist(x0.new, prob = TRUE, xlab = "", main = "Posterior Predictive P(y_0|data)")
lines(density(x0.new), # density plot
      lwd = 2, # thickness of line
      col = "red")

hist(x0, prob = TRUE, xlab = "", main = "True y")
lines(density(x0), # density plot
      lwd = 2, # thickness of line
      col = "black")

hist(x0_p, prob = TRUE, xlab = "", main = "Prior Predicitve Density")
lines(density(x0_p), # density plot
      lwd = 2, # thickness of line
      col = "blue")

plot(density(x0.new), main = "Comparisons", col = "red", lwd = 2) # plots the results
lines(density(x0), # density plot
      lwd = 2, # thickness of line
      col = "black")
lines(density(x0_p), # density plot
      lwd = 2, # thickness of line
      col = "blue")


# Posterior Predictive of new y0
x1.new= rep(NA, M)
for (i in 1:M){
  x1.new[i]  = rnorm(1, mean = max(theta.new[i], phi.new[i]), sd = sqrt(sample1$sigma2[B+i]))
}




# Prior Predictive of y
x1_p = rep(NA, n1)
for (i in 1:n1){
  x1_p[i]  = rnorm(1, mean = max(prior_theta.1[i], prior_phi[i]), sd = sqrt(sigma2))
}

par(mfrow = c(2,2))
hist(x1.new, prob = TRUE, xlab = "", main = "Posterior Predictive P(y_0|data)")
lines(density(x1.new), # density plot
      lwd = 2, # thickness of line
      col = "red")

hist(x1, prob = TRUE, xlab = "", main = "True y")
lines(density(x1), # density plot
      lwd = 2, # thickness of line
      col = "black")

hist(x1_p, prob = TRUE, xlab = "", main = "Prior Predicitve Density")
lines(density(x1_p), # density plot
      lwd = 2, # thickness of line
      col = "blue")

plot(density(x1.new), main = "Comparisons", col = "red", lwd = 2) # plots the results
lines(density(x1), # density plot
      lwd = 2, # thickness of line
      col = "black")
lines(density(x1_p), # density plot
      lwd = 2, # thickness of line
      col = "blue")
lines(density(x0.new), main = "Comparisons", col = "green", lwd = 2) # plots the results


F0 = matrix(NA, nrow = M, ncol = 10000)
for (i in 1:M){
  F0[i, ] = pnorm(sort(x0.new), mean = theta.new[i], sd = sqrt(sample1$sigma2[B+i]))
}
F1 = matrix(NA, nrow = M, ncol = 10000)
for (i in 1:M){
  F1[i, ] = pnorm(sort(x1.new), mean = max(theta.new[i], phi.new[i]), sd = sqrt(sample1$sigma2[B+i]))
}

plot(1-colMeans(F0), 1-colMeans(F1), type = "l")
lines(1-apply(F0, 2, function(x) quantile(x, 0.975)), 1-colMeans(F1), type = "l")
lines(c(0, 0.2, 0.4, 0.6, 0.8, 1), c(0, 0.2, 0.4, 0.6, 0.8, 1))




par(mfrow = c(1,2))
hist(x0, probability = TRUE)
lines(density(x0)$x, MPT_density_mean, col = "red")


library(nimble)
data(airquality)

yd = na.omit(airquality[, 1:4])
n = dim(yd)[1]

m = colMeans(yd)
s2 = 1/cov(yd[2:4])
prec_w = 1/cov(yd)



y1 = x0
y0 = x1
code <- nimbleCode({
  alpha.G ~ dgamma(5, 1)
  alpha.H ~ dgamma(5, 1)
  sigma2 ~ dinvgamma(3, 30)
  mu.H ~ dnorm(0, var = 100)
  mu.G ~ dnorm(0, var = 100)
  tao2.H ~ dinvgamma(3, 10)
  tao2.G ~ dinvgamma(3, 10)
  
  z.H[1:(n0+n1)] ~ dCRP(alpha.H, size = n0+n1)
  z.G[1:n1] ~ dCRP(alpha.G, size = n1)


  for(i in 1:M0) {
    thetatilde[i] ~ dnorm(mu.H, var = tao2.H)
  }
  for(i in 1:M1) {
    phitilde[i] ~ dnorm(mu.G, var = tao2.G)
  }
  for(i in 1:n0) {
    y0[i] ~ dnorm(thetatilde[z.H[i]], var = sigma2)  
  }
  for(i in 1:n1) {
    y1[i] ~ dnorm(max(thetatilde[z.H[i+n0]], phitilde[z.G[i]]), var = sigma2)  
  }
})

set.seed(1)
a_alpha.H = a_alpha.G = 5
b_alpha.H = b_alpha.G = 1
a_sigma2 = 3
b_sigma2 = 30
a_mu.H = a_mu.G = 0
b_mu.H = b_mu.G = 10
a_tao2.H = a_tao2.G = 3
b_tao2.H = b_tao2.G = 10
constants <- list(n0 = 500, n1 = 500, M0 = 50, M1 = 50)
data <- list(y0 = y0, y1 = y1)
inits <- list(thetatilde = rnorm(constants$M0, 1, 10), 
              phitilde = rnorm(constants$M1, 1, 10), 
              z.H = sample(1:10, size = constants$n0 + constants$n1, replace = TRUE),
              z.G = sample(1:10, size = constants$n1,  replace = TRUE),
              alpha.H  = rgamma(1, shape = a_alpha.H, rate = b_alpha.H),
              alpha.G  = rgamma(1, shape = a_alpha.G, rate = b_alpha.G),
              sigma2 = 1/rgamma(1, shape = a_sigma2, rate = b_sigma2),
              mu.H = rnorm(1, a_mu.H, b_mu.H),
              mu.G = rnorm(1, a_mu.G, b_mu.G),
              tao2.H = 1/rgamma(1, shape = a_tao2.H, rate = b_tao2.H),
              tao2.G = 1/rgamma(1, shape = a_tao2.G, rate = b_tao2.G))
model <- nimbleModel(code, constants, data, inits)


code <- nimbleCode({
  z[1:N] ~ dCRP(alpha, size = N)
  alpha ~ dgamma(1, 1)
  for(i in 1:M) {
    thetatilde[i] ~ dnorm(0, var = 100)
    s2tilde[i] ~ dinvgamma(1, 1)
  }
  for(i in 1:N) {
    y[i] ~ dnorm(thetatilde[z[i]], var = s2tilde[z[i]])  
  }
})

set.seed(1)
constants <- list(N = 100, M = 50)
data <- list(y = c(rnorm(50, -5, sqrt(3)), rnorm(50, 5, sqrt(4))))
inits <- list(thetatilde = rnorm(constants$M, 0, 10), 
              s2tilde = rinvgamma(constants$M, 1, 1), 
              z = sample(1:10, size = constants$N, replace = TRUE),
              alpha  = 1)
model <- nimbleModel(code, constants, data, inits)




