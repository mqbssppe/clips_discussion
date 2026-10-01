# specifications for figures
mymargin <- c(4, 4, 0.75, 0.75)
pch = 16
cex = 0.25
type = 'l'
#-------------------------------------------------
# required packages
#-------------------------------------------------
library("telescope")
library("label.switching")
library("mvtnorm")
library("RColorBrewer")
library("mclust")
library("coda")
#-------------------------------------------------
# gibbs sampler and related functions
#-------------------------------------------------
source('gibbs_sampler_student_mix.R')
#-------------------------------------------------
# simulate "synthetic dataset 1"
#-------------------------------------------------
K = 3
p = 4
n <- 100
mu.true <- c(-3, 0, 3)
p.true <- c(1,1,1)
s2.true <- c(1,30,1)

# simulate data
x <- numeric(n)
#set.seed(my_seed)
set.seed(3)
q.true <- rchisq(n, df = p)
z.true <- sample(K, n, prob = p.true, replace= TRUE)
for(i in 1:n){
	k = z.true[i]
	x[i] = mu.true[k] + sqrt(s2.true[k]/q.true[i]) * rnorm(1)
}

par(mar = mymargin  - c(2,0,0,0))
hist(x, 20, main = "")

# ground-truth values set equal to the estimates arising from the true data partition:
p.true2 <- mu.true2 <- s2.true2 <- numeric(K)
for(k in 1:K){
	ind <- which(z.true == k)
	p.true2[k] <- length(ind)/n
	mu.true2[k] <- mean(x[ind])
	s2.true2[k] <- (p-2) * var(x[ind])/p
}

# Run Gibbs sampler, save every 10th iteration
gibbs_sim_data <- gibbs_sampler_student_mix(x, m = 21000, thin = 10, K = 3)
# burn-in period
burn <- 100
# simulated means
mu_sim <- gibbs_sim_data$mu[-(1:burn),]
# simulated mixing proportions
w_sim <- gibbs_sim_data$w[-(1:burn),]
# simulated squared scale parameters
s2_sim <- gibbs_sim_data$s2[-(1:burn),]
# simulated allocations
z_sim <- gibbs_sim_data$z[-(1:burn),]
# values of the logarithm of the (non-normalized) posterior
log_posterior <- gibbs_sim_data$logPosterior[-(1:burn)]
# classification probabilities per MCMC iteration (for STEPHENS method)
p_matrix <- gibbs_sim_data$classification_probs[-(1:burn),,]
#-------------------------------------------------
# post-processing using label.switching package
#-------------------------------------------------
M <- length(log_posterior)
J <- 3  #   # three different types of parameters:
# j=1 corresponds to simulated means
# j=2 corresponds to simulated variances
# j=3 corresponds to simulated weights

mcmc.pars <- array(data = NA, dim = c(M, K, J))
mcmc.pars[, , 1] <- mu_sim
mcmc.pars[, , 2] <- s2_sim
mcmc.pars[, , 3] <- w_sim
zmap <- z_sim[which.max(log_posterior),]

set <- c("ECR", "ECR-ITERATIVE-1", "STEPHENS")
ls <- label.switching(method = set, zpivot = zmap, z = z_sim, K = K, mcmc = mcmc.pars, p = p_matrix)


colors <- brewer.pal(K, name = "Set1")  # #define color pallete     
pal <- colorRampPalette(colors)

# plot the simulated means (label switching phenomenon)

par(mar = mymargin+c(0,0,1,0))
plot(mcmc.pars[, 2, 1], type = type, lty = 2, col = pal(K)[2], ylab = "means", pch = pch, cex = cex,
  xlab = "iteration", ylim = c(-4,4),main = "Raw MCMC output")
points(mcmc.pars[, 3, 1], type = type, lty = 3, col = pal(K)[3], ylab = "means", pch = pch, cex = cex,
  xlab = "iteration")
points(mcmc.pars[, 1, 1], type = type, lty = 1, col = pal(K)[1], ylab = "means", pch = pch, cex = cex,
  xlab = "iteration")
abline(h = mu.true2, lty = 2)

# Reordered outputs
# (a) ECR
par(mar = mymargin+ c(0,0,1,0))
matplot(permute.mcmc(mcmc.pars, ls$permutations$ECR)$output[, c(3,2,1), 1], type = type , pch = pch, cex = cex, 
  xlab = "iteration", main = "ECR", ylab = "means",  col = pal(K)[c(3,2,1)],ylim = c(-4, 4), lty = c(3,2,1))
  abline(h = mu.true2, lty = 2)
# (b) ECR-ITERATIVE-1
par(mar = mymargin)
matplot(permute.mcmc(mcmc.pars, ls$permutations$`ECR-ITERATIVE-1`)$output[, , 1], type = type , pch = pch, cex = cex, 
  xlab = "iteration", main = "ECR-ITERATIVE-1", ylab = "means",  col = pal(K),ylim = range(x))
    abline(h = mu.true2, lty = 2)
# (c) STEPHENS
par(mar = mymargin+ c(0,0,1,0))
matplot(permute.mcmc(mcmc.pars, ls$permutations$STEPHENS)$output[, c(3,2,1), 1],type = type , pch = pch, cex = cex, 
  xlab = "iteration", main = "STEPHENS", ylab = "means",  col = pal(K)[c(3,2,1)],ylim =  c(-4, 4), lty = c(3,2,1))
    abline(h = mu.true, lty = 2)
#-------------------------------------------------
# Post-processing using CliPS
#-------------------------------------------------
theta_functional <- array(data = NA, dim = c(M, 1, K))
Mu <- array(data = NA, dim = c(M, 2, K))
for(k in 1:K){
	theta_functional[, ,k] <- mcmc.pars[,k,1]
	Mu[, ,k] <- mcmc.pars[,k,1:2]
}
Eta <- as.matrix(mcmc.pars[,,3])
Nk <- sapply(1:K, function(k) rowSums(z_sim == k))
K_vals <- rep(K,M)
Kplus <- apply(Nk, 1, function(y)sum(y>0))

clips_func <- identifyMixture(Func = theta_functional, Mu = Mu, Eta = Eta, S = z_sim, centers = K)
clips_func$non_perm_rate
# CliPS output for means
par(mar = mymargin + c(0,0,1,0))
matplot(clips_func$Mu[,1,], ylim = c(-4, 4),type = type , pch = pch, cex = cex, 
  xlab = "retained iteration", main = "CliPS", ylab = "means", lty = c(2,3,1), col = pal(K)[c(2,3,1)])
    abline(h = mu.true2, lty = 2)


# clustered data according to ECR
zhat_ecr <- ls$clusters['ECR',]
par(mar = mymargin + c(0,0,1,0))	
plot(x, zhat_ecr, yaxt = "n", ylab = "cluster",  main = "ECR", ylim = c(0.5,3.5), col = pal(K)[zhat_ecr])
axis(side = 2, at = 1:3, labels = 1:3)

# clustered data according to CliPS
zhat_clips <- apply(clips_func$S, 2, function(y) {
            uy <- unique(y)
            uy[which.max(tabulate(match(y, uy)))]
        })
	par(mar = mymargin + c(0,0,1,0))
	plot(x, c(2,3,1)[zhat_clips], yaxt = "n", ylab = "cluster",  ylim = c(0.5,3.5), col = pal(K)[c(2,3,1)][zhat_clips], main = "CliPS")
	axis(side = 2, at = 1:3, labels = 1:3)

# adjusted Rand index
adjustedRandIndex(zhat_ecr, z.true)
adjustedRandIndex(zhat_clips, z.true)
# estimates according to ECR
theta_ecr <- permute.mcmc(mcmc.pars, ls$permutations$ECR)$output
colMeans(theta_ecr[,,1])
colMeans(theta_ecr[,,2])
colMeans(theta_ecr[,,3])
# estimates according to CliPS
colMeans(clips_func$Mu[,1,])[c(3,1,2)]
colMeans(clips_func$Eta)[c(3,1,2)]
colMeans(clips_func$Mu[,2,])[c(3,1,2)]
# plot the raw values of means versus the log-posterior density
set.seed(1)
ind <- sample(M, 500, replace = FALSE)
par(mar = mymargin)
plot(mu_sim[ind,1], log_posterior[ind], xlim = range(x),ylim = c(-245, -234), col = 'gray80', cex = 0.5, yaxt = 'n', xlab = bquote(mu), ylab = 'log-posterior')
points(mu_sim[ind,2], log_posterior[ind], col = 'gray80', cex = 0.5)
points(mu_sim[ind,3], log_posterior[ind], col = 'gray80', cex = 0.5)
axis(2, las = 1)
#add ground-truth values
abline(v = mu.true2, lty = 2)

# point process represenation on (\mu,\sigma^2) space
par(mar = mymargin + c(0,0.25,0,0))
plot(mcmc.pars[ind,1,1], mcmc.pars[ind,1,2], log = 'y', xlim = range(x), ylim = c(0.05,5), col = 'gray80', xlab = bquote(mu), ylab = bquote(sigma^2), yaxt = 'n', cex = 0.5)
axis(2, las = 2)
#add ground-truth values
points(mcmc.pars[ind,2,1], mcmc.pars[ind,2,2], col = 'gray80', cex = 0.5)
points(mcmc.pars[ind,3,1], mcmc.pars[ind,3,2], col = 'gray80', cex = 0.5)
points(mu.true2, s2.true2, pch = 16, cex = 3 * p.true2)

