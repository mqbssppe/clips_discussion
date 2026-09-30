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
library("bmixture")
library("mvtnorm")
library("RColorBrewer")
library("coda")
#-------------------------------------------------
# gibbs sampler and related functions
#-------------------------------------------------
source('gibbs_sampler_student_mix.R')
#-------------------------------------------------
# load galaxy dataset from bmixture package
#-------------------------------------------------
data( galaxy )
x <- galaxy[,1]
par(mar = mymargin  - c(2,0,0,0))
hist(x, 20, main = "")
#-------------------------------------------------
# specify number of components
#-------------------------------------------------
K <- 3
# Run Gibbs sampler, save every 10th iteration
set.seed(456)
gibbs_galaxy_data <- gibbs_sampler_student_mix(x, m = 21000, thin = 10, K = K, df = 4)
# burn-in period
burn <- 100
# simulated means
mu_sim <- gibbs_galaxy_data$mu[-(1:burn),]
# simulated mixing proportions
w_sim <- gibbs_galaxy_data$w[-(1:burn),]
# simulated squared scale parameters
s2_sim <- gibbs_galaxy_data$s2[-(1:burn),]
# simulated allocations
z_sim <- gibbs_galaxy_data$z[-(1:burn),]
# values of the logarithm of the (non-normalized) posterior
log_posterior <- gibbs_galaxy_data$logPosterior[-(1:burn)]
# classification probabilities per MCMC iteration (for STEPHENS method)
p_matrix <- gibbs_galaxy_data$classification_probs[-(1:burn),,]
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

par(mar = mymargin + c(0,0,1,0))
matplot(mcmc.pars[, , 1], main = "Raw MCMC output", type = type, col = pal(K), ylab = "means", pch = pch, cex = cex,
  xlab = "iteration",  ylim = range(x)+ c(-1,2))

sub_ind <- 1:M


par(mar = mymargin)
plot(mu_sim[sub_ind,1], log_posterior[sub_ind], ylim = c(-240,-225), xlim = range(x), xlab = bquote(mu[1]), ylab = 'log-posterior', col = pal(K)[1], main = "Raw MCMC output (Galaxy)", pch = pch, cex = cex)

# Reordered outputs
# (a) ECR
par(mar = mymargin + c(0,0,1,0))
matplot(permute.mcmc(mcmc.pars, ls$permutations$ECR)$output[, , 1], type = type , pch = pch, cex = cex, 
  xlab = "iteration", main = "ECR", ylab = "means",  col = pal(K),ylim = range(x)+ c(-1,2))
# (b) ECR-ITERATIVE-1
par(mar = mymargin)
matplot(permute.mcmc(mcmc.pars, ls$permutations$`ECR-ITERATIVE-1`)$output[, , 1], type = type , pch = pch, cex = cex, 
  xlab = "iteration", main = "ECR-ITERATIVE-1 (Galaxy)", ylab = "means", col = pal(K),ylim = range(x))
# (c) ECR-ITERATIVE-1
par(mar = mymargin+ c(0,0,1,0))
matplot(permute.mcmc(mcmc.pars, ls$permutations$STEPHENS)$output[, , 1],type = type , pch = pch, cex = cex, 
  xlab = "iteration", main = "STEPHENS", ylab = "means",  col = pal(K),ylim = range(x)+ c(-1,2))


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

clips_func <- identifyMixture(Func = theta_functional, Mu = Mu, Eta = Eta, S = z_sim, centers = K)
clips_func$non_perm_rate
# CliPS output for means
par(mar = mymargin + c(0,0,1,0))
matplot(clips_func$Mu[,1,], ylim = range(x)+ c(-1,2),type = type , pch = pch, cex = cex, lty = c(3,1,2),
  xlab = "retained iteration", main = "CliPS", ylab = "means", col = pal(K))




#-------------------------------------------------
# k-means on ECR output 
#-------------------------------------------------

theta_ecr <- permute.mcmc(mcmc.pars, ls$permutations$ECR)$output
set.seed(1)
km <- kmeans(theta_ecr[,,1], centers = 2, nstart = 100)
m1 <- sum(km$cluster == 1)
m2 <- sum(km$cluster == 2)
ind1 <- which(km$cluster == 1)
ind2 <- which(km$cluster == 2)

# posterior estimates of main and minor modes
# means
mean1 <- summary(as.mcmc(theta_ecr[ind1,,1]))
mean2 <- summary(as.mcmc(theta_ecr[ind2,,1]))
mean1$statistics[,1:2]
mean2$statistics[,1:2]
# mixing proportions
p1 <- summary(as.mcmc(theta_ecr[ind1,,3]))
p2 <- summary(as.mcmc(theta_ecr[ind2,,3]))
p1$statistics[,1:2]
p2$statistics[,1:2]
# squared scale parameter
s21 <- summary(as.mcmc(theta_ecr[ind1,,2]))
s22 <- summary(as.mcmc(theta_ecr[ind2,,2]))
s21$statistics[,1:2]
s22$statistics[,1:2]

z1 <- z_sim[ind1, ]
znew1 <- z1
perms1 <- ls$permutations$ECR[ind1, ]
for (iter in 1:m1) {
    znew1[iter, ] <- order(perms1[iter, ])[z1[iter, ]]
}
zhat1 <- apply(znew1, 2, function(y) {
            uy <- unique(y)
            uy[which.max(tabulate(match(y, uy)))]
        })

z2 <- z_sim[ind2, ]
znew2 <- z2
perms2 <- ls$permutations$ECR[ind2, ]
for (iter in 1:m2) {
    znew2[iter, ] <- order(perms2[iter, ])[z2[iter, ]]
}
zhat2 <- apply(znew2, 2, function(y) {
            uy <- unique(y)
            uy[which.max(tabulate(match(y, uy)))]
        })

# clustering of observations for ECR's main mode 
	par(mar = mymargin + c(0,0,1,0))
	plot(x, zhat1, yaxt = "n", ylab = "cluster", main = "main mode", ylim = c(0.5,3.5), col = pal(K)[zhat1])
	axis(side = 2, at = 1:3, labels = 1:3)

# clustering of observations for ECR's minor mode
	par(mar = mymargin + c(0,0,1,0))
	plot(x, zhat2, yaxt = "n", ylab = "cluster", main = "minor mode", ylim = c(0.5,3.5), col = pal(K)[zhat2])
	axis(side = 2, at = 1:3, labels = 1:3)

# plot the raw values of means versus the log-posterior density

set.seed(1)
ind <- sample(M, 500, replace = FALSE)
par(mar = mymargin)
plot(mu_sim[ind,1], log_posterior[ind], xlim = range(x),ylim = c(-235, -225), col = 'gray80', cex = 0.5, yaxt = 'n', xlab = bquote(mu), ylab = 'log-posterior')
points(mu_sim[ind,2], log_posterior[ind], col = 'gray80', cex = 0.5)
points(mu_sim[ind,3], log_posterior[ind], col = 'gray80', cex = 0.5)
axis(2, las = 1)

# point process represenation on (\mu,\sigma^2) space
par(mar = mymargin + c(0,0.25,0,0))
plot(mcmc.pars[ind,1,1], mcmc.pars[ind,1,2], xlim = range(x),log = 'y', ylim = c(0.05,20), col = 'gray80', xlab = bquote(mu), ylab = bquote(sigma^2), yaxt = 'n', cex = 0.5)
axis(2, las = 2)
points(mcmc.pars[ind,2,1], mcmc.pars[ind,2,2], col = 'gray80', cex = 0.5)
points(mcmc.pars[ind,3,1], mcmc.pars[ind,3,2], col = 'gray80', cex = 0.5)

# Filter out the draws corresponding to minor mode and then apply CliPS

filter_out <- which(apply(mcmc.pars,1,max) > 26)

mcmc.pars_filtered <- mcmc.pars[-filter_out,,]
z_sim_filtered <- z_sim[-filter_out,]
M <- dim(mcmc.pars_filtered)[1]
theta_functional <- array(data = NA, dim = c(M, 1, K))
Mu <- array(data = NA, dim = c(M, 2, K))
for(k in 1:K){
	theta_functional[, ,k] <- mcmc.pars_filtered[,k,1]
	Mu[, ,k] <- mcmc.pars_filtered[,k,1:2]
}
Eta <- as.matrix(mcmc.pars_filtered[,,3])
Nk <- sapply(1:K, function(k) rowSums(z_sim_filtered == k))
K_vals <- rep(K,M)
Kplus <- apply(Nk, 1, function(y)sum(y>0))

clips_func <- identifyMixture(Func = theta_functional, Mu = Mu, Eta = Eta, S = z_sim_filtered, centers = K)
clips_func$non_perm_rate

# CliPS output for Filtered MCMC draws
par(mar = mymargin+ c(0,0,1,0))
matplot(clips_func$Mu[,1,], type = type , pch = pch, cex = cex, 
  xlab = "retained iteration", main = "CliPS", ylab = "means", lty = 1, col = pal(K))


colMeans(clips_func$Mu[,1,])
colMeans(clips_func$Mu[,2,])
colMeans(clips_func$Eta)


