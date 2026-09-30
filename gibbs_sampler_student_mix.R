# define functions

myDirichlet <- function (alpha) {
    k <- length(alpha)
    theta <- rgamma(k, shape = alpha, rate = 1)
    return(theta/sum(theta))
}

mixLoglikelihood <- function (y, mu, s2, pi, p){
    g <- length(pi)
    n <- length(y)
    logLike <- rep(0, n)
    index <- 1:g
    epsilon <- exp(-720)
    thresh <- -745
    nn <- 0
    logpi <- log(pi)
    ef <- matrix(logpi, nrow = n, ncol = g, byrow = T)
    for (k in 1:g) {
        ef[, k] <- ef[, k] +  dmvt(as.matrix(y), delta = as.matrix(mu[k]), sigma = as.matrix(s2[k]), df = p, log = TRUE)
    }
    efmax <- apply(ef, 1, max)
    wf_k <- ef
    ef <- ef - efmax
    logLike <- efmax + log(rowSums(exp(ef)))
    return(list(ll = sum(logLike, na.rm = TRUE), log_weighted_density_k = wf_k))
}


gibbs_sampler_student_mix <- function(x, m, thin, K, df = 4, delta = 1){

	n <- length(x)
	n_parameters <- 3 * K - 1
	# prior parameters
	R <- diff(range(x))
	xi <- min(x) + R/2
	kappa <- 1/R^2
	alpha <- 2
	g <- 0.2
	h <- 100*g/(alpha * R^2)
	#delta <- 1

	
	mSave <- floor(m/thin)
	w_sim <- matrix(NA, mSave,K)
	mu_sim <- matrix(NA, mSave, K)
	s2_sim <- matrix(NA, mSave, K)
	z_sim <- matrix(NA, mSave, n)
	llValues <- numeric(mSave)
	# initialization
	w <- myDirichlet(rep(delta,K))
	b <- rgamma(1, shape = g, rate = h)
	mu <- rnorm(K, xi, sqrt(1/kappa))
	s2 <- rgamma(K, shape = alpha, rate = b)
	z <- sample(K, n, prob = w, replace = TRUE)
	p_matrix <- array(data = NA, dim = c(mSave, n, K))
	logL <- mixLoglikelihood(x, mu = mu, s2 = s2, pi = w, p = df)
	cluster_size <- numeric(K)
	qsim <- numeric(n)
	s <- 0
	maxLogL <- logL$ll
	bic <- -2*maxLogL + n_parameters * log(n)
	if(thin == 1){
		s <- s + 1
		llValues[s] <- logL$ll
		w_sim[s,] <- w
		mu_sim[s,] <- mu
		s2_sim[s,] <- s2
		z_sim[s,] <- z
	}
	for(iter in 2:m){
		# allocation variables
		probs <- logL$log_weighted_density_k
		probs_log <- probs
		probs <- array(t(apply(probs, 1, function(tmp) {
		    return(exp(tmp - max(tmp)))
		})), dim = c(n, K))
		z <- apply(probs, 1, function(tmp) {
		    if (anyNA(tmp)) {
			tmp <- rep(1, K)
		    }
		    return(sample(K, 1, prob = tmp))
		})
		for (k in 1:K) {
		    index <- which(z == k)
		    cluster_size[k] <- length(index)
		}
		# mixing proportions
		w <- myDirichlet(delta + cluster_size)
		# latent q-variables + means + variances
		for(k in 1:K){
			ind <- which(z == k)
			mu_k <- xi
			var_k <- 1/kappa
			x_minus_mu <- 0
			if(cluster_size[k] > 0){
				x_minus_mu <- (x[ind] - mu[ k])^2
				qsim[ind] <- rgamma(cluster_size[k], shape = 0.5 * (df + 1), rate = 0.5 * (df + x_minus_mu/s2[k]))
				var_k <- (kappa + sum(qsim[ind])/s2[k])^{-1}		
				mu_k <- (kappa*xi + sum(qsim[ind] * x[ind])/s2[k]) * var_k
			}
			mu[k] <- rnorm(1, mu_k, sqrt(var_k))
			s2[k] <- 1/rgamma(1, alpha + 0.5 * cluster_size[k], b + 0.5 * sum(qsim[ind] * x_minus_mu))
		}
		# beta 
		b <- rgamma(1, g + K * alpha, h + sum(1/s2))
		# logL
		logL <- mixLoglikelihood(x, mu = mu, s2 = s2, pi = w, p = df)

		if(logL$ll > maxLogL){
			maxLogL <- logL$ll
			bic <- -2*maxLogL + n_parameters * log(n)
		}

		if(iter %% thin == 0){
			s <- s + 1
			llValues[s] <- logL$ll + sum(dnorm(mu, xi, sqrt(kappa^{-1}), log = TRUE)) + # prior mu
						sum(dgamma(1/s2, shape = alpha, rate = b, log = TRUE)) + # prior s2
						dgamma(b, shape = g, rate = h, log = TRUE) + 		#prior b
						sum((delta - 1) * log(w))				# prior w
			w_sim[s,] <- w
			mu_sim[s,] <- mu
			s2_sim[s,] <- s2
			z_sim[s,] <- z
			p_matrix[s,,] <- probs
		}
		if(iter %% 1000 == 0){
			cat(paste0("iteration: ", iter),'\r')
		}
		
	}	
	cat('\n')
	cat(paste0("Gibbs sampling finished.", "\n"))




	results <- vector("list", length = 7)
	
	results[[1]] <- llValues
	results[[2]] <- w_sim
	results[[3]] <- mu_sim
	results[[4]] <- s2_sim
	results[[5]] <- z_sim
	results[[6]] <- p_matrix
	results[[7]] <- bic
	names(results) <- c("logPosterior", "w", "mu", "s2", "z", "classification_probs", "BIC")
	return(results)

}



