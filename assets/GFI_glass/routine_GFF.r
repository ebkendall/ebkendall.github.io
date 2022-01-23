
source('routine_GFF_helper.r')

log_GFF_routine = function( Y_u, Y_s, Y_a, WID, steps, burnin, n_post){
	
	# Estimate parameters for the alternative source model ----------------------
	par = NULL
	par_index = NULL

	unique_WID = unique(WID)
	n = length(unique_WID)
	N = nrow(Y_a)
	p = ncol(Y_a)
	rep_index = rep( unique_WID, table(WID))

	Y_a_means = matrix( 0, n, p)
	for(i in unique_WID)  Y_a_means[i,] = colMeans(Y_a[ WID==i,])	
	Y_a_means_stacked = Y_a_means[ rep_index,]
	
	mu_a = colMeans(Y_a)

	df = 5
	BBt = ((df-2)/df) * (t(Y_a_means) - mu_a) %*% t(t(Y_a_means) - mu_a) / (n - 1)
	B = chol(BBt)
							 							 
	CCt = t(Y_a - Y_a_means_stacked) %*% (Y_a - Y_a_means_stacked) / (N - 1)
	C = chol(CCt)
							 
	count = length(par)
	par_index[['mu_a']] = (count +1):(count +length(mu_a))
	par = c( par, mu_a)

	count = length(par)
	par_index[['vec_B']] = (count +1):(count +length(c(B)))
	par = c( par, c(B))

	count = length(par)
	par_index[['vec_C']] = (count +1):(count +length(c(C)))
	par = c( par, c(C))

	a_out = mcmc_routine( par, par_index, cbind(WID,Y_a), steps, burnin, 'a')
	# ---------------------------------------------------------------------------


	# Estimate parameters for the specific source model -------------------------
	par = NULL
	par_index = NULL

	m = nrow(Y_s)
	mu_s = colMeans(Y_s)
	
	AAt = (t(Y_s) - mu_s) %*% t(t(Y_s) - mu_s) / (m - 1)
	A = chol(AAt)

	count = length(par)
	par_index[['mu_s']] = (count +1):(count +length(mu_s))
	par = c( par, mu_s)

	count = length(par)
	par_index[['vec_A']] = (count +1):(count +length(c(A)))
	par = c( par, c(A))

	s_out = mcmc_routine( par, par_index, Y_s, steps, burnin, 's')
	# ---------------------------------------------------------------------------


	# Estimate the GFF ----------------------------------------------------------
	index_post = (steps - burnin - n_post + 1):(steps - burnin)
	
	chain_s = s_out$chain[index_post,]
	par_index_s = s_out$par_index
	chain_a = a_out$chain[index_post,]
	par_index_a = a_out$par_index
	
	GFF = GFF_subroutine( Y_u, chain_s, par_index_s, chain_a, par_index_a)
	# ---------------------------------------------------------------------------
	
	return( list( GFF=GFF, s_out=s_out, a_out=a_out) )
}