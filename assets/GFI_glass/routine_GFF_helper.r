library(Matrix, quietly=T)
library(MASS, quietly=T)
library(bayesSurv, quietly=T)
library(LaplacesDemon, quietly=T)
library(data.table, quietly=T)
library(mvtnorm, quietly=T)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
mcmc_routine = function( par, par_index, data, steps, burnin, source, W=NULL){
	
	# Initialize various components needed for the MCMC algorithm ---------------
	if(source == 's'){
		
		Y = data
		
		group = list( par_index$mu_s, par_index$vec_A)
		# Proposal covariance matrices for each parameter group
		pcov = list( diag( length(group[[1]])), 
		             diag( length(group[[2]])))
		
		log_dens_prev = log_target_s( par, par_index, Y)
		
	} else if(source == 'a'){
		
		WID = data[,'WID']
		Y = data[,-1]
		
		group = list( par_index$mu_a, par_index$vec_B, par_index$vec_C)
		# Proposal covariance matrices for each parameter group		
		pcov = list( diag( length(group[[1]])),
		             diag( length(group[[2]])), 
		             diag( length(group[[3]])))
		
		W = do.call( 'rbind', build_W_i( unique(WID), par, par_index, Y, WID))
		log_dens_prev = log_target_a( par, par_index, W, Y)
	}
	# ---------------------------------------------------------------------------

	
	# Begin the MCMC algorithm --------------------------------------------------
	num_groups = length(group)
	pscale = rep( 1e-16, num_groups)
	accept = rep( 0, num_groups)
	chain = matrix( 0, steps, length(par))
	for(ttt in 1:steps){
		for(j in 1:num_groups){
			
			ind_j = group[[j]]
			
			# Propose an update, one parameter at a time
			proposal = par
			proposal[ind_j] = mvrnorm( n=1, mu=par[ind_j], Sigma=pcov[[j]]*pscale[j])

			# Evaluate the log density for the proposal -----------------------------
			if(source == 's'){
				log_dens = log_target_s( proposal, par_index, Y)
			} else if(source == 'a'){
				log_dens = log_target_a( proposal, par_index, W, Y)
			}
			# -----------------------------------------------------------------------
			
			if( log_dens!=-Inf & log_dens - log_dens_prev > log(runif(1,0,1)) ){
				
				log_dens_prev = log_dens
				par[ind_j] = proposal[ind_j]
				accept[j] = accept[j] +1
			} 
			chain[ttt,ind_j] = par[ind_j]
			
			
			# Proposal tuning scheme ------------------------------------------------
			if(ttt < burnin){
				# During the burnin period, update the proposal covariance in each step 
				# to capture the relationships within the parameters vectors for each 
				# transition.  This helps with mixing.
				if(ttt == 100)  pscale[j] = 1
					
				if(100 <= ttt & ttt <= 2000){  
					temp_chain = chain[1:ttt,ind_j]
					pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
					
				} else if(2000 < ttt){  
					temp_chain = chain[(ttt-2000):ttt,ind_j]
					pcov[[j]] = cov(temp_chain[ !duplicated(temp_chain),, drop=F])
				}
				if( sum( is.na(pcov[[j]]) ) > 0)  pcov[[j]] = diag( length(ind_j) )

				# Tune the proposal covariance for each transition to achieve 
				# reasonable acceptance ratios.
				if(ttt %% 30 == 0){ 
					if(ttt %% 480 == 0){  
						accept[j] = 0  

					} else if( accept[j] / (ttt %% 480) < .4 ){ 
						pscale[j] = (.75^2)*pscale[j] 
		
					} else if( accept[j] / (ttt %% 480) > .5 ){ 
						pscale[j] = (1.25^2)*pscale[j]
					} 
				}
			}
			# -----------------------------------------------------------------------
		}
		# Restart the acceptance ratio at burnin.
		if(ttt == burnin)  accept = rep( 0, num_groups)

		if(ttt%%1==0)  cat('GFF --->',ttt,'\n')
	}
	# ---------------------------------------------------------------------------
	
	return(list( chain=chain[burnin:steps,], accept=accept/(steps-burnin), 
	             pscale=pscale, par_index=par_index))
}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
log_target_s = function( par, par_index, Y){
	
	m = nrow(Y)
	p = ncol(Y)
	
	mu_s = par[par_index$mu_s]
	A = matrix( par[par_index$vec_A], ncol=p)
	
	U = t(t(Y) - mu_s)
	S_s = t(U) %*% U
	
	block_11 = diag(p)
	block_12 = diag(p) %x% matrix( colMeans(U), nrow=1)
	
	block_22 = diag(p) %x% S_s / m
	
	mat = rbind( cbind(    block_11, block_12),
		           cbind( t(block_12), block_22))
	
	eigen_vals = eigen( mat, symmetric=T, only.values=T)$values
	AAt = A %*% t(A)
	eigen_vals_AAt = eigen( AAt, symmetric=T, only.values=T)$values

	if( min(eigen_vals) > 1e-16 & min(eigen_vals_AAt) > 1e-16){

		log_Jacobian = sum(log(eigen_vals))/2
		
		inv_AAt = chol2inv(chol(AAt))
		log_dens = - sum(diag(S_s %*% inv_AAt))/2 
		
		log_det_AAt = sum(log(eigen_vals_AAt))
		
		return( log_Jacobian + log_dens - (m+p)*log_det_AAt/2 )
	
	} else{ return(-Inf) }
}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
log_target_a = function( par, par_index, W, Y){
	
	N = nrow(Y)
	p = ncol(Y)
		
	mu_a = par[par_index$mu_a]
	B = matrix( par[par_index$vec_B], ncol=p)
	C = matrix( par[par_index$vec_C], ncol=p)
	
	Q = t(t(Y - W %*% t(B)) - mu_a)
	S_a = t(Q) %*% Q
	
	block_11 = diag(p)
	block_12 = diag(p) %x% matrix( colMeans(W), nrow=1)
	block_13 = diag(p) %x% matrix( colMeans(Q), nrow=1)

	block_22 = diag(p) %x% (t(W) %*% W) / N
	block_23 = diag(p) %x% (t(W) %*% Q) / N

	block_33 = diag(p) %x% S_a / N

	mat = rbind( cbind(    block_11,    block_12, block_13),
			         cbind( t(block_12),    block_22, block_23),
				       cbind( t(block_13), t(block_23), block_33))

	eigen_vals = eigen( mat, symmetric=T, only.values=T  )$values

	if( min(eigen_vals) > 1e-16 ){
		
		inv_CCt = chol2inv(chol(C %*% t(C)))
		log_det_CCt = sum(log(eigen( C %*%t(C), symmetric=T, only.values=T)$values))
		 
		log_Jacobian = sum(log(eigen_vals))/2
		log_dens = - sum(diag(S_a %*% inv_CCt))/2 

		value = log_Jacobian + log_dens - (N+p)*log_det_CCt/2

	} else{ value = -Inf }
	
	return(value)
}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
build_W_i = function( i, par, par_index, Y, WID){
	
 	p = ncol(Y)

 	mu_a = par[par_index$mu_a]
 	B = matrix( par[par_index$vec_B], ncol=p)
 	C = matrix( par[par_index$vec_C], ncol=p)
	
	inv_CCt = chol2inv(chol(C %*% t(C)))
	inv_Sigma = t(B) %*% inv_CCt %*% B
	Sigma = chol2inv(chol(inv_Sigma))
	
	m_i = sum(WID==i)
	t = Sigma %*% t(B) %*% inv_CCt %*% (colMeans(Y[ WID==i,]) - mu_a)
	#t = rMVNorm( n=1, mean=t_hat, Q=inv_Sigma*m_i, param='standard') 
	W_i = matrix( t, m_i, p, byrow=T)

	return(W_i)
}
build_W_i = Vectorize( build_W_i, vectorize.args='i', SIMPLIFY=F)
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Function to compute GFF
# -----------------------------------------------------------------------------
GFF_subroutine = function( Y_u, chain_s, par_index_s, chain_a, par_index_a){
	
	m_u = nrow(Y_u)
	p = ncol(Y_u)
	vec_Y_u = c(t(Y_u))
	
	# Evaluate the expectation in the numerator of the GFF ----------------------
	log_f_s = rep( 0, nrow(chain_s))
	for(k in 1:nrow(chain_s)){
	
		mu_s = chain_s[k,][par_index_s$mu_s]
		A = matrix( chain_s[k,][par_index_s$vec_A], ncol=p)
	
		vec_mean = rep( 1, m_u) %x% mu_s
		vec_cov = diag(m_u) %x% (A %*% t(A))
	
		log_f_s[k] = dmvnorm( x=vec_Y_u, mean=vec_mean, sigma=vec_cov, log=T)
	
		print(k)
	}
	tau_s = max(log_f_s)
	log_E_s = tau_s + log(mean( exp(log_f_s - tau_s) ))
	# ---------------------------------------------------------------------------

	# Evaluate the expectation in the denominator of the GFF --------------------
	num_IS = 300
	df = 5
	log_f_a = rep( 0, nrow(chain_a)*num_IS)
	for(k in 1:nrow(chain_a)){
	
		mu_a = chain_a[k,][par_index_a$mu_a]
		B = matrix( chain_a[k,][par_index_a$vec_B], ncol=p)
		C = matrix( chain_a[k,][par_index_a$vec_C], ncol=p)
	
		vec_cov = diag(m_u) %x% (C %*% t(C))

		for(l in 1:num_IS){
			t = rnorm(n=p) /sqrt( rchisq( n=1, df=df) /df )
			vec_mean = rep( 1, m_u) %x% (mu_a + B %*% t)
		
			ind = (k-1) * num_IS + l
			log_f_a[ind] = dmvnorm( x=vec_Y_u, mean=vec_mean, sigma=vec_cov, log=T)
		}
		print(k)
	}
	tau_a = max(log_f_a)
	log_E_a = tau_a + log(mean( exp(log_f_a - tau_a) ))
	# ---------------------------------------------------------------------------
	
	return(exp(log_E_s - log_E_a))
}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------




