def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

import numpy as np
import pandas as pd
import scipy.special as special
import scipy.stats as stats
import scipy.optimize as opt
import scipy.linalg as linalg
import scipy.sparse as sparse
import scipy.sparse.linalg as splinalg
from sklearn.linear_model import ElasticNetCV
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import TimeSeriesSplit



# -----------------------------------------------------------------------------
# Function to propose graphs during the MCMC.
# -----------------------------------------------------------------------------
def proposal_graph( AddVar, G_old, p, weights):

	# Generate a new proposal model using weights.
	if AddVar == 0 or len(G_old) == 0: # Add an active component
		VarNotUsed = np.array( list(set(np.arange(0,p**2)) - set(G_old)) )
		if len(VarNotUsed) == 0:
			G = G_old
			qratio = 1
			
		else:
			prob = weights[VarNotUsed] / sum(weights[VarNotUsed])
			newvar = np.random.choice(VarNotUsed, size=1, p=prob)
			if len(G_old) == 0:
				G = newvar
			else:
				G = np.sort( np.concatenate( (G_old, newvar) ) )
				
			qratio =  (((1 / weights[newvar]) / sum(1 / weights[G])) 
			                        / (weights[newvar] / sum(weights[VarNotUsed])))[0]
	

	elif AddVar == 1: # Remove an active component
		if len(G_old) == 1:
			G = G_old
			qratio = 1
			
		else:
			invprob = ( 1 / weights[G_old] ) / sum( 1 / weights[G_old] )
			discardvar = np.random.choice(G_old, size=1, p=invprob)
			G = np.sort( list( set(G_old) - set(discardvar) ) )
			newVarNotUsed = np.array( list(set(np.arange(0,p**2)) - set(G)) )
			qratio = ((weights[discardvar] / sum(weights[newVarNotUsed])) 
			         / ((1 / weights[discardvar]) / sum(1 / weights[G_old])))[0]
			
	else: # Switch an active component with an inactive component
		# Add a new variable
		VarNotUsed = np.array( list(set(np.arange(0,p**2)) - set(G_old)) )
		if len(VarNotUsed) == 0:
			G = G_old
			qratio = 1
			
		else:
			prob = weights[VarNotUsed] / sum(weights[VarNotUsed])
			newvar = np.random.choice(VarNotUsed, size=1, p=prob)
			#Discard a current variable
			invprob = ( 1 / weights[G_old] ) / sum( 1 / weights[G_old] )
			discardvar = np.random.choice(G_old, size=1, p=invprob)
			
			G = np.sort( list( set(G_old) - set(discardvar) ) )
			if len(G) == 0:
				G = newvar
			else:
				G = np.sort( np.concatenate( (G, newvar) ) )
			
			newVarNotUsed = np.array( list(set(np.arange(0,p**2)) - set(G)) )
		
			denom = ((weights[newvar]/sum(weights[VarNotUsed])) * 
			                        ((1/weights[discardvar]) / sum(1/weights[G_old])))
			num = ((weights[discardvar]/sum(weights[newVarNotUsed])) * 
			                                ((1/weights[newvar]) / sum(1/weights[G])))
		
			qratio = ( num / denom )[0]	
			
	return np.array([ np.array( G, dtype=np.int32), qratio], dtype=object)
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Function to compute the VAR(1) model Jacobian term.
# -----------------------------------------------------------------------------
def JacobComp( scrB_s, scrB_A, Y, scrZ_G, G, n, p, vecA):
	
	# Jacobain columns related to the derivatives of the sigma_i.
	errors = Y - scrZ_G @ vecA.reshape((len(G),1))
	scrB = sparse.hstack([ scrB_A[:,G], scrB_s @sparse.kron( np.eye(p), errors) ])
	
	# Create the Theta matrix.
	hold = np.zeros(p**2)
	hold[G] = vecA
	A = hold.reshape((p,p), order='F')
	Theta = sparse.csr_matrix((n*p,0))
	for l in range(n-1,-1,-1):
		# Concatenate blocks of matrix powers of A, from top down.
		hold = np.eye(p)
		for k in range(l):
			hold = np.concatenate((hold, np.linalg.matrix_power(A,k+1)), axis=0)
		# Concatenate zeros above the main diagonal, and then stack the current 
		# column of blocks onto the right.
		hold = sparse.csr_matrix(hold)
		Theta = sparse.hstack([ Theta, 
		                sparse.vstack([ sparse.csr_matrix(((n-1-l)*p,p)), hold ]) ])
	
	# Create script D tilde.
	scrD = Theta @ scrB * n**(-1/2)
	#Jacobian = linalg.det( (np.transpose(scrD) @ scrD).toarray() )**0.5 
	#which is equivalent to:
	Jacobian = abs(np.prod( linalg.svd( scrD.toarray(), compute_uv=False) ))
	
	return Jacobian
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Function to estimate the log expected value of h x Jacobian.
# -----------------------------------------------------------------------------
def log_avg_h( scrB_s, scrB_A, Y, G, scrZ_G, r_G, m_G, n, p, N, po, d,vecA_hat):
	
	maxIter = 500 
	p_G = len(G)
	l = 0

	# If min(m_G) >= d, then the h-function might not be zero.
	if np.prod(m_G >= d) == 1:
		
		K = p_G - 1
		for importance_sample in range(0,N):  
			
			# Sample inverse script W.
			inv_sigma_sq = [ 1/stats.invgamma.rvs( a=(n-sum(r_G[:,j]))/2, 
			                                       scale=m_G[j]/2) for j in range(p) ]
			sampled_invW = sparse.kron( sparse.eye(n), sparse.diags(inv_sigma_sq))

			g0 = (np.transpose(scrZ_G) @ sampled_invW @ scrZ_G).toarray()
			# Sample vectorized A_G.
			sampled_vecA = np.random.multivariate_normal( vecA_hat, 
			                                         linalg.inv(g0), size=1).flatten()
			
			# Compute the spectral norm of the sampled A_g.
			temp = np.zeros(p**2)
			temp[G] = sampled_vecA
			sampled_A = temp.reshape((p,p), order='F')
			norm_A = max(abs( linalg.svd( sampled_A, compute_uv=False) ))
			
			# Additional constraint in the h-function that the GF A_g are stable.
			if norm_A < 1:
				
				Lambda_G = np.sum(g0.diagonal())
				r_G_max = max(np.sum( r_G, axis=0))
				
				asympt_comp = max( 1, n**.51 *p**2 *(np.log(np.log(n))*len(G)/2 -po))
				epsilon = Lambda_G * asympt_comp
			
				if epsilon > 0:
					# -------------------------------------------------------------------
					# Initialize: Set the coefficient estimate with smallest magnitude
					# equal to 0.
					b = np.array(sampled_vecA)
					b[ np.where(abs(b) == np.min(abs(b)))[0][0] ] = 0
					# Max eigenvalue, squared.
					L = linalg.eigh(g0)[0][-1]**2
					g1 = g0 @ g0
					objG = np.transpose(sampled_vecA - b) @ g1 @(sampled_vecA - b)/2
					objGprev = objG + 1
				
					if K > 0:
						Iter = 0
						while (objGprev-objG > 10**(-4) and objG>=epsilon and Iter<maxIter):
					
							objGprev = objG
							c = b - g1 @ (b - sampled_vecA) / L 
							Kth_Largest = np.partition( abs(c), p_G-K)[p_G-K]
							b = np.zeros(p_G)
							b[abs(c) >= Kth_Largest] = c[abs(c) >= Kth_Largest]
		
							objG = (np.transpose(sampled_vecA-b) @ g1 @ (sampled_vecA-b) /2)
							Iter = Iter + 1	
					
					if objGprev - objG < 0:
						objG = objGprev
					#print('objG = ',objG,' vs epsilon = ',epsilon)
					if objG >= epsilon:
						l =l +JacobComp( scrB_s, scrB_A, Y, scrZ_G, G, n, p, sampled_vecA)
					# -------------------------------------------------------------------
				else:
					l = l + JacobComp( scrB_s, scrB_A, Y, scrZ_G, G, n, p, sampled_vecA)
	
	return -np.inf if l==0 else np.log(l/N)
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Function to estimate the unnormalized mass function of a given graph G.
# -----------------------------------------------------------------------------
def log_mass_fun( scrY, scrX, Y, scrZ, scrB_s, scrB_A, G, n, p, N, po, d):
	
	m_G, vecA_hat = comp_m_G( scrY, scrX, Y, scrZ, G, p)
	
	graph = np.zeros(p**2)
	graph[G] = 1
	r_G = (graph.reshape((p,p), order='F') > 0)
	arg = (n - np.sum( r_G, axis=0))/2
	term1 = len(G)*np.log(2*np.pi)/2 + sum( -arg*np.log(m_G/2) + 
	                                                        special.gammaln(arg) )
	
	term2 = 0
	for j in range(p):
		hold = scrX[r_G[:,j],]
		if hold.size > 0:
			#hold = hold @ np.transpose(hold)
			#term2 = term2 + (np.linalg.slogdet(hold)[0]/2 if hold.shape[0] > 1 
			#else np.asscalar(np.log(hold)/2))
			term2 = term2 + np.sum(np.log(linalg.svd( hold, compute_uv=False)**2 /2))
	
	return log_avg_h( scrB_s, scrB_A, Y, G, scrZ[:,G], r_G, m_G, n, p, N, po, 
	                                                  d, vecA_hat) + term1 - term2



def comp_m_G( scrY, scrX, Y, scrZ, G, p):
	
	scrZ_G = scrZ[:,G]
	vecA_hat = ( splinalg.inv( np.transpose(scrZ_G) @ scrZ_G ) @ 
	                                          np.transpose(scrZ_G) @ Y ).flatten()
	temp = np.zeros(p**2)
	temp[G] = vecA_hat
	A_hat = temp.reshape((p,p), order='F') 
	error = scrY - A_hat @ scrX
	m_G = np.diag( error @ np.transpose(error) )
	
	return [ m_G, vecA_hat]
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# The MCMC routine.
# -----------------------------------------------------------------------------
def EAS_VAR( scrY, scrX, steps, burnin, po=None, d=None, N=None, weights=None):
	
	p, n = scrY.shape
	N = ( 200 if N is None else N )
	
	# By defaut, assume the least restrictive sparsity for po, given fixed n 
	# and p. Note that this is not compatible with asymptotic considerations, 
	# but suffices for finite samples.
	po = ( min( p**2, n) if po is None else po ) 
	
	Y = ( scrY[:,].flatten('F') ).reshape((n*p,1))
	vecX = ( scrX[:,].flatten('F') ).reshape((n*p,1))
	scrZ = sparse.csr_matrix(np.kron( np.transpose(scrX), np.eye(p)))

	# Jacobian columns corresponding only to nonzero A.
	Ip = np.eye(p)
	scrB_s = sparse.csr_matrix((n*p,0))
	scrB_A = sparse.csr_matrix((n*p,0))
	for k in range(p):
		for l in range(p):
			coef_mat = linalg.block_diag(*([ Ip[:,k].reshape((p,1)) @ 
			                                              Ip[:,l].reshape((1,p)) ]*n))
			scrB_A = sparse.hstack([ scrB_A, sparse.csr_matrix(coef_mat@vecX)],
			                                                             format='csr')
			if l==k:
				scrB_s = sparse.hstack([ scrB_s, 
				                            sparse.csr_matrix(coef_mat) ], format='csr')
	
	# If d or weights are not provided, compute default values using elastic net
	if d is None or weights is None:
		tscv = TimeSeriesSplit(n_splits=3) 
		preprocess = ElasticNetCV(l1_ratio=[.1, .5, .7, .9, .95, .99, 1], cv=tscv, 
		                                         verbose=False, fit_intercept=False)
		#preprocess = RidgeCV(cv=tscv, fit_intercept=False)
		preprocess.fit( scrZ, Y.flatten())
		G_enet = np.where((preprocess.coef_)**2 > 0)[0]
		G_enet = ( G_enet if len(G_enet)>0 else [0] )
		
		if d is None:
			m_G_enet = comp_m_G( scrY, scrX, Y, scrZ, G_enet, p)[0]
			# Set d equal to 10 percent of the minimum m_G component for the lasso
			# selected graph
			d = min(m_G_enet) * .1
	
		if weights is None:
			min_coef_sq = min((preprocess.coef_[G_enet])**2)
			weights = (preprocess.coef_)**2 + (min_coef_sq/10 if min_coef_sq>0 else 1)

	# Starting graph
	G = [ np.where(weights == np.max(weights))[0][0] ]
	logf =log_mass_fun( scrY, scrX, Y, scrZ, scrB_s, scrB_A, G, n, p, N, po,d)
	
	
	# Begin MCMC ----------------------------------------------------------------
	AcceptRatio = 0
	for t in range(0,steps):
		
		# Propose a new graph
		AddVar = np.random.randint(0,3)
		proposed_G, qratio = proposal_graph( AddVar, G, p, weights)

		# Compute/estimate the log mass of the proposed graph
		newlogf = log_mass_fun( scrY, scrX, Y, scrZ, scrB_s, scrB_A, proposed_G,
		                                                             n, p, N, po, d)
		
		# Determine whether to accept or reject the proposed graph 
		if np.log(np.random.uniform(0,1)) < newlogf - logf + np.log(qratio):
			G = proposed_G
			logf = newlogf
			AcceptRatio = AcceptRatio + 1
				
				
		# Record the Markov Chain after the burn-in period ------------------------
		if t == burnin:
			graph = np.zeros( p**2, dtype=bool)*1
			graph[G] = 1
			postSample = graph.reshape((1,p**2))
			postProbs = [1]
			chain = np.zeros( ( steps-burnin-1, p**2), dtype=bool)*1
			
		elif t > burnin:
			graph = np.zeros( p**2, dtype=bool)*1
			graph[G] = 1
			chain[t-burnin-1,:] = graph
			for k in range(postSample.shape[0]):
				if np.array_equal( graph, postSample[k,:]) == True:
					postProbs[k] = postProbs[k] + 1
					break
				
			if (k == (postSample.shape[0]-1) and 
			                        np.array_equal( graph, postSample[k,:]) == False):
				postSample = np.append( postSample, graph.reshape((1,p**2)), axis=0)
				postProbs.append(1)
		
		#print(G)
		print('---> ', t, flush=True)
		
	postProbs = np.array(postProbs) / sum(postProbs)
	AcceptRatio = AcceptRatio / steps

	return pd.Series([ chain , postSample , postProbs , AcceptRatio , d], 
			   index=['chain','postSample','postProbs','AcceptRatio','d'])
# -----------------------------------------------------------------------------



