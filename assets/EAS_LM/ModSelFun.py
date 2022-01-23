import numpy as np
import pandas as pd
import scipy.special as special
import scipy.stats as stats
import scipy.optimize as opt
import scipy.linalg as spyla
from sklearn.linear_model import ElasticNetCV



#----------------------------------------------------------------------------------------------------------------
# The scheme for proposing index sets M.
#----------------------------------------------------------------------------------------------------------------
def proposal(AddVar,m_old,p,PropWeights):

	# Add a covariate.
	if AddVar == 0 or len(m_old) == 0:
		VarNotUsed = np.array( list(set(np.arange(p)) - set(m_old)) )
		if len(VarNotUsed) == 0:
			m = m_old
			qratio = 1
			
		else:
			prob = PropWeights[VarNotUsed] / sum(PropWeights[VarNotUsed])
			newvar = np.random.choice(VarNotUsed, size=1, p=prob)
			m = ( newvar if len(m_old) == 0 else np.sort( np.append(m_old,newvar) ) )
				
			qratio =  ( ((1 / PropWeights[newvar]) / sum(1 / PropWeights[m])) 
							         / (PropWeights[newvar] / sum(PropWeights[VarNotUsed])) )[0]
	
	# Remove a covariate.
	elif AddVar == 1:
		if len(m_old) == 1:
			m = m_old
			qratio = 1
			
		else:
			invprob = ( 1 / PropWeights[m_old] ) / sum( 1 / PropWeights[m_old] )
			discardvar = np.random.choice(m_old, size=1, p=invprob)
			m = np.sort( list( set(m_old) - set(discardvar) ) )
			newVarNotUsed = np.array( list(set(np.arange(p)) - set(m)) )
			qratio = ( (PropWeights[discardvar] / sum(PropWeights[newVarNotUsed])) 
									/ ((1 / PropWeights[discardvar]) / sum(1 / PropWeights[m_old])) )[0]	
	
	# Exchange covariates.	
	else:
		VarNotUsed = np.array( list(set(np.arange(p)) - set(m_old)) )
		if len(VarNotUsed) == 0:
			m = m_old
			qratio = 1
			
		else:
			#Add a new variable.
			prob = PropWeights[VarNotUsed] / sum(PropWeights[VarNotUsed])
			newvar = np.random.choice(VarNotUsed, size=1, p=prob)
		
			#Discard a current variable
			invprob = ( 1 / PropWeights[m_old] ) / sum( 1 / PropWeights[m_old] )
			discardvar = np.random.choice(m_old, size=1, p=invprob)

			m = np.sort( list( set(m_old) - set(discardvar) ) )
			m = ( newvar if len(m) == 0 else np.sort( np.append(m,newvar) ) )

			newVarNotUsed = np.array( list(set(np.arange(p)) - set(m)) )
		
			denom = ( (PropWeights[newvar] / sum(PropWeights[VarNotUsed])) 
							* ((1 / PropWeights[discardvar]) / sum(1 / PropWeights[m_old])) )
			num = ( (PropWeights[discardvar] / sum(PropWeights[newVarNotUsed])) 
							* ((1 / PropWeights[newvar]) / sum(1 / PropWeights[m])) )
		
			qratio = ( num / denom )[0]	
			
	return np.array([m,qratio], dtype=object)
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
# The L_0 minimization procedure for estimating E(h(beta_M)).
#----------------------------------------------------------------------------------------------------------------
def model_info(y,X,g2,L,n,p,m,N,po):
			  
	maxIter = 500
	pm = len(m)
	K = pm - 1
	Xm = X[:,m]
	IXm = np.linalg.inv( np.transpose(Xm) @ Xm )
	Hm = Xm @ IXm @ np.transpose(Xm)
	RSSm = ( np.transpose(y) @ ( np.diag(np.ones(n)) - Hm ) @ y )
	sigmaHatm = np.sqrt( RSSm / (n - pm) )
	XtransXm = np.transpose(X) @ Xm
	Lambdam = np.sum( np.diag(np.transpose(X) @ Hm @ X) )
	nu = n - pm
	betaHatm = IXm @ np.transpose(Xm)@ y
	Sigma = IXm * RSSm / nu
	
	epsilon = ( n**(0.51) /9 + pm*(np.log(p*np.pi)**1.1) /9 - po )*Lambdam *sigmaHatm**2
	
	#--------------------------------------------------------------------------------------------	
	#Iterates through all of the sampled multivariate t vectors to check if a valid smaller model exists.
	if epsilon > 0:
		
		Sample = ( np.random.multivariate_normal( np.zeros(pm), Sigma, size=N) 
							* np.sqrt( nu / np.random.chisquare(nu) ) + betaHatm )
		l = 0
		for jj in range(0,N):  
			sampled_beta = Sample[jj,:]

			g1 = XtransXm @ sampled_beta
			#Set the coefficient estimate with smallest magnitude equal to 0.
			sampled_beta[np.where(abs(sampled_beta) == np.min(abs(sampled_beta)))] = 0  
			bm = np.zeros(p)
			bm[m] = sampled_beta
			iterations = 0
			objG = np.transpose(g1 - g2 @ bm) @ (g1 - g2 @ bm) / 2
			objGprev = objG + 1

			if K > 0:
				while ( objGprev - objG > pow(10,-4) ) and ( objG >= epsilon ) and ( iterations < maxIter ):
					objGprev = objG

					c = bm - ( g2 @ (g2 @ bm - g1) / L )
					Kth_Largest = np.partition(abs(c), p-K)[p-K]
					bm = np.zeros(p)
					bm[abs(c) >= Kth_Largest] = c[abs(c) >= Kth_Largest]

					objG = np.transpose(g1 - g2 @ bm) @ (g1 - g2 @ bm) / 2

					iterations = iterations + 1	

			if objGprev - objG < 0:
				objG = objGprev

			if objG >= epsilon:
				l = l + 1
			
	else:
		
		l = N
	#-----------------------------------------------------------------------------------------------
	
	
	Avg_h = l / N
	
	if Avg_h == 0:
		logAvg_h = -np.inf
	else:
		logAvg_h = np.log(Avg_h)
		
		
	logf = 0.5*pm*np.log(np.pi) + special.gammaln(0.5*(n-pm)) + ( 0.5 * (1 + pm - n) )*np.log(RSSm) + logAvg_h
	
	return np.array([Avg_h,logf], dtype=object)
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
# The main algorithm.
#----------------------------------------------------------------------------------------------------------------
def admissible_subsets(y,X,N=None,steps,burnin,po,PropWeights=None):
	
	N = ( 100 if N is None else N )
	
	X = X / np.sqrt(np.diag(np.transpose(X) @ X))
	n = X.shape[0]
	p = X.shape[1]
	
	#-----------------------------------
	#Choose the starting model with cross-validation over the elastic net.
	if PropWeights is None:
		preprocess = ElasticNetCV(l1_ratio=[.1,.5,.7,.9,.95,.99,1], cv=10, fit_intercept=False, max_iter=100000)
		preprocess.fit( X, y)
		PropWeights = (preprocess.coef_)**2 + 1/n**2	
	
	m = [ np.where(PropWeights == np.max(PropWeights))[0][0] ]
	#-----------------------------------
	
	g2 = np.transpose(X) @ X
	L = np.max( pow( np.linalg.svd(g2)[1] ,2) )
	
	Info = model_info(y,X,g2,L,n,p,m,N,po)
	oldlogf = Info[1]
	
	#-----------------------------------

	chain = np.zeros((steps-burnin,p),dtype=bool)*1
	AcceptRatio = 0
	for t in range(0,steps):
		
		m_old = m
		
		#Propose a new model --------------------------------------------------
		AddVar = np.random.randint(0,3)
		newModel = proposal(AddVar,m_old,p,PropWeights)
		m = np.array(newModel[0],dtype=np.int32)
		qratio = newModel[1]
		
		#Gather the necessary information to calculate the acceptance ratio ---
		if np.array_equal( m, m_old) == False:			
			Info = model_info(y,X,g2,L,n,p,m,N,po)
			newlogf = Info[1]
				
		else:
			newlogf = oldlogf

		logrho = newlogf - oldlogf + np.log(qratio)
		
		
		#Determine whether to accept or reject the proposed model ---------------
		if np.log(np.random.uniform(0,1)) < logrho :
			oldlogf = newlogf
			if t > burnin:
				AcceptRatio = AcceptRatio + 1
			
		else:
			m = m_old
			
		#Record the Markov Chain after the burn-in period -----------------------
		if t == burnin:
			model = np.zeros(p,dtype=bool)*1
			model[m] = 1
			postSample = model.reshape((1,p))
			postProbs = [1]
			
		elif t > burnin:
			model = np.zeros(p,dtype=bool)*1
			model[m] = 1
			chain[t-burnin-1,:] = model
			for k in range(postSample.shape[0]):
				if np.array_equal( model, postSample[k,:]) == True:
					postProbs[k] = postProbs[k] + 1
					break
					
			if k == (postSample.shape[0]-1) and np.array_equal( model, postSample[k,:]) == False:
				postSample = np.append( postSample, model.reshape((1,p)), axis=0)
				postProbs.append(1)
				
		#print(m, ' ---> ', po)
		#print('---> ',t)
		
	postProbs = np.array(postProbs) / sum(postProbs)
	AcceptRatio = AcceptRatio / steps
	
	return pd.Series([ chain , postSample , postProbs , AcceptRatio], 
			   index=['chain','postSample','postProbs','AcceptRatio'])
	

#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------------
# A default cross-validation procedure for selecting po.
#----------------------------------------------------------------------------------------------------------------
def admissible_subsets_cv(y,X,N=None,steps=None,burnin=None,grid=None,num_folds=None):
	
	N = ( 30 if N is None else N )
	steps = ( 200 if steps is None else steps )
	burnin = ( 100 if burnin is None else burnin )
	grid = ( [1,2,3,4,5,6,7,8,9,10] if grid is None else grid )
	num_folds = ( 10 if num_folds is None else num_folds )
	
	n = X.shape[0]
	
	preprocess = ElasticNetCV(l1_ratio=[.1,.5,.7,.9,.95,.99,1], cv=10, fit_intercept=False, max_iter=1000000)
	preprocess.fit( X, y)
	tmp_wght = (preprocess.coef_)**2

	bins = np.arange( n, step=int(n/num_folds))
	fold_index = [ np.arange( bins[k], bins[k+1]) for k in np.arange(len(bins)-1) ]
	fold_index.append( np.arange(fold_index[-1][-1]+1,n) if fold_index[-1][-1]!=(n-1) else fold_index[-1] )
	
	rss = np.zeros((len(grid),num_folds))
	for j in range(num_folds):

		ind1 = fold_index[j]
		ind2 = np.sort(np.array( list(set(np.arange(n)) - set(ind1)) ))

		for k, p_grid in enumerate(grid):

			Output = admissible_subsets( y[ind2], X[ind2,], N, steps, burnin, po=p_grid, PropWeights=tmp_wght)
			m = Output['postSample'][np.where(Output['postProbs'] ==max(Output['postProbs']))[0][0]].astype(bool)
			beta_hat = ( np.linalg.inv( np.transpose(X[ind2,][:,m]) @ X[ind2,][:,m]) @ 
									  		 		     np.transpose(X[ind2,][:,m]) @ y[ind2] )

			rss[k,j] = np.sum( (y[ind1] - X[ind1,][:,m] @ beta_hat)**2 )


	rss = np.sum( rss, axis=1)
	cv_rmse = [ n*np.log( rss[k] ) + p_grid*np.log(n) for k, p_grid in enumerate(grid) ]

	po = grid[ np.where( cv_rmse == np.min(cv_rmse) )[0][0] ]

	return po
#----------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------











