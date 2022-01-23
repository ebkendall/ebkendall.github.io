# -----------------------------------------------------------------------------
# Replicate the MCMC example in the accompanying slides
# Jonathan Williams - jpwill@live.unc.edu - williams.jonathan1@mayo.edu
# -----------------------------------------------------------------------------

# Setting this seed should generate exactly what is shown in the slides.
set.seed(27584)

# Generate a data set with true value lambda = 2.5.
n=1000
x = rexp(n=n,rate=2.5)

# Exponential likelihood function of the data.
fn = function(x,lambda){ prod(lambda * exp(-lambda * x)) }  

# Diffuse hyper-parameters.
a = 0.001
b = 0.001

# Function to compute the Metropolis-Hastings ratio.
MHRatio = function(lambda_new,lamba_curr,a,b){
	numerator = fn(x,lambda_new) * dgamma(lambda_new,shape = a,rate = b) 
	denominator = fn(x,lamba_curr) * dgamma(lamba_curr,shape = a,rate= b) 

	return( numerator / denominator )
}

# Number of steps to run the MCMC algorithm for.
steps = 10000

# Create objects to store the MCMC output.
proposal = rep(1,steps)
ratios = rep(0,steps)
coinflip = rep(0,steps)
trace = rep(0,steps)
acceptRatio = 0



# Begin the algorithm at some arbitrary value of lambda, say lambda = 1.
trace[1] = 1

# The MCMC algorithm ----------------------------------------------------------
for(ttt in 2:steps){ 
 
	proposal[ttt] = trace[ttt-1] + rnorm(n=1, mean=0, sd=0.5)

	ratios[ttt] = MHRatio(proposal[ttt],trace[ttt-1],a,b) 
	coinflip[ttt] = runif(1,0,1)

	if( coinflip[ttt] < min( ratios[ttt], 1) ){
	
		trace[ttt] = proposal[ttt]
		acceptRatio = acceptRatio +1

	} else{ trace[ttt] = trace[ttt-1] }
 
}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# Create and save the plots shown in the slides.
# -----------------------------------------------------------------------------
for(k in 1:10){
	pdf(paste('MCMCExamplePlots',k,'.pdf',sep=''))
	plot( proposal[1:k], xlab='step', ylab='lamba', xlim=range(1,10), 
	      ylim=range(0.5,3), type='p', pch=19, col='orangered', xaxt='n')
	lines( trace[1:k], xlab='step', ylab='lamba', xlim=range(1,10), type='o', 
	       lwd=3, pch=19, col='steelblue2', xaxt='n')
	axis( 1, at=seq(1,10,by=1))
	dev.off()
}

pdf('MCMCExamplePlots.pdf')
par(mfrow=c(1, 2))
plot( trace,xlab='step',ylab='lamba',col='steelblue2')
hist( trace,breaks=sqrt(steps),freq=FALSE,xlab='lambda',main=NA,xlim=range(2,3),
      col='steelblue2')
posterior = dgamma(seq(from=2,to=3,by=0.01), shape=a+n, rate=b+sum(x))
lines(seq(from=2,to=3,by=0.01) ,posterior,type='l',lwd=3)
dev.off()

print(trace[1:10])
print(proposal[1:10])
print(ratios[1:10])
print(coinflip[1:10])

# [1] 1 1      2.0667 2.2337 2.7115 2.7115 2.3685 2.5939 2.5939 2.5939
# [1] 1 0.7613 2.0667 2.2337 2.7115 2.9276 2.3685 2.5939 2.0695 2.0542
# [1] 0 1e-78 4e+134 3e+05 1964 0.0005 0.2142 21 5e-10  1e-10
# [1] 0 0.2788 0.5027 0.3707 0.2875 0.1298 0.1653 0.0457 0.8348 0.3117

