### metropolis hastings basic - Non_conjugate with one param - ###

#Let's talk about a models that don't have nice, clean posterior distributions.

#Suppose we have values(data) that represent the percentage change:`y(company_i)` in total personnel from last year to this year
#for,
#we'll say, 10 companies `n=10` coming from a particular industry. We're going to assume for now, that these are independent
#measurements from a **normal** with a known `variance = 1`, but an unknown mean `??` N(?, 1)

#Suppose we decide that our prior believes about `??` are better reflected using a **standard t-distribution** with `df=1`. 
#This particular prior distribution(`t(0,1,1)`: Cauchy Distribution) has heavier(fatter) tails than the conjugate normal 
#distribution, which can better accommodate the possibility of extreme values for `??`

#likelihood:  y ~ N(??, 1)
#prior: ?? ~ t(0,1,1)

#g(??) = ??(likelihood)*prior

#In computation, because our `g(??)` distribution includes likelihoods, which are the product of many numbers that are potentially
#small, our `g(??)` might evaluate to such a small number that the computer treats it effectively as a zero. 
#To avoid this problem,
#we're going to work on the **logarithmic scale** which will be more numerically stable. 

lg = function(mu, n, ybar){
  mu2 = mu^2
  n * (ybar * mu - mu2 / 2.0) - log(1.0 + mu2)
}
# n: sample size
# mu: what we want....
# ybar: sample mean


# function to execute the Random-Walk Metropolis-Hastings sampler with normal proposal distribution.
MH = function(n, ybar, n_iter, mu_init, cand_sd) {
  
  ## step 1, initialize
  mu_out = numeric(n_iter)                # pool
  accpt = 0                               # counting
  mu_now = mu_init                        # ***starting value***
  lg_now = lg(mu=mu_now, n=n, ybar=ybar)  #initial "g( )"....but alpha = g(??*)/g(??) ..that's all we need.
  
  ## step 2, iterate
  for (i in 1:n_iter) {
    ## step 2a
    mu_cand = rnorm(n=1, mean=mu_now, sd=cand_sd) # draw a candidate # ??* 
    
    ## step 2b
    lg_cand = lg(mu=mu_cand, n=n, ybar=ybar)      # "g(??*)"
    lalpha = lg_cand - lg_now                     # "g(??*)/g(??)"
    alpha = exp(lalpha)
    
    ## step 2c
    u = runif(1) # draw a uniform variable which will be less than alpha with probability min(1, alpha) # 0 < u < 1
    if (u < alpha) { 
      mu_now = mu_cand       # then accept the candidate
      accpt = accpt + 1      # keep track of acceptance
      lg_now = lg_cand
    }
    
    ## collect results
    mu_out[i] = mu_now 
  }
  
  ## return a list of output
  list(mu=mu_out, acc_rate=accpt/n_iter)
}


### Let's execute it!
# our data
y = c(1.2, 1.4, -0.5, 0.3, 0.9, 2.3, 1.0, 0.1, 1.3, 1.9)
ybar = mean(y)
n = length(y)
# hist
hist(y, freq = F, xlim = c(-1, 3))

# add individual data point!!(let them show up on x axis..) so..in Y-axis: rep(0,n)
points(y, rep(0, n))     #samples
points(ybar, 0, pch=19)  # sample mean

# let's plot our prior distribution of the mean 
curve(dt(x, df=1), lty=2, add=T)
# As you can see, there's a little bit of a discrepancy between our prior for mu, and what the data says mu should be...
# The prior has most of its probability mass in this region here, near 0. 
# However, the data suggest that the mean is up here near 1. 
# We expect that the posterior distribution from mu will have a mean as a compromise somewhere between 0 and 1. 

# now posterior sampling!!!!
set.seed(43)
post = MH(n, ybar, 1e3, 0, 3)
str(post)
library("coda")
traceplot(as.mcmc(post$mu))
# This last plot is called a trace plot. It shows the history of the chain and provides basic feedback about whether the chain 
#has reached its stationary distribution. It appears our proposal step size was too large (acceptance rate below 23%). 


# Let's try another.
post = MH(n, ybar, 1e3, 0.0, cand_sd=0.05)
post
post$acc_rate
traceplot(as.mcmc(post$mu))

#Oops, the acceptance rate is too high (above 50%). Let's try something in between.
post = MH(n, ybar, 1e3, 0.0, cand_sd=0.9)
post$acc_rate
post
traceplot(as.mcmc(post$mu)) # satisfied....

#Just for fun, let's see what happens if we initialize the chain at some far-off value.
post = MH(n, ybar, 1e3, 30.0, cand_sd=0.9)
post
post$acc_rate
traceplot(as.mcmc(post$mu))

#It took awhile to find the stationary distribution, but it looks like we succeeded! If we discard the first 100 or so values, 
#it appears like the rest of the samples come from the stationary distribution, our posterior distribution! 


#Let's plot the posterior density against the prior to see how the data updated our belief about ??
post$mu_keep = post$mu[-c(1:100)] # discard the first 200 samples
# plot density estimate of the posterior
plot(density(post$mu_keep, adjust=2.0), main="", xlim=c(-1.0, 3.0), xlab=expression(mu))

curve(dt(x=x, df=1), lty=2, add=TRUE) # prior for mu

points(ybar, 0, pch=19) # sample mean

# approximation to the true posterior in blue
curve(0.017*exp(lg(mu=x, n=n, ybar=ybar)), from=-1.0, to=3.0, add=TRUE, col="blue") 

#These results are encouraging, but they are preliminary. We still need to investigate more formally whether our Markov chain 
#has converged to the stationary distribution.

#Obtaining posterior samples using the Metropolis-Hastings algorithm can be time-consuming and require some fine-tuning, 
#as we've just seen. The good news is that we can rely on software to do most of the work for us.

# JAGS.....



