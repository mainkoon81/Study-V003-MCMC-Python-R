# Study-V003-MCMC-Python-R

> Random_Variable
<img src="https://user-images.githubusercontent.com/31917400/47187126-d62a8c80-d32a-11e8-90e2-a361686a2981.png" />

> Distribution families
<img src="https://user-images.githubusercontent.com/31917400/47054048-be76cb00-d1a7-11e8-8fbd-db042248333e.png" />

> Estimating the AVG & VAR
<img src="https://user-images.githubusercontent.com/31917400/47054380-5923d980-d1a9-11e8-8226-a2a5b1991117.png" />

> Likelihood
<img src="https://user-images.githubusercontent.com/31917400/47054381-5b863380-d1a9-11e8-98f0-08cd8524d5f0.png" />

> Estimating Maximum Likelihood 
<img src="https://user-images.githubusercontent.com/31917400/47054860-3d6e0280-d1ac-11e8-948e-da71c188bb54.png" />

> In bayesian approach, prediction is a `weighted AVG of output` of our model for all possible values of parameters while in frequentist approach, prediction is finding the best-fitted values of `parameters`. 

# Intro to Monte-Carlo
Monte-Carlo methods are methods for `generating random variables` directly or indirectly from a target distribution. 
 - **Applications**: 
   - 1> hypothesis testing 
   - 2> Bayesian computation
 - Monte Carlo is not strictly Bayesian nor strictly frequentist(Neyman-Pearson hypothesis testing), but it's **Parametric** sampling or bootstrapping (target distribution part of a family of distributions) as a **statistical simulation**, which helps us understand reality without needing to run multiple experiments or costly calculations.
 - It also helps us understand probabilities of **extreme events**. 
 - The use of Monte-Carlo methods to **calculate p-values** has become popular because:
   - Many test statistics do not have a standard asymptotic distribution. Even if a standard asymptotic distribution does exist, it may not be reliable in realistic sample sizes. 
   - In contrast, Monte-Carlo methods can be used to obtain an **Empirical p-value** that approximates the exact p-value without relying on asymptotic distributional theory or exhaustive enumeration. 

__Example 1> hypothesis testing (in R):__ Contingency table with too small sample_size
```
study = matrix(data = c(21, 2, 15, 3), nrow = 2, ncol = 2, byrow = TRUE,
               dimnames = list(c("surgery", "radiation"), c("controlled", "not controlled")))
```
<img src="https://user-images.githubusercontent.com/31917400/47085511-bf8f1300-d20e-11e8-91f8-d8b517cc484a.png" />

 - The contingency table above shows results of patients who underwent cancer treatment and either saw their cancer controlled or not. The two treatments are "surgery" and "radiation". 
 - The question we are interested in is **"Q. Is there a difference between treatment and controll of cancer?""**. Of course, a Chi-squared test would usually be used for this type of analyses.
   - There are two ways that the Chi-squared test is used:
     - 1> test the **Goodness of fit** of the `theoretical distribution` to the `observations`.
     - 2> test for **independence** between different factors(row, col)??

 - **A disadvantage of the Chi-squared test** is that it requires a `sufficient sample size` in order for the chi-square approximation to be valid. When `cell_counts` are low (below 5), asymptotic properties do not hold well. Therefore, a simple Chi-squred test may report an **invalid p-value** which would increase a **Type-I, II error** rate. 
## Here, Monte-Carlo Method can solve this issue. 
   - **Simulating contingency tables**: This function below takes some characteristics of contingency table and generates lots of contingency tables, creating a distribution of those tables. Then we can use it to calculate the Chi-squred statistics. 
   - Here, `r2dtable(n, r, c)` refers **Random 2-way Tables with Given Marginals**  where 
     - `n`: giving the **number of tables** to be drawn.
     - `r`: giving the **row totals**, to be coerced to integer. Must sum to the same as c.
     - `c`: giving the **column totals**, to be coerced to integer.
```
simulateChisq <- function(B, E, sr, sc){
    results = numeric(B)
    for(i in 1:B){
        dataa = unlist(r2dtable(1, sr, sc))   
        M = matrix(dataa, ncol = length(sc), nrow = length(sr))
        Chi_val = sum(sort( (M - E)^2 / E, decreasing = TRUE))
        results[i] = Chi_val
    }
    return(results)
}

ChisqTest <- function(data, Simulations){
    x = data                                           ## data
    B = Simulations                                    ## number of simulations to generate
    n <- sum(x)                                        ## total number of observations
    sr <- rowSums(x)                                   ## sum of rows
    sc <- colSums(x)                                   ## sum of cols
    E <- outer(sr, sc, "*")/n                          ## ORDER MATTERS
    dimnames(E) <- dimnames(study)
    tmp <- simulateChisq(B, E, sr, sc)                 ## simulated data
    Stat <- sum(sort((x - E)^2/E, decreasing = TRUE))  ## chi^2 statistic
    pval <- (1 + sum(tmp >=  Stat))/(B + 1)            ## MC p-value
    rawPVal = pchisq(q = Stat, df = 2, lower.tail = FALSE)
    
    out = list(PearsonStat = Stat, MonteCarloPVal = pval, rawPVal = rawPVal)
    return(out)
}
```
<img src="https://user-images.githubusercontent.com/31917400/47089124-fa497900-d217-11e8-9fee-8f89523dfa9a.png" />

__Example 2> Bayesian Computation (in R)__
### Without Bayesian
> This is the traditional way of Inferencing on a single proportion as a population parameter without Bayesian.
 - Let's say, a simple random sample of 1,028 US adults in March 2013 found that 56% support nuclear arms reduction. Damn + 6% !!! **"Q. Does this provide convincing evidence that a majority of Americans supported nuclear_arms_reduction at the 5% significance level?"**
 - Using a **Pearson-frequentist perspective**, we might simply do the following:
   - the number of US people supporting nuclear_arms_reduction(X) follows ~ Bin(n, p), and when `n` goes to infinity, our X also follows ~ N(np, npq). We call this `Normal Approximation for Binomial Distribution`.  
   - the proportion 'p' of them (X/n) follows ~ N(p, pq/n)
<img src="https://user-images.githubusercontent.com/31917400/47097156-b743d180-d228-11e8-975d-36c5e808ec3a.png" />

 - Under the Null-hypothesis, `p0 = q0`, so `np0 = nq0`, then `1028*0.5 = 514 > 10` which is the **mean** number of people supporting nuclear_arms_reduction.
 - Based on the normal model, the test statistic can be computed as the Z-score of the point estimate: `Z = (p_obv - p0)/SE`.  
 - SE can be computed: `SE = sqrt(p0*q0/n)` and the Null-hypothesis `p0 = 0.5` is used again here because this is a hypothesis test for a single proportion `sqrt(0.5*0.5/1028) = 0.016`, so our Z is `(0.56-0.5)/0.016 = 3.75`.
 - p-value is `1 - pnorm(q=3.75) = 8.841729e-05`. We can then look up the upper tail area, the p-value, and see that it is less than 0.001. With a `p-value < 0.05`, we reject the null hypothesis and conclude that the **poll provides evidence that a majority (greater than 50%) of Americans supported the nuclear_arms_reduction**. The 95% CI for `p_obv` would be `0.56 + c(-1,1)*1.96*0.016`. So on AVG, around from 52% to 59% of US people support the nuclear_arms_reduction. 
<img src="https://user-images.githubusercontent.com/31917400/47958986-c0c28d00-dfcf-11e8-8708-13315760a4c7.jpg" />

### With Bayesian
> Another perspective on this problem(Inferencing on a single proportion as a population parameter) is that of a Bayesian. 
<img src="https://user-images.githubusercontent.com/31917400/47188651-d463c780-d330-11e8-85b3-5e5d5d34acc9.png" />

 - Look, we have a likelihood which is Binomial. 
 - `Beta` is a conjugate prior for the `Binomial` likelihood. That's why we choose Beta as our prior...then what the posterior will be? 
 - If the posterior distributions `p(θ|x)` are in the same **distribution family** as the prior distribution `p(θ)`:
   - the prior and posterior are then called **conjugate distributions**, 
   - the `**likelihood function**` is usually well-determined from a statement of the data-generating process.
 - Let's see our prior. Prior expresses one's beliefs about this quantity before some evidence is taken into account. Here, the prior could be the probability distribution representing the relative proportions of advocaters who will support nuclear_arms_reduction. 
   - we chose Beta(1,1) as our prior and this is equivalent to Unif(0,1) and this is a non-informative prior, which means we don't have any prior information to add to this model. 

__Q. So..for Binomial Likelihood, why choose "Beta" as a prior?__ how to elicit prior distribution?
 - Theoretically, a prior(the form of the conjugate prior can generally be determined by) is a **CDF for the parameter θ distribution**. In practice, based on likelihood we have, we choose a conjugate prior from a conjugate family that's sufficiently flexible such that a member of the family will represent our prior beliefs(of course in general, if one has enough data, the information in the data will overwhelm this invasion of prior). And **then any reasonable choice of prior will lead to approximately the same posterior**.
   - However, Notice!! there are somethings that can go wrong. In the Bayesian context, events with `P(θ)=0` will have `P(θ|y)=0`. And events with `P(θ)=1` will have `P(θ|x)=1`. Thus a good bayesian will not assign probability of `0` or `1` to any event that has already occurred or already known not to occur. 
   
# Scenario_01. No-data? "Estimate data points" (with respect to θ : **Prior Predictive Distribution** for X)
<img src="https://user-images.githubusercontent.com/31917400/47260255-aa84df00-d4af-11e8-9d2c-eee68bd26b2c.png" />

   - Before observe data points, we compute a prior predictive interval (such that 95% of new observations are expected to fall into it). It's an interval for the `data points` rather than an interval for parameter we've been looking at. Prior predictive intervals are useful because they reveal the `consequences of the θ` at the data (observation) level. See, our predictive distribution of `data points` is **marginal:** `P(x) = S P(θ,x)dθ = S P(x|θ)P(θ)dθ `. 
   - To find this data point intervals, we first work with `prior` **before we observe any data**. 
   - For example, Bin(n,θ): 
     - Flip a coin 'n' times and count the number of 'H' we see. This, of course, will depend on the coin itself. "What's the probability that it shows up 'H's?" which is referring a `θ` distribution. `X` for the number of 'H's(success), as `X` being the sum of y: `X = SUM(y...)` and as we go from '1 to n' of y which is each individual coin flip(HEAD: y=1, TAIL: y=0)...but now set this aside for a while.  
     - Let's start. As for the prior(θ,θ,θ,θ,...), if we first assume that **all possible coins are equally likely**(all same θ), then `p(θ) = 1 where {0 <= θ <= 1}`, which means the probability of θ: `p(θ)` will flat...over the interval from θ=0 to θ=1. We first assume our prior is `p(θ) = 1`. 
     - So now go back to `X` and ask "what's our **predictive distribution of X** (for the number of 'H'. of course, `X` can take possible values 0, 1, 2,..up to n). The **marginal**: `P(X) = S P(X|θ)P(θ)dθ = S P(θ,X)dθ` so if n=10, we have 
     <img src="https://user-images.githubusercontent.com/31917400/47243997-d2127380-d3eb-11e8-87e0-717f9f022b50.png" />

     - Because we're interested in X now, it's important that we distinguish between a binomial density and a Bernoulli density. So here we just care about the total count rather than the exact ordering which would be Bernoulli's. But **for most of the analyses we're doing, where we're interested in θ rather than x, the binomial and the Bernoulli are interchangeable** because Binomial distribution is a distribution of sum of i.i.d. Bernoulli random variables and their likelihoods are equivalent, thus we will get the same posterior at the end. But here we care about x for a predicted distribution so we do need to specify that we're looking at **binomial** because we're looking H-counts. 
     - If we `simplify` the pdf of our data points distribution, first recall that we can write `n! = gamma(n + 1)`, then this model look like a Beta density. `Z ~ Beta(a,b)`:
     <img src="https://user-images.githubusercontent.com/31917400/47245353-d725f180-d3f0-11e8-8695-59a63dead938.png" />
     
     - Because it's a Beta density, we know all densities integrate up to 1. Thus we see that **if we start with a uniform prior, we then end up with a discrete uniform predictive density for X**. If all possible coins or all possible probabilities(θ) are equally likely, then all possible X outcomes are equally likely. That's why when we choose Beta(1,1) as our prior and this is equivalent to Unif(0,1) and this is a non-informative prior. 
     <img src="https://user-images.githubusercontent.com/31917400/47245743-73042d00-d3f2-11e8-83d4-ccded86edcb5.png" />
### **Hey! we just found the form of the prior function!! which is `sth x Beta(a,b)` !!! 

> Note: **What is Beta(α,β)?** it represents a distribution of probabilities(Random Variable is probability values) that is, it represents all the possible values of a probability. 
 - Imagine we have a baseball player, and we predict what his season-long batting average `θ` will be. You might say we can just use his batting average so far- but this will be a very poor measure at the start of a season! Given our batting average problem, which can be represented with a `binomial(X: Sum of successes)`, the best way to represent these prior expectations is with the `Beta(X: probability)` - it's saying, before we've seen the player take his first swing, what we roughly expect his batting average to be. The domain of the Beta distribution is [0, 1], because X value is a probability. We may expect that the player's season-long batting average will be most likely around `0.27`, but that it could reasonably range from `0.21 to 0.35`. This can be represented with a Beta with parameters `α=81` and `β=219`(where α=No.success, β=No.failure and the mean of Beta is `α/(α+β)=0.27`) 
<img src="https://user-images.githubusercontent.com/31917400/47266266-efe6f200-d52b-11e8-96e4-3ad67d0d988d.png" />
 
### Next,
 - __prior -> posterior__
   - When our prior for a Bernoulli likelihood(such as `p(θ) = 1`) is a `uniform`, we get a beta posterior with hyper-parameter:  
   <img src="https://user-images.githubusercontent.com/31917400/47256574-970b5100-d47a-11e8-8e66-d182c48ac514.png" />
   - In fact, the uniform distribution is a `Beta(1,1)`, and any beta distribution is conjugate for the Bernoulli distribution.
   - In posterior, the hyper-parameters are transformed: `a + sum(x), b + n - sum(x)` 
   
# Scenario_02. Don't go with a flat prior! We have some data-point (with respect to θ: **Posterior Predictive Distribution** for X)
<img src="https://user-images.githubusercontent.com/31917400/47260902-f2a9fe80-d4bb-11e8-8c80-6944cbcdf2cf.png" />
 
   - What about after we've observed data? Suppose we observe, after one flip, we got a 'H' the first time. We want to ask, what's our **predicted distribution for the second flip(H or T)**, given that we saw a 'H' on the first flip? 
   - `P(y2|y1) = S P(θ|y1, y2)dθ = S P(y2|θ,y1)P(θ|y1)dθ = S P(y2|θ)P(θ|y1)dθ`: using posterior distribution instead of prior, and `'y1 and y2' is independent` so conditional is meaningless..so we take y1 out, then..
   <img src="https://user-images.githubusercontent.com/31917400/47248777-bb2c4b00-d404-11e8-9b5d-2c67ff7f3d24.png" />
   
   - We can see here, that the posterior is a combination of the information in the prior and the information in the data. In this case, our prior is like having two data points, one 'H' and one 'T'. Saying we have a uniform prior for θ, is actually equivalent to saying we have observed one 'H' and one 'T'. And then, when we do go ahead and observe one head, it's like we now have seen two heads and one tail, and so our posterior predictive distribution for the second flip, says, if we have two heads and one tail, then we have a probability of two-thirds of getting another head, and a probability of one-third of getting a tail. 

### posterior mean & sample size
 - The effective sample size(`a+b`) gives you an idea of how much data you would need to make sure that you're prior doesn't have much influence on your posterior. If `a+b` is small compared to `sample_size: n` (non-informative prior), then the posterior will largely just be driven by the data `X`. If `a+b` is large relative to `sample_size: n` (informative prior), then your posterior will be largely driven by the prior `θ`. 
 - posterior_mean is `(prior_weight*prior_mean) + (data_weight*data_weight)`
 <img src="https://user-images.githubusercontent.com/31917400/47261034-aca26a00-d4be-11e8-897a-a9270ac5d2cb.png" />
 <img src="https://user-images.githubusercontent.com/31917400/47261080-c09a9b80-d4bf-11e8-84ef-0f8910747894.png" />

### Find posterior
 - When a family of `conjugate priors` exists, choosing a prior from that family simplifies calculation of the posterior distribution.
 - `Parameters` of prior distributions are a kind of `hyperparameter`. For example, if one uses `Beta(a,b)` to model the distribution of the parameter `p` of Bernoulli, then:
   - `p` is a parameter of the underlying system (Bernoulli), and
   - `a` and `b` are parameters of the prior distribution (Beta); hence hyperparameters
 - Sometimes hyper-parameters themselves in prior have `hyper distributions` expressing beliefs about their values in the posterior. A Bayesian model with more than one level of prior like this is called a `hierarchical Bayes model`.  

table of conjugate distribution
<img src="https://user-images.githubusercontent.com/31917400/47190665-91a6ed00-d33a-11e8-8f51-c3ab391a4871.png" />

 - find the posterior
<img src="https://user-images.githubusercontent.com/31917400/47119888-9c8f4e00-d264-11e8-9846-a03b7cb95e4d.png" />

 - the Bayesian Data model is `y|θ ~ Bin(n,θ)` thus `θ ~ Beta(a,b)`
 - the resulting posterior is then `θ|y ~ Beta(y+a, n-y+b)`. **We can now simulate the posterior distribution**, to choose `θ` !
 - Did you find the posterior? then build a Credible Interval. 
   - In Confidence Interval, the true value(population_parameter) is not a random variable. It is a fixed but unknown quantity. In contrast, our estimate is a random variable as it depends on our data x. Thus, we get different estimates each time, repeating our study. 
   - In Credible Intervals, we assume that the true value(population parameter θ) is a random variable. Thus, we capture the uncertainty about the true parameter value by a **imposing a prior distribution** on the true parameter vector. <img src="https://user-images.githubusercontent.com/31917400/47217688-0ad92b00-d3a1-11e8-9e2e-9efd544efc06.png" />

   - Using bayes theorem, we construct the posterior distribution for the parameter vector by blending the prior and the data(likelihood) we have, then arrive at **a point estimate** using the posterior distribution(use the mean of the posterior for example). However, the true parameter vector is a random variable, we also want to know the extent of uncertainty we have in our point estimate. Thus, we construct a 95% credible interval such that the following holds: `P( l(θ)<=θ<=u(θ) ) = 0.95` 

```
N = 10^4
set.seed(123)
x = rbeta(n = N, shape1 = 576 + 1, shape2 = 1028 - 576 + 1)
d = density(x)
hist(x = x, probability = TRUE, main = "Beta Posterior Distribution",
     xlab = expression(theta), ylab = "Density",
     ylim = c(0,40), col = "gray", border = "white")
lines(x = d$x , y = d$y, type = "l", col = 2) # add chart
abline(v = median(x), lty = 3, col = "3")

print("Median: ")
print(quantile(x = x, probs = c(0.025, 0.5, 0.975)))
```
<img src="https://user-images.githubusercontent.com/31917400/47120410-7e2a5200-d266-11e8-9613-75a2fe8e6995.png" />

### other prior
 - __Normal__
   - If you are collecting normal data to make inferences about `μ`, but you don't know the true `σ^2` of the data, three options available to you are:
     - 1. Fix σ^2 at your best guess.
     - 2. Estimate σ^2 from the data and fix it at this value.
     - 3. Specify a prior for σ^2 and μ to estimate them jointly
   - Options 1 and 2 allow you to use the methods "pretending you know the true value of σ^2 and this leads to a simpler posterior calculation for μ.
   - Option 3 allows you to more honestly reflect your uncertainty in σ^2, thereby protecting against overly confident inferences. 
 

### Two catagories of prior
 - __Informative prior__
   - An informative prior is a prior that is not dominated by the likelihood and that has an impact on the posterior distribution. If a prior distribution dominates the likelihood, it is clearly an informative prior.
   - A reasonable approach is to make the prior a `normal distribution` with expected value equal to the given mean value, with variance equal to the given variance. 
   - pre-existing evidence which has already been taken into account is part of the prior and, as more evidence accumulates, the posterior is determined largely by the evidence rather than any original assumption, provided that the original assumption admitted the possibility of what the evidence is suggesting.
   
 - __Non-informative prior__
   - Roughly speaking, a prior distribution is non-informative if the prior is "flat". 
   - Non-informative priors can express "objective" information such as "the variable is positive" or "the variable is less than some limit". The simplest rule for determining a non-informative prior is the principle of indifference, which assigns equal probabilities to all possibilities, thus the use of an uninformative prior typically yields results which are not too different from conventional statistical analysis.
   - However, Non-informative priors are useful when 'stronger' priors would unjustifiably favor some hypotheses in a way that's inconsistent with your actual (lack of) knowledge/beliefs. 
   - Or Priors can also be chosen according to some principle, such as symmetry or maximizing entropy given constraints; examples are `Jeffreys' prior` for the Bernoulli random variable.
   - Jeffreys prior
     - Choosing a uniform prior depends upon the particular parameterization. For example, thinking about normal distribution. Suppose 
     <img src="https://user-images.githubusercontent.com/31917400/47270910-5e4aa500-d56a-11e8-8831-e8c25df1f9ca.png" />
     
       - Say, I use a prior which is uniform `σ^2*f(σ^2)=1`. Suppose somebody just want to put a uniform prior on `σ^2` itself `f(σ^2)=1`. These are both uniform on certain parameterizations, but they are different priors. And so when we compute the posteriors, we will get different posteriors. The key thing is that **uniform priors are not invariant with respect to transformation**. 
     - One attempt to round to this is the Jeffreys prior. Jeffreys prior is to find this `prior` proportional to the `sqrt(fisher information)`. Fisher Information is `E[{d(log(P(x|θ))/dθ}^2]`, and we square-root it. In most cases, this will be an improper prior.
     - For the example of normal data, `σ^2*f(σ^2) ∝ 1` and `f(μ) ∝ 1`
     - For the example of Bernoulli data, `Beta(0.5,0.5)`
     - Jeffreys’ prior is locally uniform and hence noninformative . It provides an automated scheme for finding a noninformative prior for any parametric model. 
     - Jeffreys’ prior is invariant with respect to one-to-one transformations or the change of scale.
     - Jeffreys’ prior has some shortcomings: the prior is improper for many models, which leads to improper posterior.
     
---------------------------------------------------------------------------------------------------------   
## Before saying Monte-Carlo
### Hierarchical form of the Conjugate model
<img src="https://user-images.githubusercontent.com/31917400/47272952-c4442600-d584-11e8-947b-0c128e83f9d8.png" />
 
 - The hierarchical representations above show how you could hypothetically simulate data from this model. 
   - You start with the variables that don't have any dependence on any other variables. You would simulate those
   - And then given those draws, you would simulate from the distributions for these other variables further down the chain.
   - This is also how you might simulate from a **prior predictive distribution** for X. 
 - The posterior distribution will be derived from the components of the hierarchical structure.
 - If we add more layers, then surprisingly...(Here `μ is unknown`, and `σ^2 is known`, then the `conjugate prior from μ` would be a `normal distribution`)
<img src="https://user-images.githubusercontent.com/31917400/47273757-fe67f480-d591-11e8-8731-9e273af547c2.png" />

 - If we can recognize this standard form as being proportional to a common distribution, then our work is done, and we know what our posterior distribution is and we can do **direct simulation from a posterior distribution**. However, if we do not use conjugate priors or if the models are more complicated, then the posterior distribution will not have a standard form that we can recognize.

### Non-Conjugate model !!!!!!!!!!!!!!!!!!!!!!
 - Suppose we're now going to estimate `μ` and `σ^2`, because they're both `unknown` (If sigma squared were known, the conjugate prior from `μ` would be a `normal distribution`. And if `μ` were known, the conjugate prior we could choose for `σ^2` would be an `inverse gamma`). 
 - In the more general case that we have here(both unknown), the posterior distribution does not appear as a distribution that we can simulate or integrate. We are unable to integrate it to obtain important quantities, such as the posterior mean or probability intervals. However, the **computational methods** invented in the 1950's revolutionized this field. We do have the ability to simulate from this challenging posterior distributions 

### Monte-Carlo methods will be helpful for generating samples from difficult to sample target distributions.
How? by generating random number from target distributions through **transformation methods**??
 - Monte Carlo simulation uses random sampling and statistical modeling to mimic the operations of complex systems. 
   - Monte Carlo models a system as a series of PDF.
   - Monte Carlo repeatedly samples from the PDF.
   - Monte Carlo computes the statistics of interest.
   
__Integration & Error__
 - Monte Carlo estimation refers to **simulating hypothetical draws** from a probability distribution in order to calculate important quantities such as mean, variance, the probability of some event, etc(all of these calculations involve `integration`, which can be very difficult to compute analytically). So we will find the area of the swimming pool by throwing balls. 
 - But how good is an approximation by Monte-Carlo sampling? Again, we can turn to the central limit theorem, which tells us that the variance of our estimate is controlled in part by our sample size. Increase the number of samples simulated. 
<img src="https://user-images.githubusercontent.com/31917400/47397900-600e9700-d729-11e8-93e9-ada0d603129b.png" />

---------------------------------------------------------------------------------------------------------------
# Markov chain Monte Carlo (MCMC)
## 1. Markov Chain
 - Let's say we can use a chain rule to calculate the probability of the entire sequence. 
 <img src="https://user-images.githubusercontent.com/31917400/47972179-99cc8f80-e091-11e8-9044-9561f334188b.jpg" />

 - Markov Chain simplifies this expression by using `Markov_assumption`. In this assumption, given the entire past history, *the probability distribution for the RandomVariable at the next time step only depends on the **current variable**. 
 <img src="https://user-images.githubusercontent.com/31917400/47972176-9933f900-e091-11e8-99c5-505818e9e11c.jpg" />

 - For all `t = 2, 3,..., n`. Under this assumption, the first expression can be rewrite as:
 <img src="https://user-images.githubusercontent.com/31917400/47972177-9933f900-e091-11e8-9f64-c6f0488f64df.jpg" />

 - which is much simpler than the original. It consists of an initial distribution for the first variable `p(X1)` and and n−1 other transition probabilities. We usually make one more assumption: that *these transition probabilities do not change with time*. Hence, the transition from time `t` to time `t+1` depends only on the value of Xt. 

### > Discrete Markov Chain
   - Suppose you have a secret number between {`1`, `2`, `3`, `4`, `5`}. We will call it your **initial number** at step 1. 
   - Now for each time step, your secret number will change by:
     - Assume that the coin is fair, so that with each flip, the probability of `H` and `T` are both `0.5`.
     - Flip a coin.
     - If the coin turns up `H`, then increase(+) your secret number by one. for example...`3->4`, `4->5`, `5->1`
     - If the coin turns up `T`, then decrease(-) your secret number by one. so...`3->2`, `2->1`, `1->5`
     - Repeat n times, and record the evolving history of your secret number.
   - Before the experiment, we can think of the sequence of secret numbers as a sequence of random variables, each taking on a value in {`1`, `2`, `3`, `4`, `5`}.
   - Suppose your secret number is **currently `4`** and that the history of your secret numbers is `(2, 1, 2, 3)`. What is the probability that on the next step, your secret number will be 5? 
   - What about the other four possibilities? Because of the rules of this game, the probability of the next transition will depend only on the fact that your **current number is `4`**. `The numbers further back in your history are irrelevant`, so this is a Markov chain.
   - If we assume that transition probabilities `0.5` do not change with time, then there are a total of 5^2=25 potential transition probabilities. Potential transition probabilities would be from State 1 to State 2, or from State 1 to State 3,... and so forth. These **transition probabilities** can be arranged into a matrix. 
<img src="https://user-images.githubusercontent.com/31917400/47991109-1d619d00-e0e1-11e8-9a62-d4346ff60abb.jpg" />

```
Q = matrix(c(0.0, 0.5, 0.0, 0.0, 0.5,  #for state_1
             0.5, 0.0, 0.5, 0.0, 0.0,  #for state_2
             0.0, 0.5, 0.0, 0.5, 0.0,  #for state_3
             0.0, 0.0, 0.5, 0.0, 0.5,  #for state_4
             0.5, 0.0, 0.0, 0.5, 0.0), #for state_5
           nrow=5, byrow=TRUE)
Q %*% Q
(Q %*% Q)[1,3]
```

<img src="https://user-images.githubusercontent.com/31917400/47990508-94963180-e0df-11e8-9da3-a0837dadfe79.jpg" />
Therefore, if your secret number is currently 1, the probability that the number will be 3 two steps from now is `0.25`.

### > Stationary Distribution of the Markov chain (discrete)
   - Suppose we want to know the probability distribution of the your secret number in the **distant future**. 
   <img src="https://user-images.githubusercontent.com/31917400/47998647-05e1de80-e0f8-11e8-83b6-f8ac0e74869e.jpg" />
   
   - Let’s calculate this for a few different values of h
   ```
   Q5 = Q %*% Q %*% Q %*% Q %*% Q   # h=5 steps in the future
   
   Q10 = Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q %*% Q   # h=10 steps in the future
   
   Q30 = Q
   for (i in 2:30) {
     Q30 = Q30 %*% Q
   }                  # h=30 steps in the future

   ```
   <img src="https://user-images.githubusercontent.com/31917400/47999841-a4bc0a00-e0fb-11e8-90bb-bebb89e6287d.jpg" />
   
   - Notice that as the future horizon gets more distant, the transition distributions appear to `converge`. 
     - The state you are currently in becomes less important in determining the more distant future. 
     - If we let `h` get really large, and take it to the limit, all the rows of the long-range transition matrix will become equal to `(0.2, 0.2, 0.2, 0.2, 0.2)`. That is, if you run the Markov chain for a very long time, the probability that you will end up in any particular state becomes **the same** for each of the states. This is what is called the **stationary distribution** of the Markov chain.
     - The stationary distribution of a chain is the `initial state distribution` for which **performing a transition will not change the probability of ending up in any given state**. 
     - `π*p = π` stationary*transition=stationary
     - `c(0.2, 0.2, 0.2, 0.2, 0.2) %*% Q` will give `[1,] 0.2  0.2  0.2  0.2  0.2`. One consequence of this property is that once a chain reaches its stationary distribution, the stationary distribution will remain the **distribution of the states**(transition matrix) thereafter.
   - Let's demonstrate the stationary distribution by simulating a long chain from this example.
   ```
   n = 5000
   x = numeric(n)
   x[1] = 1 # fix the state as 1 for time 1
   for (i in 2:n) {
     x[i] = sample.int(5, size=1, prob=Q[x[i-1],]) 
     # draw the next state from the intergers 1 to 5 with probabilities from the transition matrix Q, based on the previous value of X.
   }
   
   table(x) / n
   ```
   - it will give `1: 0.1996, 2: 0.2020, 3: 0.1980, 4: 0.1994, 5: 0.2010`. The overall distribution of the visits to the states is approximately equal to the stationary distribution. 

> As we have just seen, if you simulate a Markov chain for many iterations, the samples can be used as a Monte Carlo sample **from the stationary distribution**(each state -sample- has the same probability - iid?). This is exactly how we are going to `use Markov chains for Bayesian inference`. In order to simulate from a complicated posterior distribution, we will set up and run a Markov chain **whose stationary distribution is the posterior distribution**(each state -sample- follows the same probability distribution). 
- It is important to note that the stationary distribution doesn’t always exist for any given Markov chain. The Markov chain must have certain properties, which we won’t discuss here. However, the Markov chain algorithms we’ll use in future lessons for Monte Carlo estimation are guaranteed to produce stationary distributions.

### > Continuous Markov Chain
<img src="https://user-images.githubusercontent.com/31917400/47972178-9933f900-e091-11e8-9b5d-8654b21cea47.jpg" /> That is, the probability distribution for the next state is **Normal** with **variance = 1** and **mean = the current state**. This is often referred to as a “random walk.” Clearly, it is a Markov chain because the transition to the next state Xt+1 only depends on the current state Xt.
   ```
   n = 100
   x = numeric(n)
   
   for (i in 2:n) {
     x[i] = rnorm(1, mean=x[i-1], sd=1.0)
   }
   
   plot.ts(x)
   ```
<img src="https://user-images.githubusercontent.com/31917400/47972356-52df9980-e093-11e8-8732-e94ba0e2a862.jpg" />

The continuous "random walk" example here does not have a stationary distribution.

### > Stationary Distribution of the Markov chain (continuous)
However, we can modify it so that it does have a stationary distribution. Let the transition distribution be
<img src="https://user-images.githubusercontent.com/31917400/48024145-bae7bb80-e137-11e8-8286-379299ba9c1a.jpg" /> That is, the probability distribution for the next state is Normal with **variance** `1` and **mean** equal to `ϕ * current_mean`. As long as ϕ is between −1 and 1, then the **stationary distribution** will exist for this model. 

 - Let’s simulate this chain for ϕ=−0.6.
```
n = 1500
x = numeric(n)
phi = -0.6

for (i in 2:n) {
  x[i] = rnorm(1, mean=phi*x[i-1], sd=1.0)
}

plot.ts(x)
```
<img src="https://user-images.githubusercontent.com/31917400/48024547-ddc69f80-e138-11e8-90c8-9c97906bb35b.jpg" /> The theoretical stationary distribution for this chain is normal with **mean** `0` and **variance** `1 / (1 − ϕ^2)`, which in our example approximately equals `1.562`. Let’s look at a histogram of our chain and compare that with the theoretical stationary distribution. Oh, my..they are in sync!
```
hist(x, freq=FALSE)  # 'mean': ϕ * current_mean, 'var': 1
curve(dnorm(x, mean=0.0, sd=sqrt(1.0/(1.0-phi^2))), col="red", add=TRUE) # 'mean': 0, 'var': 1/(1-ϕ^2)
legend("topright", legend="theoretical stationary\ndistribution", col="red", lty=1, bty="n")
```
<img src="https://user-images.githubusercontent.com/31917400/48024783-680f0380-e139-11e8-86d9-ff224089c34e.jpg" /> Can you see that? 
 - Our MC chain(histogram): Normal with **variance** `1`, **mean** `ϕ * current_mean` (current_mean always change..)
 - perhaps posterior(target distribution): Normal with **mean** `0`, **var** `1 / (1 - ϕ^2)`
 
It appears that the chain has reached the stationary distribution. Therefore, we could treat this simulation from the chain like a Monte Carlo sample from the stationary distribution, a normal with mean `0` and variance `1.562`. Because most posterior distributions we will look at are continuous, our Monte Carlo simulations with Markov chains will be similar to this example.

## 2. Metropolis Hastings
Metropolis_Hastings algorithm allows us to sample from a **generic probability distribution**(target distribution), even if we don't know the `normalizing constant`(the bottom marginal stuff -the data-probability distribution- in Bayes theorem) because perhaps it is difficult to integrate. To do this, we construct and sample from a `Markov chain` whose **stationary distribution** is the target distribution that we're looking for. 
 - It consists of picking an arbitrary starting value and then iteratively accepting or rejecting candidate samples drawn from another distribution, one that is easy to sample. 
 - Let's say we want to produce samples from a target distribution called `P(θ)`, but we only know it up to a `normalizing constant`(denominator) or `proportionality`(nominator) so what we have is `P(θ) ∝ g(θ)`.
 <img src="https://user-images.githubusercontent.com/31917400/48061890-4dc83a80-e1b8-11e8-9d8d-0e7359918875.jpg" />

 - However, you still may want to have `q( )` have a larger variance than `P( )`, and see some rejection of candidates to be as an assurance that `q( )` is covering the space well. 
   - A high acceptance rate for random walk Metropolis-Hastings samplers is not a good thing. If the random walk is taking too small of steps, it will accept candidates often, but will take a very long time to fully explore the posterior distribution. 
   - If the random walk is taking too large of steps, many of its proposals will have low probability and the acceptance rate will be low. That will cause us to waste many of the draws. 
   - Ideally, a random walk sampler should accept somewhere between `23% to 50%` of the candidates proposed. 
### Example for Metropolis Hastings (discrete MarkovChain)
 - You are given a single coin, but don't know if the coin is fair(0.5, 0.5) or loaded(0.7, 0.3). People said, 6 out of 10 get a loaded coin. You flipped the coin 5 times and got H,T,H,T,T. Given this result, what's the posterior probability `P(θ|x)` you got messed(you were flipping a loaded coin)?  
<img src="https://user-images.githubusercontent.com/31917400/48072702-0d2aea00-e1d5-11e8-9f7e-d7b9acca60c1.JPG" />

 - Above is a simple example. What if we had a more complicated problem, where we couldn't work this all out in closed form? We'll know the likelihood and the prior, but we may not be able to get this `normalizing constant`. Can we instead do this by simulation? 
 - Let's use Metropolis Hastings in MCMC. 
   - Set up a Markov chain whose `equilibrium distribution` has this posterior distribution.
   - Consider a Markov chain with two states: `θ=fair` and `θ=loaded`. And we'll allow the chain to move between those two states, with certain **transition probabilities**. 
 - We set this up using the Metropolis–Hastings algorithm.
<img src="https://user-images.githubusercontent.com/31917400/48143888-b3dfbb00-e2a7-11e8-88c1-d0950b84f3fa.JPG" />

### Example for Metropolis Hastings (continuous MarkovChain)
Let's talk about a models that don't have nice, clean posterior distributions. 
 - Example of a one parameter model that is not conjugate:
   - Suppose we have values(data) that represent the percentage `μ` change in total personnel from last year to this year for, we'll say, ten companies `n=10` coming from a particular industry. We're going to assume for now, that these are independent measurements from a **normal** with a known `variance = 1`, but an unknown mean `μ`. In this case, the unknown mean could represent growth for this particular industry. It's the average of the growth of all the different companies. The **small variance** between the `companies and **percentage growth**` might be appropriate if the industry is stable. 
   - Suppose we decide that our prior believes about `μ` are better reflected using a **standard t-distribution** with `df=1`. This particular prior distribution(`t(0,1,1)`: Cauchy Distribution) has heavier(fatter) tails than the conjugate normal distribution, which can better accommodate the possibility of extreme values for `μ`. It is centered on `0` so that in our prior, there is a 50% chance that the growth is positive and a 50% chance that the growth is negative.  
   - We used a **normal-distribution likelihood** with `known variance` and **t-distribution prior** on the `unknown mean`.
<img src="https://user-images.githubusercontent.com/31917400/48153863-26a86080-e2bf-11e8-9f95-33e80fc867e2.jpg" />

### Because this model is not conjugate, the posterior distribution does not have a standard form that we can easily sample. 
 - To get posterior samples, we're going to need to setup a `Markov chain`, who's stationary distribution is the posterior distribution we want. <img src="https://user-images.githubusercontent.com/31917400/48151892-6caef580-e2ba-11e8-8cee-3566487bcc61.jpg" />















































