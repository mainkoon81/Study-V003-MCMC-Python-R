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

## Intro to Monte-Carlo
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
> This is the traditional way of Inferencing on a single proportion as a population parameter without Bayesian.
 - Let's say, a simple random sample of 1,028 US adults in March 2013 found that 56% support nuclear arms reduction. Damn + 6% !!! **"Q. Does this provide convincing evidence that a majority of Americans supported nuclear_arms_reduction at the 5% significance level?"**
 - Using a **Pearson-frequentist perspective**, we might simply do the following:
   - the number of US people supporting nuclear_arms_reduction follows ~ Bin(n, p), and follows ~ N(np, npq)
   - the proportion of them follows ~ N(p, pq/n)
<img src="https://user-images.githubusercontent.com/31917400/47097156-b743d180-d228-11e8-975d-36c5e808ec3a.png" />

 - Under the Null-hypothesis, `p0 = q0`, so `np0 = nq0`, then `1028*0.5 = 514 > 10` which is the **mean** number of people supporting nuclear_arms_reduction.
 - Based on the normal model, the test statistic can be computed as the Z-score of the point estimate: `Z = (p_obv - p0)/SE`.  
 - SE can be computed: `SE = sqrt(p0*q0/n)` and the Null-hypothesis `p0 = 0.5` is used again here because this is a hypothesis test for a single proportion `sqrt(0.5*0.5/1028) = 0.016`, so our Z is `(0.56-0.5)/0.016 = 3.75`.
 - p-value is `1 - pnorm(q=3.75) = 8.841729e-05`. We can then look up the upper tail area, the p-value, and see that it is less than 0.001. With a `p-value < 0.05`, we reject the null hypothesis and conclude that the **poll provides evidence that a majority (greater than 50%) of Americans supported the nuclear_arms_reduction**. The 95% CI for `p_obv` would be `0.56 + c(-1,1)*1.96*0.016`. So on AVG, around from 52% to 59% of US people support the nuclear_arms_reduction. 

> Another perspective on this problem(Inferencing on a single proportion as a population parameter) is that of a Bayesian. 
<img src="https://user-images.githubusercontent.com/31917400/47188651-d463c780-d330-11e8-85b3-5e5d5d34acc9.png" />

 - Look, we have a likelihood which is Binomial. 
 - Beta is conjugate to the Binomial distribution, i.e. `Beta` is a conjugate prior for the `Binomial` likelihood. That's why we choose Beta as our prior...then what the posterior will be?  
 - If the posterior distributions `p(θ|x)` are in the same **distribution family** as the prior distribution `p(θ)`:
   - the prior and posterior are then called **conjugate distributions**, 
   - the `**likelihood function**` is usually well-determined from a statement of the data-generating process.
 - Let's see our prior. Prior expresses one's beliefs about this quantity before some evidence is taken into account. Here, the prior could be the probability distribution representing the relative proportions of advocaters who will support nuclear_arms_reduction. 
   - we chose Beta(1,1) as our prior and this is equivalent to Unif(0,1) and this is a non-informative prior, which means we don't have any prior information to add to this model.
   - __Non-informative prior__
     - Non-informative prior expresses vague or general information about a variable.
     - Uninformative priors can express "objective" information such as "the variable is positive" or "the variable is less than some limit". The simplest and oldest rule for determining a non-informative prior is the principle of indifference, which assigns equal probabilities to all possibilities. In parameter estimation problems, the use of an uninformative prior typically yields results which are not too different from conventional statistical analysis, as the likelihood function often yields more information than the uninformative prior.
     - However, the non-informative prior can be created to reflect a balance among outcomes when no information is available. 
     - Or Priors can also be chosen according to some principle, such as symmetry or maximizing entropy given constraints; examples are `Jeffreys' prior` for the Bernoulli random variable. 
     
   - __Informative prior__
     - An informative prior expresses specific, definite information about a variable.
     - A reasonable approach is to make the prior a `normal distribution` with expected value equal to the given mean value, with variance equal to the given variance. 
     - pre-existing evidence which has already been taken into account is part of the prior and, as more evidence accumulates, the posterior is determined largely by the evidence rather than any original assumption, provided that the original assumption admitted the possibility of what the evidence is suggesting.
        
   - When a family of `conjugate priors` exists, choosing a prior from that family simplifies calculation of the posterior distribution.
   - `Parameters` of prior distributions are a kind of `hyperparameter`. For example, if one uses `Beta(a,b)` to model the distribution of the parameter `p` of Bernoulli, then:
     - `p` is a parameter of the underlying system (Bernoulli), and
     - `a` and `b` are parameters of the prior distribution (Beta); hence hyperparameters
   - Hyperparameters themselves may have `hyper-prior distributions` expressing beliefs about their values. A Bayesian model with more than one level of prior like this is called a `hierarchical Bayes model`.  

table of conjugate distribution
<img src="https://user-images.githubusercontent.com/31917400/47116021-53390180-d258-11e8-98c7-fa14a36415fe.png" />

 - find the posterior
<img src="https://user-images.githubusercontent.com/31917400/47119888-9c8f4e00-d264-11e8-9846-a03b7cb95e4d.png" />

 - the Bayesian Data model is `y|θ ~ Bin(n,θ)` thus `θ ~ Beta(a,b)`
 - the resulting posterior is then `θ|y ~ Beta(y+a, n-y+b)`. **We can now simulate the posterior distribution**, to choose `θ` !
 - Did you find the posterior? then build a Credible Interval. 
   - In Confidence Interval, the true value(population_parameter) is not a random variable. It is a fixed but unknown quantity. In contrast, our estimate is a random variable as it depends on our data x. Thus, we get different estimates each time, repeating our study. 
   - In Credible Intervals, we assume that the true value(population parameter) is a random variable. Thus, we capture the uncertainty about the true parameter value by a **imposing a prior distribution** on the true parameter vector. Using bayes theorem, we construct the posterior distribution for the parameter vector by blending the prior and the data(likelihood) we have, then arrive at **a point estimate** using the posterior distribution(use the mean of the posterior for example). However, the true parameter vector is a random variable, we also want to know the extent of uncertainty we have in our point estimate. Thus, we construct a 95% credible interval such that the following holds: `P( l(θ)<=θ<=u(θ) ) = 0.95` 
```
N = 10^4
set.seed(123)
x = rbeta(n = N, shape1 = 576 + 1, shape2 = 1028 - 576 + 1)
d = density(x)
hist(x = x, probability = TRUE, main = "Beta Posterior Distribution",
     xlab = expression(theta), ylab = "Density",
     ylim = c(0,40), col = "gray", border = "white")
lines(x = d$x , y = d$y, type = "l", col = 2)
abline(v = median(x), lty = 3, col = "3")

print("Median: ")
print(quantile(x = x, probs = c(0.025, 0.5, 0.975)))
```
<img src="https://user-images.githubusercontent.com/31917400/47120410-7e2a5200-d266-11e8-9613-75a2fe8e6995.png" />

Here is the thing. This example focused on **direct simulation from a posterior distribution**. However, there are some posteriors that will not be as easily identifiable. 
### Monte-Carlo methods will be helpful for generating samples from difficult to sample target distributions.
How? by generating random number from target distributions through **transformation methods**??
 - Monte Carlo simulation uses random sampling and statistical modeling to mimic the operations of complex systems. 
   - Monte Carlo models a system as a series of probability density functions.
   - Monte Carlo repeatedly samples from the PDF.
   - Monte Carlo computes the statistics of interest.

------------------------------------------------------------------------------------------------------------
## Generating Random Variables






















































































