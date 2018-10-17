# Study-V003-MCMC-Python-R

> Random_Variable
<img src="https://user-images.githubusercontent.com/31917400/47053287-20cdcc80-d1a4-11e8-98cb-11d63c2108c6.png" />

> Distribution families
<img src="https://user-images.githubusercontent.com/31917400/47054048-be76cb00-d1a7-11e8-8fbd-db042248333e.png" />

> Estimating the AVG & VAR
<img src="https://user-images.githubusercontent.com/31917400/47054380-5923d980-d1a9-11e8-8226-a2a5b1991117.png" />

> Likelihood
<img src="https://user-images.githubusercontent.com/31917400/47054381-5b863380-d1a9-11e8-98f0-08cd8524d5f0.png" />

> Estimating Maximum Likelihood 
<img src="https://user-images.githubusercontent.com/31917400/47054860-3d6e0280-d1ac-11e8-948e-da71c188bb54.png" />

## Intro to Monte-Carlo
Monte-Carlo methods are methods for `generating random variables` directly or indirectly from a target distribution. 
 - **Applications**: 
   - 1> hypothesis testing 
   - 2> Bayesian computation
 - Monte Carlo is not strictly Bayesian nor strictly frequentist(Neyman-Pearson hypothesis testing), but it's **Parametric** sampling or bootstrapping (target distribution part of a family of distributions) as a **statistical simulation**, which helps us understand reality without needing to run multiple experiments or costly calculations.
 - It also helps us understand probabilities of **extreme events**. 
 - The use of Monte-Carlo methods to **calculate p-values** has become popular because:
   - many test statistics do not have a standard asymptotic distribution
   - even if a standard asymptotic distribution does exist, it may not be reliable in realistic sample sizes
   - In contrast, Monte Carlo methods can be used to obtain an **empirical p-value** that approximates the exact p-value without relying on asymptotic distributional theory or exhaustive enumeration. 

__Example 1> hypothesis testing (in R)__
> (a)Contingency table with too small sample_size
```
study = matrix(data = c(21, 2, 15, 3), nrow = 2, ncol = 2, byrow = TRUE,
               dimnames = list(c("surgery", "radiation"), c("controlled", "not controlled")))
```
<img src="https://user-images.githubusercontent.com/31917400/47085511-bf8f1300-d20e-11e8-91f8-d8b517cc484a.png" />

 - The contingency table above shows results of patients who underwent cancer **treatment** and either saw their cancer **controlled or not**. The two treatments are "surgery" and "radiation". 
 - The question we are interested in is **"Is there a difference between treatment and controll of cancer?""**. Of course, a Chi-squared test would usually be used for this type of analyses.
   - There are two ways that the Chi-squared test is used:
     - 1> test the **Goodness of fit** of the `theoretical distribution` to the `observations`.
     - 2> test for **independence** between different factors(row, col)??

 - **A disadvantage of the Chi-squared test** is that it requires a `sufficient sample size` in order for the chi-square approximation to be valid. When `cell_counts` are low (below 5), asymptotic properties do not hold well. Therefore, a simple Chi-squred test may report an **invalid p-value** which would increase a **Type-I, II error** rate. 
 - Monte-Carlo Method can solve this issue. 
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

 > (b)Inference on a single proportion
 - A simple random sample of 1,028 US adults in March 2013 found that 56% support nuclear arms reduction. Does this provide **convincing evidence** that a majority of Americans supported nuclear arms reduction at the 5% significance level?
 - Using a **Pearson-frequentist perspective**, we might simply do the following:
   - the number of US people supporting nuclear arms follows ~ B(n, p), and follows ~ N(np, npq)
   - the proportion of them follows ~ N(p, pq/n)
<img src="https://user-images.githubusercontent.com/31917400/47097156-b743d180-d228-11e8-975d-36c5e808ec3a.png" />













__Example 2> Bayesian Computation (in R)__























































































