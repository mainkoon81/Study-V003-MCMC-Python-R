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
 - Example 1> hypothesis testing (in R)
   - The contingency table below shows results of patients who underwent cancer **treatment** and either saw their cancer **controlled or not**. The two treatments are "surgery" and "radiation". 
   - The question we are interested in is **"Is there a difference between treatment and controll of cancer?""**
```
study = matrix(data = c(21, 2, 15, 3), nrow = 2, ncol = 2, byrow = TRUE,
               dimnames = list(c("surgery", "radiation"), c("controlled", "not controlled")))
```
<img src="https://user-images.githubusercontent.com/31917400/47085511-bf8f1300-d20e-11e8-91f8-d8b517cc484a.png" />



























































































