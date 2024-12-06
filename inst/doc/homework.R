## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
set.seed(0)
GF <- function(f,a,b){
  x <- sort(runif(1000 ,min = a,max = b))
  y <- f(x)
  plot(x,y,type = "l")
}

f1 <- function(x) x^2
GF(f1 ,-1,1)

## -----------------------------------------------------------------------------
A <- matrix(c(1,2,3,4),2,2)
B <- matrix(c(0,1,1,0),2,2)
A%*%B

## -----------------------------------------------------------------------------
set.seed(0)
N <- 1000
x <- runif(N)
y <- runif(N)
z <- x^2 + y^2
cat("pi =",sum(z<1)/N*4)

par(pty="s")
plot(x,y,col=ifelse((z) < 1, "red","black"),pch=16)
lx <-seq(from = 0, to = 1,length=N)
lines(lx ,sqrt(1-lx^2))

## -----------------------------------------------------------------------------
set.seed(0)
rRayleigh <- function(n,sigma)
{
  x <- runif(n)
  r <- sqrt(-2*sigma^2*log(1-x))
  return(r)
}
hist(rRayleigh(1000,1),main = "σ=1")
hist(rRayleigh(1000,2),main = "σ=2")
hist(rRayleigh(1000,3),main = "σ=3")
hist(rRayleigh(1000,4),main = "σ=4")


## -----------------------------------------------------------------------------
fun <- function(p){
  set.seed(1)
  x1 <- rnorm(1000)
  x2 <- rnorm(1000, 3, 1)
  k <- sample(0:1, size = 1000, replace = TRUE, prob = c(p,1-p))
  x <- (1-k)*x1 + (k)*x2
  hist(x, 
       probability = TRUE,
       main = bquote('p='*.(p))
       )
  lines(density(x), lwd = 2, col = "red")
}

fun(0.75)
fun(0.1)
fun(0.3)
fun(0.5)
fun(0.7)
fun(0.9)

## -----------------------------------------------------------------------------
set.seed(0)
rPoi_Gamma <- function(n,t,l,a,b)
{
  Nt <- rpois(n,t*l)
  x <- c()
  for (i in 1:n) {
    y <- rgamma(Nt[i],shape=a,rate=b)
    x[i] <- sum(y)
  }
  return(x)
}

t <- 10
l <- c(1,2,3) #lambda of poisson distribution
a <- c(1,2,3) #alpha of gamma distribution
b <- c(3,2,1) #beta of gamma distribution
tm <- l*t*a/b #theoretical mean
tv <- l*t*(a/b^2+(a/b)^2) #theoretical variance
rm <- rv <- c() #real mean and variance
for (i in 1:3) {
  x <- rPoi_Gamma(1000,t,l[i],a[i],b[i])
  rm[i] <- mean(x)
  rv[i] <- var(x)
}
rbind(tm,rm,tv,rv)

## -----------------------------------------------------------------------------
set.seed(0)
N <- 1e5
cdfbeta <- function(x)
{
  s <- rbeta(N,3,3)
  p <- mean(s <= x)
  return(p)
}
x1 <- x2 <- numeric(9)
for (i in 1:9) {
  x1[i] <- cdfbeta(0.1*i)
  x2[i] <- pbeta(0.1*i,3,3)
}
result <-rbind(x1,x2)
row.names(result) <- c("MC_estimation","pbeta")
colnames(result) <- 1:9*0.1
result


## -----------------------------------------------------------------------------
rRay_av <- function(n,sigma)
{
  u <- runif(n/2)
  x1 <- sqrt(-2*sigma^2*log(1-u))
  x2 <- sqrt(-2*sigma^2*log(u))
  return(c(x1,x2))
}

## -----------------------------------------------------------------------------
set.seed(0)
u <- runif(N)
mean(sqrt(log(1-u)*log(u)))

## -----------------------------------------------------------------------------
 # Quick sort algorithm:
 quickSort <- function(arr) {
  # Pick a number at random.
   mid <- sample(arr, 1)

   # Place-holders for left and right values.
   left <- c()
   right <- c()

  # Move all the smaller values to the left, bigger values to the right.
   lapply(arr[arr != mid], function(d) {
    if (d < mid) {
       left <<- c(left, d)
     }
     else {
      right <<- c(right, d)
     }
   })

   if (length(left) > 1) {
     left <- quickSort(left)
   }

  if (length(right) > 1) {
    right <- quickSort(right)   }

   # Finally, return the sorted values.
  c(left, mid, right)
 }

## -----------------------------------------------------------------------------
set.seed(0)
N <- 100
n <- c(1,2,4,6,8)*1e2 # 4次幂时间过长，改为3次
a <- numeric(5)
for (j in 1:5) {
  t <- numeric(N)
  for (i in 1:N) {
    x <- sample(1:n[j])
    st <- Sys.time()
    nop <- quickSort(x)
    et <- Sys.time()
    t[i] <- et-st
  }
  a[j] <- mean(t)
}

## -----------------------------------------------------------------------------
x <- n*log(n)
m <- lm(a~x)
m
c <- m$coefficients
x1 <- 1:8e2
plot(n,a)
lines(c[2]*x1*log(x1)+c[1])

## -----------------------------------------------------------------------------
sk <- function(x){ #计算样本偏度
  xbar <- mean(x)
  m3 <- mean((x-xbar)^3)
  m2 <- mean((x-xbar)^2)
  return(m3/m2^1.5)
}

## -----------------------------------------------------------------------------
set.seed(0)
n <- m <- 1000
x <- replicate(n,sk(rnorm(m)))
quantile(x, probs = c(0.025,0.05,0.95,0.975))#估计分位数
qnorm(c(0.025,0.05,0.95,0.975), mean = 0, sd = sqrt(6/n))#大样本近似

## -----------------------------------------------------------------------------
library(MASS)
#二元正态
alpha <- 0.05 #置信水平
rho <- 0.8 #相关系数
size <- 5:30
m <- 100 #蒙特卡洛重复次数
f <- function(rho, size)
{
  set.seed(0)
  pearson <- pearson.sd <- kendall <- kendall.sd <- spearman <- spearman.sd <- c()
  for(n in size){
    result <- replicate(m, {
      S <- matrix(c(1, rho, rho, 1), nrow = 2)
      x <- mvrnorm(n = n, Sigma = S, mu = c(0, 0))
      c(cor.test(x[,1], x[,2], method = "pearson")$p.value < alpha,
      cor.test(x[,1], x[,2], method = "kendall")$p.value < alpha,
      cor.test(x[,1], x[,2], method = "spearman")$p.value < alpha)
    })
    pearson[n] <- mean(result[1,])
    #pearson.sd[n] <- sd(result[1,])
    kendall[n] <- mean(result[2,])
    #kendall.sd[n] <- sd(result[2,])
    spearman[n] <- mean(result[3,])
    #spearman.sd[n] <- sd(result[3,])
  }
  plot(pearson, xlim = c(size[1],size[length(size)]), ylim = c(0,1),
       main = bquote(rho == .(rho)),
       xlab = "样本量",
       ylab = "检验功效",
       col = "red", pch = 16)
  points(kendall, col = "blue", pch = 16)
  points(spearman, col = "green", pch = 16)
  lines(pearson, col = "red")
  lines(kendall, col = "blue")
  lines(spearman, col = "green")
  legend("bottomright", 
         col = c("red", "blue", "green"), 
         lwd = rep(1,3), 
         pch = rep(16,3),
         legend = c("pearson", "kendall", "spearman"))
}

## -----------------------------------------------------------------------------
f(0.8, 5:30)

## -----------------------------------------------------------------------------
f(0.5, 5:60)

## -----------------------------------------------------------------------------
f(0.3, 20:120)

## -----------------------------------------------------------------------------
#对Y取倒数
rho <- 0.8
size <- 5:100
set.seed(0)
pearson <- pearson.sd <- kendall <- kendall.sd <- spearman <- spearman.sd <- c()
for(n in size){
  result <- replicate(m, {
    S <- matrix(c(1, rho, rho, 1), nrow = 2)
    x <- mvrnorm(n = n, Sigma = S, mu = c(0, 0))
    x[,2] <- 1/x[,2]
    c(cor.test(x[,1], x[,2], method = "pearson")$p.value < alpha,
    cor.test(x[,1], x[,2], method = "kendall")$p.value < alpha,
    cor.test(x[,1], x[,2], method = "spearman")$p.value < alpha)
  })
  pearson[n] <- mean(result[1,])
  #pearson.sd[n] <- sd(result[1,])
  kendall[n] <- mean(result[2,])
  #kendall.sd[n] <- sd(result[2,])
  spearman[n] <- mean(result[3,])
  #spearman.sd[n] <- sd(result[3,])
}
plot(pearson, xlim = c(size[1],size[length(size)]), ylim = c(0,1),
     main = "",
     xlab = "样本量",
     ylab = "检验功效",
     col = "red", pch = 16)
points(kendall, col = "blue", pch = 16)
points(spearman, col = "green", pch = 16)
lines(pearson, col = "red")
lines(kendall, col = "blue")
lines(spearman, col = "green")
legend("bottomright", 
       col = c("red", "blue", "green"), 
       lwd = rep(1,3), 
       pch = rep(16,3),
       legend = c("pearson", "kendall", "spearman"))


## -----------------------------------------------------------------------------
library(boot)
set.seed(0)
m <- 10000
n_null <- 950
n_alt <- 50
alpha <- 0.1

fwer_bonferroni <- numeric(m)
fdr_bonferroni <- numeric(m)
tpr_bonferroni <- numeric(m)
fwer_bh <- numeric(m)
fdr_bh <- numeric(m)
tpr_bh <- numeric(m)

for (i in 1:m) {
  # Generate p-values
  p_null <- runif(n_null)  # Null hypothesis p-values
  p_alt <- rbeta(n_alt, 0.1, 1)  # Alternative hypothesis p-values
  p_values <- c(p_null, p_alt)  # Combine all p-values

  # Bonferroni correction
  p_bonferroni <- p.adjust(p_values, method = "bonferroni")
  reject_bonferroni <- p_bonferroni < alpha
  
  # FWER and TPR calculation for Bonferroni
  fwer_bonferroni[i] <- mean(reject_bonferroni[1:n_null])  # Type I error rate
  tpr_bonferroni[i] <- sum(reject_bonferroni[(n_null + 1):(n_null + n_alt)]) / n_alt  # True Positive Rate

  # FDR calculation for Bonferroni
  num_rejected_bonferroni <- sum(reject_bonferroni)
  num_false_rejections_bonferroni <- sum(reject_bonferroni[1:n_null])
  if (num_rejected_bonferroni > 0) {
    fdr_bonferroni[i] <- num_false_rejections_bonferroni / num_rejected_bonferroni
  } else {
    fdr_bonferroni[i] <- 0
  }

  # B-H correction
  p_bh <- p.adjust(p_values, method = "BH")
  reject_bh <- p_bh < alpha

  fwer_bh[i] <- mean(reject_bh[1:n_null])  # Type I error rate
  num_rejected_bh <- sum(reject_bh)  # Total rejections
  num_true_rejections_bh <- sum(reject_bh[(n_null + 1):(n_null + n_alt)])  # True positives
  num_false_rejections_bh <- sum(reject_bh[1:n_null])  # False positives

  # FDR calculation for B-H
  if (num_rejected_bh > 0) {
    fdr_bh[i] <- num_false_rejections_bh / num_rejected_bh
  } else {
    fdr_bh[i] <- 0
  }

  # TPR calculation for B-H
  tpr_bh[i] <- num_true_rejections_bh / n_alt
}

mean_fwer_bonferroni <- mean(fwer_bonferroni)
mean_fdr_bonferroni <- mean(fdr_bonferroni)
mean_fdr_bh <- mean(fdr_bh)
mean_fwer_bh <- mean(fwer_bh)
mean_tpr_bonferroni <- mean(tpr_bonferroni)
mean_tpr_bh <- mean(tpr_bh)


# Create results matrix
results <- matrix(c(mean_fwer_bonferroni, mean_fdr_bonferroni, mean_tpr_bonferroni,mean_fwer_bh, mean_fdr_bh, mean_tpr_bh), nrow = 3)
colnames(results) <- c("Bonferroni", "B-H")
rownames(results) <- c("FWER", "FDR", "TPR")

results


## -----------------------------------------------------------------------------
set.seed(0)

data <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)

# Function to calculate the MLE of lambda
mle_lambda <- function(data, indices) {
  sample_data <- data[indices]  # Resample with replacement
  n <- length(sample_data)
  return(n / sum(sample_data))
}

# Calculate the MLE of lambda from the original data
original_lambda <- mle_lambda(data, 1:length(data))

# Perform bootstrap
n_boot <- 10000  # Number of bootstrap samples
boot_results <- boot(data, mle_lambda, R = n_boot)

# Bias and standard error
bias <- mean(boot_results$t) - original_lambda
se <- sd(boot_results$t)

# Summary results
cat("MLE of λ:", original_lambda, "\n")
cat("Bootstrap Bias:", bias, "\n")
cat("Bootstrap Standard Error:", se, "\n")


## -----------------------------------------------------------------------------
set.seed(0)
# Mean time between failures
mean_time_between_failures <- 1 / original_lambda

# Calculate the bootstrap estimates for mean time between failures
boot_mean_times <- 1 / boot_results$t

# Standard Normal Method CI
z <- qnorm(0.975)  
se <- sd(boot_mean_times)
ci_normal <- mean_time_between_failures + c(-z * se, z * se)

# Basic Method CI
ci_basic <- mean_time_between_failures + c(mean(boot_mean_times) - mean(boot_mean_times), mean(boot_mean_times) - mean(boot_mean_times))

# Percentile Method CI
ci_percentile <- quantile(boot_mean_times, c(0.025, 0.975))

# BCa Method CI
ci_bca <- boot.ci(boot_results, type = "bca")

# Extract BCa CI limits
bca_limits <- ci_bca$bca[4:5]  

# Combine results
results <- data.frame(
  Method = c("Standard Normal", "Basic", "Percentile", "BCa"),
  Lower_CI = c(ci_normal[1], ci_basic[1], ci_percentile[1], bca_limits[1]),
  Upper_CI = c(ci_normal[2], ci_basic[2], ci_percentile[2], bca_limits[2])
)

# Print results
print(results)

## -----------------------------------------------------------------------------
library(bootstrap)
scores <- scor
cov_matrix <- cov(scores)
eigen_values <- eigen(cov_matrix)$values
theta_hat <- eigen_values[1] / sum(eigen_values)

n <- nrow(scores)
jackknife_estimates <- numeric(n)

for (i in 1:n) {
  jackknife_sample <- scores[-i, ]
  jackknife_cov_matrix <- cov(jackknife_sample)
  jackknife_eigen_values <- eigen(jackknife_cov_matrix)$values
  jackknife_estimates[i] <- jackknife_eigen_values[1] / sum(jackknife_eigen_values)
}

bias <- (n - 1) * (mean(jackknife_estimates) - theta_hat)
jackknife_se <- sqrt((n - 1) / n * sum((jackknife_estimates - mean(jackknife_estimates))^2))

cat("Jackknife Bias: ", bias, "\n")
cat("Jackknife Standard Error: ", jackknife_se, "\n")

## -----------------------------------------------------------------------------
library(DAAG)
data(ironslag)

chemical <- ironslag$chemical
magnetic <- ironslag$magnetic

n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- e5 <- numeric(n)

for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]

  J1 <- lm(y ~ x)
  yhat1 <- predict(J1, newdata = data.frame(x = chemical[k]))
  e1[k] <- magnetic[k] - yhat1

  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- predict(J2, newdata = data.frame(x = chemical[k]))
  e2[k] <- magnetic[k] - yhat2

  J3 <- lm(log(y) ~ x)
  logyhat3 <- predict(J3, newdata = data.frame(x = chemical[k]))
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3

  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  yhat4 <- predict(J4, newdata = data.frame(x = chemical[k]))
  e4[k] <- magnetic[k] - yhat4
}

mse <- c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
print(mse)

adj_r2 <- c(
  summary(lm(magnetic ~ chemical))$adj.r.squared,
  summary(lm(magnetic ~ chemical + I(chemical^2)))$adj.r.squared,
  summary(lm(log(magnetic) ~ chemical))$adj.r.squared,
  summary(lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3)))$adj.r.squared
)

print(adj_r2)

best_model_mse <- which.min(mse)
best_model_adj_r2 <- which.max(adj_r2)

cat("MSE best model", best_model_mse, "\n")
cat("R-squared best model", best_model_adj_r2, "\n")

## -----------------------------------------------------------------------------
x <- c(158, 171, 193, 199, 230, 243, 248, 248, 250, 267, 271, 316, 327, 329)
y <- c(141, 148, 169, 181, 203, 213, 229, 244, 257, 260, 271, 309)

cramer_von_mises <- function(x, y) {
  n <- length(x)
  m <- length(y)
  combined <- c(x, y)
  ranks <- rank(combined)
  Rx <- sum(ranks[1:n])^2 / (n^2)
  Ry <- sum(ranks[(n + 1):(n + m)])^2 / (m^2)
  return((Rx + Ry) / (n + m) - (n^2 + m^2) / ((n + m)^2))
}

observed_stat <- cramer_von_mises(x, y)

set.seed(123)
n_permutations <- 1000
permutation_stats <- numeric(n_permutations)

combined_data <- c(x, y)
n <- length(x)
m <- length(y)

for (i in 1:n_permutations) {
  permuted_data <- sample(combined_data)
  perm_x <- permuted_data[1:n]
  perm_y <- permuted_data[(n + 1):(n + m)]
  permutation_stats[i] <- cramer_von_mises(perm_x, perm_y)
}

p_value <- mean(permutation_stats >= observed_stat)

cat("Observed Cramér-von Mises statistic:", observed_stat, "\n")
cat("P-value from permutation test:", p_value, "\n")


## -----------------------------------------------------------------------------
set.seed(0)
x <- rnorm(30)
y <- rnorm(30)

spearman_stat <- cor(x, y, method = "spearman")

n_permutations <- 1000
permutation_stats <- numeric(n_permutations)

for (i in 1:n_permutations) {
  permuted_y <- sample(y)
  permutation_stats[i] <- cor(x, permuted_y, method = "spearman")
}

p_value_perm <- mean(abs(permutation_stats) >= abs(spearman_stat))

cor_test_result <- cor.test(x, y, method = "spearman")
p_value_cor_test <- cor_test_result$p.value

cat("Observed Spearman statistic:", spearman_stat, "\n")
cat("P-value from permutation test:", p_value_perm, "\n")
cat("P-value from cor.test:", p_value_cor_test, "\n")


## -----------------------------------------------------------------------------
# Load necessary libraries
library(coda)

# Metropolis-Hastings sampler for the standard Cauchy distribution
metropolis_hastings_cauchy <- function(n, burn_in) {
  samples <- numeric(n)
  samples[1] <- rnorm(1)  # Initial value

  for (i in 2:n) {
    # Propose a new sample
    y <- rnorm(1, samples[i - 1], 1)  # Normal proposal
    r <- dcauchy(y) / dcauchy(samples[i - 1])  # Acceptance ratio
    if (runif(1) < min(1, r)) {
      samples[i] <- y  # Accept the new sample
    } else {
      samples[i] <- samples[i - 1]  # Reject and keep the old sample
    }
  }

  return(samples[-(1:burn_in)])  # Discard burn-in samples
}

# Run multiple chains for convergence check
set.seed(0)
n_samples <- 5000
burn_in <- 1000
num_chains <- 4
chains_cauchy <- lapply(1:num_chains, function(i) metropolis_hastings_cauchy(n_samples, burn_in))

# Check convergence using Gelman-Rubin statistic
R_cauchy <- function(chains) {
  m <- length(chains)  # Number of chains
  n <- sapply(chains, length)  # Length of each chain
  theta_hat <- sapply(chains, mean)  # Mean of each chain
  B <- var(theta_hat) * n[1] / (m - 1)  # Between-chain variance
  W <- mean(sapply(chains, var))  # Within-chain variance
  var_theta_hat <- (1 + 1 / m) * W + B / m  # Total variance estimate
  return(sqrt(var_theta_hat / W))  # Gelman-Rubin statistic
}

# Run until convergence
R_value <- R_cauchy(chains_cauchy)
timer <- 0
while (R_value >= 1.2) {
  timer <- timer + 1
  chains_cauchy <- lapply(1:num_chains, function(i) metropolis_hastings_cauchy(n_samples, burn_in))
  R_value <- R_cauchy(chains_cauchy)
  if(timer >= 100) break
}

# Calculate deciles of the samples
combined_samples <- unlist(chains_cauchy)
sample_deciles <- quantile(combined_samples, probs = seq(0.1, 0.9, by = 0.1))
theoretical_deciles <- qcauchy(seq(0.1, 0.9, by = 0.1))

# Print results
print(sample_deciles)
print(theoretical_deciles)

## -----------------------------------------------------------------------------
# Gibbs sampler for the bivariate density
gibbs_sampler_bivariate <- function(n, a, b) {
  x_samples <- numeric(n)
  y_samples <- numeric(n)

  # Initial values
  x_samples[1] <- rbinom(1, n, 0.5)  # Initial value for x
  y_samples[1] <- runif(1)            # Initial value for y

  for (i in 2:n) {
    # Sample from the conditional distributions
    x_samples[i] <- rbinom(1, n, y_samples[i - 1])
    y_samples[i] <- rbeta(1, x_samples[i] + a, n - x_samples[i] + b)
  }

  return(list(x = x_samples, y = y_samples))
}

# Run multiple chains for convergence check
set.seed(0)
n_samples_gibbs <- 1000
a <- 2
b <- 3
num_chains <- 4
chains_gibbs <- lapply(1:num_chains, function(i) gibbs_sampler_bivariate(n_samples_gibbs, a, b))

# Check convergence using Gelman-Rubin statistic
R_gibbs <- function(chains) {
  m <- length(chains)  # Number of chains
  n <- sapply(chains, length)  # Length of each chain
  theta_hat <- sapply(chains, mean)  # Mean of each chain
  B <- var(theta_hat) * n[1] / (m - 1)  # Between-chain variance
  W <- mean(sapply(chains, var))  # Within-chain variance
  var_theta_hat <- (1 + 1 / m) * W + B / m  # Total variance estimate
  return(sqrt(var_theta_hat / W))  # Gelman-Rubin statistic
}

# Run until convergence
R_x_value <- R_gibbs(lapply(chains_gibbs, `[[`, "x"))
R_y_value <- R_gibbs(lapply(chains_gibbs, `[[`, "y"))

timer <- 0
while (R_x_value >= 1.2 || R_y_value >= 1.2) {
  timer <- timer + 1
  chains_gibbs <- lapply(1:num_chains, function(i) gibbs_sampler_bivariate(n_samples_gibbs, a, b))
  R_x_value <- R_gibbs(lapply(chains_gibbs, `[[`, "x"))
  R_y_value <- R_gibbs(lapply(chains_gibbs, `[[`, "y"))
  if(timer >= 100) break
}
plot(chains_gibbs[[num_chains]]$x,chains_gibbs[[num_chains]]$y)


## -----------------------------------------------------------------------------
kth_term <- function(k, d, a) {
  norm_a <- sqrt(sum(a^2)) 
  
  gamma_d_half <- lgamma((d + 1) / 2)
  gamma_k_3_half <- lgamma(k + 1.5)
  gamma_k_d_half_plus_1 <- lgamma(k + (d / 2) + 1)
  
  coefficient <- ((-1)^k) / (factorial(k) * (2^k))
  power_term <- (norm_a^(2 * k + 2)) / ((2 * k + 1) * (2 * k + 2))
  
  gamma_ratio_log <- gamma_d_half + gamma_k_3_half - gamma_k_d_half_plus_1
  
  term <- coefficient * power_term * exp(gamma_ratio_log)
  return(term)
}


## -----------------------------------------------------------------------------
sum_series <- function(d, a, tol = 1e-5, max_k = 100) {
  sum <- 0
  k <- 0
  term <- kth_term(k, d, a)
  
  while (abs(term) > tol && k < max_k) {
    sum <- sum + term
    k <- k + 1
    term <- kth_term(k, d, a)
  }
  
  return(sum)
}

## -----------------------------------------------------------------------------
a <- c(1, 2)
d <- 1e4

result <- sum_series(d, a)
result


## -----------------------------------------------------------------------------
library(gmp)

integral <- function(c, k) {
  integrand <- function(u) {
    (1 + u^2 / k) ^ (-(k + 1) / 2)
  }
  return(integrate(integrand, 0, c)$value)
}

equation <- function(a, k) {
  c_k <- sqrt(a^2 * k / (k + 1 - a^2))
  c_k_minus_1 <- sqrt(a^2 * (k - 1) / (k - a^2))
  
  lhs <- (2 * gamma(k / 2) / sqrt(pi * (k - 1) * gamma((k - 1) / 2))) * integral(c_k_minus_1, k - 1)
  
  rhs <- (2 * gamma((k + 1) / 2) / (sqrt(pi * k) * gamma(k / 2))) * integral(c_k, k)
  
  return(lhs - rhs)
}

solve_for_a <- function(k) {
  result <- uniroot(equation, c(0.1, 1.5), k = k)
  return(result$root)
}

k <- 3
a_solution <- solve_for_a(k)
print(paste("Solution for a when k =", k, "is a =", a_solution))


## -----------------------------------------------------------------------------

Y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
tau <- 1

delta <- ifelse(Y < tau, 1, 0)

lambda <- mean(Y) 
tolerance <- 1e-6
max_iter <- 1000
diff <- 1
iter <- 0

# E-M algorithm
while (diff > tolerance && iter < max_iter) {
  iter <- iter + 1
  lambda_old <- lambda
  
  T_imputed <- ifelse(delta == 1, Y, tau + lambda)
  
  lambda <- mean(T_imputed)

  diff <- abs(lambda - lambda_old)
}

cat("Estimated lambda using E-M algorithm:", lambda, "\n")

## -----------------------------------------------------------------------------

delta <- ifelse(Y < tau, 1, 0)


n1 <- sum(delta)


S <- sum(delta * Y + (1 - delta) * tau)

lambda_mle <- S / n1

cat("Estimated lambda (MLE) considering censoring:", lambda_mle, "\n")

## -----------------------------------------------------------------------------

initialize_tableau <- function() {
  matrix(c(
    2, 1, 1, 1, 0, 2,   
    1, -1, 3, 0, 1, 3,  
    -4, -2, -9, 0, 0, 0
  ), nrow = 3, byrow = TRUE)
}


simplex_pivot <- function(tableau, pivot_row, pivot_col) {

  tableau[pivot_row, ] <- tableau[pivot_row, ] / tableau[pivot_row, pivot_col]

  for (i in seq_len(nrow(tableau))) {
    if (i != pivot_row) {
      tableau[i, ] <- tableau[i, ] - tableau[i, pivot_col] * tableau[pivot_row, ]
    }
  }
  return(tableau)
}


is_optimal <- function(tableau) {
  all(tableau[nrow(tableau), -ncol(tableau)] >= 0)
}


get_pivot <- function(tableau) {

  pivot_col <- which.min(tableau[nrow(tableau), -ncol(tableau)])
  

  ratios <- tableau[-nrow(tableau), ncol(tableau)] / tableau[-nrow(tableau), pivot_col]
  ratios[ratios <= 0] <- Inf 
  pivot_row <- which.min(ratios)
  
  return(list(pivot_row = pivot_row, pivot_col = pivot_col))
}

simplex <- function() {
  tableau <- initialize_tableau()
  
  while (!is_optimal(tableau)) {
    pivot <- get_pivot(tableau)
    tableau <- simplex_pivot(tableau, pivot$pivot_row, pivot$pivot_col)
  }
  
  
  solution <- rep(0, ncol(tableau) - 1)  # 初始化解
  for (i in seq_len(nrow(tableau) - 1)) {
    basic_var_index <- which(tableau[i, 1:(ncol(tableau) - 1)] == 1)
    if (length(basic_var_index) == 1) {
      solution[basic_var_index] <- tableau[i, ncol(tableau)]
    }
  }
  
  objective_value <- tableau[nrow(tableau), ncol(tableau)] * -1  # 原问题是最小化
  list(solution = solution, objective_value = objective_value)
}

result <- simplex()
cat("最优解:\n")
print(result$solution)
cat("最小目标值:\n")
print(result$objective_value)


## -----------------------------------------------------------------------------
data(mtcars)

formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

models_for <- list() 
for (i in seq_along(formulas)) {
  models_for[[i]] <- lm(formulas[[i]], data = mtcars)
}


cat("Models fitted using for loop:\n")
print(models_for)


models_lapply <- lapply(formulas, function(formula) lm(formula, data = mtcars))

cat("\nModels fitted using lapply():\n")
print(models_lapply)


## -----------------------------------------------------------------------------
data(mtcars)

set.seed(0)  
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})


models_for1 <- list()  
for (i in seq_along(bootstraps)) {
  models_for1[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}


cat("Models fitted using for loop:\n")
print(models_for1)


fit_model <- function(data) lm(mpg ~ disp, data = data)  # Define a named function
models_lapply1 <- lapply(bootstraps, fit_model)


cat("\nModels fitted using lapply():\n")
print(models_lapply1)


## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

sapply(models_for, rsq)
sapply(models_lapply, rsq)
sapply(models_for1, rsq)
sapply(models_lapply1, rsq)


## -----------------------------------------------------------------------------
set.seed(0)
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

sapply(trials, function(test) test$p.value)


## -----------------------------------------------------------------------------
parallel_apply <- function(fun, ..., FUN.VALUE) {
  inputs <- list(...)
  results <- Map(fun, inputs[[1]], inputs[[2]])
  vapply(results, identity, FUN.VALUE = FUN.VALUE)
}

#Example
set.seed(0)
list1 <- 1:5
list2 <- c(2, 4, 6, 8, 10)

sum_fun <- function(x, y) x + y

result <- parallel_apply(sum_fun, list1, list2, FUN.VALUE = numeric(1))
print(result)


## -----------------------------------------------------------------------------
fast_chisq_test <- function(x, y) {
  
  observed <- table(x, y)
  
  row_totals <- margin.table(observed, 1)
  col_totals <- margin.table(observed, 2)
  grand_total <- sum(observed)
  
  expected <- outer(row_totals, col_totals, FUN = "*") / grand_total
  chi_square_stat <- sum((observed - expected)^2 / expected)
  return(chi_square_stat)
}

# test
set.seed(1)
x <- sample(1:5, 1e6, replace = TRUE)
y <- sample(1:5, 1e6, replace = TRUE)

fast_chisq_time <- system.time({
  result_fast <- fast_chisq_test(x, y)
})
cat("fast_chisq_test 执行时间: ", fast_chisq_time[3], "秒\n")

chisq_time <- system.time({
  result_chisq <- chisq.test(x, y)
})
cat("chisq.test 执行时间: ", chisq_time[3], "秒\n")


## -----------------------------------------------------------------------------
set.seed(0)
fast_table <- function(x, y) {
  ux <- seq(min(x), max(x)) 
  uy <- seq(min(y), max(y))
  
  xi <- match(x, ux)
  yi <- match(y, uy)
  
  result <- matrix(0L, nrow = length(ux), ncol = length(uy))
  index <- (xi - 1L) * length(uy) + yi
  counts <- tabulate(index, nbins = length(ux) * length(uy))
  result[] <- counts
  dimnames(result) <- list(as.character(ux), as.character(uy))
  return(result)
}



x <- sample(1:5, 1e6, replace = TRUE)
y <- sample(1:5, 1e6, replace = TRUE)
fast_table_time <- system.time({
  result_fast <- fast_table(x, y)
})
cat("fast_table 执行时间: ", fast_table_time[3], "秒\n")

table_time <- system.time({
  result_table <- table(x, y)
})
cat("table 执行时间: ", table_time[3], "秒\n")

## -----------------------------------------------------------------------------
fast_chisq_test <- function(x, y) {
  
  observed <- fast_table(x, y)
  
  row_totals <- margin.table(observed, 1)
  col_totals <- margin.table(observed, 2)
  grand_total <- sum(observed)
  
  expected <- outer(row_totals, col_totals, FUN = "*") / grand_total
  chi_square_stat <- sum((observed - expected)^2 / expected)
  return(chi_square_stat)
}

# test
set.seed(1)
x <- sample(1:5, 1e6, replace = TRUE)
y <- sample(1:5, 1e6, replace = TRUE)

fast_chisq_time <- system.time({
  result_fast <- fast_chisq_test(x, y)
})
cat("fast_chisq_test 执行时间: ", fast_chisq_time[3], "秒\n")

chisq_time <- system.time({
  result_chisq <- chisq.test(x, y)
})
cat("chisq.test 执行时间: ", chisq_time[3], "秒\n")

## -----------------------------------------------------------------------------
# R-function
gibbs_sampler_r <- function(n, a, b, num_iter, n_init = 0) {
  x_chain <- numeric(num_iter)
  y_chain <- numeric(num_iter)
  
  x <- n_init
  y <- 0.5

  for (i in 1:num_iter) {
    x <- rbinom(1, n, y)
    y <- rbeta(1, x + a, n - x + b)
    
    x_chain[i] <- x
    y_chain[i] <- y
  }
  
  list(x = x_chain, y = y_chain)
}

## -----------------------------------------------------------------------------
library(Rcpp)
sourceCpp("gibbs_sampler.cpp")

n <- 10
a <- 2
b <- 2
num_iter <- 10000

set.seed(0)
result_rcpp <- gibbs_sampler_rcpp(n, a, b, num_iter)
result_r <- gibbs_sampler_r(n, a, b, num_iter)

# QQ Plot comparison for x and y
par(mfrow = c(1, 2))
qqplot(result_rcpp$x, result_r$x, main = "QQ Plot for x", xlab = "Rcpp x", ylab = "R x")
abline(0, 1, col = "red")

qqplot(result_rcpp$y, result_r$y, main = "QQ Plot for y", xlab = "Rcpp y", ylab = "R y")
abline(0, 1, col = "red")


## -----------------------------------------------------------------------------
library(microbenchmark)

ts <- microbenchmark(
  Rcpp = gibbs_sampler_rcpp(n, a, b, num_iter),
  R = gibbs_sampler_r(n, a, b, num_iter),
  times = 10
)

summary(ts)[,c(1,3,5,6)]


