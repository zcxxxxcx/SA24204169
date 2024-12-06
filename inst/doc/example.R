## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
# Load the data from the .rda file
load(system.file("data", "data.rda", package = "SA24204169"))
d <- dim(data)[2] #dim of the data
w0 <- rep(1/d,d) # Initial weight
n <- 120 # Window
delta <- 1e-3 # Robustness parameter
alpha <- -0.1 # Adjusted expected return

result <- SA24204169::optimize_weights(data, w0, n, delta, alpha)
result
hist(result$returns, main = "Distribution of earning rates", xlab = "Returns")

## -----------------------------------------------------------------------------
set.seed(0)
r1 <- result$returns # returns we get above
r2 <- rnorm(length(r1)) # simulated data for convenience
SA24204169::pvalue(r1,r2)


