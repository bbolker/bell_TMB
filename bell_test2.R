library(Rcpp)
library(lamW)      ## for lambert W
library(VGAMdata)  ## for dbell()
library(multicool) ## for Bell()
library(glmmTMB)
library(bellreg)

z <- 1.0
all.equal(lambertW(z), lamW::lambertW0(1.0)) ## TRUE



# Load the functions from the C++ file
Rcpp::sourceCpp("bell_functions.cpp")

## copied from the web (https://www.statisticshowto.com/bells-numbers-bell-triangle/)
bell_nums <- c(1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975)
## Test the bell function against these and existing implementation
as.integer(sapply(1:8, Bell))  ## uh-oh
multicool::Bell(1:8)

maxint <- 2147483647
for (i in 1:20) {
    if (Bell(i) > maxint) break
}
print(i)

## Test the dbell() function against the VGAMdata version
d1 <- VGAMdata::dbell(0:15, shape = 1.5)
d2 <- dbell(0:15, shape = 1.5)
all.equal(d1[1:3], d2[1:3])
## things start to go wrong at x = 4
d1[4]
d2[4]

rbell_samples <- rbell(100, 1.5)  # Generate 100 random samples from Bell distribution with shape = 1.5

quantiles <- qbell(c(0.1, 0.5, 0.9), 1.5)  # Compute quantiles for p = 0.1, 0.5, 0.9 with theta = 1.5
print(quantiles)
##

data(faults)

fit <- bellreg(nf ~ lroll, data = faults)

#  Prepare the data for C++ function
object <- list(
  y = faults$nf,    # Dependent variable
  mu = fitted(fit)       # Fitted values from the model
)

## BMB: whatever this is is slow ...
# Step 4: Compute quantile residuals using the custom C++ function
if (FALSE) {
    system.time(cpp_resid <- qresiduals(object))
}

plot(cpp_resid, main = "C++ Implementation Residuals", ylab = "Residuals", xlab = "Index")
cpp_vmu(fitted(fit))

x <- c(1, 2, 3, 4, 5)  # Example data
theta <- 1.5  # Example theta

# Calculate the log-likelihood
log_likelihood <- log_likelihood_bell(x, theta)
print(log_likelihood)

devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")

remotes::install_github("glmmTMB/glmmTMB/glmmTMB", ref = "bell")
library(glmmTMB)
dd <- data.frame(x = 1:5)
fit1 <- glmmTMB(x ~ 1, data = dd, family = list(family = "bell", link = "log"))
## mean value
c0 <- fixef(fit1)$cond[1]
exp(c0)*exp(exp(c0))  ## 3

set.seed(101)
dd <- data.frame(x = bellreg::rbell(100, theta = 2))
glmmTMB(x ~ 1, data = dd, family = list(family = "bell", link = "log"), verbose = TRUE)

Rcpp::sourceCpp("bell_functions.cpp")
b1 <- sapply(1:20, numbers::bell)
b2 <- multicool::Bell(1:20)
b3 <- Bell(1:20)


## need microbenchmark instead of rbenchmark because our second version is too fast ...
testfun <- function(n) {
    ## rbenchmark::benchmark(
    microbenchmark::microbenchmark(
        numbers = numbers::bell(n),
        multicool = multicool::Bell(n),
        ours = Bell(n),
        ours2 = Bell2(n),
        ## replications = 50
        times = 50
        ) |> summary()
}
single_fun <- function(n) {
    data.frame(numbers = numbers::bell(n),
               multicool = multicool::Bell(n),
               ours2 = Bell2(n)
               )
}

ss <- function(n) setNames(n, n)
bb <- purrr::map_dfr(ss(10:25), testfun, .id = "n")
vv <- purrr::map_dfr(ss(10:30), single_fun, .id = "n")

save("bb", "vv", file = "bell_bench.rda")           
library(ggplot2); theme_set(theme_bw())
ggplot(bb, aes(as.numeric(n),
               ##elapsed,
               median,
               colour = expr
               ## test
               )) + geom_point() + geom_smooth() +
    scale_x_log10() + scale_y_log10()


