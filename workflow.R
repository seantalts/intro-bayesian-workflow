library(rstan)
library(ggplot2)
util = new.env()
source("stan_utility.R", local = util)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#############
## Scope out your problem
#############
df = read.csv("data.csv")
head(df)

#############
## Model!
#############
## 1. Create a new Stan file with your proposed model. We'll start with pooled.stan

#############
## Check the model with fake data! (skip if short on time)
## 1. Create a new Stan program (or use R, or some other language) and
##    generate fake data according to the data-generating process described by
##    your model.
#############

#############
## Fit the model
#############
model = stan_model("../pooled.stan")
fit = sampling(model, list(log_radon=df$log_radon,
                           basement = df$basement,
                           N=nrow(df)))
#############
## Check diagnostics
#############
util$check_all_diagnostics(fit)

#############
## Check & graph fit
#############
print(fit)
sample = extract(fit)

# This line just plots our data by basement indicator
ggplot(df, aes(basement, log_radon)) + geom_count() +
  # This next snippet will select every 10th posterior sample and plot the 
  # fit lines corresponding to that sample. It's a good way to get a sense of
  # the uncertainty / distribution of the posterior
  sapply(seq(1, nrow(sample$alpha), 10), function(i) {
    geom_abline(intercept = sample$alpha[i],
                slope = sample$basement_effect[i], alpha=0.03)
  })
#############
## Check PPCs
#############
# Plot a single house's posterior predictive distribution (house 919 here)
ppcs = as.data.frame(sample$yppc)
ggplot(ppcs) + geom_density(aes(V919)) + geom_vline(xintercept=df$log_radon[919])

# Plot an entire replication's distribution against the actual data's distribution
rep = 2000
ppcs = data.frame(log_radon = df$log_radon,
                  model = sample$yppc[rep,])
library(reshape2) # I wish I didn't have this dependency here
ppcs = reshape2::melt(ppcs)
ggplot(ppcs, aes(x=value, linetype=variable)) + geom_density(alpha=0.2)

# Shinystan!
#launch_shinystan(fit)