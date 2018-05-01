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
## 2. Follow the rest of this workflow with your fake data!
#############
# To use Stan to generate fake data, we use the `generated quantities` block
# and the 'Fixed_param' algorithm. The following snippet will generate 
# 200 fake data points.
# gen.model = stan_model("gen.stan")
# gen.data = sampling(data=d, algorithm='Fixed_param', warmup=0, iter=200)
# gendf = as.data.frame(gen.data)
# dataset.1 = gendf[1,]

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

## This snippet makes the unpooled and hierarchical fit checking graph
## This will fail if you try to use it on the pooled model with a single alpha,
## so I'm starting with it commented out.
# estimates = apply(sample$a, 2, mean)
# sd = apply(sample$a, 2, sd)
# fit.df = data.frame(estimates, 1:85, 
#                          estimates + 2*sd,
#                          estimates - 2*sd,
#                          sd)
# names(fit.df) = c("a", "county", "upper", "lower", "sd")
# fit.df = fit.df[order(fit.df$a),]
# 
# ggplot(fit.df, aes(x=1:85, y=a)) + geom_pointrange(aes(ymin = lower, ymax = upper)) +
#   ylim(-1.1, 3.5)


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
library(reshape2)
ppcs = reshape2::melt(ppcs)
ggplot(ppcs, aes(x=value, linetype=variable)) + geom_density(alpha=0.2)

# Shinystan!
#launch_shinystan(fit)