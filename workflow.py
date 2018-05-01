import pystan
import pandas as pd
from ggplot import *
import stan_utility as util

#############
## Scope out your problem
#############
df = pd.read_csv("data.csv")
print(df.head())

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
model = util.compile_model("../pooled.stan")
fit = model.sampling(dict(log_radon=df.log_radon,
                          basement = df.basement,
                          N=df.shape[0]))
# #############
# ## Check diagnostics
# #############
util.check_all_diagnostics(fit)

# #############
# ## Check & graph fit
# #############
print(fit)
sample = fit.extract()

# # This line just plots our data by basement indicator
g = ggplot(df, aes(x='basement', y='log_radon')) + geom_point()

# This next snippet will select every 10th posterior sample and plot the
# fit lines corresponding to that sample. It's a good way to get a sense of
# the uncertainty / distribution of the posterior
for i in range(0, len(sample['alpha']), 100):
    g += geom_abline(intercept = sample['alpha'][i],
              slope = sample['basement_effect'][i], alpha=0.3)

g += xlim(-0.1, 1.1)
g += ylim(-2.5, 4.1)
g.show()
#
# ## This snippet makes the unpooled and hierarchical fit checking graph
# ## This will fail if you try to use it on the pooled model with a single alpha,
# ## so I'm starting with it commented out.
# # estimates = apply(sample$a, 2, mean)

# estimates = od["a"].mean(0)
# sd = od["a"].std(0)
# fitdf = pd.DataFrame(dict(a = estimates, county=range(0, 85),
#                          upper = estimates + 2*sd,
#                          lower = estimates - 2*sd,
#                          sd = sd))
# fitdf = fitdf.sort_values(by='a')
# fitdf["index"] = np.arange(0, 85)
#
# (ggplot(fitdf, aes(x='index'))
#  + geom_point(aes(y='a'))
#  + geom_ribbon(aes(ymin='lower', ymax='upper'), alpha=0.2)
#  + ylim(-1.1, 3.5))


#############
## Check PPCs
#############
# Plot a single house's posterior predictive distribution (house 919 here)
#ppcs = pd.DataFrame(sample['yppc'])
#ppcs.columns = map(str, ppcs.columns)
#(ggplot(ppcs, aes('918'))
# + geom_density()
# + geom_vline(x=df.log_radon[918]))
#
## Plot an entire replication's distribution against the actual data's distribution
#rep = 2000
#ppcs = pd.DataFrame(dict(measured = df.log_radon,
#                         predicted = sample['yppc'][rep,:]))
#ppcs = pd.melt(ppcs)
#ggplot(ppcs, aes(x='value', linetype='variable')) + geom_density()
