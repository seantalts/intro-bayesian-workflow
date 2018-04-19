library(rstan)
library(ggplot2)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

df = read.csv("data.csv")
head(df)

## Pooled
pooled.model = stan_model("pooled.stan")
pooled.fit = sampling(pooled.model, list(y=df$log.radon,
                           x = df$floor.measure,
                           N=nrow(df)))
print(pooled.fit)

pooled.sample = extract(pooled.fit)
b = mean(pooled.sample$beta[,1])
m = mean(pooled.sample$beta[,2])

ggplot(df, aes(floor.measure, log.radon)) + geom_count() +
  geom_abline(intercept = b, slope = m, linetype=2)



## Unpooled
unpooled.model = stan_model("unpooled.stan")
unpooled.fit = sampling(unpooled.model, list(y=df$log.radon,
                                             x = df$floor.measure,
                                             county = df$county,
                                             N=nrow(df)))
print(unpooled.fit, probs = c(0.025, 0.5, 0.975))

unpooled.samples = extract(unpooled.fit)
unpooled.estimates = apply(unpooled.samples$a, 2, mean)
unpooled.se = apply(unpooled.samples$a, 2, sd)
unpooled.df = data.frame(unpooled.estimates, 1:85, 
                         unpooled.estimates + unpooled.se,
                         unpooled.estimates - unpooled.se)
names(unpooled.df) = c("a", "county", "upper", "lower")
unpooled.df = unpooled.df[order(unpooled.df$a),]
ggplot(unpooled.df, aes(x=1:85, y=a)) + geom_pointrange(aes(ymin = lower, ymax = upper))

## Partial pooling for alpha
partial.model = stan_model("partial.stan")
partial.fit = sampling(partial.model, list(y=df$log.radon,
                                           x = df$floor.measure,
                                           J=length(unique(df$county)),
                                           county = df$county,
                                           uranium = df$uranium,
                                           N=nrow(df)))
samples = extract(partial.fit)
ggplot(df, aes(floor.measure, log.radon)) + geom_count() +
  sapply(1:85, function(i){
    geom_abline(intercept = mean(samples$a[,i]),
                slope = mean(samples$beta),
                alpha=0.3)
  })

## Partial pooling for alpha and beta
partial.model = stan_model("partial.stan")
partial.fit = sampling(partial.model, list(y=df$log.radon,
                                           x = df$floor.measure,
                                           J=length(unique(df$county)),
                                           county = df$county,
                                           uranium = df$uranium,
                                           N=nrow(df)))
samples = extract(partial.fit)
ggplot(df, aes(floor.measure, log.radon)) + geom_count() +
  sapply(1:85, function(i){
    geom_abline(intercept = mean(samples$a[,i]),
                slope = mean(samples$beta[,i]),
                alpha=0.3)
  })
#launch_shinystan(partial.fit)