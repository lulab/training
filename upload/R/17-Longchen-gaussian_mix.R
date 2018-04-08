a <- rbinom(1000, 1, 0.7)
z <- a*rnorm(1000,10,4) + (1 - a)*rnorm(1000,30,5)
hist(z,100)


require(mixtools)
gaussian_mix <- normalmixEM(z, k=2)

fit <- gev(z)
test <- ks.test("pgev", fit$par.ests[1], fit$par.ests[3], fit$par.ests[2])
param <- data.frame(c('xi', 'mu', 'sigma', 'pvalue'), c(fit$par.ests[1], fit$par.ests[3], 
                                                        fit$par.ests[2], test$p.value))


x <- seq(min(z)-1,max(z),0.5)
y <- gaussian_mix$lambda[1]*dnorm(x, gaussian_mix$mu[1], gaussian_mix$sigma[1]) + 
  gaussian_mix$lambda[2]*dnorm(x, gaussian_mix$mu[2], gaussian_mix$sigma[2])

par(new=TRUE)
plot(x,y,'l',col='blue', xlab = '', ylab = '', axes = FALSE)