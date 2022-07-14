#HMSC Single Species Distribution Modelling Example

####Generate simulated data####

#generate data where we simulate variation in a single covariate x on n=50
#sampling units. you then construct the linear predictor assuming that beta1
#is teh intercept and beta 2 the slope. the linear predictor is then used to 
#construct 3 response variables which conform the the assumptions of 
#y1 normal model, y2 probit model and y3 lognormal poisson model

set.seed(1)
n = 50
x = rnorm(n)
beta1 = 0
beta2 = 1
L = beta1 + beta2*x
y1 = L + rnorm(n, sd = 1)
y2 = 1* ((L + rnorm(n, sd = 1)) > 0)
y3 = rpois (n = n, lambda = exp(L + rnorm(n, sd = 1)))
par(mfrow = c(1,3))
plot(x,y1, main = "Normal")
plot(x,y2, main = "probit")
plot(x,y3, main = "Lognormal Poisson")

####Fitting Model and Examining Parameter Estimates####

#analyse the normally distributed data

df = data.frame(x, y1)
m.lm = lm(y1 ~ x, data = df)
summary(m.lm)

#conduct the analogous analyses with Hmsc

library(Hmsc)

Y = as.matrix(y1)
XData = data.frame(x = x)
m.normal = Hmsc(Y=Y, XData = XData, XFormula = ~x)

#fitting an HMSC model with Bayesian inference

nChains = 2
thin = 5
samples = 1000
transient = 500*thin
verbose = 500*thin

m.normal = sampleMcmc(m.normal, thin=thin, samples=samples,
                      transient=transient, nChains=nChains,
                      verbose=verbose)

#extract the posterior distribution from the model object and convert it to a
#format understood by the coda package

mpost = convertToCodaObject(m.normal)

#use the summary function to inspect the estimates of intercept and slope
#which are referred to as Beta in HMSC

summary(mpost$Beta)

#Rsquared can then be assessed with the function evaluateModelFit once the 
#posterior distribution of the predicted values is computed using 
#computePredictedValues

preds = computePredictedValues(m.normal)
MF = evaluateModelFit(hM = m.normal, predY = preds)
MF$R2

####Checking MCMC Convergence Diagnostics####

#showing trace plots of the bete-parameters

plot(mpost$Beta)

effectiveSize(mpost$Beta)

gelman.diag(mpost$Beta,multivariate = FALSE)$psrf

####Check the Assumptions of the Linear Model####

#diagnostic plot to check assumptions of linear model

nres.lm = rstandard(m.lm)
preds.lm = fitted.values(m.lm)
par(mfrow = c(1,2))
hist(nres.lm, las = 1)
plot(preds.lm, nres.lm, las = 1)
abline(a = 0,b = 0)

#to generate the diagnostic plots for the linear model we just fitted with HMSC,
#we first summarise the posterior distribution of predicted values into the
#posterior mean and then extract and standardize the residuals

preds.mean = apply(preds, FUN = mean, MARGIN = 1)
nres = scale(y1-preds.mean)
par(mfrow = c(1,2))
hist(nres)
plot(preds.mean, nres)
abline(a = 0, b = 0)

####Fitting Generalised Linear Models####

#can fit the model with arguments under distr depending on the model chosen

m.normal = Hmsc(Y=Y, XData = XData, XFormula = ~x,
                distr = "normal")

#fit the probit model for presence-absence data

Y = as.matrix(y2)
m.probit = Hmsc(Y = Y, XData = XData, XFormula = ~x,
                distr = "probit")
#distr = "probit" or "poisson" or "lognormal poisson"

#obtain posterior samples

verbose = 0

m.probit = sampleMcmc(m.normal, thin=thin, samples=samples,
                      transient=transient, nChains=nChains,
                      verbose=verbose)

#evaluate MCMC convergence

mpost = convertToCodaObject(m.probit)
effectiveSize(mpost$Beta)

gelman.diag(mpost$Beta,multivariate = FALSE)$psrf

#look at the parameter estimates and evaluate the models explanatory power

round(summary(mpost$Beta)$quantiles, 2)

preds = computePredictedValues(m.probit)
evaluateModelFit(hM = m.probit, predY = preds)

#fitting a lognormal Poisson model to response variable y3 and evaluating
#MCMC convergence

Y = as.matrix(y3)
m.lognormal.poisson=Hmsc(Y = Y, XData = XData, XFormula = ~x,
                         distr = "lognormal poisson")
m.lognormal.poisson = sampleMcmc(m.lognormal.poisson,
                                 thin = thin, samples = samples,
                                 transient = transient, nChains = nChains,
                                 verbose = verbose)
mpost = convertToCodaObject(m.lognormal.poisson)
effectiveSize(mpost$Beta)

gelman.diag(mpost$Beta,multivariate = FALSE)$psrf

#look at parameter estimates and evaluate model fit

round(summary(mpost$Beta)$quantiles,2)

preds = computePredictedValues(m.lognormal.poisson, expected = FALSE)
evaluateModelFit(hM = m.lognormal.poisson, predY = preds)

####Predicting New Sampling Units####

par(mfrow = c(1,3))
for (i in 1:3){
  m = switch (i, m.normal, m.probit, m.lognormal.poisson)
  Gradient = constructGradient(m, focalVariable = "x")
  predY = predict(m, Gradient = Gradient, expected = TRUE)
  plotGradient(m, Gradient, pred = predY, measure = "Y",
               index = 1, showData = TRUE, main = c("Normal",
                                                    "Probit","Lognormal Poisson")[ i])
}
