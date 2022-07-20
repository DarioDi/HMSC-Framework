#HMSC Joint Species Distribution Model - Simulated Case Study

library(Hmsc)
library(ape)
library(MASS)
library(beepr)

#### Simulating Species Niches ####

#generate a niche variation. construct a random phylogeny for 100 species
#vcv function in the ape package turns the phylogenetic tree into a phylogenetic
#correlation matrix C

ns = 100
phy = rcoal(n = ns, tip.label =
              sprintf('sp_%.3d',1:ns), br = "coalescent")
C = vcv(phy, model = "Brownian", corr = TRUE)

#Community A trait values are sampled independently for each species from the
# standard normal distribution. 
# Community B trait values are sampled from the multivariate normal distribution
# for which the variance-covariance matrix equals the phylogenetic correlation
# matrix C

Tr.A = cbind(rep(1,ns), rnorm(ns))
Tr.B = cbind(rep(1,ns), mvrnorm(n = 1, mu = rep(0, ns), Sigma = C))

# define the T matrix that describes the link between species traits and niches
# and use the matrix product %*% to compute the expected values of species niches

gamma = cbind(c(-2,2), c(-1,1))
mu.A = gamma %*% t(Tr.A)
mu.B = gamma %*% t(Tr.B)

# generate species niches assuming the Community A residual variation is 
# phylogenetically fully structured and that in Community B the residual 
# variation is fully independent among species

V2 = diag(2)
beta.A = matrix(mvrnorm(n=1,mu=as.vector(mu.A),
                        Sigma = kronecker(C, V2)), ncol = ns)
beta.B = matrix(mvrnorm(n=1,mu=as.vector(mu.B),
                        Sigma = kronecker(diag(ns), V2)), ncol = ns)

#### Simulating Species Data ####

# we consider a single environmental covariate x and use the standard normal
# distribution to simulation variation in x over n = 50 sampling units

n = 50
X = cbind(rep(1, n), rnorm(n))

# the intercept was included in X matrix so we can compute the linear predictors

L.A = X %*% beta.A
L.B = X %*% beta.B

# now convert the linear predictors into community data matrix Y using occurrence
# data (probit model)

Y.A = 1*((L.A+matrix(rnorm(n*ns), ncol = ns)) > 0)
Y.B = 1*((L.B+matrix(rnorm(n*ns), ncol = ns)) > 0)

#### Exploring Raw Data ####

# before fitting to HMSC we explore the raw data. since we are using species
# occurrences we may wish to look at species richness (row sums) and species 
# prevalences (row means)

S.A = rowSums(Y.A)
P.A = colMeans(Y.A)
S.B = rowSums(Y.B)
P.B = colMeans(Y.B)

par(mfrow = c(2,2))
hist(S.A)
hist(P.A)
hist(S.B)
hist(P.B)

#### Fitting ab HNSC Model for the Community A with Phylogenetically Structured Species Niches ####

#start by formatting the data so they are suitable for input into HMSC

community = "A"
Y = switch(community, "A" = Y.A, "B" = Y.B)
colnames(Y) = phy$tip.label
Tr = switch(community, "A" = Tr.A, "B" = Tr.B)
TrData = data.frame(trait = Tr[,2])
rownames(TrData) = phy$tip.label
XData = data.frame(x = X[,2])

# had to put row names in TrData to match Y

# define the HMSC model. model species occurrences as a linear function to the 
# environmental variable x and species niches as a linear function of the
# trait covarite. since data is species occurences, we use the probit model

m = Hmsc(Y=Y, XData = XData, XFormula = ~x, TrData = TrData,
         TrFormula = ~trait, phyloTree = phy, distr = "probit")

# perform the model fitting

nChains = 2 #how many chains to sample
thin = 5 #how much thinning to apply
samples = 1000 #how many samples to obtain per chain
transient = 500*thin #length of transient to include
verbose = 500*thin #how frequently we see the progress of the MCMC sampling

m = sampleMcmc(m, thin = thin, samples = samples,
               transient = transient, nChains = nChains, verbose = verbose)
beep(sound = 1, expr = NULL)

# check MCMC convergence

mpost = convertToCodaObject(m)
# plot(mpost$Rho)

effectiveSize(mpost$Rho)


gelman.diag(mpost$Rho, multivariate = FALSE,
            autoburnin = FALSE)$psrf



#### Explanatory and Predictive Powers of the HMSC Model ####

preds = computePredictedValues(m)
MF = evaluateModelFit(hM = m, predY = preds)

partition = createPartition(m, nfolds = 2)
preds = computePredictedValues(m, partition = partition)
MFCV = evaluateModelFit(hM = m, predY = preds)
beep(sound = 1, expr = NULL)

par(mfrow = c(1,2))
plot(MF$TjurR2 ~ MF$AUC)
plot(MFCV$TjurR2 ~ MFCV$AUC)

dev.off()

#### Examining Parameter Estimates ####

# we use plotBeta to visualise the estimated species niches

postBeta = getPostEstimate(m, parName = "Beta")
plotBeta(m, post = postBeta, param = "Sign", plotTree = TRUE,
         supportLevel = 0.95, split = 0.4, spNamesNumbers = c(F, F))

# when the species are estimated to be positive, it means they are responding
# to covariate x

# we use plotGamma to visualise how species niches are estimed to depend on
# species traits

postGamma = getPostEstimate(m, parName = "Gamma")
plotGamma(m, post = postGamma, param = "Sign", supportLevel = 0.95)

# as the covariate is positive, this indicates that species with a high trait 
# value respond positively to the environmental covariate.
# if the parameter estimate is negative, this indicates that species with a 
# high trait value have a low baseline occurrence probability 

# we next look at the parameter estimate for phylogenetic signal parameter p

summary(mpost$Rho)$quantiles

# the posterior distribution reveals strong evidence for a phylogenetic signal
# as the 95% credible interval of p is positive and very close to 1

#### Repeating the Analyses for the Community B Where Species Niches Are Structured by Their Traits ####

# repeat the analyses now for community B

community = "B"

m = sampleMcmc(m, thin = thin, samples = samples,
               transient = transient, nChains = nChains, verbose = verbose)
beep(sound = 1, expr = NULL)

mpost = convertToCodaObject(m)

postBeta = getPostEstimate(m, parName = "Beta")
plotBeta(m, post = postBeta, param = "Sign", plotTree = TRUE,
         supportLevel = 0.95, split = 0.4, spNamesNumbers = c(F, F))

postGamma = getPostEstimate(m, parName = "Gamma")
plotGamma(m, post = postGamma, param = "Sign", supportLevel = 0.95)

summary(mpost$Rho)$quantiles











