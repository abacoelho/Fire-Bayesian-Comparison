library(readxl)
library(ggplot2)
library(reshape2)
library(BEST)
library(rjags)
library(latex2exp)
setwd("~/MS Stats/Research/Wildfire/data")

# Read in data and subset the training data
data <- read_excel("fire modeling data.xlsx")
data = data[complete.cases(data), ]
data = data[data$DISCOVER_YEAR > 2009,]
data$LGHTNG = as.factor(data$LGHTNG)

# Creating data visualizations
ggplot() + geom_density(data = data, aes(x=TMAX, group=LGHTNG, fill=LGHTNG),alpha=0.3, adjust=2) + 
  xlab("Max Temperature") +
  ylab("Density")

ggplot() + geom_density(data = data, aes(x=AWND, group=LGHTNG, fill=LGHTNG),alpha=0.3, adjust=2) + 
  xlab("Average Windspeed") +
  ylab("Density")

ggplot() + geom_density(data = data, aes(x=HWY_DSTNC, group=LGHTNG, fill=LGHTNG),alpha=0.3, adjust=2) + 
  xlab("Miles") +
  ylab("Density") +
  ggtitle('Distance to Highway Density Plot of True Values') +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))


# Estimating the differences in distributions using BEST
# Lightning vs Human-cased ignited fires
preform_BEST <- function(var_name){
  pooled = as.vector(data[[var_name]])
  Y1 = as.vector(data[data$LGHTNG == 1,][[var_name]])
  Y2 = as.vector(data[data$LGHTNG == 0,][[var_name]])
  
  best_out = BESTmcmc(Y1, Y2, numSavedSteps = 1000,
                       thinSteps = 10, burnInSteps = 100)
  
  plotAll(best_out, credMass=0.8, ROPEm=c(-0.1,0.1),
          ROPEeff=c(-0.2,0.2), compValm=0.5)
  best_out
}

best.tmax = preform_BEST('TMAX')
best.awnd = preform_BEST('AWND')
#best.hwy  = preform_BEST('HWY_DSTNC') <- not a t-distribution


# Modifying BEST for log-normal likelihood distributions

Y1 = as.vector(data[data$LGHTNG == 1,][['HWY_DSTNC']])
Y2 = as.vector(data[data$LGHTNG == 0,][['HWY_DSTNC']])
pooled = c(Y1, Y2)
x <- c(rep(1, length(Y1)), rep(2, length(Y2))) # indicator variable

model.str <- 'model {
  for (i in 1:Ntotal) {
      # log-normal likelihood:
      y[i] ~ dlnorm(mu[x[i]], tau[x[i]])
  }
  for ( j in 1:2 ) {
      # uninformed priors
      mu[j] ~ dnorm( muM, muP )
      tau[j] <- 1/pow( sigma[j], 2 )
      sigma[j] ~ dunif( sigmaLow, sigmaHigh )
  }
}'

cpd.model <- jags.model(textConnection(model.str),
                        data=list(y=pooled, x=x,
                                  muM=mean(pooled),
                                  muP=1/(1000 * sd(pooled))^2,
                                  sigmaLow=sd(pooled) / 1000,
                                  sigmaHigh=sd(pooled) * 1000,
                                  Ntotal=length(pooled)))

update(cpd.model, 1000)
chain <- coda.samples(model = cpd.model, n.iter = 5000,
                      variable.names = c('mu', 'sigma'))
rchain <- as.matrix(chain)

# Plot creditable values from the marginal posterior for each parameter
par(mfrow=c(3,2))
hist(rchain[, 'mu[1]'], freq=FALSE,
     main = 'Lightning Location', yaxt='n',
     ylab = '', xlab = TeX('$\\mu_1$'))
hist(rchain[, 'mu[2]'], freq=FALSE,
     main = 'Human Location', yaxt='n',
     ylab = '', xlab = TeX('$\\mu_2$'))

hist(rchain[, 'sigma[1]'], freq=FALSE,
     main = 'Lightning Scale', yaxt='n',
     ylab = '', xlab = TeX('$\\sigma_1$'))
hist(rchain[, 'sigma[2]'], freq=FALSE,
     main = 'Human Scale', yaxt='n',
     ylab = '', xlab = TeX('$\\sigma_2$'))

hist(rchain[, 'mu[1]'] - rchain[, 'mu[2]'], freq=FALSE,
     main = 'Difference of Locations', yaxt='n',
     ylab = '', xlab = TeX('$\\mu_1 - \\mu_2$'))
hist(rchain[, 'sigma[1]'] - rchain[, 'sigma[2]'], freq=FALSE,
     main = 'Difference of Scales', yaxt='n',
     ylab = '', xlab = TeX('$\\sigma_1 - \\sigma_2$'))


# Plot the predicted posterior using mean values of parameters
x<-seq(from=0,to=+40,length.out=1000)
plot(x, dlnorm(x, meanlog = mean(rchain[, 'mu[2]']),
               sdlog = mean(rchain[, 'sigma[2]'])),
     type = 'l', xlab = 'Y', ylab = 'P(Y)', 
     main = "Predicted PDF of Distance from Highway", col = 'red')
lines(x, dlnorm(x, meanlog = mean(rchain[, 'mu[1]']), 
                sdlog = mean(rchain[, 'sigma[1]'])), 
      col = 'blue') 
legend("topright",c("LGHTNG=0","LGHTNG=1"),col=c("red", "blue"),lty=1)   
