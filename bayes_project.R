# Import necessary libraries
library(rjags)
library(coda)
library(ggplot2)
library(dplyr)

# Import data
my.data <- datasets::Puromycin
N <- nrow(my.data)
my.data.jags <- list(
  Y = my.data$rate,
  x = my.data$conc,
  state = as.integer(my.data$state),
  N = N
)


# QUESTION 1
cat("model {

  # LIKELIHOOD
  for (i in 1:N) {
      Y[i] ~ dnorm(mu[i], tau)
      mu[i] <- (Vmax[state[i]] * x[i]) / (K[state[i]] + x[i])
  }
  
  # PRIORS
  for (j in 1:2) {
    Vmax[j] ~ dunif(0, 1000)
  }
  
  for (j in 1:2) {
    K[j] ~ dunif(0, 10)
  }
  
  tau ~ dgamma(0.001, 0.001)
  sigma2 <- 1 / tau
  sigma <- 1 / sqrt(tau)
}",file="enzyme_model.txt")


# QUESTION 2
my.inits <- function()list(
  Vmax = c(200, 100),
  K = c(0.1, 0.1),
  tau = 1
)

params <- c("Vmax", "K", "sigma")

jags <- jags.model(file="enzyme_model.txt",
                   data=my.data.jags,
                   inits=my.inits,
                   n.chains=3)

update(jags,5000)

enzyme.sim <- coda.samples(
  model = jags,
  variable.names = params,
  n.iter = 10000,
  thin=1)

# check trace plots
traceplot(enzyme.sim)

# convergence test
gelman.diag(enzyme.sim)


# QUESTION 3
# Full posterior summary
summary(enzyme.sim)

# Clean combined table
sim.summary <- summary(enzyme.sim)
cbind(sim.summary$statistics[, c("Mean","SD")],
      sim.summary$quantiles[, c("2.5%","50%","97.5%")])

# Posterior density plots
par(mar = c(3,3,2,1))
plot(enzyme.sim)


# QUESTION 4
sim.matrix <- as.matrix(enzyme.sim)

# Compute ratios sample by sample
Vmax.ratio <- sim.matrix[, "Vmax[1]"] / sim.matrix[, "Vmax[2]"]
K.ratio <- sim.matrix[, "K[1]"] / sim.matrix[, "K[2]"]

# Posterior summaries
# Vmax summary
cat("Vmax Ratio (Treated / Untreated)\n")
cat("Mean:", mean(Vmax.ratio),"\n")
cat("SD:", sd(Vmax.ratio), "\n")
cat("Median:", median(Vmax.ratio), "\n")
cat("95% CI:", quantile(Vmax.ratio, c(0.025, 0.975)), "\n\n")

# K summary
cat(" K Ratio (Treated / Untreated)\n")
cat("Mean:", mean(K.ratio), "\n")
cat("SD:", sd(K.ratio), "\n")
cat("Median:", median(K.ratio), "\n")
cat("95% CI:", quantile(K.ratio, c(0.025, 0.975)), "\n")

# Density plots
par(mfrow=c(1,2), mar=c(4,4,3,1))

hist(Vmax.ratio, breaks=50, probability=TRUE,
     main="Posterior: Vmax Treated/Untreated",
     xlab="Ratio", col="steelblue", border="white")
abline(v=1, col="red",   lwd=2, lty=2)
abline(v=quantile(Vmax.ratio, c(0.025,0.975)), col="black", lwd=1.5, lty=3)

hist(K.ratio, breaks=50, probability=TRUE,
     main="Posterior: K Treated/Untreated",
     xlab="Ratio", col="darkorange", border="white")
abline(v=1, col="red",   lwd=2, lty=2)
abline(v=quantile(K.ratio, c(0.025,0.975)), col="black", lwd=1.5, lty=3)


# QUESTION 5
# Concentration grid
conc.seq  <- seq(0.1, 1, length.out = 200)
n.samples <- nrow(sim.matrix)

# Storage matrices - one row per posterior sample, one column per concentration
curve.untreated <- matrix(NA, nrow=n.samples, ncol=length(conc.seq))
curve.treated <- matrix(NA, nrow=n.samples, ncol=length(conc.seq))

# Plug each posterior sample into the Michaelis-Menten equation
for (s in 1:n.samples) {
  Vmax1 <- sim.matrix[s, "Vmax[1]"]
  Vmax2 <- sim.matrix[s, "Vmax[2]"]
  K1 <- sim.matrix[s, "K[1]"]
  K2 <- sim.matrix[s, "K[2]"]
  
  curve.treated[s, ] <- (Vmax1 * conc.seq) / (K1 + conc.seq)
  curve.untreated[s, ] <- (Vmax2 * conc.seq) / (K2 + conc.seq)
}

# Posterior mean and 95% credible bands at each concentration
mean.treated <- apply(curve.treated,   2, mean)
lo.treated <- apply(curve.treated,   2, quantile, 0.025)
hi.treated <- apply(curve.treated,   2, quantile, 0.975)

mean.untreated <- apply(curve.untreated, 2, mean)
lo.untreated <- apply(curve.untreated, 2, quantile, 0.025)
hi.untreated <- apply(curve.untreated, 2, quantile, 0.975)

# Plot
par(mar=c(4,4,3,2))

plot(conc.seq, mean.treated, type="n",
     ylim=c(80, 250), xlim=c(0.1, 1),
     xlab="Concentration", ylab="Reaction Rate",
     main="Michaelis-Menten Fitted Curves with 95% Credible Bands")

# Uncertainty bands
polygon(c(conc.seq, rev(conc.seq)),
        c(hi.treated, rev(lo.treated)),
        col=adjustcolor("steelblue", alpha=0.2), border=NA)

polygon(c(conc.seq, rev(conc.seq)),
        c(hi.untreated, rev(lo.untreated)),
        col=adjustcolor("darkorange", alpha=0.2), border=NA)

# Posterior mean curves
lines(conc.seq, mean.treated, col="steelblue", lwd=2.5)
lines(conc.seq, mean.untreated, col="darkorange", lwd=2.5)

# Observed data points
points(my.data$x[my.data$state==1], my.data$Y[my.data$state==1],
       col="steelblue", pch=16, cex=1.2)
points(my.data$x[my.data$state==2], my.data$Y[my.data$state==2],
       col="darkorange", pch=16, cex=1.2)

legend("bottomright",
       legend = c("Treated", "Untreated"),
       col = c("steelblue", "darkorange"),
       lwd=2, pch=16, bty="n")


# QUESTION 6
x.star <- 0.5

# Plug x = 0.5 into Michaelis-Menten for every posterior sample
rate.treated <- (sim.matrix[, "Vmax[1]"] * x.star) /
  (sim.matrix[, "K[1]"] + x.star)

rate.untreated <- (sim.matrix[, "Vmax[2]"] * x.star) /
  (sim.matrix[, "K[2]"] + x.star)

# Posterior summaries
cat(" Treated: Rate at x = 0.5 \n")
cat("Mean:", mean(rate.treated),"\n")
cat("SD:", sd(rate.treated),"\n")
cat("95% CI:", quantile(rate.treated, c(0.025, 0.975)),"\n\n")

cat("Untreated: Rate at x = 0.5 \n")
cat("Mean:", mean(rate.untreated),"\n")
cat("SD:", sd(rate.untreated),"\n")
cat("95% CI:", quantile(rate.untreated, c(0.025, 0.975)),"\n")

# Combined density plot
par(mar=c(4,4,3,1))

# Set up plot window using treated range
plot(density(rate.treated),
     main = "Posterior Reaction Rate at x = 0.5",
     xlab = "Reaction Rate",
     col = "steelblue",
     lwd = 2.5,
     xlim = range(c(rate.treated, rate.untreated)))

# Overlay untreated
lines(density(rate.untreated),
      col = "darkorange",
      lwd = 2.5)

# Vertical mean lines
abline(v=mean(rate.treated), col="steelblue",  lty=2)
abline(v=mean(rate.untreated), col="darkorange",  lty=2)

legend("topleft",
       legend = c("Treated", "Untreated"),
       col = c("steelblue", "darkorange"),
       lwd= 2, bty = "n")


# QUESTION 7
# Take the mean of the simulation matrix 
prob_Vmax <- mean(sim.matrix[,"Vmax[1]"] > sim.matrix[,"Vmax[2]"])
prob_K <- mean(sim.matrix[,"K[1]"] < sim.matrix[,"K[2]"])

cat("P(Vmax_treated>Vmax_untreated | data) =",prob_Vmax,"\n")
cat("P(K_treated<K_untreated | data) =",prob_K)


# QUESTION 8
# Calculate posteriors
k_1 <- sim.matrix[,"K[1]"]
treated_x <- k_1 * 4

k_2 <- sim.matrix[,"K[2]"]
untreated_x <- k_2 * 4

# Create result table and plots
df <- data.frame(
  value = c(treated_x, untreated_x),
  group = c(
    rep("treated", length(treated_x)),
    rep("untreated", length(untreated_x))
  )
)

df_summary <- df %>% group_by(group) %>% summarise(
  mean = mean(value),
  sd = sd(value),
)

df_summary

ggplot(data = df) + geom_density(aes(x = value, fill = group),alpha = 0.7) + labs(
  title = "Posterior density: Concetration for 80% of Vmax",
  x = "Concentration",
  y = "Density"
)