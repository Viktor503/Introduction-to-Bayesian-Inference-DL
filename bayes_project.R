# Import necessary libraries
library(rjags)
library(coda)
library(ggplot2)
library(dplyr)
library(kableExtra)
library(knitr)
theme_set(theme_minimal())

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

cat(readLines("enzyme_model.txt"), sep="\n")

# QUESTION 2
my.inits <- list(
  list(
    Vmax = c(200, 100),
    K = c(0.1, 0.1),
    tau  = 1,
    .RNG.name = "base::Wichmann-Hill",
    .RNG.seed = 1
  ),
  list(
    Vmax = c(150, 80),
    K = c(0.2, 0.15),
    tau = 2,
    .RNG.name = "base::Marsaglia-Multicarry",
    .RNG.seed = 2
  ),
  list(
    Vmax = c(250, 120),
    K = c(0.05, 0.08),
    tau = 0.5,
    .RNG.name = "base::Super-Duper",
    .RNG.seed = 3
  )
)

params <- c("Vmax", "K", "sigma")

jags <- jags.model(
  file = "enzyme_model.txt",
  data = my.data.jags,
  inits = my.inits,
  n.chains = 3
)

# check trace plots
param.names <- c("K (treated)", "K (untreated)", 
                 "Vmax (treated)", "Vmax (untreated)", 
                 "Sigma")

par(mfrow = c(2, 3))

for (i in 1:5) {
  traceplot(enzyme.sim[, i],
            main = paste("Trace plot:", param.names[i]),
            xlab = "Iteration",
            ylab = "Parameter value")
}

# convergence test
# Gelman-Rubin convergence diagnostics
gel <- gelman.diag(enzyme.sim)

gel.table <- data.frame(
  Parameter = c("K (treated)", "K (untreated)", 
                "Vmax (treated)", "Vmax (untreated)", 
                "Sigma"),
  `Point Estimate` = gel$psrf[, 1],
  `Upper C.I` = gel$psrf[, 2],
  check.names = FALSE
)

gel.table[, c("Point Estimate", "Upper C.I")] <- 
  round(gel.table[, c("Point Estimate", "Upper C.I")], 2)

kable(gel.table,
      caption = "Gelman-Rubin convergence diagnostics",
      align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE)


# QUESTION 3
# Full posterior summary
sim.summary <- summary(enzyme.sim)

posterior.table <- cbind(
  sim.summary$statistics[, c("Mean", "SD")],
  sim.summary$quantiles[, c("2.5%", "50%", "97.5%")]
)

# Round the values
posterior.table <- round(posterior.table, 4)

# Rename rows
rownames(posterior.table) <- c("K (treated)", "K (untreated)", 
                               "Vmax (treated)", "Vmax (untreated)", 
                               "Sigma")

kable(posterior.table,
      caption = "Posterior summaries for all model parameters",
      col.names = c("Mean", "SD", "2.5%", "Median", "97.5%"),
      align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE) %>%
  add_header_above(c(" " = 1, 
                     "Point estimates" = 2, 
                     "95% Credible interval" = 3))


param.names <- c("K (treated)", "K (untreated)", 
                 "Vmax (treated)", "Vmax (untreated)", 
                 "Sigma")

param.units <- c("Concentration", "Concentration",
                 "Reaction rate", "Reaction rate",
                 "Reaction rate")

par(mfrow = c(2, 3),
    cex.main = 1)

for (i in 1:5) {
  densplot(enzyme.sim[, i],
           main = paste("Posterior density:", param.names[i]),
           xlab = param.units[i],
           ylab = "Posterior density")
}

# QUESTION 4
# Extract posterior samples as a matrix
samples <- as.matrix(enzyme.sim)

# Compute ratios sample by sample
vmax.ratio <- samples[, "Vmax[1]"] / samples[, "Vmax[2]"]
k.ratio <- samples[, "K[1]"] / samples[, "K[2]"]

# Posterior summaries
ratio.summary <- rbind(
  c(mean(vmax.ratio), sd(vmax.ratio), 
    quantile(vmax.ratio, c(0.025, 0.5, 0.975))),
  c(mean(k.ratio), sd(k.ratio), 
    quantile(k.ratio, c(0.025, 0.5, 0.975)))
)

rownames(ratio.summary) <- c("Vmax (treated / untreated)", 
                             "K (treated / untreated)")
colnames(ratio.summary) <- c("Mean", "SD", "2.5%", "Median", "97.5%")

kable(round(ratio.summary, 4),
      caption = "Posterior summaries for parameter ratios",
      align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE) %>%
  add_header_above(c(" " = 1,
                     "Point estimates" = 2,
                     "95% Credible interval" = 3))


# QUESTION 5
# Concentration sequence
x.seq <- seq(0.1, 1, length.out = 100)

# Extract samples
samples <- as.matrix(enzyme.sim)
n.samp  <- nrow(samples)

# Compute posterior curves for each sample and concentration
curve.treated  <- matrix(NA, nrow = n.samp, ncol = length(x.seq))
curve.untreated <- matrix(NA, nrow = n.samp, ncol = length(x.seq))

for (j in seq_along(x.seq)) {
  x <- x.seq[j]
  curve.treated[, j]   <- (samples[, "Vmax[1]"] * x) / 
    (samples[, "K[1]"] + x)
  curve.untreated[, j] <- (samples[, "Vmax[2]"] * x) / 
    (samples[, "K[2]"] + x)
}

# Summarise: mean and 95% credible band at each concentration
summarise.curve <- function(mat) {
  data.frame(
    x = x.seq,
    mean = apply(mat, 2, mean),
    lower = apply(mat, 2, quantile, 0.025),
    upper = apply(mat, 2, quantile, 0.975)
  )
}

df.treated <- summarise.curve(curve.treated)
df.untreated <- summarise.curve(curve.untreated)

df.treated$Group <- "Treated"
df.untreated$Group <- "Untreated"

df.all <- rbind(df.treated, df.untreated)

ggplot(df.all, aes(x = x, colour = Group, fill = Group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2, colour = NA) +
  geom_line(aes(y = mean), linewidth = 1) +
  labs(title = "Fitted Michaelis-Menten curves with 95% credible bands",
       x = "Substrate concentration",
       y = "Reaction rate",
       colour = "Group",
       fill = "Group") +
  theme_minimal() +
  theme(legend.position = "top")


# QUESTION 6
# Extract posterior samples
samples <- as.matrix(enzyme.sim)

# Compute posterior reaction rate at x = 0.5 for both groups
x <- 0.5

# Plug x = 0.5 into Michaelis-Menten for every posterior sample
rate.treated <- (samples[, "Vmax[1]"] * x) /
  (samples[, "K[1]"] + x)

rate.untreated <- (samples[, "Vmax[2]"] * x) /
  (samples[, "K[2]"] + x)

# Posterior summaries
rate.summary <- rbind(
  c(mean(rate.treated),   sd(rate.treated),
    quantile(rate.treated,   c(0.025, 0.5, 0.975))),
  c(mean(rate.untreated), sd(rate.untreated),
    quantile(rate.untreated, c(0.025, 0.5, 0.975)))
)

rownames(rate.summary) <- c("Rate at x=0.5 (treated)", 
                            "Rate at x=0.5 (untreated)")
colnames(rate.summary) <- c("Mean", "SD", "2.5%", "Median", "97.5%")

kable(round(rate.summary, 4),
      caption = "Posterior reaction rate at concentration x = 0.5",
      align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE) %>%
  add_header_above(c(" " = 1,
                     "Point estimates" = 2,
                     "95% Credible interval" = 3))


# Combine into data frame for plotting
rate.df <- data.frame(
  Rate  = c(rate.treated, rate.untreated),
  Group = rep(c("Treated", "Untreated"), each = nrow(samples))
)

# Density plot
ggplot(rate.df, aes(x = Rate, fill = Group, colour = Group)) +
  geom_density(alpha = 0.7) +
  labs(title = "Posterior reaction rate at x = 0.5",
       x = "Reaction rate",
       y = "Posterior density") +
  theme_minimal() +
  theme(legend.position = "top")


# QUESTION 7
# Take the mean of the simulation matrix 
prob_Vmax <- mean(samples[,"Vmax[1]"] > samples[,"Vmax[2]"])
prob_K <- mean(samples[,"K[1]"] < samples[,"K[2]"])

cat("P(Vmax treated > Vmax untreated | data) =", round(prob_Vmax, 4), "\n",
    "P(K treated < K untreated | data) =", round(prob_K,    4), "\n")


# QUESTION 8
# Calculate posteriors
k_1 <- samples[,"K[1]"]
treated_x <- k_1 * 4

k_2 <- samples[,"K[2]"]
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
  Mean = mean(value),
  SD = sd(value),
  `2.5%` = quantile(value, 0.025),
  Median = quantile(value, 0.50),
  `97.5%`= quantile(value, 0.975)
) %>%
  mutate(across(where(is.numeric), ~ round(., 4)),
         group = c("Treated", "Untreated"))

kable(df_summary,
      caption = "Posterior summaries for concentration at 80\\% of Vmax",
      col.names = c("Group", "Mean", "SD", "2.5%", "Median", "97.5%"),
      align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE) %>%
  add_header_above(c(" " = 1,
                     "Point estimates" = 2,
                     "95% Credible interval" = 3))

ggplot(data = df) + geom_density(aes(x = value, fill = group),alpha = 0.7) + labs(
  title = "Posterior summaries for concentration at 80% of Vmax",
  x = "Concentration",
  y = "Density"
)