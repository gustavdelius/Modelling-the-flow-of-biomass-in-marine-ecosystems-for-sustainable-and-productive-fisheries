library(mizer)
library(mizerExperimental)
library(tidyverse)

params <- newSingleSpeciesParams(lambda = 2.05)#lambda dictates the resource/background spectrum
plotSpectra(params, power = 0)#power is essentially just saying whther the y axis will be number density or biomass
####### Exercise 1 ########
params_2 <- newSingleSpeciesParams(lambda=2.1)
plotSpectra(params,power = 0)
##########################
n <- initialN(params)#gives us the actual numerical values of the number density at each size (binned)
dimnames(n)

n[1,61]#returns average number density,not number of fish

numbers <- n * dw(params)#actually contains the data for the number of fish in each size class

numbers[1, 61]#returns number of individuals in the size class between 1-1.12 grams
#Exercise2############
dimnames(n)
n[1,80:87]#2d array
sum(numbers[1,81:86])
#######################

plotSpectra(params)#plots against biomass density

plotSpectra(params, power = 2)#plots against biomass density against log weight

#plotCDF(params, power = 0)#doesn't work for some reason

biomass_density <- n * w(params)#w returns the weights at the start of eaach size

biomass <- biomass_density * dw(params)

# Initialise an array with the right dimensions
cumulative_biomass <- biomass
# Calculate the cumulative sum of all biomasses in previous bins
cumulative_biomass[] <- cumsum(biomass)
# Normalise this so that it is given as a percentage of the total biomass
cdf <- cumulative_biomass / cumulative_biomass[1, 101] * 100
# Melt the array to a data frame and then plot
p_biomass_cdf <- ggplot(melt(cdf), aes(x = w, y = value)) +
  geom_line() + 
  labs(x = "Weight [g]",
       y = "% of total biomass")
p_biomass_cdf

plotCDF(params, log_x = FALSE)

growth_rate <- getEGrowth(params)
growth_rate[1, 61]
growth_rate[1, 61] / 365
plot(growth_rate, log_x = TRUE, log_y = TRUE)

w_small_fish <- w(params)[w(params) <= 10]
g_small_fish <- growth_rate[w(params) <= 10]

lm(log(g_small_fish) ~ log(w_small_fish))

############ Exercise 3 #############
mort_rate <- getMort(params)
w_small_fish <- w(params)[w(params) <= 10]
m_small_fish <- mort_rate[w(params) <= 10]
lm(log(m_small_fish) ~ log(w_small_fish))
#####################################

plot(n, log_x = TRUE, log_y = TRUE)

n_small_fish <- n[w(params) <= 10]
lm(log(n_small_fish) ~ log(w_small_fish))#alternate way to find the allemetry coefficient

plotSpectra(params, wlim = c(10, NA))#w axis is 10 and above

p <- plot(setReproduction(params), log_x = FALSE)
p

species_params(params)
select(species_params(params), w_mat, w_mat25, w_max, m)
p + geom_vline(xintercept = species_params(params)$w_mat, lty = 2) +
  geom_vline(xintercept = species_params(params)$w_mat25, lty = 2, col = "grey")
#adds line at 25 and 50 percent maturity rate

params_changed_maturity <- params

given_species_params(params_changed_maturity)$w_mat <- 40
given_species_params(params_changed_maturity)$w_mat25 <- 30
select(species_params(params_changed_maturity), w_mat, w_mat25, w_max, m)

pc <- plot(psi(params_changed_maturity), log_x = FALSE) +
  geom_vline(xintercept = species_params(params_changed_maturity)$w_mat, 
             lty = 2) +
  geom_vline(xintercept = species_params(params_changed_maturity)$w_mat25, 
             lty = 2, col = "grey")
pc

addPlot(p, psi(params_changed_maturity))#for comparison

p <- plot(getEGrowth(params), log_x = FALSE)
addPlot(p, getEGrowth(params_changed_maturity))

params_changed_maturity <- steadySingleSpecies(params_changed_maturity)

plotSpectra2(params, name1 = "Early maturity",
             params_changed_maturity, name2 = "Late maturity",
             power = 2, resource = FALSE, wlim = c(10, NA))