library(tidyr)
library(dplyr)
library(mice)
library(ggplot2)
library(plm)
library(splm)
library(sf)
library(spdep)
library(sf)
library(viridis)
library(jtools)


theme_set(theme_minimal())

d_wide <- read.csv("inkar/data.csv")

# Data pre-processing
col_eng <- c("code_num", "spatial_unit", "year", "employed", "unemployment_rate",
             "land_price", "employees_academic_degree", "employees_no_qualification",
             "hospital_bed", "prop_foreigners", "waste_amount", "gdp")

d <- d_wide %>% 
  slice(-1) %>%
  select(- Aggregat) %>%
  pivot_longer(
    cols = - c(Kennziffer, Raumeinheit),
    names_to = c(".value", "year"),
    names_pattern = "(.+?)\\.?(\\d*)$") %>%
  replace(.=="", "0") %>%
  setNames(col_eng) %>%
  mutate(year = recode(year,
                       "0" = "2015",
                       "1" = "2016",
                       "2" = "2017",
                       "3" = "2018",
                       "4" = "2019",
                       "5" = "2020",
                       "6" = "2021"),
         d_employed = 100 * (log(employed) - log(dplyr::lag(employed))))%>%
  filter(year != "2015")
  

# Impute data
imputed <- mice(d, method="cart")
d <- complete(imputed)
anyNA(d)

# transform some variables
d <- d %>% mutate(ln_gdp = log(gdp),
                  sc_land_price = land_price / 1000,
                  sc_waste_amount = waste_amount / 10000,
                  year = as.factor(year),
                  spatial_unit = as.factor(spatial_unit)) %>% as_tibble()

# Visualisation
d %>%
  select("d_employed", "unemployment_rate", "sc_land_price", 
         "prop_foreigners", "sc_waste_amount", "ln_gdp", "hospital_bed") %>%
  tidyr::gather(key = "variable", value = "value") %>%
  ggplot(aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~variable, scale = "free") +
  labs(x = "", y = "", title = "", subtitle = "") +
  theme(legend.position = "none")


## with shape data
shapefile <- st_read("vg2500_krs/vg2500_krs.shp")
shapefile$ARS <- as.numeric(shapefile$ARS)

d_geo <- merge(shapefile, d, by.x = "ARS", by.y = "code_num", all.y = TRUE)

d_geo %>%
    filter(year == "2021") %>%
    mutate(d_employed_q = cut(d_employed, breaks = quantile(d_employed, probs = c(0, 0.20, 0.5, 0.8, 1)), 
                          include.lowest = TRUE)) %>%
    ggplot() +
    geom_sf(aes(fill = d_employed_q)) +
    scale_fill_viridis_d(option = "magma", direction = -1)

d_geo %>%
    filter(year == "2016") %>%
    mutate(d_employed_q = cut(d_employed, breaks = quantile(d_employed, probs = c(0, 0.20, 0.5, 0.8, 1)), 
                              include.lowest = TRUE)) %>%
    ggplot() +
    geom_sf(aes(fill = d_employed_q)) +
    scale_fill_viridis_d(option = "magma", direction = -1)


d %>%
  filter(year == "2021") %>%
  select("d_employed") %>%
  ggplot(aes(x = d_employed)) +
  geom_density(alpha = 0.5)

d1 <- d %>%
  filter(year == "2016") %>%
  select("d_employed")


# delete row in shapefile which is not recorded in main data
setdiff(unique(shapefile$ARS), unique(d$code_num))
shapefile <- shapefile %>% filter(ARS != 16056)
  
# panel data with group specific effect
d <- d %>% select(- code_num)

f <- d_employed ~ ln_gdp + unemployment_rate + sc_land_price + prop_foreigners +
  sc_waste_amount + hospital_bed
  
fit1 <- plm(f, data = d, model = "random")
plot_summs(fit1, robost  = TRUE)

# panel data with spacially correlated group specific effect

## creating weight matrix
krs_shp_mtrx <- sf::st_read("vg2500_krs")  %>% filter(ARS != 16056)
krs_mtrx <- poly2nb(krs_shp_mtrx)
W <- nb2listw(krs_mtrx, style="W")

# LM test for.\lambda = 0.
test1 <- bsktest(x = f, data = d, listw = W, test = "CLMlambda")
# adding special error term
fit2 <- spml(f, data = d, listw = W, model = "random", lag = FALSE, effect = c("individual"))

fit2_f <- spml(f, data = d, listw = W, model = "within", lag = FALSE, effect = c("individual"))



fit2_vis <- fit1
fit2_vis$coefficients <- fit2$coefficients; fit2_vis$vcov <- fit2$vcov

plot_summs(fit2_vis, robost  = TRUE)


# adding autoregressive term
fit3 <- spml(f, data = d, listw = W, model = "random", lag = TRUE)
fit3_vis <- fit1
fit3_vis$coefficients <- fit3$coefficients; fit3_vis$vcov <- fit3$vcov

plot_summs(fit3_vis, robost  = TRUE)

plot_lambda_conf_int <- function(lambda, vcov_lambda, crit){
  df <- tibble(var_name = "lag term",
               est = lambda,
               lower = lambda - crit * sqrt(vcov_lambda),
               upper = lambda + crit * sqrt(vcov_lambda))
  ggplot(df, aes(x = est, y = var_name)) +
    geom_point() +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, colour= "blue") +
    geom_vline(aes(xintercept = 0), colour = "black", linetype="dotted")
}

plot_lambda_conf_int(fit3$arcoef, fit3$vcov.arcoef, 1.96)


# impact measure
impac <- impacts(fit3, listw = W, time = length(unique(d$year)))
summary(impac)

  



