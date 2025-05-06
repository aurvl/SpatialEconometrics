rm(list=ls())
source("functions.R")

library(sf)
library(plm)
library(splm)
library(sp)
library(spdep)
library(dplyr)
library(ggplot2)
library(lmtest)
library(modelsummary)
library(plm)
library(stargazer)
library(tidyr)


data <- read.csv('./data.csv')
columns <- c('country', 'employment_15_64', 'year')

data[data == ""] <- NA
col_confu <- c('employment_15_64', 'tech_employment', 'HRST', 'regional_GDP', 
               'pop_density')

for (col in col_confu) {
   data[[col]] <- as.numeric(gsub("[^0-9.]", "", data[[col]]))
}

shp <- st_read("./data/shp/NUTS_RG_20M_2024_3035.shp", quiet = TRUE)
shp <- shp[shp$LEVL_CODE == 2, ]
bbox_europe <- st_bbox(c(xmin = -5, xmax = 50, ymin = 0, ymax = 70), 
                       crs = st_crs(4326))
bbox_europe <- st_transform(bbox_europe, 3035)
shp <- st_crop(shp, bbox_europe)

pdata <- data %>%
   select(-internet_usage) %>%
   filter(country %in% c('FR', 'ES', 'BE', 'IT', 'NL')) %>%
   mutate(regional_GDP = regional_GDP/1000, 
          log_pop_density = log(pop_density)) %>%
   drop_na() %>%
   group_by(id) %>%
   filter(all(2012:2021 %in% year)) %>%
   filter(year >= 2012 & year <= 2021) %>%
   ungroup() %>% 
   select(id, year, everything())

sdata <- shp %>% 
   select(NAME_LATN, geometry) %>%
   rename(geo = NAME_LATN) %>%
   left_join(pdata, by = "geo")

codes <- pdata$id
shp_filter <- shp[shp$NUTS_ID %in% codes,]

# Plot map
map <- ggplot(shp_filter, aes(fill = CNTR_CODE)) + geom_sf() + theme_light() + 
   theme(axis.text = element_blank())
map

head(pdata)
glimpse(pdata)
attach(pdata)


vars <- data.frame(pdata$female_unemp, pdata$tech_employment, pdata$HRST, 
                   pdata$regional_GDP, pdata$higher_edu, pdata$pop_density)

stargazer(vars, type = "text", title = "Descriptive stats")


## Data Visualization
femp23 <- data %>%
   filter(year == 2021)

weak_femp_reg <- femp23 %>%
   arrange(desc(female_unemp)) %>%
   slice_head(n = 10) %>%
   pull(geo)

reg_no_na <- data %>%
   group_by(geo) %>%
   filter(all(!is.na(female_unemp))) %>%
   pull(geo) %>% unique()

evol_f_unempl <- data %>%
   filter(geo %in% weak_femp_reg & geo %in% reg_no_na)

# Plot 1
ggplot(evol_f_unempl, aes(x = year, y = female_unemp, color= geo)) +
   geom_line(size = .8) +
   labs(y = 'Unemploy rate (%)', x = 'Time') +
   theme_light() +
   theme()


sdata_13 <- sdata[sdata$year == 2012, ]
sdata_23 <- sdata[sdata$year == 2021, ]

# Plot 2
ggplot() +
   geom_sf(data = sdata_13, aes(fill = female_unemp), color = "white") +
   scale_fill_viridis_c(option = "magma", name = "Unemployment rate (%)") +
   theme_minimal() + theme(axis.text = element_blank())

ggplot() +
   geom_sf(data = sdata_23, aes(fill = female_unemp), color = "white") +
   scale_fill_viridis_c(option = "magma", name = "Unemployment rate (%)") +
   theme_minimal() + theme(axis.text = element_blank())


# Plot 3
ggplot() +
   geom_sf(data = sdata_23, aes(fill = regional_GDP), color = "white") +
   scale_fill_viridis_c(option = "plasma", name = 'Annual Revenu (in K€)') +
   theme_minimal() + theme(axis.text = element_blank())


# Spatial models estimation
formula <- female_unemp ~ tech_employment + HRST + regional_GDP + 
   higher_edu + pop_density

pdata <- pdata.frame(pdata, index = c('id', 'year'))
pdim(pdata)


## Spatial Weighted matrix (W)
nb_queen <- poly2nb(shp_filter, queen = TRUE)
W <- nb2listw(nb_queen, style = "W", zero.policy = TRUE)  # W =Normalisation row-standardized
# zoro-policy : droping no neighboors
no_neighbors <- which(card(nb_queen) == 0)
filtered <- shp_filter[-no_neighbors, ]
new_codes <- filtered$NUTS_ID
pdata <- pdata[pdata$id %in% new_codes, ]

nb_queen <- poly2nb(filtered, queen = TRUE)
W <- spdep::nb2listw(nb_queen, style = "W", zero.policy = TRUE)
attr(W, "region.id") <- filtered$NUTS_ID
W

# Plot 4
plot(st_geometry(filtered), border="gray")
plot(nb_queen, st_coordinates(st_centroid(filtered)), add=TRUE, col="blue", pch=20)


## Non-spatial models
pdata <- pdata.frame(pdata, index = c('id', 'year'))

ols_model <- plm(formula, data = pdata, model = "pooling")
fe_ind <- plm(formula, data = pdata, model = "within", effect = "individual")
fe_tim <- plm(formula, data = pdata, model = "within", effect = "time")
fe_woway <- plm(formula, data = pdata, model = "within", effect = "twoways")
re_model <- plm(formula, data = pdata, model = "random")

models <- list("Pooled (OLS)" = ols_model, "Individual Fixed effects" = fe_ind, 
               "Time Fixed effects" = fe_tim, "Twoway Fixed effects" = fe_woway,
               "Random effects" = re_model)

modelsummary(models, stars = TRUE, gof_map = "all", statistic = "std.error", 
             output = "latex")


### Panel data Tests
# BP test
cat("====== BP Test : OLS vs FE/RE ======")
plmtest(ols_model, type = "bp")

# Test F
cat("\n====== F Test : FE vs OLS ======")
pFtest(fe_ind, ols_model)

# Twoway vs Indiv
cat("\n===== LR Test : Twoway vs Indiv =====\n")
lrtest(fe_tim, fe_woway)

# Hausman test (plm)
cat("\n====== Hausman Test : FE vs RE ======")
phtest(fe_woway, re_model)


### Spatial Hausman Tests
cat("\n====== Spatial Hausman Test : FE vs RE - SEM ======")
sphtest(formula, data =pdata, listw = W, spatial.model = "error", method = "ML")

cat("\n====== Spatial Hausman Test : FE vs RE - SAR ======") 
sphtest(formula, data =pdata, listw = W, spatial.model = "lag", method = "ML")



### Robust Spatial LM Tests
cat('OLS ============================\n')
slmtest(ols_model, W, test = "rlml")
slmtest(ols_model, W, test = "rlme")

cat('\n\nTwoway FE ======================\n')
slmtest(fe_woway, W, test = "rlml")
slmtest(fe_woway, W, test = "rlme")

cat('\n\nIndividual FE ==================\n')
slmtest(fe_ind, W, test = "rlml")
slmtest(fe_ind, W, test = "rlme")

cat('\n\nTime FE ========================\n')
slmtest(fe_tim, W, test = "rlml")
slmtest(fe_tim, W, test = "rlme")

cat('\n\nRE =============================\n')
slmtest(re_model, W, test = "rlml")
slmtest(re_model, W, test = "rlme")


## Spatial econometrics models
pool_model <- spml(formula, data = pdata, listw = W, lag = FALSE,
                   model = "pooling")

sem_model <- spml(formula, data = pdata, listw = W, spatial.error = "b",
                  model = "within", effect = "twoway")

sar_model <- spml(formula, data = pdata, listw = W, lag = TRUE, model = "within",
                  effect = "twoway", spatial.error ="none")

sac_model <- spgm(formula, data = pdata, listw = W, listw2 = W, model = "within",
                  lag = TRUE, spatial.error = TRUE, moments = "fullweights")

models <- list(
   "Pooled OLS" = pool_model,
   "SEM FE" = sem_model,
   "SAR FE" = sar_model,
   "SAC FE" = sac_model
)
table <- combine_models(models)
mod_sumsum(table, format = "latex")


# Results
## SAR Fixed effects
summary(sar_model)

## Direct, indirect, totel effects
T <- length(unique(pdata$year))
impacts_sar <- impacts(sar_model, listw = W, time = T)
summary(impacts_sar)

# ==============================================================================

# Calculer la surface de chaque région (en m²)
filtered$surface <- st_area(filtered)
# Convertir en km² si besoin
filtered$surface_km2 <- as.numeric(filtered$surface) / 1e6

coords <- st_coordinates(st_centroid(filtered))

# Matrice des distances euclidiennes
library(fields)
D <- rdist(coords)  # matrice n x n

# Surface des j (colonnes)
S <- matrix(rep(filtered$surface_km2, each = nrow(D)), nrow = nrow(D), byrow = FALSE)
# Éviter les divisions par 0 (diagonale)
D[D == 0] <- NA

W_custom <- S / D
diag(W_custom) <- 0  # pas d'auto-influence

# Remplacer les NA par 0 avant normalisation
W_custom[is.na(W_custom)] <- 0

# Row standardisation (somme des lignes = 1)
row_sums <- rowSums(W_custom)
row_sums[row_sums == 0] <- 1  # éviter division par zéro si une ligne est toute à 0

W_custom_norm <- W_custom / row_sums

# Créer la matrice de voisinage
rownames(W_custom_norm) <- filtered$NUTS_ID
colnames(W_custom_norm) <- filtered$NUTS_ID
W_listw <- mat2listw(W_custom_norm, style = "W", zero.policy = TRUE)
W_listw

sar_model_custom <- spml(formula, data = pdata, listw = W_listw, lag = TRUE, 
                         model = "within", effect = "twoway", spatial.error ="none")
summary(sar_model_custom)
summary(sar_model)

l1 <- sar_model_custom$logLik
l2 <-  sar_model$logLik
LR <- 2 * (as.numeric(l1) - as.numeric(l2))
LR
# Degrés de liberté = nb paramètres ajoutés
df <- 1

# p-value
pval <- 1 - pchisq(LR, df)

# Résultat
cat("LR test statistic =", LR, "\n")
cat("Degrees of freedom =", df, "\n")
cat("p-value =", pval, "\n")

res1 <- residuals(sar_model)         # modèle avec W "classique"
res2 <- residuals(sar_model_custom)  # modèle avec W pondéré par surface

# Test de Moran's I sur les résidus
# Test avec ton modèle custom et la matrice pondérée par surface
moran_resid_panel(sar_model_custom, pdata, W_listw)

# Ou avec la matrice de contiguïté classique
moran_resid_panel(sar_model, pdata, W)

plot_residuals_map(pdata, filtered, sar_model_custom, title = "SAR Résidus (W = S/d)")
plot_residuals_map(pdata, filtered, sar_model, title = "SAR Résidus (W = contiguïté)")

models <- spml_metrics(sar_model, sar_model_custom, model_names = c("SAR W classique", "SAR W = S/d"))

# Ensuite :
models[[1]]$AIC
models[[2]]$logLik_value

# ==============================================================================
W
W_iden <- build_surface_weighted_W(filtered, f_surface = identity)
W_log      <- build_surface_weighted_W(filtered, f_surface = log)
W_sqrt     <- build_surface_weighted_W(filtered, f_surface = sqrt)

sar_W <- spml(formula, data = pdata, listw = W, lag = TRUE, model = "within",
              effect = "twoway", spatial.error = "none")

sar_iden <- spml(formula, data = pdata, listw = W_iden, lag = TRUE, model = "within",
                     effect = "twoway", spatial.error = "none")

sar_log <- spml(formula, data = pdata, listw = W_log, lag = TRUE, model = "within",
                effect = "twoway", spatial.error = "none")

sar_sqrt <- spml(formula, data = pdata, listw = W_sqrt, lag = TRUE, model = "within",
                 effect = "twoway", spatial.error = "none")

models <- spml_metrics(
   sar_W,
   sar_iden,
   sar_log,
   sar_sqrt,
   model_names = c("Contiguity W", "S/d", "log(S)/d", "sqrt(S)/d")
)

plot_model_residuals_comparison(
   models = list(sar_W, sar_iden, sar_log, sar_sqrt),
   model_names = c("Contiguity W", "S/d", "log(S)/d", "sqrt(S)/d"),
   pdata = pdata
)

plot_residuals_maps_grid(
   models = list(sar_W, sar_iden, sar_log, sar_sqrt),
   model_names = c("Contiguity W", "S/d", "log(S)/d", "sqrt(S)/d"),
   pdata = pdata,
   shp = filtered
)

evaluate_spatial_models(
   models = list(sar_W, sar_iden, sar_log, sar_sqrt),
   model_names = c("Contiguity W", "S/d", "log(S)/d", "sqrt(S)/d"),
   pdata = pdata,
   listws = list(W, W_iden, W_log, W_sqrt),
   shp = filtered,
   T = length(unique(pdata$year)),
   varnb = 2
)

compare_predictive_metrics(
   models = list(sar_W, sar_iden, sar_log, sar_sqrt),
   model_names = c("Contiguity W", "S/d", "log(S)/d", "sqrt(S)/d"),
   pdata = pdata
)
