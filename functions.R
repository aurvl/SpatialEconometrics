library(kableExtra)

showNA <- function(data) {
   miss <- colSums(is.na(data))
   for (col in names(miss)){
      n <- 20 - nchar(col)
      cat(col, strrep(" ", n), miss[col], "\n")
   }
}

extract_splm <- function(model) {
   # Extraire les coefficients principaux
   coefs <- model$coefficients
   
   # Vérifier si les coefficients existent
   if (is.null(coefs) || length(coefs) == 0) {
      stop("Pas de coefficients trouvés dans le modèle.")
   }
   
   # Vérifier si les erreurs standard existent
   if (!is.null(model$vcov)) {
      std_errors <- sqrt(diag(model$vcov))
   } else {
      std_errors <- rep(NA, length(coefs))
   }
   
   # Calcul des p-values
   t_values <- coefs / std_errors
   p_values <- 2 * (1 - pnorm(abs(t_values)))  # Approximation des p-values
   
   # Ajouter des étoiles de significativité
   stars <- ifelse(p_values < 0.001, "***",
                   ifelse(p_values < 0.01, "**",
                          ifelse(p_values < 0.05, "*",
                                 ifelse(p_values < 0.1, ".", ""))))
   
   # Fusionner les coefficients et les étoiles
   coef_labels <- paste0(round(coefs, 4), " ", stars)
   
   # Vérifier si un coefficient d'autocorrélation rho existe
   if (!is.null(model$arcoef)) {
      rho_df <- data.frame(
         Term = "rho",
         Estimate = paste0(round(model$arcoef, 4), " "),
         Std.Error = "",  
         t.value = "",
         p.value = ""
      )
   } else {
      rho_df <- NULL
   }
   
   # Construire le tableau des coefficients
   coef_df <- data.frame(
      Term = names(coefs),
      Estimate = coef_labels,
      Std.Error = ifelse(is.na(std_errors), "", round(std_errors, 4)),
      t.value = ifelse(is.na(t_values), "", round(t_values, 4)),
      p.value = ifelse(is.na(p_values), "", round(p_values, 4))
   )
   
   # Fusionner avec rho si disponible
   if (!is.null(rho_df)) {
      coef_df <- rbind(coef_df, rho_df)
   }
   
   # **Extraire logLik, AIC et BIC correctement**
   logLik_value <- if (!is.null(model$logLik)) round(model$logLik, 4) else ""
   AIC_value <- if (!is.null(model$logLik)) round(-2 * model$logLik + 2 * length(coefs), 4) else ""
   BIC_value <- if (!is.null(model$logLik)) round(-2 * model$logLik + log(length(residuals(model))) * length(coefs), 4) else ""
   nobs_value <- if (!is.null(residuals(model))) length(residuals(model)) else ""
   
   stats_df <- data.frame(
      Term = c("Log-Likelihood", "AIC", "BIC", "Observations"),
      Estimate = c(logLik_value, AIC_value, BIC_value, nobs_value),
      Std.Error = "",  
      t.value = "",
      p.value = ""
   )
   
   # Ajouter les statistiques globales en bas du tableau
   coef_df <- rbind(coef_df, stats_df)
   
   return(coef_df)
}


combine_models <- function(models) {
   # Extraire les coefficients et stats globales de chaque modèle
   model_results <- lapply(names(models), function(model_name) {
      df <- extract_splm(models[[model_name]])
      df <- df %>%
         select(Term, Estimate) %>%
         rename(!!model_name := Estimate)  # Renommer la colonne Estimate par le nom du modèle
   })
   
   # Fusionner tous les modèles sur la colonne "Term"
   final_table <- Reduce(function(x, y) full_join(x, y, by = "Term"), model_results)
   
   return(final_table)
}

mod_sumsum <- function(final_table, format = "html") {
   kable(final_table, format = format) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"))
}

moran_resid_panel <- function(model, pdata, listw) {
   # Extraction des résidus
   res <- residuals(model)
   
   # Vérification : même taille que le panel ?
   if (nrow(pdata) != length(res)) {
      stop("Longueur des résidus ≠ nombre de lignes dans les données panel.")
   }
   
   # Ajout des résidus aux données
   pdata$residuals <- res
   
   # Moyenne des résidus par région
   library(dplyr)
   resid_region <- pdata %>%
      group_by(id) %>%
      summarise(mean_resid = mean(residuals, na.rm = TRUE)) %>%
      ungroup()
   
   # Identifiants dans listw
   W_regions <- attr(listw, "region.id")
   
   # Message diagnostic
   missing <- setdiff(W_regions, resid_region$id)
   if (length(missing) > 0) {
      message("Régions manquantes dans les résidus : ", paste(missing, collapse = ", "))
      resid_region <- resid_region[resid_region$id %in% W_regions, ]
   }
   
   # Réordonner selon l’ordre de listw
   resid_region <- resid_region[match(W_regions, resid_region$id), ]
   
   # Vérification finale
   if (any(is.na(resid_region$mean_resid))) {
      stop("Résidus manquants après réorganisation. Vérifie les correspondances entre id et listw.")
   }
   
   # Test de Moran
   library(spdep)
   moran.test(resid_region$mean_resid, listw)
}


plot_residuals_map <- function(pdata, shp, model, title = "Residuals Map") {
   library(dplyr)
   library(ggplot2)
   library(sf)
   
   # Ajouter les résidus
   res <- residuals(model)
   pdata$residuals <- res
   
   # Moyenne par région
   resid_region <- pdata %>%
      group_by(id) %>%
      summarise(mean_resid = mean(residuals, na.rm = TRUE)) %>%
      ungroup()
   
   # Fusion avec shapefile
   shp_map <- shp %>%
      left_join(resid_region, by = c("NUTS_ID" = "id"))
   
   # Carte
   ggplot(shp_map) +
      geom_sf(aes(fill = mean_resid), color = "white") +
      scale_fill_viridis_c(option = "plasma", name = "Residuals") +
      labs(title = title) +
      theme_minimal() +
      theme(axis.text = element_blank())
}


build_surface_weighted_W <- function(shp, f_surface = identity) {
   library(sf)
   library(fields)
   library(spdep)
   
   # Vérification des colonnes
   if (!"geometry" %in% colnames(shp)) {
      stop("Le shapefile ne contient pas de colonne 'geometry'")
   }
   if (!"NUTS_ID" %in% colnames(shp)) {
      stop("Le shapefile doit contenir une colonne 'NUTS_ID'")
   }
   
   # 1. Surface en km²
   shp$surface_km2 <- as.numeric(st_area(shp)) / 1e6
   
   # 2. Application de la transformation demandée (log, sqrt, etc.)
   transformed_surface <- f_surface(shp$surface_km2)
   
   # Vérification que pas de NA ou -Inf
   if (any(!is.finite(transformed_surface))) {
      stop("Transformation de surface a produit des valeurs non finies (NA, Inf).")
   }
   
   # 3. Coordonnées des centroïdes
   coords <- st_coordinates(st_centroid(st_geometry(shp)))
   
   # 4. Matrice des distances
   D <- rdist(coords)
   D[D == 0] <- NA
   
   # 5. Matrice des surfaces S_j (colonnes)
   S <- matrix(rep(transformed_surface, each = nrow(D)), nrow = nrow(D), byrow = FALSE)
   
   # 6. Calcul de W = S_j / d_ij
   W <- S / D
   diag(W) <- 0
   W[is.na(W)] <- 0
   
   # 7. Normalisation par ligne
   row_sums <- rowSums(W)
   row_sums[row_sums == 0] <- 1
   W_norm <- W / row_sums
   
   # 8. Attribution des noms
   rownames(W_norm) <- shp$NUTS_ID
   colnames(W_norm) <- shp$NUTS_ID
   
   # 9. Conversion en listw
   listw_custom <- mat2listw(W_norm, style = "W")
   
   return(listw_custom)
}

spml_metrics <- function(..., model_names = NULL) {
   models <- list(...)
   
   if (is.null(model_names)) {
      model_names <- paste0("Model_", seq_along(models))
   }
   
   results <- data.frame(
      Model = character(),
      LogLik = numeric(),
      AIC = numeric(),
      BIC = numeric(),
      N = numeric(),
      K = numeric(),
      stringsAsFactors = FALSE
   )
   
   for (i in seq_along(models)) {
      model <- models[[i]]
      logL <- model$logLik
      k <- length(model$coefficients) + 1
      n <- length(model$residuals)
      AIC <- -2 * logL + 2 * k
      BIC <- -2 * logL + log(n) * k
      
      # Enrichir l'objet modèle avec les métriques
      model$logLik_value <- logL
      model$AIC <- AIC
      model$BIC <- BIC
      model$n <- n
      model$k <- k
      models[[i]] <- model
      
      # Ajouter à la table
      results[i, ] <- c(model_names[i], round(logL, 4), round(AIC, 4), round(BIC, 4), n, k)
   }
   
   print(results, row.names = FALSE)
   invisible(models)
}

plot_model_residuals_comparison <- function(models, model_names, pdata) {
   library(dplyr)
   library(tidyr)
   library(ggplot2)
   
   resid_df <- data.frame()
   
   for (i in seq_along(models)) {
      model <- models[[i]]
      pdata_tmp <- pdata
      pdata_tmp$residuals <- residuals(model)
      
      mean_res <- pdata_tmp %>%
         group_by(id) %>%
         summarise(mean_residual = mean(residuals, na.rm = TRUE)) %>%
         mutate(model = model_names[i])
      
      resid_df <- bind_rows(resid_df, mean_res)
   }
   
   # Plot
   ggplot(resid_df, aes(x = reorder(id, mean_residual), y = mean_residual, fill = model)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, size = 6),
            legend.position = "bottom")
}

plot_residuals_maps_grid <- function(models, model_names, pdata, shp) {
   library(dplyr)
   library(ggplot2)
   library(sf)
   
   all_maps <- list()
   
   for (i in seq_along(models)) {
      model <- models[[i]]
      pdata_tmp <- pdata
      pdata_tmp$residuals <- residuals(model)
      
      mean_resid <- pdata_tmp %>%
         group_by(id) %>%
         summarise(mean_resid = mean(residuals, na.rm = TRUE)) %>%
         mutate(model = model_names[i])
      
      all_maps[[i]] <- mean_resid
   }
   
   resid_df <- bind_rows(all_maps)
   
   # Fusion avec shapefile
   shp_combined <- shp %>%
      left_join(resid_df, by = c("NUTS_ID" = "id"))
   
   # Carte facettée
   ggplot(shp_combined) +
      geom_sf(aes(fill = mean_resid), color = "white") +
      scale_fill_viridis_c(option = "plasma", name = "Résidu moyen") +
      facet_wrap(~ model) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            strip.text = element_text(size = 12, face = "bold"))
}

extract_residual_metrics <- function(model, pdata, listw) {
   res <- residuals(model)
   pdata$residuals <- res
   
   resid_region <- pdata %>%
      group_by(id) %>%
      summarise(
         mean_resid = mean(residuals, na.rm = TRUE),
         sd_resid = sd(residuals, na.rm = TRUE),
         abs_mean = mean(abs(residuals), na.rm = TRUE),
         var_resid = var(residuals, na.rm = TRUE)
      )
   
   mean_abs <- mean(resid_region$abs_mean, na.rm = TRUE)
   var_resid <- mean(resid_region$var_resid, na.rm = TRUE)
   
   # Moran's I
   region_ids <- attr(listw, "region.id")
   resid_vector <- resid_region$mean_resid[match(region_ids, resid_region$id)]
   I <- moran.test(resid_vector, listw)$estimate[1]
   
   return(data.frame(
      abs_resid_mean = round(mean_abs, 4),
      var_resid_mean = round(var_resid, 4),
      moran_I = round(I, 4)
   ))
}

extract_spatial_impacts <- function(model, listw, time, varnb) {
   imp <- summary(impacts(model, listw = listw, time = time, R = 100))
   data.frame(
      Direct = round(imp$res$direct[varnb], 4),
      Indirect = round(imp$res$indirect[varnb], 4),
      Total = round(imp$res$total[varnb], 4)
   )
}

plot_standardized_residuals_maps <- function(models, model_names, pdata, shp) {
   all_maps <- list()
   
   for (i in seq_along(models)) {
      model <- models[[i]]
      pdata_tmp <- pdata
      pdata_tmp$residuals <- residuals(model)
      
      std_resid <- pdata_tmp %>%
         group_by(id) %>%
         summarise(std_resid = mean(scale(residuals), na.rm = TRUE)) %>%
         mutate(model = model_names[i])
      
      all_maps[[i]] <- std_resid
   }
   
   resid_df <- bind_rows(all_maps)
   shp_combined <- shp %>% left_join(resid_df, by = c("NUTS_ID" = "id"))
   
   p <- ggplot(shp_combined) +
      geom_sf(aes(fill = std_resid), color = "white") +
      scale_fill_viridis_c(option = "plasma", name = "Résidu standardisé") +
      facet_wrap(~ model) +
      theme_minimal() +
      theme(axis.text = element_blank(), strip.text = element_text(size = 12, face = "bold"))
   print(p)
}


evaluate_spatial_models <- function(models, model_names, pdata, listws, shp, T, varnb = 1) {
   message("▶️ Comparaison graphique des résidus standardisés...")
   plot_standardized_residuals_maps(models, model_names, pdata, shp)
   
   message("\n📊 Comparaison tabulaire des résidus et effets spatiaux :")
   
   # Résidus et impacts
   table <- data.frame()
   for (i in seq_along(models)) {
      metrics <- extract_residual_metrics(models[[i]], pdata, listws[[i]])
      impacts <- extract_spatial_impacts(models[[i]], listws[[i]], T, varnb)
      
      row <- cbind(Model = model_names[i], metrics, impacts)
      table <- rbind(table, row)
   }
   
   print(table)
   invisible(table)
}

compare_predictive_metrics <- function(models, model_names, pdata) {
   library(dplyr)
   
   results <- data.frame(
      Model = character(),
      MSE = numeric(),
      RMSE = numeric(),
      MAE = numeric(),
      stringsAsFactors = FALSE
   )
   
   for (i in seq_along(models)) {
      model <- models[[i]]
      y_true <- model$model[[1]]  # variable dépendante
      y_hat <- fitted(model)      # valeurs prédites par le modèle
      
      mse <- mean((y_true - y_hat)^2, na.rm = TRUE)
      rmse <- sqrt(mse)
      mae <- mean(abs(y_true - y_hat), na.rm = TRUE)
      
      results[i, ] <- c(model_names[i], round(mse, 4), round(rmse, 4), round(mae, 4))
   }
   
   print(results, row.names = FALSE)
   invisible(results)
}
