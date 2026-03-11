##########################################################################
# This script contains the recursive algorithm to generate a tree        #
# from the potential surface.  Make sure to run covariates_notdr.r first #
# to grab the euclid_vec_sc2                                             #
##########################################################################

library(mvtnorm)
library(gganimate)
library(rnaturalearth)

source("./covariates_nodt.r")

# Disable s2 processing to avoid self-intersection issues
sf_use_s2(FALSE)

# Load Africa map
africa <- ne_countries(continent = "Africa", returnclass = "sf") #shape file of Africa

# Load natural water bodies
africa_lakes <- ne_download(scale = 10, category = "physical", 
                            type = "lakes", returnclass = "sf")
africa_rivers <- ne_download(scale = 10, category = "physical", 
                             type = "rivers_lake_centerlines", 
                             returnclass = "sf")


# Ensure geometries are valid
africa_lakes <- st_make_valid(africa_lakes)
africa_rivers <- st_make_valid(africa_rivers)

# Remove invalid geometries if any still exist
africa_lakes <- africa_lakes[st_is_valid(africa_lakes), ]
africa_rivers <- africa_rivers[st_is_valid(africa_rivers), ]

# Perform spatial intersection
africa_lakes <- st_intersection(africa_lakes, africa)
africa_rivers <- st_intersection(africa_rivers, africa)

# Convert polygons to individual coordinate points
africa_coords <- africa %>%
  st_cast("MULTIPOLYGON") %>%
  st_cast("POLYGON") %>%
  st_cast("LINESTRING") %>%
  st_cast("POINT") %>%
  st_coordinates() %>%
  as.data.frame()


# --- load the betameans for potential surface --- #

betaMeans = read.csv('./betameans_r.csv')

# --- Load Krieging data --- #

etaDF = read.csv("./etaDF.csv")

# --- load the OOS set --- #

grid_land_df2 = read.csv("./for_krieg.csv")

grid_land_df2_copy = grid_land_df2

grid_land_df2$pred1 <- with(grid_land_df2,
    long_sc * betaMeans$betaMeans_r[1] +
    lat_sc * betaMeans$betaMeans_r[2] +
    elevation * betaMeans$betaMeans_r[3] +
    euclid_vec_sc2 * betaMeans$betaMeans_r[4] +
    (long_sc^2) * betaMeans$betaMeans_r[5] +
    (lat_sc^2) * betaMeans$betaMeans_r[6] +
    (long_sc * lat_sc) * betaMeans$betaMeans_r[7] +
    cat_5000 * betaMeans$betaMeans_r[8]
    )

grid_land_df2$pred2 <- with(grid_land_df2,
    long_sc * betaMeans$betaMeans_r[1] +
    lat_sc * betaMeans$betaMeans_r[2] +
    elevation * betaMeans$betaMeans_r[3] +
    euclid_vec_sc2 * betaMeans$betaMeans_r[4] +
    (long_sc^2) * betaMeans$betaMeans_r[5] +
    (lat_sc^2) * betaMeans$betaMeans_r[6] +
    (long_sc * lat_sc) * betaMeans$betaMeans_r[7] +
    cat_4_5000 * betaMeans$betaMeans_r[8]
    )

grid_land_df2$pred3 <- with(grid_land_df2,
    long_sc * betaMeans$betaMeans_r[1] +
    lat_sc * betaMeans$betaMeans_r[2] +
    elevation * betaMeans$betaMeans_r[3] +
    euclid_vec_sc2 * betaMeans$betaMeans_r[4] +
    (long_sc^2) * betaMeans$betaMeans_r[5] +
    (lat_sc^2) * betaMeans$betaMeans_r[6] +
    (long_sc * lat_sc) * betaMeans$betaMeans_r[7] +
    cat_3200 * betaMeans$betaMeans_r[8]
    )


grid_land_df2$Epi1 = etaDF$z + grid_land_df2$pred1
grid_land_df2$Epi2 = etaDF$z + grid_land_df2$pred2
grid_land_df2$Epi3 = etaDF$z + grid_land_df2$pred3

grid_land_df2_copy$Epi3 = with(grid_land_df2_copy,
  cat_3200 * betaMeans$betaMeans_r[8]
)

# --- Data Preparation --- #



sites = as.matrix(
  time_data[, c("longitude", "latitude")]
)


r0_idx = 9
diff_matrix <- diag(2)

# Combine your prediction vectors into a list
# pred1: 5000BP, pred2: 4200BP, pred3: 3000BP
pred_list <- list(
  grid_land_df2$Epi1, 
  grid_land_df2$Epi2, 
  grid_land_df2$Epi3
)


local_grad_2d <- function(x, z, y, x0, z0, h, degree = 2, k = NULL) {
    n <- length(y)
    u <- x - x0
    v <- z - z0
    d <- sqrt(u^2 + v^2)
    
    # 1. Ensure we have enough neighbors for a Quadratic fit (needs at least 6)
    # If the radius h is too small, we force it to take the 15 nearest neighbors
    idx <- which(d <= 3 * h)
    if (length(idx) < 15) {
        idx <- order(d)[1:min(20, n)] 
    }
    
    u_loc <- u[idx]; v_loc <- v[idx]; y_loc <- y[idx]; d_loc <- d[idx]
    w <- exp(-(d_loc^2) / (2 * h^2))
    sw <- sqrt(w)
    
    # 2. Build Design Matrix
    # We use tryCatch to handle cases where the Quadratic math fails
    grad_out <- tryCatch({
        if (degree == 2) {
            X <- cbind(1, u_loc, v_loc, u_loc^2, u_loc * v_loc, v_loc^2)
            beta_hat <- qr.solve(X * sw, y_loc * sw)
            c(dx = beta_hat[2], dz = beta_hat[3])
        } else {
            X <- cbind(1, u_loc, v_loc)
            beta_hat <- qr.solve(X * sw, y_loc * sw)
            c(dx = beta_hat[2], dz = beta_hat[3])
        }
    }, error = function(e) {
        # FALLBACK: If Quadratic fails, try a Simple Linear fit
        X_simple <- cbind(1, u_loc, v_loc)
        beta_simple <- qr.solve(X_simple * sw, y_loc * sw)
        return(c(dx = beta_simple[2], dz = beta_simple[3]))
    })
    
    return(list(grad = grad_out))
}



# 1. Define the 415 Era Labels (The "Who")
site_eras <- ifelse(time_data$x5000bp == 1, "5000BP",
               ifelse(time_data$x4_5000bp == 1, "4000BP", 
                 ifelse(time_data$x3200bp == 1, "3200BP", NA)))



library(mvtnorm)

build_migration_tree <- function(S, r0_idx, pred_list, site_eras, grid_coords, dt_g, Q, rho, p_min, K_max) {
  
  # 1. Setup Data Structures
  S_clean <- unname(matrix(as.numeric(as.matrix(S[, 1:2])), ncol = 2))
  N <- nrow(S_clean)
  
  # Initialize the collector for ALL probabilities seen in each era
  era_probs_collector <- list("5000BP" = c(), "4000BP" = c(), "3200BP" = c())
  
  # Root node
  root_node <- list(idx = r0_idx, loc = S_clean[r0_idx, ], depth = 0, 
                    parent_idx = NULL, era = site_eras[r0_idx], prob_table = NULL)
  
  V_global <- c(r0_idx)
  All_Nodes <- list(root_node)
  era_labels <- c("5000BP", "4000BP", "3200BP")
  
  # 2. Main Era Loop
  for (e_idx in 1:length(era_labels)) {
    current_era <- era_labels[e_idx]
    current_potentials <- pred_list[[e_idx]] 
    
    target_pool <- which(site_eras == current_era & !(1:N %in% V_global))
    Frontier <- All_Nodes
    
    k_era <- 0
    while (length(target_pool) > 0 && k_era < K_max) {
      F_next <- list()
      for (v in Frontier) {
        U <- intersect(target_pool, setdiff(1:N, V_global))
        if (length(U) == 0) next
        
        # 3. Calculate Gradient and Target
        grad_H <- local_grad_2d(grid_coords$longitude, grid_coords$latitude, 
                                current_potentials, v$loc[1], v$loc[2], h = 2)
        
        m <- -grad_H$grad * dt_g
        target_mu <- v$loc + m
        
        # 4. Probability Calculations
        S_subset <- S_clean[U, , drop = FALSE]
        diffs_U <- S_subset - rep(target_mu, each = nrow(S_subset))
        
        # dmvnorm gives raw weights (density)
        weights <- dmvnorm(diffs_U, mean = c(0, 0), sigma = Q * dt_g)
        
        if (sum(weights, na.rm = TRUE) == 0) next 
        
        # Normalize weights to get individual probabilities
        probs <- weights / sum(weights)

        # STORE: All individual probabilities found in this decision
        era_probs_collector[[current_era]] <- c(era_probs_collector[[current_era]], probs)

        # 5. Create Enhanced Report Table (Warning-Free)
        current_prob_report <- data.frame(
          site_idx = U,
          prob     = probs,
          m_lon    = unname(m[1]),      
          m_lat    = unname(m[2]),      
          target_x = unname(target_mu[1]), 
          target_y = unname(target_mu[2]),
          row.names = NULL 
        )
        current_prob_report <- current_prob_report[order(current_prob_report$prob, decreasing = TRUE), ]

        # 6. Branching Logic
        is_split <- FALSE
        if (length(U) >= 2) {
          j1 <- current_prob_report$site_idx[1]; p_j1 <- current_prob_report$prob[1]
          j2 <- current_prob_report$site_idx[2]; p_j2 <- current_prob_report$prob[2]
          
          if ((p_j1 / p_j2 <= rho) && (p_j2 >= p_min)) {
            is_split <- TRUE
            for (idx_child in c(j1, j2)) {
              if (!(idx_child %in% V_global)) {
                child <- list(idx = idx_child, loc = S_clean[idx_child,], 
                              depth = v$depth + 1, parent_idx = v$idx, 
                              era = current_era, prob_table = current_prob_report)
                F_next <- c(F_next, list(child))
                V_global <- c(V_global, idx_child)
                target_pool <- setdiff(target_pool, idx_child)
              }
            }
          }
        }
        
        if (!is_split) {
          U_curr <- intersect(target_pool, setdiff(1:N, V_global))
          if (length(U_curr) > 0) {
            # Use the calculated probs for sampling
            curr_probs <- probs[U %in% U_curr]
            if (sum(curr_probs, na.rm = TRUE) > 0) {
              J <- if(length(U_curr) == 1) U_curr else sample(U_curr, 1, prob = curr_probs/sum(curr_probs))
              
              child <- list(idx = J, loc = S_clean[J,], depth = v$depth + 1, 
                            parent_idx = v$idx, era = current_era,
                            prob_table = current_prob_report)
              F_next <- c(F_next, list(child))
              V_global <- c(V_global, J)
              target_pool <- setdiff(target_pool, J)
            }
          }
        }
      }
      if (length(F_next) == 0) break 
      All_Nodes <- c(All_Nodes, F_next)
      Frontier <- F_next
      k_era <- k_era + 1
    }
  }
  
  return(list(tree = All_Nodes, histograms = era_probs_collector))
}

output <- build_migration_tree(
  S = sites, 
  r0_idx = 9, 
  pred_list = pred_list, 
  site_eras = site_eras, 
  grid_coords = grid_land_df2, 
  dt_g = 1, Q = diag(2), rho = 1.05, p_min = 0.01, K_max = 1000
)

# Visualization Script
par(mfrow = c(1, 3)) # One plot per era

for (era in names(output$histograms)) {
  # We use log10 to see the spread of probabilities clearly
  # Adding a tiny constant to avoid log(0)
  log_probs <- log10(output$histograms[[era]] + 1e-10)
  
  hist(log_probs, 
       main = paste("Prob. Distribution:", era),
       xlab = "log10(Probability)", 
       col = "darkseagreen", 
       border = "white",
       breaks = 30)
}

  hist(log10(output$histograms[["4000BP"]] + 1e-10), 
       main = paste("Prob. Distribution:", era),
       xlab = "log10(Probability)", 
       col = "darkseagreen", 
       border = "white",
       xlim = c(-20, 10),
       breaks = 30)





# --- 5. Plotting --- #
plot_migration_tree <- function(tree_nodes, all_sites) {
    all_sites_df <- as.data.frame(all_sites); colnames(all_sites_df) <- c("x", "y")
    node_lookup <- setNames(tree_nodes, sapply(tree_nodes, `[[`, "idx"))

    segments <- do.call(rbind, lapply(tree_nodes, function(node) {
        if (is.null(node$parent_idx)) return(NULL)
        parent <- node_lookup[[as.character(node$parent_idx)]]
        data.frame(x = parent$loc[1], y = parent$loc[2], xend = node$loc[1], 
                   yend = node$loc[2], depth = node$depth, era = node$era,
                   id = paste(parent$idx, node$idx, sep = "_"))
    }))

    ggplot() +
        geom_sf(data = africa, fill = "gray95", color = "gray80") +
        geom_sf(data = africa_lakes, fill = "blue", color = "blue") +
        geom_sf(data = africa_rivers, color = "blue") +
        geom_contour_filled(
        data = grid_land_df2_copy, aes(x = longitude, y = latitude, z = Epi3), 
        alpha = 0.4, bins = 10) +
        geom_point(data = all_sites_df, aes(x, y), color = "black", alpha = 0.2, size = 0.5) +
        geom_segment(data = segments, aes(x=x, y=y, xend=xend, yend=yend, color=era),
                     arrow = arrow(length = unit(0.12, "cm")), linewidth = 0.7) +
        scale_color_brewer(palette = "Set1", name = "Climate Era") +
        labs(title = "Migration Tree: Temporal Transitions by Site Count",
             subtitle = paste("Total Sites Reached:", length(tree_nodes)),
             x = "Longitude", y = "Latitude") +
        theme_minimal()
  ggsave("./plot_tree.pdf", height = 8, width = 6)
}

plot_migration_tree(migration_tree, sites)


animate_migration_tree <- function(tree_nodes, all_sites) {
    # 1. Prepare data
    all_sites_df <- as.data.frame(all_sites)
    colnames(all_sites_df) <- c("x", "y")
    node_lookup <- setNames(tree_nodes, sapply(tree_nodes, `[[`, "idx"))

    segments <- do.call(rbind, lapply(tree_nodes, function(node) {
        if (is.null(node$parent_idx)) return(NULL)
        parent <- node_lookup[[as.character(node$parent_idx)]]
        data.frame(
            x = parent$loc[1], y = parent$loc[2], 
            xend = node$loc[1], yend = node$loc[2], 
            depth = node$depth, 
            era = node$era,
            id = paste(parent$idx, node$idx, sep = "_")
        )
    }))

    # 2. Build the ggplot object
    anim_plot <- ggplot() +
        geom_sf(data = africa, fill = "gray95", color = "gray80") +
        geom_point(data = all_sites_df, aes(x, y), color = "black", alpha = 0.1, size = 0.5) +
        # The key for animation: we map 'group' to 'id' to keep segments distinct
        geom_segment(data = segments, 
                     aes(x=x, y=y, xend=xend, yend=yend, color=era, group=id),
                     arrow = arrow(length = unit(0.12, "cm")), 
                     linewidth = 0.8) +
        scale_color_brewer(palette = "Set1", name = "Climate Era") +
        labs(title = "Migration Tree Growth: {closest_state} Steps",
             subtitle = "Era: {current_frame}", # This is a placeholder, custom labels below
             x = "Longitude", y = "Latitude") +
        theme_minimal() +
        # 3. Animation Settings
        transition_states(depth, transition_length = 2, state_length = 1) +
        shadow_mark(past = TRUE, future = FALSE) + # This keeps old arrows on screen
        enter_fade()

    return(anim_plot)
}

# --- Execution and Saving --- #

# Generate the animation object
my_animation <- animate_migration_tree(migration_tree, sites)

# Render the animation (you can adjust fps and duration for speed)
rendered_anim <- animate(my_animation, 
                         nframes = 200, 
                         fps = 5, 
                         width = 600, 
                         height = 800, 
                         renderer = gifski_renderer())

# Save as GIF
anim_save("./migration_tree_animation.gif", animation = rendered_anim)

# --- Simulation Study --- #

run_sim_summary <- function(rho_val) {
  # Run the simulation
  migration_tree <- build_migration_tree(
    S = sites, r0_idx = 9, pred_list = pred_list, 
    site_eras = site_eras, grid_coords = grid_land_df2, 
    dt_g = 1, Q = diag(2), rho = rho_val, p_min = 0.01, K_max = 1000
  )
  
  # Return metrics
  n_nodes <- length(migration_tree)
  return(data.frame(rho = rho_val, nodes = n_nodes, edges = n_nodes - 1))
}


rho_values <- c(1.0, 1.05, 1.1, 1.2)
results_list <- lapply(rho_values, run_sim_summary)
results_df <- do.call(rbind, results_list)

# View the impact of rho on tree size
print(results_df)

n_sims <- 1000

# Create a sequence of rho values to test
rho_values <- 4#seq(2, 4.0, length.out = n_sims)

sim_results <- data.frame(
  Run = 1, #:n_sims, 
  Rho = rho_values,
  Nodes = 0, 
  Max_Depth = 0, 
  Avg_Children = 0
)

for (i in 1:n_sims) {
  # Use the specific rho for this run
  current_rho <- rho_values[i]
  
  current_tree <- build_migration_tree(
    S = sites, r0_idx = 9, pred_list = pred_list, 
    site_eras = site_eras, grid_coords = grid_land_df2, 
    dt_g = 1, Q = diag(2), rho = current_rho, p_min = 0.01, K_max = 1000
  )
    
  # 1. Total Nodes
  sim_results$Nodes[i] <- length(current_tree)
  
  # 2. Max Depth
  sim_results$Max_Depth[i] <- max(sapply(current_tree, function(x) x$depth))
  
  # 3. Branching Factor (Average children per non-leaf node)
  parents <- sapply(current_tree, function(x) x$parent_idx)
  parents <- unlist(parents) # Remove NULLs
  
  if (length(parents) > 0) {
    parent_counts <- table(parents)
    sim_results$Avg_Children[i] <- mean(parent_counts)
  } else {
    sim_results$Avg_Children[i] <- 0
  }
  
  # Progress tracker every 100 runs
  if (i %% 100 == 0) cat("Completed Run:", i, "| Rho:", round(current_rho, 2), "\n")
}

sim_results
