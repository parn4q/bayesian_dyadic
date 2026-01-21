########################################################################
# Plotting the predictions of new data: Potential Surface and Krieging #
########################################################################

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(scales)
library(elevatr) 
library(ggpubr)

# --- Set working Directories if needed --- #

#setwd("E:\\Likun\\Bayesian Dyadic") 
#setwd("/Volumes/Lexar/Likun/Bayesian Dyadic") 

# --- load the betameans for potential surface --- #

betaMeans = read.csv('./betameans_r.csv')

# --- Load Krieging data --- #

etaDF = read.csv("./etaDF.csv")

# --- load the OOS set --- #

grid_land_df2 = read.csv("./for_krieg.csv")


grid_land_df2_copy = grid_land_df2



# --- load the Africa Data --- #


# Disable s2 processing to avoid self-intersection issues
sf_use_s2(FALSE)

# Load Africa map
africa <- ne_countries(continent = "Africa", returnclass = "sf") #shape file of Africa

# ggplot() +
#   geom_sf(data = africa, fill = "lightgray", color = "black") +
#   geom_point(data = grid_land_df2, aes(x = longitude, y = latitude), alpha = 0.1) +
#   ggtitle("Freshwater Bodies in Africa for each Language") +
#   theme_minimal()

##################################################
##### Creating the potential surface    ##########
##################################################

# across ALL datasets you will plot
z_min <- min(, na.rm = TRUE)
z_max <- max(df_all$z, na.rm = TRUE)

breaks <- seq(z_min, z_max, length.out = 10)


# --- Longitude and Latitude --- #

grid_land_df2$pred <- with(grid_land_df2,
  #long_sc * betaMeans$betaMeans_r[1] #+
  #lat_sc * betaMeans$betaMeans_r[2] #
)

# Plot with shape overlay
plot1 = ggplot() +
    geom_sf(data = africa, fill = "lightgray", color = "black") +
    geom_sf(data = africa_lakes, fill = "blue", color = "blue") +
    geom_sf(data = africa_rivers, color = "blue") +
    #geom_point(data = grid_land_df2, aes(x = longitude, y = latitude), color = "black", size = 0.8, alpha = 0.5) +
    geom_contour_filled(data = grid_land_df2, 
        aes(x = longitude, y = latitude, z = pred), 
            alpha = 0.7, bins = 24) +
    labs(title = "Latitude Effect on the Potential Surface",
        x = "Longitude", y = "Latitude",
        fill = "Prediction") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "top")
plot1

ggsave("lat_eff.pdf", width = 10, height = 8)  # wider PDF


# --- Longitude^2 and Latitude^2 and Interaction --- #

grid_land_df2$pred <- with(grid_land_df2,
  #(long_sc^2) * betaMeans$betaMeans_r[5] #+
   #(lat_sc^2) * betaMeans$betaMeans_r[6] #+
  (long_sc * lat_sc) * betaMeans$betaMeans_r[7]
)

# Plot with shape overlay
plot1 = ggplot() +
    geom_sf(data = africa, fill = "lightgray", color = "black") +
    geom_sf(data = africa_lakes, fill = "blue", color = "blue") +
    geom_sf(data = africa_rivers, color = "blue") +
    #geom_point(data = grid_land_df2, aes(x = longitude, y = latitude), color = "black", size = 0.8, alpha = 0.5) +
    geom_contour_filled(data = grid_land_df2, 
        aes(x = longitude, y = latitude, z = pred), 
            alpha = 0.7, bins = 24) +
    labs(title = "Interaction Between Longitude and Latitude Effect on the Potential Surface",
        x = "Longitude", y = "Latitude",
        fill = "Prediction") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "top")
plot1

ggsave("lon_lat_int_eff.pdf", width = 10, height = 8)

# --- Elevation and Distance to Fresh Water --- #

grid_land_df2$pred <- with(grid_land_df2,
   #elevation * betaMeans$betaMeans_r[3] #+
   euclid_vec_sc2 * betaMeans$betaMeans_r[4] #+
)

# Plot with shape overlay
plot1 = ggplot() +
    geom_sf(data = africa, fill = "lightgray", color = "black") +
    geom_sf(data = africa_lakes, fill = "blue", color = "blue") +
    geom_sf(data = africa_rivers, color = "blue") +
    #geom_point(data = grid_land_df2, aes(x = longitude, y = latitude), color = "black", size = 0.8, alpha = 0.5) +
    geom_contour_filled(data = grid_land_df2, 
        aes(x = longitude, y = latitude, z = pred), 
            alpha = 0.7, bins = 24) +
    labs(title = "Distance to Fresh Water Effect on the Potential Surface",
        x = "Longitude", y = "Latitude",
        fill = "Prediction") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "top")
plot1

ggsave("diff_to_fresh_eff.pdf", width = 10, height = 8)

# --- Geo Regions --- #

grid_land_df2_copy$pred1 <- with(grid_land_df2,
  cat_5000 #* betaMeans$betaMeans_r[8]
)

grid_land_df2_copy$pred2 <- with(grid_land_df2,
  cat_4_5000 #* betaMeans$betaMeans_r[8]
)

grid_land_df2_copy$pred3 <- with(grid_land_df2,
  cat_3200 #* betaMeans$betaMeans_r[8]
)

###
### Kriging
###


plot2 = ggplot() +
    geom_sf(data = africa, fill = "lightgray", color = "black") +
    geom_sf(data = africa_lakes, fill = "blue", color = "blue") +
    geom_sf(data = africa_rivers, color = "blue") +
    geom_point(data = grid_land_df2, aes(x = longitude, y = latitude), color = "black", size = 0.8, alpha = 0.6) +
    geom_contour_filled(data = grid_land_df2, aes(x = longitude, y = latitude, z = etaDF$z), alpha = 0.7, bins = 10) +
    labs(title = expression(eta ~ " Surface") ,
        x = "Longitude", y = "Latitude",
        fill = "Prediction") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

plot2


# --- Standard Deviation --- #

Epi = etaDF$z + grid_land_df2$pred

togdf = data.frame(x = grid_land_df2$long_sc,
                   y = grid_land_df2$lat_sc,
                   z = Epi)


plot3 = ggplot() +
    geom_sf(data = africa, fill = "lightgray", color = "black") +
    geom_sf(data = africa_lakes, fill = "blue", color = "blue") +
    geom_sf(data = africa_rivers, color = "blue") +
    #geom_point(data = grid_land_df2, aes(x = longitude, y = latitude), color = "black", size = 0.8, alpha = 0.6) +
    geom_contour_filled(
            data = grid_land_df2, 
            aes(x = longitude, y = latitude, z = togdf$z), 
            alpha = 0.7, bins = 24) +
    labs(x = "", y = "", title = "Potential Surface (All Times)", fill = "ρ") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) 

plot3

ggsave("./ps_pot_alltime.png", plot3, width = 8, height = 6)


# --- Create the PS plots for 5000, 4200-5000, 3200 BP --- #


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


p1 = ggplot() +
    geom_sf(data = africa, fill = "lightgray", color = "black") +
    geom_sf(data = africa_lakes, fill = "blue", color = "blue") +
    geom_sf(data = africa_rivers, color = "blue") +
    #geom_point(data = grid_land_df2, aes(x = longitude, y = latitude), color = "black", size = 0.8, alpha = 0.5) +
    geom_contour_filled(
    data = grid_land_df2, 
            aes(x = longitude, y = latitude, z =Epi1), 
                                            alpha = 0.7, bins = 24) +
            labs(title = "Potential Surface Plot 5000BP+",
        x = "Longitude", y = "Latitude",
                                                    fill = "Prediction") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

grid_land_df2$Epi2 = etaDF$z + grid_land_df2$pred2


p2 = ggplot() +
geom_sf(data = africa, fill = "lightgray", color = "black") +
geom_sf(data = africa_lakes, fill = "blue", color = "blue") +
geom_sf(data = africa_rivers, color = "blue") +
#geom_point(data = grid_land_df2, aes(x = longitude, y = latitude), color = "black", size = 0.8, alpha = 0.5) +
geom_contour_filled(
data = grid_land_df2, 
aes(x = longitude, y = latitude, z = Epi2), 
alpha = 0.7, bins = 24) +
labs(title = "Potential Surface Plot 4-5000BP",
x = "Longitude", y = "Latitude",
fill = "Prediction") +
theme_minimal() +
theme(plot.title = element_text(hjust = 0.5))

grid_land_df2$Epi3 = etaDF$z + grid_land_df2$pred3


p3 = ggplot() +
geom_sf(data = africa, fill = "lightgray", color = "black") +
geom_sf(data = africa_lakes, fill = "blue", color = "blue") +
geom_sf(data = africa_rivers, color = "blue") +
#geom_point(data = grid_land_df2, aes(x = longitude, y = latitude), color = "black", size = 0.8, alpha = 0.5) +
geom_contour_filled(
data = grid_land_df2, 
aes(x = longitude, y = latitude, z = Epi3), 
alpha = 0.7, bins = 24) +
labs(title = "Potential Surface Plot 3200BP",
x = "Longitude", y = "Latitude",
fill = "Prediction") +
theme_minimal() +
theme(plot.title = element_text(hjust = 0.5))


ggarrange(p1,p2,p3, nrow = 1, common.legend = T)

ggsave("./overall_ps.pdf",  width = 10, height = 6)






# --- Create Isolated Geographic Region Covariate plots --- #

# 5000BP

grid_land_df2_copy$Epi1 = grid_land_df2_copy$pred1 #+ etaDF$z


p1 = ggplot() +
    geom_sf(data = africa, fill = "lightgray", color = "black") +
    geom_sf(data = africa_lakes, fill = "blue", color = "blue") +
    geom_sf(data = africa_rivers, color = "blue") +
    #geom_point(data = grid_land_df2, aes(x = longitude, y = latitude), color = "black", size = 0.8, alpha = 0.5) +
    geom_contour_filled(
            data = grid_land_df2_copy, 
            aes(x = longitude, y = latitude, z =Epi1), 
            alpha = 0.7, bins = 24) +
    labs(title = "Isolated Geographic Region Covariate Potential Surface Plot 5000BP+",
        x = "Longitude", y = "Latitude",
        fill = "Prediction") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 7))


grid_land_df2_copy$Epi2 = grid_land_df2_copy$pred2 # + etaDF$z 


p2 = ggplot() +
    geom_sf(data = africa, fill = "lightgray", color = "black") +
    geom_sf(data = africa_lakes, fill = "blue", color = "blue") +
    geom_sf(data = africa_rivers, color = "blue") +
    #geom_point(data = grid_land_df2, aes(x = longitude, y = latitude), color = "black", size = 0.8, alpha = 0.5) +
    geom_contour_filled(
        data = grid_land_df2_copy, 
        aes(x = longitude, y = latitude, z = Epi2), 
        alpha = 0.7, bins = 24) +
    labs(title = "Isolated Geographic Region Covariate Potential Surface Plot 4-5000BP",
        x = "Longitude", y = "Latitude",
        fill = "Prediction") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 7))

grid_land_df2_copy$Epi3 = grid_land_df2_copy$pred3 #+ etaDF$z


p3 = ggplot() +
    geom_sf(data = africa, fill = "lightgray", color = "black") +
    geom_sf(data = africa_lakes, fill = "blue", color = "blue") +
    geom_sf(data = africa_rivers, color = "blue") +
    #geom_point(data = grid_land_df2, aes(x = longitude, y = latitude), color = "black", size = 0.8, alpha = 0.5) +
    geom_contour_filled(
        data = grid_land_df2_copy, 
        aes(x = longitude, y = latitude, z = Epi3), 
        alpha = 0.7, bins = 24) +
    labs(title = "Isolated Geographic Region Covariate Potential Surface Plot 3200BP",
        x = "Longitude", y = "Latitude",
        fill = "Prediction") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 7))


ggarrange(p1,p2,p3, nrow = 1, common.legend = T)


ggsave("./iso_geo_region.pdf",  width = 10, height = 6)
