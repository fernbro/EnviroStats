# STAT 574E Fall 2025
# Independent Research Project
# Fern Bromley

# Multivariate spatial analysis for plant communities

library(adespatial)
library(ade4)
library(tidyverse)
library(ggspatial)
library(vegan)
library(spdep)
library(spmodel)
library(sf)

# read in sampling locations
# references as points but are actually 1x1 m2 plots

plot_xy_raw <- st_read("irp/MuseumFirePlots.shp") %>% 
  transmute(plot = as.numeric(title), 
            severity = case_when(severity == "High" ~ "H",
                                 severity == "Low" ~ "L",
                                 severity == "Unburned" ~ "U"), 
            geometry = geometry)
  

plot_xy_utm <- st_transform(plot_xy_raw, crs = "epsg:32612")
plot_xy <- st_zm(plot_xy_utm, drop = T)


# quickly visualize spatial distribution of plots:

ggplot()+
  geom_sf(data = plot_xy, aes(color = severity))+
  theme_light()+
  annotation_scale()

# now read in the data associated with these points

com_data <- read_csv("irp/CoverSep2024.csv")
com_data[is.na(com_data)] <- 0
plot_xy <- filter(plot_xy, plot %in% com_data$plot) # some plots were damaged and dont have data

# let's first do some ordination so that we can look at this scores across space

# we wisconsin standardize (double standardized; row then column) data; common for plant ecology

com <- com_data %>%
  column_to_rownames(var = "plot") %>%  
  select(ARLU:VETH) %>% # only select species columns for the standardization
  wisconsin()

# define distance-based neighbors:
plot_nbs_dist <- as.matrix(st_is_within_distance(plot_xy, plot_xy, 100)) # 100 meters
plotw <- mat2listw(plot_nbs_dist, style = "W")

# NMDS

com_ord <- metaMDS(com, distance = "bray", k = 2, trymax = 50)

# check stress (measure of how well the multivar dissimilarity fits with ordination distance):
com_ord$stress < 0.2
stressplot(com_ord)


# extract NMDS scores for each plot, assign to severities
com_nmds <- as_tibble(scores(com_ord, display = "sites"), 
                      rownames = "sites") %>% 
  transmute(plot = sites, NMDS1 = NMDS1, NMDS2 = NMDS2) %>% 
  mutate(severity = case_when(plot <= 20 ~ "U",
                              plot > 20 & plot <= 40 ~ "L", 
                              plot > 40 ~ "H")) %>% 
  mutate(plot = as.numeric(plot))

plot_nmds <- full_join(select(com_nmds, plot, NMDS1, NMDS2), plot_xy)

ggplot()+
  geom_sf(data = plot_nmds, size = 3,
          aes(geometry = geometry, shape = severity,
              fill = NMDS1))+
  scale_fill_gradient2(low="darkgreen", high = "magenta")+
  scale_shape_manual(values = c(21, 23, 24))+
  theme_minimal()+
  ggtitle("NMDS Axis 1")

ggplot()+
  geom_sf(data = plot_nmds, size = 3,
          aes(geometry = geometry, shape = severity,
              fill = NMDS2))+
  scale_fill_gradient2(low="darkgreen", high = "magenta")+
  scale_shape_manual(values = c(21, 23, 24))+
  theme_minimal()+
  ggtitle("NMDS Axis 2")

# # plot NMDS scores for each monitoring plot, with 95% "CI" ellipses
ggplot(com_nmds, aes(x = NMDS1, y = NMDS2))+
  geom_point(aes(color = severity))+
  stat_ellipse(aes(color = severity))+
  theme_light()


# to a PCA now:

# turn the community data to a matrix object
com_mat <- as.matrix(com)

# run a PCA
# dudi.pca from the ade4 package
# dudi stands for duality diagram, just another name for a "reduced space ordination"
com_pca <- dudi.pca(com_mat, scale = FALSE, scannf = FALSE, nf = 2)
summary(com_pca)
scores <- com_pca$li
scores$plot <- com_data$plot
plot_scores <- full_join(scores, plot_xy)

ggplot()+
  geom_sf(data = plot_scores, size = 3, 
          aes(geometry = geometry, shape = severity,
              fill = Axis1))+
  scale_fill_gradient2(low="darkgreen", high = "magenta")+
  scale_shape_manual(values = c(21, 23, 24))+
  theme_minimal()+
  ggtitle("PCA Axis 1")

ggplot()+
  geom_sf(data = plot_scores, size = 3, 
          aes(geometry = geometry, shape = severity,
              fill = Axis2))+
  scale_fill_gradient2(low="darkgreen", high = "magenta")+
  scale_shape_manual(values = c(21, 23, 24))+
  theme_minimal()+
  ggtitle("PCA Axis 2")

# multispati: (PCA x Moran's I)

multi <- multispati(com_pca, plotw, scannf = F)
summary(multi)

plots <- rownames(multi$li)
multi$li$plot <- as.numeric(plots)
multi_scores <- full_join(multi$li, plot_xy)

ggplot()+
  geom_sf(data = multi_scores, size = 3,
          aes(geometry = geometry, shape = severity,
              fill = CS1))+
  scale_fill_gradient2(low="darkgreen", high = "magenta")+
  scale_shape_manual(values = c(21, 23, 24))+
  theme_minimal()+
  ggtitle("CS1")

ggplot()+
  geom_sf(data = multi_scores, size = 3,
          aes(geometry = geometry, shape = severity,
              fill = CS2))+
  scale_fill_gradient2(low="darkgreen", high = "magenta")+
  scale_shape_manual(values = c(21, 23, 24))+
  theme_minimal()+
  ggtitle("CS2")

# comparative plots for Axis1 (PCA) and CS1 (MULTISPATI-PCA):

color_breaks <- c(-1, 0, 1)
colors <- c("darkgreen", "lightgreen", "white", "pink", "magenta")

ggplot()+
  geom_sf(data = plot_scores, size = 3, alpha = 0.9,
          aes(geometry = geometry, shape = severity,
              fill = Axis1))+
  scale_fill_gradientn(
    limits  = c(-1, 1),
    colours = colors,
    values  = c(0, scales::rescale(color_breaks, from = range(multi_scores$CS1)), 1),
  )+
  scale_shape_manual(values = c(21, 23, 24))+
  theme_minimal(base_size = 20)+
  ggtitle("PCA Axis 1")

ggplot()+
  geom_sf(data = multi_scores, size = 3, alpha = 0.9,
          aes(geometry = geometry, shape = severity,
              fill = CS1))+
  scale_fill_gradientn(
    limits  = c(-1, 1),
    colours = colors,
    values  = c(0, scales::rescale(color_breaks, from = range(multi_scores$CS1)), 1),
  )+
  scale_shape_manual(values = c(21, 23, 24))+
  theme_minimal(base_size = 20)+
  ggtitle("CS2")

# # picking up, let's merge our plot locations to the 2024 community data

# join the community data with just the plot and geometry columns from that shapefile, since it has a bunch of info i don't care about
com_sf <- inner_join(com_data, plot_xy)

esv_test <- spmodel::esv(MUWR~severity, com_sf)

# let's look at a multivariate empirical variogram..

# just need the x, y locations of the sample locations
plot_loc <- st_zm(com_sf$geometry, drop = T) %>% 
  st_coordinates() %>% 
  data.frame()

# first we'll remove any linear effect of fire severity
# and order the factors too..

com_data$severity <- factor(com_data$severity, levels = c("U", "L", "H"))
com_resid <- data.frame(residuals(lm(com_mat ~ com_data$severity)))


com_resid$xcoord <-  plot_loc$X ; com_resid$ycoord <-  plot_loc$Y
com_resid$severity <- com_data$severity

esv_test <- esv(MUST~severity, data = com_resid, xcoord = "xcoord")
esv_test <- esv(PSMA~severity, data = com_resid, xcoord = "xcoord", bins = 10) ; plot(esv_test)

# using permutations, we can test the significance of Moran's I using our different (univariate and multivariate) community composition metrics

morI_pca <- moran.randtest(com_pca$li, plotw, alter = "two-sided")
morI_pca

nmds_scores <- select(com_nmds, NMDS1, NMDS2)
morI_nmds <- moran.randtest(nmds_scores, plotw, alter = "two-sided")
morI_nmds

morI_multi <- moran.randtest(com_mat, plotw, alter = "two-sided")
morI_multi # separately tests every variable in the community
morI_multi_np <- moranNP.randtest(com_mat, plotw, alter = "two-sided") # separately tests positive & negitive eigenvectors
morI_multi_np
summary(morI_multi_np)

# other exploration:

# compute the multivariate variogram

varg <- variogmultiv(com_resid, plot_loc, nclass = 20)
# the multivariate variogram is equal to the sum of univariate variograms for all the variables
plot(varg$d, varg$var, pch = 20)





# varg$d is the centers of distance classes, 
# varg$var is the empirical semivariance

# compare to a semivariogram of ordination scores

com_nmds_sf <- inner_join(com_nmds, com_sf) %>% 
  st_as_sf()

# univariate ESV with NMDS scores:
esv1 <- esv(NMDS1 ~ severity, data = com_nmds_sf, bins = 20) 
esv2 <- esv(NMDS2 ~ severity, data = com_nmds_sf, bins = 20)
plot(esv1)
plot(esv2)

# univariate ESV with PCA scores:
plot_scores <- st_as_sf(plot_scores)
esv3 <- esv(Axis1 ~ severity, data = plot_scores, bins = 20)
esv4 <- esv(Axis2 ~ severity, data = plot_scores, bins = 20)
plot(esv3)
plot(esv4)

