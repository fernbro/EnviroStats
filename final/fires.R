library(tidyverse)
library(sf)

ignitions_raw <- st_read("final/data/National_USFS_Fire_Occurrence_Point_(Feature_Layer).shp")

ignitions <- ignitions_raw %>% 
  filter(FIREYEAR == 2020)

# download state of arizona

az <- tigris::states() %>% 
  filter(NAME == "Arizona")

# st_write(select(az, geometry), "final/data/AZ.shp", append = F)

st_crs(az) == st_crs(ignitions) # true! yay

az_ignitions <- st_intersection(ignitions, az) %>% 
  transmute(date = as.POSIXct(DISCOVERYD),
            cause = STATCAUSE,
            id = UNIQFIREID) %>% 
  unique()

fires <- az_ignitions %>% 
  select(geometry, id) %>% 
  mutate(fire_name = id) %>% 
  select(-id)

st_write(fires, "final/data/AZ_2020_Fires.shp", append = F)

fires <- st_transform(fires, "epsg:32612")

# extract full year time series of SPEI, interpolate temporally, and match
# yes this is a time varying factor but ultimately the time aspect will get removed from the data
# we are just using time to match the drought to the fire location

min(az_ignitions$date);max(az_ignitions$date)

# extracted from GEE:

# spei_raw <- read_csv('final/data/IgnitionPoints_SPEI.csv')

# 2020 mean SPEI raster:
library(terra)
library(maptools)
spei_unproj <- rast('final/data/SPEI1y_2020_AnnualMean.tif')
spei <- project(spei_unproj, "epsg:32612")
spei <- terra::as.data.frame(spei, xy = T)
spei <- as.im(spei)
# spei <- terra::as.data.frame(spei, xy = T)

# spei <- spei_raw %>%
#   mutate(month = as.numeric(str_sub(`system:index`, 5, 6)),
#          day = as.numeric(str_sub(`system:index`, 7, 8)),
#          date = make_date(year = 2020, month = month, day = day)) %>% 
#   filter()

# spei_dates <- spei %>% 
#   select(fire_name, date, spei1y) %>% 
#   group_by(fire_name) %>% 
#   arrange(date) %>% 
#   complete(date=seq(min(date), max(date), by = "1 day")) %>% 
#   ungroup()

# spei_interp <- spei_dates %>% 
#   group_by(fire_name) %>% 
#   arrange(date, .by_group = T) %>% 
#   mutate(
#     spei_int = ifelse(
#       is.na(spei1y), # only interpolate where evi is NA
#       approx(x = date, y = spei1y, xout = date, rule = 1)$y, # rule=1 does not allow extrapolation
#       spei1y
#     )
#   ) %>%
#   ungroup() %>% 
#   select(-spei1y)

# fire_spei <- az_ignitions %>% 
#   mutate(fire_name = id) %>% 
#   inner_join(spei_interp)
# fire_spei <- st_as_sf(fire_spei)
# fire_spei <- st_transform(fire_spei, "epsg:32612")
# 
# st_write(fire_spei, "final/Fires_SPEI.shp", append = T)

library(spatstat)

fire_pp <- as.ppp(fires)

spei_mod <- kppm(unmark(fire_pp) ~ spei, owin = az,
                 clusters = "LGCP")

plot(Fest(fire_pp))









