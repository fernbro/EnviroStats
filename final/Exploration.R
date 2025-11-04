library(tidyverse)
library(sf)

pipo <- read_csv("final/PIPO_GBIF_10312025.csv") %>% 
  filter(!is.na(decimalLongitude))

pipo_sf <- st_as_sf(pipo, coords = c("decimalLongitude", "decimalLatitude"),
                    crs = "proj+latlong")

# clip to North America bbox

bbox_coords <- c(18.65, -140.66, 55.58, -50.31)
names(bbox_coords) = c("ymin","xmin","ymax","xmax")
bbox_na <- st_as_sfc(st_bbox(bbox_coords))

pipo_na <- pipo_sf[st_intersection(pipo_sf, bbox_na),]

ggplot()+
  geom_sf(data = pipo_na, alpha = 0.3, size = 0.5)+
  theme_minimal(base_size = 26)+
  labs(title = "Ponderosa pine occurrences in North America")
