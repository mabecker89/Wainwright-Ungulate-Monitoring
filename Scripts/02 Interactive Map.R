# ----- Fix CRS definition -> Reproject -> Leaflet map -------------------------
library(terra)
library(raster)   # for addRasterImage()
library(sf)
library(leaflet)
library(leaflet.extras)

# 0) Paths
RASTER_PATH   <- "G:/Shared drives/ABMI Camera Mammals/Projects/Wain/lc_clip_recl1.tif"
SITES_PATH    <- "Data/Wainwright_GRTS_sites_base.shp"  # trusted vector for CRS

# 1) Read data
lc <- rast(RASTER_PATH)
sites <- st_read(SITES_PATH)

# 2) Ensure raster has a FULL, CORRECT CRS (WKT) before any reprojection
#    Borrow WKT from trusted vector (same project CRS)
crs_vec <- st_crs(sites)                  # full WKT from shapefile .prj
if (is.na(crs(lc)) || grepl("^NAD_1983_Transverse_Mercator$", crs(lc))) {
  crs(lc) <- crs_vec$wkt                  # assign explicit WKT (no reprojection yet)
}


# 3) Reproject raster to WGS84 with nearest neighbor (keeps categories intact)
lc_ll <- terra::project(lc, "EPSG:4326", method = "near")

# 4) (Optional) drop Water from display
lc_ll[lc_ll == 1] <- NA

# 5) Convert to RasterLayer for leaflet::addRasterImage()
lc_r_ll <- raster::raster(lc_ll)

# 6) Prepare palette only for PRESENT values
landcover_cols <- c(
  "1"  = "#7FA6E5",  # Water
  "6"  = "#A67C52",  # Shrubland
  "7"  = "#556B2F",  # Wetland
  "8"  = "#DAA520",  # Grassland
  "9"  = "#006400",  # Coniferous
  "10" = "#228B22",  # Broadleaf
  "11" = "#2E8B57",  # Mixedwood
  "5"  = "#B38B6D",  # Exposed/Barren
  "12" = "#C2B280",  # Agriculture
  "13" = "#B0A199",  # Developed
  "14" = "#8B5A2B"   # Burned
)
landcover_labels <- c(
  "1"="Water","6"="Shrubland","7"="Wetland","8"="Grassland",
  "9"="Coniferous","10"="Broadleaf","11"="Mixedwood",
  "5"="Exposed/Barren","12"="Agriculture","13"="Developed","14"="Burned"
)
vals_present      <- sort(unique(na.omit(raster::values(lc_r_ll))))
vals_present_chr  <- as.character(vals_present)
pal <- colorFactor(unname(landcover_cols[vals_present_chr]), domain = vals_present)

# 7) Camera sites to WGS84
sites_ll <- st_transform(sites, 4326)

# 8) Camera icon
cam <- makeAwesomeIcon(icon = "camera", iconColor = "black",
                       library = "ion", markerColor = "white")

# 9) Leaflet
map <- leaflet(options = leafletOptions(preferCanvas = TRUE)) |>
  addProviderTiles(providers$Esri.WorldTopoMap, group = "Topo") |>
  addProviderTiles(providers$Esri.WorldImagery, group = "Imagery") |>
  addRasterImage(lc_r_ll, colors = pal, opacity = 0.8, project = FALSE, group = "Land cover") |>
  addAwesomeMarkers(data = sites_ll, icon = cam, popup = paste0(sites_ll$siteID, "<br>", sites_ll$stratum), group = "Camera sites") |>
  addLegend(position = "bottomright", title = "Land cover",
            colors = unname(landcover_cols[vals_present_chr]),
            labels =  landcover_labels[vals_present_chr],
            opacity = 0.8) |>
  addLayersControl(baseGroups    = c("Topo", "Imagery"),
                   overlayGroups = c("Land cover", "Camera sites"),
                   options       = layersControlOptions(collapsed = FALSE)) |>
  #addMouseCoordinates() |>
  addMeasure(primaryLengthUnit = "meters", secondaryLengthUnit = "kilometers")

map

htmlwidgets::saveWidget(map, "docs/Camera_Locations.html",
                        selfcontained = TRUE)


