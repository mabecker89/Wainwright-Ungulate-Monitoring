# ---- Packages ----
library(terra)
library(sf)
library(dplyr)
library(spsurvey)
library(units)
library(mapview)
library(stars)
library(lwgeom)

set.seed(2024)

# ---- Config ----
RASTER_PATH <- "G:/Shared drives/ABMI Camera Mammals/Projects/Wain/lc_clip_recl1.tif"
AGG_FACTOR  <- 1    # 1=10m; try 2 or 3 if polygonizing is heavy
TOTAL_SITES <- 100
MINDIS_M    <- 400  # set to 0 to disable spacing

# ---- 1) Read raster (must be projected in metres) ----
lc <- rast(RASTER_PATH)
if (is.lonlat(lc)) stop("Raster is in lon/lat. Reproject to a metric CRS first.")

# Optional: coarsen to speed polygonization
if (AGG_FACTOR > 1) {
  lc <- aggregate(lc, fact = AGG_FACTOR, fun = modal, na.rm = TRUE)
}

# ---- 2) Raster -> polygons (dissolve contiguous cells by value) ----
# SpatVector polygons with an attribute equal to raster values
polys_v <- as.polygons(lc, dissolve = TRUE, values = TRUE, na.rm = TRUE)
polys_sf <- st_as_sf(polys_v)

# Normalize value column name to 'gridcode'
val_col <- setdiff(names(polys_sf), "geometry")[1]
polys_sf <- polys_sf |>
  rename(gridcode = all_of(val_col)) |>
  st_make_valid() |>
  st_set_precision(1) |>
  st_snap_to_grid(1) |>
  filter(!st_is_empty(geometry))

# ---- 3) Map classes; drop Water; keep target strata ----
# Codes: 1 Water | 5 Exposed/Barren | 6 Shrubland | 7 Wetland | 8 Grassland
#        9 Coniferous | 10 Broadleaf | 11 Mixedwood | 12 Agriculture | 13 Developed | 14 Burned
frame0 <- polys_sf |>
  mutate(
    group_cat = case_when(
      gridcode == 1                ~ "Water",
      gridcode == 6                ~ "Shrubland",
      gridcode == 7                ~ "Wetland",
      gridcode == 8                ~ "Grassland",
      gridcode %in% c(9,10,11)     ~ "Forested",
      gridcode %in% c(5,12,13,14)  ~ "Other",
      TRUE                         ~ "Other"
    )
  ) |>
  filter(group_cat != "Water")

targets <- c("Shrubland","Wetland","Grassland","Forested")
frame_strata <- frame0 |>
  filter(group_cat %in% targets) |>
  mutate(area_m2 = as.numeric(set_units(st_area(geometry), "m^2")))

if (nrow(frame_strata) == 0) stop("No target strata present after dropping Water.")

# --- Up-weight Wetland to 8 and Forested to 20 (of 100) ----------------------
TOTAL <- 100L
WETLAND_TARGET   <- 8L
GRASSLAND_TARGET  <- 40L

# Area per stratum (you already have frame_strata with group_cat)
area_by_stratum <- frame_strata |>
  st_drop_geometry() |>
  dplyr::group_by(group_cat) |>
  dplyr::summarise(area_m2 = sum(area_m2, na.rm = TRUE), .groups = "drop")

# Sanity checks
needed <- c("Shrubland","Grassland","Forested","Wetland")
stopifnot(all(needed %in% area_by_stratum$group_cat))

# Fixed counts for Wetland & Forested
fixed_vec <- c("Wetland" = WETLAND_TARGET, "Grassland" = GRASSLAND_TARGET)

# Distribute the remainder across Shrubland + Grassland proportionally by area
other_strata <- c("Shrubland","Forested")
rem_sites <- TOTAL - sum(fixed_vec) 

props_other <- area_by_stratum |>
  dplyr::filter(group_cat %in% other_strata) |>
  dplyr::mutate(prop = area_m2 / sum(area_m2))

alloc_raw   <- props_other$prop * rem_sites
alloc_floor <- floor(alloc_raw)
remainder   <- rem_sites - sum(alloc_floor)

alloc <- alloc_floor
if (remainder > 0) {
  bump <- order(alloc_raw - alloc_floor, decreasing = TRUE)[seq_len(remainder)]
  alloc[bump] <- alloc[bump] + 1
}

# Build the final named vector in a stable order
n_base_vec <- setNames(integer(0), character(0))
n_base_vec[other_strata] <- alloc[match(other_strata, props_other$group_cat)]
n_base_vec["Grassland"]   <- GRASSLAND_TARGET
n_base_vec["Wetland"]    <- WETLAND_TARGET

# Optional: ensure totals sum exactly to 100
stopifnot(sum(n_base_vec) == TOTAL)

print(n_base_vec)
# Example: Shrubland=?? Grassland=?? Forested=20 Wetland=8

# ---- 5) GRTS draw (stratified; spacing; oversample) ----
design <- grts(
  sframe      = frame_strata,     # sf POLYGON frame
  n_base      = n_base_vec,       # named per-stratum targets
  stratum_var = "group_cat",
  seltype     = "equal",          # equal within stratum; totals pre-allocated
  mindis      = 400,         # set to 0 if you don't want spacing
  maxtry      = 10,
  n_over      = 100              # replacements available in `design$sites_over`
  #n_near    = 2                 # optional nearest-neighbor backups
)

sites_base <- design$sites_base   # your 100 base sites
sites_over <- design$sites_over   # backups if you need to swap later

# ---- 6) Quick map — show the raster (fast), not the polygons ----
# Build a categorical palette for the raster values you care about
# (Water will display too unless you reclass or mask the raster; it's only excluded from sampling)
landcover_cols <- c(
  "1"  = "#7FA6E5",  # Water (only for display)
  "6"  = "#A67C52",  # Shrubland (tan/brown)
  "7"  = "#556B2F",  # Wetland (olive)
  "8"  = "#DAA520",  # Grassland (golden)
  "9"  = "#006400",  # Coniferous -> Forested (dark green)
  "10" = "#228B22",  # Broadleaf   -> Forested (forest green)
  "11" = "#2E8B57",  # Mixedwood   -> Forested (sea green)
  "5"  = "#B38B6D",  # Exposed/Barren -> Other (beige/brown)
  "12" = "#C2B280",  # Agriculture -> Other (khaki)
  "13" = "#B0A199",  # Developed   -> Other (grey-brown)
  "14" = "#8B5A2B"   # Burned      -> Other (brown)
)

# Convert terra raster to stars for mapview
Landcover <- st_as_stars(lc)

# Map (set maxpixels high so you see everything)
mapview(Landcover, col.regions = landcover_cols, maxpixels = ncell(lc)) +
  mapview(sites_base, zcol = NULL, col.regions = "red", cex = 4, alpha = 1)

# Manual removal

sites_to_remove <- c(
  "Site-093", "Site-027",
  "Site-039", "Site-044",
  "Site-016", "Site-047",
  "Site-084", "Site-092",
  "Site-030", "Site-002",
  "Site-100", "Site-065",
  "Site-057", "Site-006",
  "Site-055"
)

unsuitable <- c(
  "Site-103", "Site-106"
  
)

# --- settings ---
MIN_PAIR_M <- 400  # set to 0 to ignore spacing during replacement

# --- prep: split keep/remove, set precision for stable distances -------------
sites_base <- st_set_precision(sites_base, 1) |> st_snap_to_grid(1)

# HARD EXCLUDE unsuitable oversamples up front
sites_over <- st_set_precision(sites_over, 1) |> st_snap_to_grid(1) |>
  dplyr::filter(!(siteID %in% unsuitable))

to_drop <- tibble(siteID = sites_to_remove)
base_keep <- anti_join(sites_base, to_drop, by = "siteID")
base_drop <- semi_join(sites_base, to_drop, by = "siteID")

# helper for spacing check against already kept
ok_spacing <- function(candidate, kept, min_pair) {
  if (min_pair <= 0 || nrow(kept) == 0) return(TRUE)
  d <- as.numeric(st_distance(candidate, kept))
  all(is.infinite(d) | d >= min_pair)
}

# we’ll build up a final set: start with the kept base sites
final_sites <- base_keep
used_over_ids <- character(0)

# iterate each removed site; choose a replacement
for (i in seq_len(nrow(base_drop))) {
  s <- base_drop[i, , drop = FALSE]
  
  # 1) primary pool: oversamples explicitly linked to this base via `replsite`
  pool <- sites_over |>
    dplyr::filter(replsite == s$siteID,
                  !(siteID %in% used_over_ids))
  
  # 2) fallback: same-stratum oversamples not yet used
  if (nrow(pool) == 0 && "stratum" %in% names(sites_over) && "stratum" %in% names(s)) {
    pool <- sites_over |>
      dplyr::filter(stratum == s$stratum,
                    !(siteID %in% used_over_ids))
  }
  
  # pick the first that satisfies spacing; if none do, take the farthest from current set
  picked <- NULL
  if (nrow(pool) > 0) {
    # try to find the first spacing-valid candidate (keep current order)
    for (j in seq_len(nrow(pool))) {
      cand <- pool[j, , drop = FALSE]
      if (ok_spacing(cand, final_sites, MIN_PAIR_M)) { picked <- cand; break }
    }
    # if spacing failed for all, choose the one maximizing min distance to final set
    if (is.null(picked)) {
      dmat <- st_distance(pool, final_sites)
      # min distance from each candidate to the current set
      mins <- apply(dmat, 1, function(x) min(as.numeric(x)))
      picked <- pool[which.max(mins), , drop = FALSE]
      message(sprintf(
        "Spacing relaxed for %s -> picked oversample %s with min distance %.1f m",
        s$siteID, picked$siteID, max(mins)))
    }
  }
  
  if (!is.null(picked)) {
    picked$siteuse <- "over" # mark as replacement
    final_sites <- bind_rows(final_sites, picked)
    used_over_ids <- c(used_over_ids, picked$siteID)
  } else {
    warning(sprintf(
      "No available oversample for %s (stratum %s). Keeping original base site.",
      s$siteID, ifelse("stratum" %in% names(s), s$stratum, "NA")))
    final_sites <- bind_rows(final_sites, s) # keep the base if no replacement exists
  }
}

# ensure we have exactly the same count as original base set
if (nrow(final_sites) != nrow(sites_base)) {
  warning(sprintf("Final count %d differs from original base %d.", nrow(final_sites), nrow(sites_base)))
}

# quick diagnostics ------------------------------------------------------------
if (nrow(final_sites) > 1) {
  D <- as.matrix(st_distance(final_sites))
  diag(D) <- Inf
  cat(sprintf("Min inter-site distance: %.1f m\n", min(D)))
}
cat("By stratum (final):\n")
print(table(final_sites$stratum))

# optional: a tidy subset of columns for export
sites_final <- final_sites %>%
  mutate(replaced = ifelse(siteID %in% used_over_ids, TRUE, FALSE)) %>%
  dplyr::select(siteID, siteuse, replsite, stratum, wgt, gridcode, group_cat, geometry)

# Map (set maxpixels high so you see everything)
mapview(Landcover, col.regions = landcover_cols, maxpixels = ncell(lc)) +
  mapview(sites_final, zcol = NULL, col.regions = "red", cex = 4, alpha = 1) +
  mapview(sites_base, zcol = NULL, col.regions = "yellow", cex = 4, alpha = 1)

# ---- 7) Save outputs ----
st_write(sites_final, "Data/Wainwright_GRTS_sites_base.shp",
         append = FALSE)
st_write(sites_over, "Data/Wainwright_GRTS_sites_over.shp")
