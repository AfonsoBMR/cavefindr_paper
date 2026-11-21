# DATA PREPARATION -----------------------------

# import las files for prediction
las1 <- lidR::readLAS( here::here("data","LIDAR_tiles","LO-141283-06-2024_v01.laz"))
las2 <- lidR::readLAS( here::here("data","LIDAR_tiles","LO-141284-06-2024_v01.laz"))

# merging both sections into single pointcloud
las <- rbind( las1, las2 )

# declare LAS catalog with all tiles
my_catalog <- lidR::readLAScatalog( here::here("data","LIDAR_tiles"))

# Alecrineiros area in km2
area <- as.numeric( sf::st_area(las) ) / 1000000

# subset corrected rois training data to exclude caves from Alecrineiros
load( here::here( "output","corrected_cave_dataset_20260810122226", "cave_rois_corrected.Rdata" ))
ext_las <- lidR::ext(las)
intersects_las <- vapply(cave_rois_corrected, function(roi) {
  !is.null(terra::intersect( lidR::ext(roi), ext_las))
}, logical(1))
cave_rois_corrected_out <- cave_rois_corrected[!intersects_las]

# attribute different id to subsetted training data, and recover radius
attr( cave_rois_corrected_out, "id" ) <- paste0( attr( cave_rois_corrected, "id" ), "_out" )
attr( cave_rois_corrected_out, "radius" ) <- attr( cave_rois_corrected, "radius" ) 

# subset real detectable caves inside alecrineiros
cave_rois_corrected_in <- cave_rois_corrected[intersects_las]
real_coords <- lapply( cave_rois_corrected_in, function(roi){
  c(mean(roi$X), mean(roi$Y))
} )
xy <- do.call(rbind, real_coords)
colnames(xy) <- c("x", "y")

# create vector from coordinates
real <- terra::vect(xy, type = "points", crs = lidR::wkt(las))

# keep only caves whose centre lies inside the study area (edge ROIs can have centres outside)
real <- terra::crop(real, ext_las)

# drop duplicated entries (closer than 1 m to an earlier one)
d <- as.matrix(terra::distance(real))
d[lower.tri(d, diag = TRUE)] <- Inf
real <- real[colSums(d < 1) == 0, ]

# buffer (hit radius) around real detectable caves 
buffer <- terra::buffer( real, 5 )

# get real (detectable and undetectable) caves inside alecrineiros
my_cave_dataset <- terra::vect(here::here("data","cave_dataset","cave_dataset.shp"))
my_cave_dataset <- terra::project(my_cave_dataset, real)
real_all <- terra::crop(my_cave_dataset, ext_las)

# non-detectable caves: known caves with no detectable cave within tol metres
tol <- 10
d <- terra::distance(real, real_all)
stopifnot(all(apply(d, 1, min) <= tol))   # every detectable cave has its original position
real_nondet <- real_all[apply(d, 2, min) > tol, ]

# all caves: corrected positions for detectable ones, original positions for the rest
real_all <- terra::vect(
  rbind(terra::crds(real), terra::crds(real_nondet)),
  type = "points", crs = terra::crs(real)
)
# corrected positions differ from the original records by 0.7–9.0 m (median 3.6 m)

# buffer (hit radius) around (detectable and undetectable) caves 
buffer_all <- terra::buffer( real_all, 5 )

hit_area <- pi*5^2



# EVALUATING MODEL METHOD --------------------------------

# training model with all study area caves except those from Alecrineiros
cavefindr::train_model( 
  corrected_rois = cave_rois_corrected_out,
  las = my_catalog,
  metrics = c(
    "depth_plane",
    "depth_mean",
    "depth_half_radius_mean",
    "depth_w_mean",
    "slope",
    "intensity",
    "NIR"
  ),
  plot = FALSE,
  var_imp = TRUE
)

# getting candidates using model method
cavefindr::get_predictions(
  method = "model", 
  las = las1,
  model = '20260810122226_out', 
  score_threshold = 0.1,
  isolated_policy = "exclude",
  remove_edges = TRUE,
  remove_waterlines = TRUE,
  inspection = FALSE
)
cavefindr::get_predictions(
  method = "model", 
  las = las2,
  model = '20260810122226_out', 
  score_threshold = 0.1,
  isolated_policy = "exclude",
  remove_edges = TRUE,
  remove_waterlines = TRUE,
  inspection = FALSE
)

# fetching output shp
predictions1 <- terra::vect(here::here("output","prediction_202611120041", "shp", "202611120041.shp"))
predictions2 <- terra::vect(here::here("output","prediction_202611120219", "shp", "202611120219.shp"))

# merging output vectors
predictions <- rbind(predictions1, predictions2)

# decomposing candidates by thresholds
thresholds <- seq( 0.1, 0.99, by = 0.01)
predictions_by_thres <- lapply( thresholds, function(threshold){
  subset(predictions, predictions$score >= threshold)
})



# EVALUATING SIMPLE METHOD --------------------------------------------

# getting candidates using simple method
cavefindr::get_predictions(
  method = "simple", 
  las = las1,
  metric_thresholds = list(depth_mean=0.7),
  radius = 10,
  isolated_policy = "exclude",
  remove_edges = TRUE,
  remove_waterlines = TRUE,
  inspection = FALSE
)
cavefindr::get_predictions(
  method = "simple", 
  las = las2,
  metric_thresholds = list(depth_mean=0.7),
  radius = 10,
  isolated_policy = "exclude",
  remove_edges = TRUE,
  remove_waterlines = TRUE,
  inspection = FALSE
)

# fetching output shp
predictions1_s <- terra::vect(here::here("output","prediction_202611171233", "shp", "202611171233.shp"))
predictions2_s <- terra::vect(here::here("output","prediction_202611171349", "shp", "202611171349.shp"))

# merging output vectors
predictions_s <- rbind(predictions1_s, predictions2_s)

# simplify depth_mean column
predictions_s$metrics <- as.numeric( sub( "depth_mean=", "", predictions_s$metrics))

# thresholds for simple method: every distinct depth value
depths <- sort(unique(predictions_s$metrics))



# FROC FIGURE ---------------------------------------------------

# FROC for one method against one set of cave buffers
# pred: SpatVector of candidates
# scores: score per candidate 
# thr: thresholds to sweep
# buf: buffered caves
# n_caves: number of caves in buf
# marks: thresholds to label with a red point on the curve
# lab: function turning a threshold into its label text
make_froc <- function(pred, scores, thr, buf, n_caves, marks = NULL, lab = format) {
  
  # candidate x cave hit matrix, computed once
  hits      <- terra::relate(pred, buf, "intersects")
  det_hit   <- rowSums(hits) > 0                                                  # candidate hits any cave
  cave_best <- apply(hits, 2, function(h) if (any(h)) max(scores[h]) else -Inf)   # best score hitting each cave
  
  # thresholds to label (only those inside the swept range), added to the sweep
  marks <- marks[marks >= min(thr) & marks <= max(thr)]
  thr   <- sort(unique(c(thr, marks)))
  
  # total detections, false detections, sensitivity/recall for each threshold
  stats <- data.frame(
    thr         = thr,
    n           = vapply(thr, function(t) sum(scores >= t), numeric(1)),
    FP          = vapply(thr, function(t) sum(scores >= t & !det_hit), numeric(1)),
    sensitivity = vapply(thr, function(t) sum(cave_best >= t) / n_caves, numeric(1))
  )
  stats$FP_per_km2 <- stats$FP / area   # false detections per km2 searched
  
  # rows to mark with red points
  marked <- stats[stats$thr %in% marks, ]
  
  # random-placement baseline
  n_rand <- seq(0, max(stats$n), length.out = 200)
  random <- data.frame(
    FP_per_km2  = n_rand / area,
    sensitivity = 1 - (1 - hit_area / (area * 1e6))^n_rand
  )
  
  # adding (0,0) origin point to data
  # the line segment that connects to (0,0) is, strictly speaking, extrapolation, but at the scale of the plot a straight versus curved line cant be told apart
  stats <- rbind(
    data.frame(thr = NA, n = 0, FP = 0, sensitivity = 0, FP_per_km2 = 0),
    stats
  )
  
  ##print(stats)
  
  # plot FROC curve: sensitivity against false detections per km2
  ggplot2::ggplot(stats, ggplot2::aes(x = FP_per_km2, y = sensitivity)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_line(data = random, linetype = "dashed", colour = "grey40", linewidth = 0.8) +
    ggplot2::geom_point(data = marked, colour = "red", size = 1.8) +
    ggplot2::geom_text(data = marked, ggplot2::aes(label = lab(thr)),
                       colour = "red", size = 3.5, hjust = -0.25, vjust = 1.4) +
    ggplot2::scale_x_continuous(limits = c(0, NA)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::labs(
      x = "False detections per km²",
      y = "Sensitivity"
    ) +
    ggplot2::theme_bw(base_size = 13)
}

# thresholds to label on each curve (adjust to your ranges)
mark_s <- c(0.8, 1, 2, 5)             # depths (m)
mark_m <- c(0.18, 0.41, 0.90, 0.99)   # model scores
lab_s  <- function(x) paste0(x, " m")

# four FROCs: (simple / model) x (detectable caves / all caves)
FROC_s     <- make_froc(predictions_s, predictions_s$metrics, depths,     buffer,     nrow(real),     mark_s, lab_s)
FROC       <- make_froc(predictions,   predictions$score,     thresholds, buffer,     nrow(real),     mark_m)
FROC_s_all <- make_froc(predictions_s, predictions_s$metrics, depths,     buffer_all, nrow(real_all), mark_s, lab_s)
FROC_all   <- make_froc(predictions,   predictions$score,     thresholds, buffer_all, nrow(real_all), mark_m)

# keep axis titles only on the outer edges
FROC_s     <- FROC_s     + ggplot2::labs(x = NULL)
FROC       <- FROC       + ggplot2::labs(x = NULL, y = NULL)
FROC_all   <- FROC_all   + ggplot2::labs(y = NULL)

# column and row titles
col_lab <- function(x) cowplot::ggdraw() +
  cowplot::draw_label(x, fontface = "bold", size = 15)
row_lab <- function(x) cowplot::ggdraw() +
  cowplot::draw_label(x, fontface = "bold", size = 15, angle = -90)

# define plot grid
panels <- cowplot::plot_grid(
  FROC_s, FROC, FROC_s_all, FROC_all,
  nrow = 2, align = "hv", axis = "tblr"
) +
  cowplot::draw_plot_label(
    c("a)", "b)", "c)", "d)"),
    x = c(0.105, 0.607, 0.105, 0.607),
    y = c(0.98, 0.98, 0.48, 0.48),
    hjust = 0, vjust = 1, size = 16
  )

header <- cowplot::plot_grid(
  NULL, col_lab("Simple method"), NULL, col_lab("Model method"), NULL,
  nrow = 1, rel_widths = c(0.1, 1, 0.08, 1, 0.1)
)

body <- cowplot::plot_grid(
  panels,
  cowplot::plot_grid(row_lab("Detectable caves"), row_lab("All known caves"), ncol = 1),
  nrow = 1, rel_widths = c(1, 0.08)
)  

fig <- cowplot::plot_grid(header, body, ncol = 1, rel_heights = c(0.06, 1))

# export
ggplot2::ggsave(
  plot = fig,
  filename = here::here("figures", "FROC.png"),
  width = 8.5, height = 8.5,
  bg = "white", dpi = 800
)



# STUDY AREA FIGURE --------------------------------------------

# get dtm
dtm_tin <- lidR::rasterize_terrain( las, res = 2, algorithm = lidR::tin() )

# get shaded dtm
dtm_prod <- terra::terrain(dtm_tin, v = c("slope", "aspect"), unit = "radians")
dtm_hillshade <- terra::shade(slope = dtm_prod$slope, aspect = dtm_prod$aspect)

# plot
hill_df <- as.data.frame(dtm_hillshade, xy = TRUE)
names(hill_df)[3] <- "hs"

to_df <- function(v) {
  xy <- as.data.frame(terra::crds(v)); names(xy) <- c("x", "y"); xy
}
caves_df <- to_df(real)
pred_df  <- to_df(predictions_by_thres[[which(thresholds == 0.9)]])

# get non-detectable caves: no detectable cave within 10 m
nondet_df <- to_df(real_nondet)

e <- as.vector(terra::ext(dtm_hillshade))   # xmin, xmax, ymin, ymax
w <- e[4] - e[3]                            # map width on screen (original y range)
h <- e[2] - e[1]                            # map height on screen (original x range)

x0   <- e[4] + 0.06 * w                     # common left edge of legend, north, scale
ytop <- -e[1]                               # top of the panel
scale_len <- 500                            # scale bar length in map units (m)

p <- ggplot2::ggplot() +
  ggplot2::geom_raster(data = hill_df, ggplot2::aes(x = y, y = -x, fill = hs)) +
  ggplot2::scale_fill_gradientn(colours = gray(0:30/30), guide = "none") +
  ggplot2::geom_point(data = nondet_df, ggplot2::aes(x = y, y = -x),
                      colour = "darkorange", shape = 16, size = 2) +
  ggplot2::geom_point(data = caves_df, ggplot2::aes(x = y, y = -x),
                      colour = "darkgreen", shape = 16, size = 2) +
  ggplot2::geom_point(data = pred_df, ggplot2::aes(x = y, y = -x),
                      colour = "purple", shape = 21, size = 4, stroke = 0.7, fill = NA) +
  
  # legend (manual, left-aligned at x0) 
  ggplot2::annotate("point", x = x0 + 0.012 * w, y = ytop - 0.10 * h,
                    colour = "darkgreen", shape = 16, size = 2) +
  ggplot2::annotate("text",  x = x0 + 0.04 * w,  y = ytop - 0.10 * h,
                    label = "Detectable cave entrance", hjust = 0, size = 5) +
  ggplot2::annotate("point", x = x0 + 0.012 * w, y = ytop - 0.18 * h,
                    colour = "darkorange", shape = 16, size = 2) +
  ggplot2::annotate("text",  x = x0 + 0.04 * w,  y = ytop - 0.18 * h,
                    label = "Non-detectable cave entrance", hjust = 0, size = 5) +
  ggplot2::annotate("point", x = x0 + 0.012 * w, y = ytop - 0.26 * h,
                    colour = "purple", shape = 21, size = 4, stroke = 0.7, fill = NA) +
  ggplot2::annotate("text",  x = x0 + 0.04 * w,  y = ytop - 0.26 * h,
                    label = "Model prediction", hjust = 0, size = 5) +
  
  # north arrow (north = right after rotation) 
  ggplot2::annotate("segment",
                    x = x0, xend = x0 + 0.08 * w,
                    y = ytop - 0.45 * h, yend = ytop - 0.45 * h,
                    arrow = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
                    linewidth = 0.6) +
  ggplot2::annotate("text", x = x0 + 0.10 * w, y = ytop - 0.45 * h,
                    label = "N", fontface = "bold", hjust = 0, size = 5) +
  
  # scale bar (starts at x0) 
  ggplot2::annotate("segment", x = x0, xend = x0 + scale_len,
                    y = ytop - 0.80 * h, yend = ytop - 0.80 * h, linewidth = 0.6) +
  ggplot2::annotate("segment", x = c(x0, x0 + scale_len), xend = c(x0, x0 + scale_len),
                    y = ytop - 0.80 * h - 0.02 * h, yend = ytop - 0.80 * h + 0.02 * h,
                    linewidth = 0.6) +
  ggplot2::annotate("text", x = x0, y = ytop - 0.72 * h,
                    label = paste(scale_len, "m"), hjust = 0, size = 4.5) +
  
  ggplot2::coord_fixed(xlim = c(e[3], e[4] + 0.65 * w),
                       ylim = c(-e[2], -e[1]),
                       expand = FALSE) +
  ggplot2::theme_void()

ggplot2::ggsave(
  plot = p,
  filename = here::here("figures", "validation.png"),
  width = 9, height = 4,
  bg = "white",
  dpi = 800
)

# store session info 
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
