#!/usr/bin/env/Rscript --vanilla
#>------------------------------------<
#
# Script name: results_to_raster
#
# Author: Manuel R. Popp
# Email: manuel.popp@wsl.ch
#
# Date Created: 2024-02-28
#
# ---------------------------
#
# Description: Summarise the metric values calculated from multiple samples.
#
# Notes: -
#
#>----------------------------------------------------------------------------<|
#> Install/load packages
packages <- c(
  "optparse", "tools", "terra", "dplyr", "tidyterra", "progress"
)

for (i in 1:NROW(packages)) {
  if (
    !require(packages[i], character.only = TRUE)
  ) {
    install.packages(
      packages[i],
      dependencies = TRUE
    )
    require(
      packages[i],
      character.only = TRUE
    )
  }
}

#>----------------------------------------------------------------------------<|
#> Options
option_list = list(
  make_option(c("-f", "--folder"), type = "character",
              help = "Folder containting the patches.shp files and metadata."
  )
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

file <- opt$folder

#>----------------------------------------------------------------------------<|
#> Settings
file_pattern_shp <- "patches_"
file_pattern_shpinfo <- "patches_"
file_pattern_metric <- "Metric_"

dir_main <- "/lud11/poppman/shared/dami"
dir_dat <- file.path(dir_main, "dat")
dir_pat <- file.path(dir_main, "dat", "pat")
dir_out <- file.path(dir_main, "out", file)

info <- readLines(
  file.path(dir_out, "job.info"), warn = FALSE, encoding = "utf-8"
  )

habitat <- strsplit(info[[1]], split = ": ")[[1]][2]
resistance <- basename(strsplit(info[[5]], split = ": ")[[1]][2])
samplename <- strsplit(info[[1]], split = ": ")[[1]][2]
r_patches <- strsplit(info[[4]], split = ": ")[[1]][2]

dir_pat <- file.path(dir_dat, "pat")
dir_spl <- file.path(dir_dat, "spl", samplename)

patches_file <- file.path(dir_dat, "pat", r_patches)
patches_raster <- terra::rast(patches_file)

#>----------------------------------------------------------------------------<|
#> Functions
create_raster <- function(x, template = patches_raster, metric = "delta-H") {
  # Get metric values
  delta_files <- list.files(
    dir_out,
    pattern = paste0("*_", x, ".txt"),
    full.names = TRUE
  )

  f_metric <- delta_files[which(startsWith(basename(delta_files), metric))]

  table <- data.frame(read.table(f_metric, header = TRUE)[, -1])
  values <- table[2:nrow(table),]
  
  f_sample <- file.path(dir_out, paste0(file_pattern_shp, x, ".shp"))
  sample <- terra::vect(f_sample, crs = terra::crs(template))
  sample <- cbind(sample, data.frame(Metric = values))

  sample %>%
    terra::rasterize(
      template,
      field = "Metric",
      background = NA,
      touches = FALSE,
      filename = paste0(tools::file_path_sans_ext(f_sample), ".tif"),
      overwrite = TRUE
    )
  
  return(paste0(tools::file_path_sans_ext(f_sample), ".tif"))
}

combine_rasters <- function(x, bs = 10) {
  n <- ceiling(length(x) / bs)
  pb <- progress_bar$new(
    format = "[:bar] :percent ETA: :eta",
    total = length(x)
  )
  
  rst_sum <- NULL
  cnt_sum <- NULL
  
  for (i in 1:n) {
    first <- (i - 1) * bs + 1
    last <- min(i * bs, length(x))
    raster_stack <- terra::rast(x[first:last])
    subset_sum <- terra::app(raster_stack, sum, na.rm = TRUE)
    terra::values(subset_sum)[which(is.na(terra::values(subset_sum)))] <- 0
    counter <- terra::app(!is.na(raster_stack), sum, na.rm = TRUE)
    
    if (is.null(rst_sum)) {
      rst_sum <- subset_sum
      cnt_sum <- counter
    } else {
      rst_sum <- (rst_sum + subset_sum)
      cnt_sum <- (cnt_sum + counter)
    }
    
    rm(raster_stack, subset_sum, counter)
    gc()
    pb$tick()
  }
  rst_average <- (rst_sum / cnt_sum)
  
  return(rst_average)
}

raster_variance <- function(x, mean_rst) {
  means <- terra::values(mean_rst)
  
  n <- ceiling(length(x))
  pb <- progress_bar$new(
    format = "[:bar] :percent ETA: :eta",
    total = length(x)
  )
  
  for (r in x) {
    rst <- terra::rast(r)
    vals <- terra::values(rst)
    diff <- vals - means
    diff2 <- diff^2
    diff2[which(is.na(diff2))] <- 0
    
    if (r == x[1]) {
      sums <- diff2
      sum_vals <- as.numeric(!is.na(diff))
    } else {
      sums <- (sums + diff2)
      sum_vals <- (sum_vals + as.numeric(!is.na(diff)))
    }
    
    rm(rst, vals)
    gc()
    pb$tick()
  }
  
  variance <- sums / (sum_vals - 1)
  
  return(
    list(
      variance[which(sum_vals > 0 & sum_vals < length(x) - 1)],
      sum_vals[which(sum_vals > 0 & sum_vals < length(x) - 1)]
      )
    )
}

#>----------------------------------------------------------------------------<|
#> Get file numbers
print("Counting files...")
N <- length(list.files(dir_out, pattern = "*.shp"))
file_numbers <- as.numeric(
  gsub(
    "[^0-9.-]", "", list.files(
      dir_out, pattern = "*.shp"
    )
  )
)

file_numbers <- file_numbers[which(!is.na(file_numbers))]
cat("Found", length(file_numbers), "files.\n")

# Extract data------------------------------------------------------------------
print("Writing rasters...")
# out_list <- c()
# pb <- progress_bar$new(
#   format = "[:bar] :percent ETA: :eta",
#   total = length(file_numbers)
# )
# 
# for (i in file_numbers) {
#    out_list <- c(out_list, create_raster(i))
#    pb$tick()
#  }

# Combine outputs---------------------------------------------------------------
print("Combining rasters...")
out_list <- list.files(dir_out, pattern = "*.tif", full.names = TRUE)
# rst_mean <- combine_rasters(out_list)
# 
# rst_mean_msk <- rst_mean %>%
#   terra::mask(mask = is.na(patches_raster), maskvalues = 1, updatevalue = NA)
# 
# terra::writeRaster(
#   rst_mean_msk,
#   filename = file.path(dir_out, "Mean_delta_H.tif"),
#   overwrite = TRUE
# )

print("Calculating statistics...")
set.seed(7)
all_patches <- patches_raster %>%
  terra::as.polygons()

centers <- terra::centroids(all_patches, inside = TRUE)
sample_points <- centers[runif(50000, min = 1, max = length(centers))]

pb <- progress_bar$new(
  format = "[:bar] :percent ETA: :eta",
  total = length(out_list)
)

results <- list()
# i <- 1
# for (r in out_list) {
#   raster <- terra::rast(r)
#   vals <- tryCatch(
#     terra::extract(raster, sample_points),
#     error = function(e){e}
#   )
#   
#   if (!inherits(vals, "error")) {
#     results[[i]] <- vals
#   }
#   
#   rm(raster)
#   rm(vals)
#   gc()
#   
#   if (i %% 100 == 0) {
#     tmp <- results
#     save(
#       tmp,
#       file = file.path(dir_dat, paste0("Value_samples_", i, ".Rdata"))
#       )
#     rm(tmp)
#     rm(results)
#     gc()
#     results <- list()
#   }
#   
#   i <- i + 1
#   pb$tick()
# }

parts <- list.files(dir_dat, pattern = "Value_samples_", full.names = TRUE)

for (p in parts) {
  load(p)
  results <- c(results, tmp)
}

result_matrix <- do.call(
  cbind, lapply(results, function(x){as.numeric(x[, 2])})
  )

nonzero <- apply(
  result_matrix, FUN = function(x){length(which(x > 0))}, MARGIN = 1
  )

resmat <- result_matrix[-which(nonzero == 5000),]

results_compact <- apply(
  resmat,
  1,
  FUN = function(x){
    as.numeric(na.omit(as.numeric(x))[1:138])
    }
  ) %>%
  t()

row.names(results_compact) <- seq(1:NROW(results_compact))

variances_list <- list()
i <- 1
for (s in seq(5, 120, 5)) {
  set.seed(s)
  random <- sample(NCOL(results_compact), size = s)
  ss <- results_compact[, random]
  print(paste("Subset rows:", NROW(ss)))
  print(paste("Subset cols:", NCOL(ss)))
  variances_list[[i]] <- apply(ss, FUN = var, MARGIN = 1)
  rm(random)
  rm(ss)
  gc()
  i <- i + 1
}

variances_mat <- do.call(c, variances_list)
variances_df <- data.frame(
  var = do.call(c, variances_list),
  n_samples = rep(seq(5, 120, 5), each = NROW(results_compact))
  )

# mod <- nls(
#   var ~ a * n_samples^(-b),
#   data = variances_df,
#   start = list(a = 0.001, b = 0.2)
#   )
# 
# summary(mod)
# 
# fit_x <- seq(5, 120, 1)
# fit_y <- predict(mod, newdata = list(n_samples = fit_x))

medians <- variances_df %>%
  dplyr::group_by(n_samples) %>%
  dplyr::summarise(medi = median(var))

means <- variances_df %>%
  dplyr::group_by(n_samples) %>%
  dplyr::summarise(mean = mean(var))

quants <- variances_df %>%
  dplyr::group_by(n_samples) %>%
  dplyr::summarise(quant = quantile(var, 0.75))

png(file = file.path(dir_out, "Sample_size.png"), width = 500, height = 380)
par(mar = c(3.1, 3.75, 0.1, 0.1), oma = c(0, 0, 0, 0), mgp = c(2, 1, 0))
plot(
  var * 10^3 ~ n_samples,
  data = variances_df,
  pch = 16, cex = 1.3, col = rgb(0, 102, 102, 80, maxColorValue = 255),
  xlab = "Number of samples that include the patch",
  ylab = expression("Variance"~(Delta~"H")%*%10^3),
  cex.lab = 1.3, las = 1
)
points(medi * 10^3 ~ n_samples, data = medians, col = "green", pch = 16)

#lines(fit_x, fit_y, col = "blue", lwd = 2)
dev.off()

# 
# stepsize <- 500
# n_samples <- seq(100, length(file_numbers), stepsize)
# pb <- progress_bar$new(
#   format = "[:bar] :percent ETA: :eta",
#   total = length(n_samples)
# )
# 
# varnce <- data.frame(x = c(), y = c())
# for (n in n_samples) {
#   print(paste0("N = ", n, " samples."))
#   sample_list <- sample(out_list, size = n)
#   rst_mean <- combine_rasters(sample_list) %>%
#     terra::mask(mask = is.na(patches_raster), maskvalues = 1, updatevalue = NA)
#   
#   #varnce <- c(varnce, raster_variance(sample_list, rst_mean))
#   out <- raster_variance(sample_list, rst_mean)
#   out_df <- data.frame(x = out[[2]], y = out[[1]])
#   
#   varnce <- rbind(varnce, out_df)
#   
#   pb$tick()
# }
# 
# save(varnce, file = "/lud11/poppman/shared/dami/fig/Samplesize_mean_var.Rsave")
# 
# png(
#   filename = paste0(
#     "/lud11/poppman/shared/dami/fig/Samplesize_mean_var.png"
#     )
#   )
# 
# plot(
#   y ~ x, data = varnce,
#   #n_samples[1:which(n_samples == n)],
#   #varnce,
#   xlab = "Number of samples",
#   ylab = "Mean variance for pixel value",
#   pch = 16,
#   col = rgb(0, 102, 102, maxColorValue = 255)
# )
# dev.off()
