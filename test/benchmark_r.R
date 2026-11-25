#!/usr/bin/env Rscript
# test/benchmark_r.R

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript benchmark_r.R <library> <filepath> <mz> <tol> <num_reps> <output_csv>")
}

target_lib <- tolower(args[1])
filepath <- args[2]
mz_value <- as.numeric(args[3])
mz_tolerance <- as.numeric(args[4])
num_reps <- as.numeric(args[5])
output_csv <- args[6]

get_rss <- function() {
  if (file.exists("/proc/self/status")) {
    lines <- readLines("/proc/self/status")
    vm_rss <- lines[grep("^VmRSS:", lines)]
    if (length(vm_rss) > 0) {
      rss_kb <- as.numeric(gsub("[^0-9]", "", vm_rss))
      return(rss_kb / 1024) # MB
    }
  }
  return(0)
}

# Ensure output directory exists
dir.create(dirname(output_csv), showWarnings = FALSE, recursive = TRUE)

if (target_lib == "maldiquant") {
  # ------------------------------------------------------------
  # Benchmarking MALDIquant
  # ------------------------------------------------------------
  message("\n[R] Loading MALDIquant and MALDIquantForeign...")
  library(MALDIquant)
  library(MALDIquantForeign)

  gc()
  rss_start <- get_rss()

  message(paste("[R] MALDIquant - Loading:", basename(filepath)))
  t_load_mq <- system.time({
    spectra <- MALDIquantForeign::importImzMl(filepath, centroided = TRUE)
  })["elapsed"]

  rss_loaded_mq <- get_rss()
  mem_load_mq <- max(0, rss_loaded_mq - rss_start)

  message(paste("[R] MALDIquant - Generating slices (repetitions:", num_reps, ")..."))
  t_slices_mq <- numeric(num_reps)
  for (i in 1:num_reps) {
    t_slices_mq[i] <- system.time({
      slices <- MALDIquant::msiSlices(spectra, center = mz_value, tolerance = mz_tolerance)
    })["elapsed"]
  }
  mean_time_mq <- mean(t_slices_mq)

  results <- data.frame(
    Library = "MALDIquant",
    LoadTime_s = t_load_mq,
    LoadMem_MB = mem_load_mq,
    AvgSliceTime_s = mean_time_mq
  )
  write.csv(results, output_csv, row.names = FALSE)
  message(paste("[R] MALDIquant results written to", output_csv))

} else if (target_lib == "cardinal") {
  # ------------------------------------------------------------
  # Benchmarking Cardinal
  # ------------------------------------------------------------
  message("\n[R] Loading Cardinal...")
  library(Cardinal)

  gc()
  rss_start_card <- get_rss()

  message(paste("[R] Cardinal - Loading:", basename(filepath)))
  t_load_card <- system.time({
    mse <- Cardinal::readMSIData(filepath)
  })["elapsed"]

  rss_loaded_card <- get_rss()
  mem_load_card <- max(0, rss_loaded_card - rss_start_card)

  message(paste("[R] Cardinal - Generating slices (repetitions:", num_reps, ")..."))
  t_slices_card <- numeric(num_reps)
  
  # Determine if we have MSImagingArrays (processed/unaligned spectra) or MSImagingExperiment
  is_arrays <- inherits(mse, "MSImagingArrays")

  for (i in 1:num_reps) {
    t_slices_card[i] <- system.time({
      if (is_arrays) {
        # Processed dataset pathway (each pixel has different m/z values)
        # Iterate over all spectra in parallel/loop and find target mass
        slice_data <- sapply(seq_along(mse), function(idx) {
          mzs <- Cardinal::mz(mse)[[idx]]
          ints <- Cardinal::intensity(mse)[[idx]]
          in_tol <- which(abs(mzs - mz_value) <= mz_tolerance)
          if (length(in_tol) > 0) {
            max(ints[in_tol])
          } else {
            0.0
          }
        })
      } else {
        # Continuous dataset pathway (shared m/z grid)
        feat_idx <- Cardinal::features(mse, mz = mz_value, tolerance = mz_tolerance, units = "mz")
        if (length(feat_idx) > 0) {
          raw_data <- Cardinal::iData(mse)[feat_idx, , drop = FALSE]
          if (length(feat_idx) > 1) {
            slice_data <- apply(raw_data, 2, max)
          } else {
            slice_data <- as.numeric(raw_data)
          }
        } else {
          slice_data <- rep(0, ncol(mse))
        }
      }
    })["elapsed"]
  }
  mean_time_card <- mean(t_slices_card)

  results <- data.frame(
    Library = "Cardinal",
    LoadTime_s = t_load_card,
    LoadMem_MB = mem_load_card,
    AvgSliceTime_s = mean_time_card
  )
  write.csv(results, output_csv, row.names = FALSE)
  message(paste("[R] Cardinal results written to", output_csv))

} else {
  stop(paste("Unknown library requested:", target_lib))
}
