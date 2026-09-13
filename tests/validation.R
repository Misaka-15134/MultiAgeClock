library(MultiAgeClock)

expect_error <- function(expression, pattern) {
  message <- tryCatch({force(expression); NULL}, error = function(e) conditionMessage(e))
  stopifnot(is.character(message), grepl(pattern, message))
}

d <- example_data()
expect_error(predict_age(d[, -3]), "Missing required columns")
d$cystatin_c[2] <- NA_real_
expect_error(predict_age(d), "cystatin_c.*rows 2")
d <- example_data()
d$cystatin_c[1] <- Inf
expect_error(predict_age(d), "non-finite")
d$cystatin_c <- "not a number"
expect_error(predict_age(d), "cystatin_c")
d <- example_data()
names(d)[3] <- names(d)[4]
expect_error(predict_age(d), "Duplicate")
expect_error(predict_age(example_data(), batch_size = 0), "positive integer")
expect_error(predict_age(example_data(), column_map = c(unknown = "sample_id")), "unknown")
expect_error(predict_age(example_data()[0, ]), "non-empty")
models <- list_models()
stopifnot(identical(models$n_features, c(17L, 4L, 47L, 2920L, 168L, 3135L)))
stopifnot(identical(model_features("clinical_k4")$feature,
                   c("chronological_age", "cystatin_c", "systolic_blood_pressure", "waist_circumference", "hba1c")))
expect_error(predict_age(example_data("clinical_k4")[, -3], "clinical_k4"), "Missing required columns")
for (id in models$model) {
  stopifnot(all(model_features(id)$feature %in% names(example_data(id))))
  stopifnot(nrow(example_data(id)) == 9L)
}

weights <- Sys.getenv("MULTIAGECLOCK_TEST_MODELS")
if (nzchar(weights)) {
  for (id in models$model) {
    d <- example_data(id)
    expected <- read.csv(system.file("extdata", paste0(id, "_expected.csv"), package = "MultiAgeClock"))
    elapsed <- system.time(result <- predict_age(d[, rev(seq_along(d))], id, weights,
                                                 id_col = "sample_id", batch_size = 4L))[["elapsed"]]
    errors <- vapply(c("ba", "raw_gap", "baa"), function(x) max(abs(result[[x]] - expected[[x]])), numeric(1))
    stopifnot(all(errors < 1e-3), identical(result$sample_id, d$sample_id))
    stopifnot(max(abs(result$raw_gap - result$ba + result$chronological_age)) < 1e-12)
    cat(id, "rows", nrow(d), "seconds", elapsed, "max_abs_errors", errors, "\n")
  }
  d <- example_data()[5, ]
  names(d)[names(d) == "chronological_age"] <- "age_years"
  result <- predict_age(d, model_dir = weights,
                       column_map = c(chronological_age = "age_years"), batch_size = 1L)
  expected <- read.csv(system.file("extdata", "clinical_k17_expected.csv", package = "MultiAgeClock"))
  stopifnot(abs(result$ba - expected$ba[5]) < 1e-3, abs(result$baa - expected$baa[5]) < 1e-3)

  folder <- tempfile("multiageclock-tables-")
  dir.create(folder)
  d <- example_data()[1:2, ]
  d$sample_id <- c("00123", "00007")
  formats <- c("csv", "tsv")
  if (requireNamespace("readxl", quietly = TRUE) && requireNamespace("writexl", quietly = TRUE)) formats <- c(formats, "xlsx")
  for (extension in formats) {
    input <- file.path(folder, paste0("input.", extension))
    output <- file.path(folder, paste0("output.", extension))
    if (extension == "xlsx") writexl::write_xlsx(d, input)
    else write.table(d, input, sep = if (extension == "csv") "," else "\t", row.names = FALSE)
    original <- readBin(input, "raw", n = file.info(input)$size)
    result <- score_file(input, output, model_dir = weights, id_col = "sample_id")
    stopifnot(identical(result$sample_id, d$sample_id), file.exists(output))
    saved <- if (extension == "xlsx") as.data.frame(readxl::read_excel(output, col_types = "text"))
      else read.table(output, header = TRUE, sep = if (extension == "csv") "," else "\t", colClasses = "character")
    stopifnot(identical(saved$sample_id, d$sample_id), nrow(saved) == 2L)
    stopifnot(identical(readBin(input, "raw", n = file.info(input)$size), original))
    expect_error(score_file(input, output, model_dir = weights), "already exists")
    unlink(c(input, output))
    cat(extension, "file roundtrip passed\n")
  }
  unlink(folder)
} else message("Model-asset numerical tests skipped: set MULTIAGECLOCK_TEST_MODELS to run them.")
cat("VALIDATION_PASSED\n")
