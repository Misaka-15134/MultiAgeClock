.model_ids <- c("clinical_k17", "clinical", "olink", "nmr", "integrated")
.release <- "v0.1.0"

.metadata <- function(model) {
  if (!is.character(model) || length(model) != 1L || is.na(model) || !model %in% .model_ids)
    stop("Unknown model. Choose: ", paste(.model_ids, collapse = ", "), call. = FALSE)
  jsonlite::read_json(system.file("models", paste0(model, ".json"), package = "MultiAgeClock"),
                      simplifyVector = FALSE)
}

.cache_dir <- function(model_dir) {
  if (is.null(model_dir)) file.path(tools::R_user_dir("MultiAgeClock", "cache"), .release) else model_dir
}

#' List the available frozen models
#' @return A data frame with model identifiers and input feature counts.
#' @examples
#' list_models()
#' @export
list_models <- function() {
  data.frame(model = .model_ids,
             label = c("Clinical K17", "Clinical full", "Olink", "NMR", "Integrated"),
             n_features = vapply(.model_ids, function(x) length(.metadata(x)$features), integer(1)),
             version = .release, row.names = NULL)
}

#' Inspect required features and measurement units
#' @param model A model identifier from [list_models()].
#' @return A data frame. Age in years is always required in addition to biomarkers.
#' @examples
#' model_features("clinical_k17")
#' @export
model_features <- function(model = "clinical_k17") {
  meta <- .metadata(model)
  records <- lapply(meta$feature_schema, function(x) {
    data.frame(feature = x$feature, label = x$label, unit = x$unit,
               ukb_field_id = if (is.null(x$ukb_field_id)) NA_integer_ else x$ukb_field_id,
               nightingale_name = if (is.null(x$nightingale_name)) NA_character_ else x$nightingale_name)
  })
  rbind(data.frame(feature = "chronological_age", label = "Chronological age", unit = "years",
                   ukb_field_id = NA_integer_, nightingale_name = NA_character_), do.call(rbind, records))
}

#' Download versioned inference weights
#' @param models One or more model identifiers. Integrated downloads its three components.
#' @param model_dir Directory to store model folders; defaults to the versioned user cache.
#' @return The model directory, invisibly. Existing files are reused.
#' @details Weights are downloaded from the release associated with this package
#'   version. Input data are never uploaded. Installation of the R package alone
#'   does not download model weights.
#' @export
download_models <- function(models = "clinical_k17", model_dir = NULL) {
  if (!length(models)) stop("Provide at least one model.", call. = FALSE)
  invisible(lapply(models, .metadata))
  previous_options <- options(timeout = max(600, getOption("timeout", 60)))
  on.exit(options(previous_options), add = TRUE)
  if ("integrated" %in% models) models <- c(setdiff(models, "integrated"), "clinical", "olink", "nmr")
  model_dir <- .cache_dir(model_dir)
  for (id in unique(models)) {
    meta <- .metadata(id)
    dest <- file.path(model_dir, id)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    for (filename in unlist(meta$seed_files)) {
      target <- file.path(dest, filename)
      if (file.exists(target)) next
      url <- paste0("https://github.com/Misaka-15134/MultiAgeClock/releases/download/", .release,
                    "/", id, "_", filename)
      partial <- tempfile(pattern = "download-", tmpdir = dest)
      tryCatch({
        status <- utils::download.file(url, partial, mode = "wb", quiet = FALSE)
        if (status != 0L) stop("Weight download failed: ", filename, call. = FALSE)
        if (!file.rename(partial, target)) stop("Cannot save weight file: ", target, call. = FALSE)
      }, finally = { if (file.exists(partial)) unlink(partial) })
    }
  }
  invisible(normalizePath(model_dir, winslash = "/", mustWork = TRUE))
}

#' Load a frozen model
#' @param model A model identifier from [list_models()].
#' @param model_dir Directory containing downloaded model folders.
#' @return A `multiage_model` object. Large tensors are read during prediction,
#'   one seed at a time to bound memory use.
#' @export
load_model <- function(model = "clinical_k17", model_dir = NULL) {
  meta <- .metadata(model)
  model_dir <- .cache_dir(model_dir)
  ids <- if (model == "integrated") unlist(meta$components) else model
  files <- unlist(lapply(ids, function(id) file.path(model_dir, id, unlist(.metadata(id)$seed_files))))
  absent <- files[!file.exists(files)]
  if (length(absent)) stop("Model weights are missing. Run download_models(\"", model,
                          "\", model_dir = ...) first. Missing: ", absent[1], call. = FALSE)
  structure(list(id = model, metadata = meta, model_dir = model_dir), class = "multiage_model")
}

#' Read synthetic example inputs
#' @param model A model identifier from [list_models()].
#' @return A data frame of nine artificial records. These are software examples,
#'   not participants and not a reference population.
#' @examples
#' head(example_data())
#' @export
example_data <- function(model = "clinical_k17") {
  .metadata(model)
  utils::read.csv(system.file("extdata", paste0(model, "_example.csv"), package = "MultiAgeClock"),
                  check.names = FALSE, stringsAsFactors = FALSE)
}
