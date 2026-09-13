.read_seed <- function(path, layout) {
  handle <- gzfile(path, "rb")
  on.exit(close(handle))
  weights <- vector("list", length(layout))
  names(weights) <- vapply(layout, `[[`, character(1), "name")
  for (i in seq_along(layout)) {
    dims <- as.integer(unlist(layout[[i]]$dim))
    n <- prod(dims)
    values <- readBin(handle, "numeric", n = n, size = 4L, endian = "little")
    if (length(values) != n || any(!is.finite(values)))
      stop("Incomplete or invalid model weight file: ", path, call. = FALSE)
    weights[[i]] <- array(values, dim = dims)
  }
  if (length(readBin(handle, "raw", n = 1L))) stop("Unexpected trailing model data: ", path, call. = FALSE)
  weights
}

.linear <- function(x, weight, bias) sweep(x %*% t(weight), 2L, as.numeric(bias), "+")

.activate <- function(x, activation) {
  if (activation == "ELU") { x[x < 0] <- expm1(x[x < 0]); return(x) }
  if (activation == "ReLU") { x[x < 0] <- 0; return(x) }
  if (activation == "LeakyReLU") { x[x < 0] <- 0.01 * x[x < 0]; return(x) }
  stop("Unsupported activation: ", activation, call. = FALSE)
}

.normalize <- function(x, w, prefix, indices = NULL) {
  pick <- function(suffix) {
    a <- as.numeric(w[[paste0(prefix, suffix)]])
    if (is.null(indices)) a else a[indices]
  }
  x <- sweep(x, 2L, pick(".running_mean"), "-")
  x <- sweep(x, 2L, pick(".weight") / sqrt(pick(".running_var") + 1e-5), "*")
  sweep(x, 2L, pick(".bias"), "+")
}

.dense_block <- function(x, w, prefix, activation, n_layers, final_linear = FALSE) {
  for (i in seq_len(n_layers) - 1L) {
    key <- paste0(prefix, ".layers.", i)
    x <- .linear(x, w[[paste0(key, ".weight")]], w[[paste0(key, ".bias")]])
    if (!(final_linear && i == n_layers - 1L)) {
      x <- .normalize(x, w, paste0(prefix, ".norms.", i))
      x <- .activate(x, activation)
    }
  }
  x
}

.head_linear <- function(x, w, prefix, head) {
  weight <- w[[paste0(prefix, ".weight")]]
  dims <- dim(weight)
  weight <- matrix(weight[head, , , drop = FALSE], nrow = dims[2], ncol = dims[3])
  .linear(x, weight, w[[paste0(prefix, ".bias")]][head, ])
}

.forward <- function(x, w, activation, head_count) {
  shared <- .dense_block(x, w, "shared", activation, 3L)
  output <- matrix(NA_real_, nrow(x), head_count + 1L)
  for (head in seq_len(head_count)) {
    bypass <- x
    for (i in 1:3) {
      bypass <- .head_linear(bypass, w, paste0("disease_tower.bypass_", i), head)
      indices <- (head - 1L) * ncol(bypass) + seq_len(ncol(bypass))
      bypass <- .normalize(bypass, w, paste0("disease_tower.norm_", i), indices)
      bypass <- .activate(bypass, activation)
    }
    hidden <- .head_linear(cbind(shared, bypass), w, "disease_tower.predictor_1", head)
    indices <- (head - 1L) * ncol(hidden) + seq_len(ncol(hidden))
    hidden <- .normalize(hidden, w, "disease_tower.predictor_norm", indices)
    hidden <- .activate(hidden, activation)
    output[, head] <- .head_linear(hidden, w, "disease_tower.predictor_2", head)
  }
  mortality <- .dense_block(x, w, "mortality_head.bypass", activation, 3L)
  output[, head_count + 1L] <- .dense_block(cbind(shared, mortality), w,
    "mortality_head.predictor", activation, 2L, final_linear = TRUE)
  output
}

.predict_scores <- function(data, meta, model_dir, batch_size, backend = "r", device = "cpu", python = NULL) {
  if (meta$model_id == "integrated") {
    out <- matrix(0, nrow(data), length(meta$endpoints))
    for (id in unlist(meta$components)) {
      scores <- .predict_scores(data, .metadata(id), model_dir, batch_size, backend, device, python)
      out <- out + sweep(scores, 2L, unlist(meta$stacking[[id]]), "*")
    }
    return(out)
  }
  x <- as.matrix(data[, unlist(meta$features), drop = FALSE])
  if (backend == "pytorch") {
    if (is.null(python)) {
      candidates <- c(Sys.getenv("RETICULATE_PYTHON"), Sys.which(c("python3", "python")))
      candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
      if (!length(candidates)) stop("Provide python = path to an environment with numpy and torch.", call. = FALSE)
      python <- candidates[[1]]
    }
    if (!is.character(python) || length(python) != 1L || is.na(python) || !file.exists(python))
      stop("python must be the path to an existing Python executable.", call. = FALSE)
    input_file <- tempfile("multiageclock-input-", fileext = ".bin")
    output_file <- tempfile("multiageclock-output-", fileext = ".bin")
    on.exit(unlink(c(input_file, output_file)), add = TRUE)
    writeBin(as.double(x), input_file, size = 8L, endian = "little")
    args <- c(system.file("python", "multiageclock_inference.py", package = "MultiAgeClock"),
              "--input", input_file, "--output", output_file,
              "--metadata", system.file("models", paste0(meta$model_id, ".json"), package = "MultiAgeClock"),
              "--weights", normalizePath(model_dir, winslash = "/", mustWork = TRUE),
              "--rows", nrow(x), "--columns", ncol(x), "--batch-size", batch_size, "--device", device)
    log <- suppressWarnings(system2(python, args = shQuote(as.character(args)), stdout = TRUE, stderr = TRUE))
    status <- attr(log, "status")
    if ((!is.null(status) && status != 0L) || !file.exists(output_file))
      stop("PyTorch inference failed. Check numpy, torch and the selected device.\n", paste(log, collapse = "\n"), call. = FALSE)
    n <- nrow(x) * length(meta$endpoints)
    if (file.info(output_file)$size != 8 * n) stop("Incomplete PyTorch output.", call. = FALSE)
    values <- readBin(output_file, "numeric", n = n, size = 8L, endian = "little")
    return(matrix(values, nrow = nrow(x), byrow = TRUE))
  }
  x <- sweep(sweep(x, 2L, unlist(meta$means), "-"), 2L, unlist(meta$sds), "/")
  out <- matrix(0, nrow(x), length(meta$endpoints))
  for (seed in unlist(meta$seed_files)) {
    weights <- .read_seed(file.path(model_dir, meta$model_id, seed), meta$tensor_layout)
    for (start in seq.int(1L, nrow(x), by = batch_size)) {
      index <- start:min(nrow(x), start + batch_size - 1L)
      out[index, ] <- out[index, , drop = FALSE] +
        .forward(x[index, , drop = FALSE], weights, meta$activation, length(meta$endpoints) - 1L)
    }
  }
  out / length(meta$seed_files)
}

.prepare_input <- function(data, meta, column_map) {
  if (!is.data.frame(data) || !nrow(data)) stop("data must be a non-empty data frame.", call. = FALSE)
  if (anyDuplicated(names(data))) stop("Duplicate input column names are not allowed.", call. = FALSE)
  required <- c("chronological_age", unlist(meta$features))
  if (!is.null(column_map)) {
    if (!is.character(column_map) || is.null(names(column_map)) ||
        anyNA(column_map) || anyNA(names(column_map)) || any(!nzchar(names(column_map))) ||
        anyDuplicated(names(column_map)) || anyDuplicated(unname(column_map)))
      stop("column_map must be a named character vector: canonical_name = input_column.", call. = FALSE)
    if (!all(names(column_map) %in% required) || !all(column_map %in% names(data)))
      stop("column_map contains unknown model features or absent input columns.", call. = FALSE)
    for (target in names(column_map)) {
      source <- column_map[[target]]
      if (target %in% names(data) && source != target)
        stop("Both canonical and mapped input columns exist for: ", target, call. = FALSE)
      names(data)[match(source, names(data))] <- target
    }
  }
  absent <- setdiff(required, names(data))
  if (length(absent)) stop("Missing required columns: ", paste(absent, collapse = ", "), call. = FALSE)
  invalid <- character()
  for (name in required) {
    original <- data[[name]]
    if (is.factor(original)) original <- as.character(original)
    if (!is.numeric(original) && !is.character(original))
      stop("Expected numeric measurements in column: ", name, call. = FALSE)
    value <- suppressWarnings(as.numeric(original))
    bad <- which(!is.finite(value))
    if (length(bad)) invalid <- c(invalid, paste0(name, " (rows ", paste(utils::head(bad, 5L), collapse = ", "), ")"))
    data[[name]] <- value
  }
  if (length(invalid)) stop("Missing or non-numeric/non-finite measurements: ",
                            paste(invalid, collapse = "; "), call. = FALSE)
  if (any(data$chronological_age <= 0)) stop("chronological_age must be positive years.", call. = FALSE)
  data
}

#' Calculate frozen biological ages
#' @param data A data frame containing complete measurements in documented model units.
#' @param model A model identifier or an object returned by [load_model()].
#' @param model_dir Directory containing model weights; ignored for a loaded model.
#' @param id_col Optional column containing sample identifiers, preserved as supplied.
#' @param column_map Optional named character vector, `canonical_name = input_column`.
#' @param batch_size Maximum records per neural-network forward pass.
#' @param backend `"r"` for the dependency-light R implementation, or `"pytorch"`
#'   to call the optional vectorized PyTorch backend in a separate local Python process.
#' @param device PyTorch device, such as `"cpu"` or `"cuda"`. The R backend uses CPU.
#' @param python Path to the Python executable with numpy and torch installed.
#'   If omitted, uses RETICULATE_PYTHON or python3/python on PATH.
#' @return A data frame in input row order with `row_id`, optional `sample_id`,
#'   `model`, `chronological_age`, `ba`, `raw_gap`, `baa`, and `age_outside_reference`.
#' @details `raw_gap` is BA minus chronological age. `baa` subtracts the expected
#'   BA from the frozen reference-age spline. No coefficients are fitted to the
#'   input table. Out-of-reference ages are flagged without clipping the result.
#'   These research models do not establish clinical diagnostic thresholds.
#' @export
predict_age <- function(data, model = "clinical_k17", model_dir = NULL, id_col = NULL,
                        column_map = NULL, batch_size = 256L, backend = c("r", "pytorch"), device = "cpu", python = NULL) {
  backend <- match.arg(backend)
  if (!is.character(device) || length(device) != 1L || is.na(device) || !nzchar(device))
    stop("device must name a compute device.", call. = FALSE)
  if (backend == "r" && device != "cpu") stop("The R backend uses device = 'cpu'.", call. = FALSE)
  if (!is.numeric(batch_size) || length(batch_size) != 1L || !is.finite(batch_size) ||
      batch_size < 1 || batch_size != floor(batch_size)) stop("batch_size must be a positive integer.", call. = FALSE)
  if (inherits(model, "multiage_model")) {
    meta <- model$metadata
    model_dir <- model$model_dir
  } else meta <- .metadata(model)
  input <- .prepare_input(data, meta, column_map)
  if (!is.null(id_col) && (!is.character(id_col) || length(id_col) != 1L || is.na(id_col) || !id_col %in% names(data)))
    stop("id_col must name one existing input column.", call. = FALSE)
  if (!inherits(model, "multiage_model")) model <- load_model(meta$model_id, model_dir)
  scores <- .predict_scores(input, meta, model$model_dir, batch_size, backend, device, python)
  c <- meta$calibration
  age <- input$chronological_age
  eta <- drop(scores %*% unlist(c$score_coefficients))
  ba <- age + c$scale * (eta - c$eta_mean)
  expected <- stats::splinefun(unlist(c$spline$knots), unlist(c$spline$values), method = "natural")(age)
  if (any(!is.finite(c(ba, expected)))) stop("Model produced non-finite outputs; check measurement units.", call. = FALSE)
  range <- unlist(c$age_range)
  out <- data.frame(row_id = seq_len(nrow(input)))
  if (!is.null(id_col)) out$sample_id <- data[[id_col]]
  out$model <- meta$model_id
  out$chronological_age <- age
  out$ba <- ba
  out$raw_gap <- ba - age
  out$baa <- ba - expected
  out$age_outside_reference <- age < range[1] | age > range[2]
  out
}
