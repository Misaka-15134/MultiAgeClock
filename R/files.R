#' Score a CSV, TSV or Excel table
#' @param input Path to a CSV, TSV or XLSX file. The file is not modified.
#' @param output Output CSV, TSV or XLSX path. Defaults to a sibling CSV file.
#' @param model A model identifier or loaded model.
#' @param model_dir Directory containing model weights.
#' @param id_col Optional sample identifier column.
#' @param column_map Optional named character vector mapping model names to input columns.
#' @param sheet Excel sheet name or index. Defaults to the first sheet.
#' @param batch_size Maximum records per forward pass.
#' @param backend `"r"` or `"pytorch"`; see [predict_age()].
#' @param device PyTorch compute device, such as `"cpu"` or `"cuda"`.
#' @param python Path to a Python executable with numpy and torch; see [predict_age()].
#' @return The predictions, invisibly; only prediction columns are written.
#' @details Text imports preserve sample identifiers such as `00123`. Excel
#'   identifiers should be stored as text in the workbook. Existing output files
#'   are never overwritten. XLSX input requires readxl; XLSX output requires writexl.
#' @export
score_file <- function(input, output = NULL, model = "clinical_k17", model_dir = NULL,
                       id_col = NULL, column_map = NULL, sheet = 1, batch_size = 256L,
                       backend = c("r", "pytorch"), device = "cpu", python = NULL) {
  backend <- match.arg(backend)
  if (!is.character(input) || length(input) != 1L || is.na(input) || !file.exists(input))
    stop("input must name an existing table file.", call. = FALSE)
  if (is.null(output)) output <- paste0(tools::file_path_sans_ext(input), "_MultiAgeClock.csv")
  if (!is.character(output) || length(output) != 1L || is.na(output) || !nzchar(output))
    stop("output must be a file path.", call. = FALSE)
  if (file.exists(output)) stop("Output already exists: ", output, call. = FALSE)
  ext <- tolower(tools::file_ext(input))
  out_ext <- tolower(tools::file_ext(output))
  if (!out_ext %in% c("csv", "tsv", "xlsx")) stop("Output format must be CSV, TSV or XLSX.", call. = FALSE)
  if (out_ext == "xlsx" && !requireNamespace("writexl", quietly = TRUE))
    stop("Install writexl for XLSX output.", call. = FALSE)
  if (ext %in% c("csv", "tsv")) {
    data <- utils::read.table(input, header = TRUE, sep = if (ext == "csv") "," else "\t",
      colClasses = "character", check.names = FALSE, comment.char = "", quote = "\"",
      fileEncoding = "UTF-8-BOM", na.strings = c("", "NA"))
  } else if (ext == "xlsx") {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl for XLSX input.", call. = FALSE)
    data <- as.data.frame(readxl::read_excel(input, sheet = sheet, col_types = "text", .name_repair = "minimal"))
  } else stop("Input format must be CSV, TSV or XLSX.", call. = FALSE)
  result <- predict_age(data, model, model_dir, id_col, column_map, batch_size, backend, device, python)
  if (out_ext == "xlsx") writexl::write_xlsx(result, output)
  else utils::write.table(result, output, sep = if (out_ext == "csv") "," else "\t",
                           row.names = FALSE, col.names = TRUE, quote = TRUE, fileEncoding = "UTF-8")
  invisible(result)
}
