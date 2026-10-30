#' Read pipeline files into a named list of data.tables
#'
#' Prepares pipeline inputs for [run_strand_filtering()]. Each element of
#' `files` can be either a file path, which is read with
#' [data.table::fread()], or a data.frame/data.table already in memory. In
#' both cases the columns are renamed to the names expected by the rest of
#' the package. In-memory inputs are copied, so the caller's objects are not
#' modified.
#'
#' @param files A named list. Each name is the pipeline name (e.g.
#'   `"Bismark"`). Each value is either a single file path (character) or a
#'   data.frame/data.table. Files and tables must have either 5 columns
#'   (`chr, pos, TRead, MRead, ML`) or 6 columns (the same five plus
#'   `strand`). Files must have no header.
#'
#' @return A named list of data.tables, with the same names as `files`.
#'   Columns are `chr, pos, TRead, MRead, ML` for 5-column inputs, or
#'   `chr, pos, TRead, MRead, ML, strand` for 6-column inputs.
#' @export
#'
#' @examples
#' \dontrun{
#' # All file paths
#' pipelines <- read_pipeline_files(list(
#'   Bismark = "Bismark.chr.txt",
#'   Biscuit = "biscuit.chr.txt"
#' ))
#'
#' # Mix of file paths and tables already loaded
#' pipelines <- read_pipeline_files(list(
#'   Biscuit    = "pipeline_test_1.txt",
#'   Bismark    = already_loaded_dt,
#'   `Bwa-meth` = "pipeline_test_4.txt"
#' ))
#' }
read_pipeline_files <- function(files) {
  if (!is.list(files) || is.null(names(files)) ||
      any(is.na(names(files))) || any(names(files) == "")) {
    stop("`files` must be a named list, with each name giving the pipeline name.")
  }

  result <- lapply(seq_along(files), function(i) {
    name <- names(files)[i]
    x <- files[[i]]

    if (is.character(x) && length(x) == 1) {
      if (!file.exists(x)) {
        stop("File not found for '", name, "': ", x)
      }
      dt <- data.table::fread(x, header = FALSE)
    } else if (is.data.frame(x)) {
      dt <- data.table::copy(data.table::as.data.table(x))
    } else {
      stop("Element '", name, "' must be a single file path or a data.frame/data.table.")
    }

    if (!(ncol(dt) %in% c(5L, 6L))) {
      stop("Input for '", name, "' has ", ncol(dt),
           " columns; expected 5 or 6 (chr, pos, TRead, MRead, ML[, strand]).")
    }

    # Keep already-correct column names, regardless of their order. Unnamed
    # files (V1, V2, ...) retain the original positional interpretation.
    if ("ml" %in% names(dt) && !"ML" %in% names(dt)) {
      data.table::setnames(dt, "ml", "ML")
    }
    if ("Strand" %in% names(dt) && !"strand" %in% names(dt)) {
      data.table::setnames(dt, "Strand", "strand")
    }
    required <- c("chr", "pos", "TRead", "MRead", "ML")
    if (all(required %in% names(dt)) &&
        (ncol(dt) == 5L || "strand" %in% names(dt))) {
      data.table::setcolorder(dt, c(required, if (ncol(dt) == 6L) "strand"))
    } else if (all(grepl("^V[0-9]+$", names(dt)))) {
      data.table::setnames(dt, c(required, if (ncol(dt) == 6L) "strand"))
    } else {
      stop("Input for '", name, "' has unexpected column names: ",
           paste(names(dt), collapse = ", "),
           ". Supply named columns or a headerless 5/6-column file.")
    }

    if (anyNA(dt$chr) || anyNA(dt$pos) || anyNA(dt$TRead) || anyNA(dt$MRead) ||
        !is.numeric(dt$pos) || !is.numeric(dt$TRead) ||
        !is.numeric(dt$MRead) || any(dt$pos < 1) ||
        any(dt$TRead < 0) || any(dt$MRead < 0) ||
        any(dt$MRead > dt$TRead)) {
      stop("Input for '", name, "' contains invalid coordinates or read counts.")
    }
    if ("strand" %in% names(dt) &&
        (anyNA(dt$strand) || any(!dt$strand %in% c("+", "-")))) {
      stop("Input for '", name, "' has strand values other than '+' and '-'.")
    }
    dt
  })

  names(result) <- names(files)
  result
}
