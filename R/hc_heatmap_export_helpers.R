.hc_export_path_with_ext <- function(path, ext) {
  ext <- base::paste0(".", base::sub("^\\.", "", base::as.character(ext[[1]])))
  base::paste0(tools::file_path_sans_ext(path), ext)
}

# Validate the internal OOXML package, not only the outer ZIP. Some openxlsx
# releases emit unused drawing relationships without the corresponding parts;
# Excel repairs those workbooks differently across platforms.
.hc_xlsx_read_zip_part <- function(path, entry) {
  con <- NULL
  tryCatch(
    {
      con <- base::unz(path, entry, open = "rb")
      chunks <- base::list()
      repeat {
        bytes <- base::readBin(con, what = "raw", n = 65536L)
        if (base::length(bytes) == 0) {
          break
        }
        chunks[[base::length(chunks) + 1L]] <- bytes
      }
      if (base::length(chunks) == 0) {
        return("")
      }
      base::rawToChar(base::do.call(base::c, chunks))
    },
    error = function(e) NA_character_,
    finally = {
      if (!base::is.null(con)) {
        base::close(con)
      }
    }
  )
}

.hc_xlsx_xml_attribute <- function(tag, attribute) {
  match <- stringi::stri_match_first_regex(
    tag,
    base::paste0("\\b", attribute, "\\s*=\\s*([\"'])(.*?)\\1"),
    opts_regex = base::list(case_insensitive = TRUE)
  )
  if (base::ncol(match) < 3L || base::is.na(match[[1L, 3L]])) {
    return(NA_character_)
  }
  match[[1L, 3L]]
}

.hc_xlsx_relationship_source <- function(rel_entry) {
  if (base::identical(rel_entry, "_rels/.rels")) {
    return("")
  }
  rel_dir <- base::dirname(rel_entry)
  source_dir <- base::dirname(rel_dir)
  source_file <- base::sub("\\.rels$", "", base::basename(rel_entry))
  base::gsub(
    "\\\\",
    "/",
    base::file.path(source_dir, source_file),
    fixed = FALSE
  )
}

.hc_xlsx_resolve_target <- function(rel_entry, target) {
  if (base::is.na(target) || !base::nzchar(target)) {
    return(NA_character_)
  }
  target <- base::sub("[?#].*$", "", target)
  target <- base::gsub("\\\\", "/", utils::URLdecode(target), fixed = FALSE)
  if (base::grepl("^[A-Za-z][A-Za-z0-9+.-]*:", target)) {
    return(NA_character_)
  }

  source <- .hc_xlsx_relationship_source(rel_entry)
  base_dir <- if (base::nzchar(source)) base::dirname(source) else ""
  if (base::startsWith(target, "/")) {
    parts <- base::strsplit(base::sub("^/+", "", target), "/", fixed = TRUE)[[1]]
  } else {
    combined <- base::paste(base::c(base_dir, target), collapse = "/")
    parts <- base::strsplit(combined, "/", fixed = TRUE)[[1]]
  }

  resolved <- base::character(0)
  for (part in parts) {
    if (!base::nzchar(part) || base::identical(part, ".")) {
      next
    }
    if (base::identical(part, "..")) {
      if (base::length(resolved) == 0) {
        return(NA_character_)
      }
      resolved <- resolved[-base::length(resolved)]
    } else {
      resolved <- base::c(resolved, part)
    }
  }
  base::paste(resolved, collapse = "/")
}

.hc_xlsx_dangling_relationships <- function(path, entries) {
  rel_entries <- entries[base::grepl("\\.rels$", entries, ignore.case = TRUE)]
  missing <- base::list()
  for (rel_entry in rel_entries) {
    xml <- .hc_xlsx_read_zip_part(path, rel_entry)
    if (base::is.na(xml)) {
      next
    }
    tags <- stringi::stri_extract_all_regex(
      xml,
      "<Relationship\\b[^>]*/\\s*>",
      opts_regex = base::list(case_insensitive = TRUE)
    )[[1]]
    tags <- tags[!base::is.na(tags)]
    for (tag in tags) {
      target_mode <- .hc_xlsx_xml_attribute(tag, "TargetMode")
      if (!base::is.na(target_mode) &&
        base::identical(base::tolower(target_mode), "external")) {
        next
      }
      target <- .hc_xlsx_xml_attribute(tag, "Target")
      resolved <- .hc_xlsx_resolve_target(rel_entry, target)
      if (base::is.na(resolved) || resolved %in% entries) {
        next
      }
      missing[[base::length(missing) + 1L]] <- base::data.frame(
        rel_entry = rel_entry,
        source_entry = .hc_xlsx_relationship_source(rel_entry),
        relationship_id = .hc_xlsx_xml_attribute(tag, "Id"),
        target = target,
        resolved_target = resolved,
        tag = tag,
        stringsAsFactors = FALSE
      )
    }
  }
  if (base::length(missing) == 0) {
    return(base::data.frame())
  }
  base::do.call(base::rbind, missing)
}

.hc_xlsx_missing_content_overrides <- function(path, entries) {
  xml <- .hc_xlsx_read_zip_part(path, "[Content_Types].xml")
  if (base::is.na(xml)) {
    return(base::data.frame())
  }
  tags <- stringi::stri_extract_all_regex(
    xml,
    "<Override\\b[^>]*/\\s*>",
    opts_regex = base::list(case_insensitive = TRUE)
  )[[1]]
  tags <- tags[!base::is.na(tags)]
  missing <- base::lapply(tags, function(tag) {
    part <- .hc_xlsx_xml_attribute(tag, "PartName")
    resolved <- base::sub("^/+", "", utils::URLdecode(part))
    if (base::is.na(resolved) || resolved %in% entries) {
      return(NULL)
    }
    base::data.frame(part = resolved, tag = tag, stringsAsFactors = FALSE)
  })
  missing <- missing[!base::vapply(missing, base::is.null, FUN.VALUE = base::logical(1))]
  if (base::length(missing) == 0) {
    return(base::data.frame())
  }
  base::do.call(base::rbind, missing)
}

.hc_xlsx_package_links_valid <- function(path, entries) {
  base::nrow(.hc_xlsx_dangling_relationships(path, entries)) == 0 &&
    base::nrow(.hc_xlsx_missing_content_overrides(path, entries)) == 0
}

.hc_xlsx_source_uses_relationship <- function(path, source_entry, relationship_id, entries) {
  if (base::is.na(source_entry) || !base::nzchar(source_entry) ||
    !(source_entry %in% entries) ||
    base::is.na(relationship_id) || !base::nzchar(relationship_id)) {
    return(FALSE)
  }
  xml <- .hc_xlsx_read_zip_part(path, source_entry)
  if (base::is.na(xml)) {
    return(TRUE)
  }
  base::grepl(base::paste0("\"", relationship_id, "\""), xml, fixed = TRUE) ||
    base::grepl(base::paste0("'", relationship_id, "'"), xml, fixed = TRUE)
}

.hc_xlsx_repair_dangling_parts <- function(path) {
  entries <- tryCatch(
    base::suppressWarnings(utils::unzip(path, list = TRUE)$Name),
    error = function(e) base::character(0)
  )
  dangling <- .hc_xlsx_dangling_relationships(path, entries)
  missing_overrides <- .hc_xlsx_missing_content_overrides(path, entries)
  if (base::nrow(dangling) == 0 && base::nrow(missing_overrides) == 0) {
    return(invisible(path))
  }

  removable <- if (base::nrow(dangling) > 0) {
    !base::vapply(
      base::seq_len(base::nrow(dangling)),
      function(i) {
        .hc_xlsx_source_uses_relationship(
          path = path,
          source_entry = dangling$source_entry[[i]],
          relationship_id = dangling$relationship_id[[i]],
          entries = entries
        )
      },
      FUN.VALUE = base::logical(1)
    )
  } else {
    base::logical(0)
  }
  if (base::nrow(dangling) > 0 && !base::all(removable)) {
    stop("The XLSX package contains a missing part that is still referenced by a worksheet.")
  }

  extract_dir <- base::tempfile("hc-xlsx-repair-")
  repaired <- base::tempfile(
    pattern = "hc-xlsx-repaired-",
    # Build the second archive on R's local temporary filesystem. Creating it
    # beside `path` can exhaust a quota-limited Docker/bind-mounted output
    # volume because the original, atomic temp, and repaired archive coexist.
    tmpdir = base::tempdir(),
    fileext = ".xlsx"
  )
  base::dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  base::on.exit(
    {
      extract_abs <- base::normalizePath(extract_dir, winslash = "/", mustWork = FALSE)
      temp_abs <- base::normalizePath(base::tempdir(), winslash = "/", mustWork = FALSE)
      if (base::startsWith(extract_abs, base::paste0(temp_abs, "/"))) {
        base::unlink(extract_dir, recursive = TRUE, force = TRUE)
      }
      if (base::file.exists(repaired)) {
        base::file.remove(repaired)
      }
    },
    add = TRUE
  )
  extracted <- base::suppressWarnings(utils::unzip(path, exdir = extract_dir))
  if (base::length(extracted) == 0) {
    stop("Could not extract the XLSX package for validation.")
  }

  if (base::nrow(dangling) > 0) {
    for (rel_entry in base::unique(dangling$rel_entry)) {
      rel_path <- base::file.path(extract_dir, rel_entry)
      xml <- base::paste(base::readLines(rel_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
      tags <- dangling$tag[dangling$rel_entry == rel_entry]
      for (tag in tags) {
        xml <- base::gsub(tag, "", xml, fixed = TRUE)
      }
      base::writeLines(xml, rel_path, useBytes = TRUE)
    }
  }
  if (base::nrow(missing_overrides) > 0) {
    types_path <- base::file.path(extract_dir, "[Content_Types].xml")
    xml <- base::paste(base::readLines(types_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    for (tag in missing_overrides$tag) {
      xml <- base::gsub(tag, "", xml, fixed = TRUE)
    }
    base::writeLines(xml, types_path, useBytes = TRUE)
  }

  files <- base::list.files(
    extract_dir,
    recursive = TRUE,
    all.files = TRUE,
    full.names = FALSE,
    include.dirs = FALSE,
    no.. = TRUE
  )
  zip::zipr(
    zipfile = repaired,
    files = files,
    include_directories = FALSE,
    root = extract_dir,
    mode = "mirror"
  )
  repaired_entries <- base::suppressWarnings(utils::unzip(repaired, list = TRUE)$Name)
  if (!.hc_xlsx_package_links_valid(repaired, repaired_entries)) {
    stop("Could not repair invalid XLSX package relationships.")
  }

  removed <- base::file.remove(path)
  if (!base::isTRUE(removed)) {
    stop("Could not remove the invalid XLSX package after creating its repaired replacement.")
  }
  moved <- base::suppressWarnings(base::file.rename(repaired, path))
  if (!base::isTRUE(moved)) {
    moved <- base::file.copy(repaired, path, overwrite = FALSE)
    if (base::isTRUE(moved)) {
      base::file.remove(repaired)
    }
  }
  if (!base::isTRUE(moved)) {
    stop("Could not replace the invalid XLSX package after repair.")
  }
  invisible(path)
}

.hc_xlsx_xml_parts_valid <- function(path, entries) {
  xml_entries <- entries[base::grepl("(\\.xml|\\.rels)$", entries, ignore.case = TRUE)]
  if (base::length(xml_entries) == 0) {
    return(FALSE)
  }

  forbidden <- base::as.raw(base::c(0:8, 11:12, 14:31))
  for (entry in xml_entries) {
    con <- NULL
    part_ok <- tryCatch(
      {
        con <- base::unz(path, entry, open = "rb")
        repeat {
          bytes <- base::readBin(con, what = "raw", n = 65536L)
          if (base::length(bytes) == 0) {
            break
          }
          if (base::any(bytes %in% forbidden)) {
            return(FALSE)
          }
        }
        TRUE
      },
      error = function(e) FALSE,
      finally = {
        if (!base::is.null(con)) {
          base::close(con)
        }
      }
    )
    if (!base::isTRUE(part_ok)) {
      return(FALSE)
    }
  }
  TRUE
}

.hc_xlsx_clean_text <- function(x) {
  if (base::is.null(x) || base::length(x) == 0) {
    return(base::as.character(x))
  }
  x <- stringi::stri_enc_toutf8(
    base::as.character(x),
    is_unknown_8bit = TRUE,
    validate = TRUE
  )
  stringi::stri_replace_all_regex(
    x,
    "[\\x{0000}-\\x{0008}\\x{000B}\\x{000C}\\x{000E}-\\x{001F}\\x{007F}-\\x{009F}\\x{FFFE}\\x{FFFF}]",
    ""
  )
}

# Excel cells are limited to 32,767 UTF-16 units. Keep chunks below that limit
# so supplementary Unicode characters cannot accidentally overflow a cell.
.hc_xlsx_split_text <- function(x, max_utf16_units = 30000L) {
  if (base::length(x) == 0 || base::is.na(x)) {
    return(x)
  }
  code_points <- base::utf8ToInt(x)
  if (base::length(code_points) == 0) {
    return("")
  }

  units <- 1L + base::as.integer(code_points > 0xffffL)
  if (base::sum(units) <= max_utf16_units) {
    return(x)
  }

  groups <- base::integer(base::length(code_points))
  group <- 1L
  used <- 0L
  for (i in base::seq_along(code_points)) {
    if (used + units[[i]] > max_utf16_units) {
      group <- group + 1L
      used <- 0L
    }
    groups[[i]] <- group
    used <- used + units[[i]]
  }
  base::vapply(
    base::split(code_points, groups),
    base::intToUtf8,
    FUN.VALUE = base::character(1),
    USE.NAMES = FALSE
  )
}

.hc_xlsx_overflow_sheet_name <- function(existing) {
  candidate <- "text_overflow"
  suffix <- 1L
  while (base::tolower(candidate) %in% base::tolower(existing)) {
    suffix <- suffix + 1L
    candidate <- base::paste0("text_overflow_", suffix)
  }
  candidate
}

.hc_xlsx_prepare_table <- function(x, sheet, overflow_sheet, overflow) {
  if (!base::is.data.frame(x) && !base::is.matrix(x)) {
    return(x)
  }

  is_matrix <- base::is.matrix(x)
  if (is_matrix) {
    x <- base::as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE)
  }

  for (j in base::seq_along(x)) {
    values <- x[[j]]
    if (base::is.factor(values)) {
      values <- base::as.character(values)
    }
    if (!base::is.character(values)) {
      next
    }

    values <- .hc_xlsx_clean_text(values)
    for (i in base::seq_along(values)) {
      if (base::is.na(values[[i]])) {
        next
      }
      chunks <- .hc_xlsx_split_text(values[[i]])
      if (base::length(chunks) <= 1L) {
        next
      }

      column <- base::names(x)[[j]]
      if (base::is.null(column) || base::is.na(column) || !base::nzchar(column)) {
        column <- base::paste0("column_", j)
      }
      overflow$rows[[base::length(overflow$rows) + 1L]] <- base::data.frame(
        source_sheet = base::rep(sheet, base::length(chunks)),
        source_row = base::rep(i, base::length(chunks)),
        source_column = base::rep(column, base::length(chunks)),
        part = base::seq_along(chunks),
        parts = base::rep(base::length(chunks), base::length(chunks)),
        text = chunks,
        stringsAsFactors = FALSE
      )
      marker <- base::paste0(
        "\n[Full value: ", overflow_sheet, " sheet, ",
        base::length(chunks), " parts]"
      )
      values[[i]] <- base::paste0(chunks[[1]], marker)
    }
    x[[j]] <- values
  }

  base::names(x) <- .hc_xlsx_clean_text(base::names(x))
  if (is_matrix) {
    x <- base::as.matrix(x)
  }
  x
}

.hc_xlsx_prepare_payload <- function(x) {
  if (!base::is.list(x) || base::is.data.frame(x)) {
    overflow <- base::new.env(parent = base::emptyenv())
    overflow$rows <- base::list()
    out <- .hc_xlsx_prepare_table(
      x,
      sheet = "Sheet 1",
      overflow_sheet = "text_overflow",
      overflow = overflow
    )
    if (base::length(overflow$rows) > 0) {
      out <- base::list("Sheet 1" = out)
      out[["text_overflow"]] <- base::do.call(
        base::rbind,
        overflow$rows
      )
    }
    return(out)
  }

  sheet_names <- base::names(x)
  if (base::is.null(sheet_names)) {
    sheet_names <- base::paste0("Sheet ", base::seq_along(x))
  } else {
    missing_names <- base::is.na(sheet_names) | !base::nzchar(sheet_names)
    sheet_names[missing_names] <- base::paste0("Sheet ", base::which(missing_names))
  }
  overflow_sheet <- .hc_xlsx_overflow_sheet_name(sheet_names)
  overflow <- base::new.env(parent = base::emptyenv())
  overflow$rows <- base::list()

  out <- base::lapply(base::seq_along(x), function(i) {
    value <- x[[i]]
    if (base::is.data.frame(value) || base::is.matrix(value)) {
      return(.hc_xlsx_prepare_table(
        value,
        sheet = sheet_names[[i]],
        overflow_sheet = overflow_sheet,
        overflow = overflow
      ))
    }
    if (base::is.character(value) || base::is.factor(value)) {
      return(.hc_xlsx_clean_text(value))
    }
    value
  })
  base::names(out) <- sheet_names

  if (base::length(overflow$rows) > 0) {
    out[[overflow_sheet]] <- base::do.call(base::rbind, overflow$rows)
  }
  out
}

# Atomic, verified file write.
#
# Sync clients (Sciebo / ownCloud / OneDrive) can grab a file mid-write when the
# output folder is a synced directory, producing a truncated or empty result
# (e.g. an unreadable PDF or a zero-row .xlsx). `producer(tmp)` writes to a
# sibling temporary file in the *same* directory; on success the temp file is
# renamed onto `final_path` (an atomic operation on the same filesystem), so a
# consumer never observes a half-written final file. If the finished file is
# missing or smaller than `min_bytes`, an error is raised so the failure is
# surfaced loudly instead of leaving a silently-broken output behind.
.hc_output_payload_valid <- function(path, min_bytes = 1) {
  if (!base::file.exists(path)) {
    return(FALSE)
  }
  size <- base::file.info(path)$size
  if (base::is.na(size) || size < min_bytes) {
    return(FALSE)
  }

  ext <- base::tolower(tools::file_ext(path))
  if (base::identical(ext, "pdf")) {
    signature <- tryCatch(
      readBin(path, what = "raw", n = 5L),
      error = function(e) raw(0)
    )
    return(base::identical(signature, charToRaw("%PDF-")))
  }
  if (base::identical(ext, "png")) {
    signature <- tryCatch(
      readBin(path, what = "raw", n = 8L),
      error = function(e) raw(0)
    )
    return(base::identical(
      signature,
      as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))
    ))
  }
  if (base::identical(ext, "xlsx")) {
    entries <- tryCatch(
      base::suppressWarnings(utils::unzip(path, list = TRUE)$Name),
      error = function(e) character(0)
    )
    return(
      "xl/workbook.xml" %in% entries &&
        base::any(base::grepl("^xl/worksheets/sheet[^/]*\\.xml$", entries)) &&
        .hc_xlsx_xml_parts_valid(path, entries) &&
        .hc_xlsx_package_links_valid(path, entries)
    )
  }

  TRUE
}

.hc_write_atomic <- function(final_path,
                             producer,
                             min_bytes = 1,
                             validator = .hc_output_payload_valid) {
  final_path <- base::as.character(final_path[[1]])
  if (base::is.na(final_path) || !base::nzchar(final_path)) {
    stop("`final_path` must be a non-empty file path.")
  }
  if (!base::is.function(producer)) {
    stop("`producer` must be a function of one argument (the temp path).")
  }
  if (!base::is.function(validator)) {
    stop("`validator` must be a function accepting `path` and `min_bytes`.")
  }
  if (!base::is.numeric(min_bytes) || base::length(min_bytes) != 1L ||
    base::is.na(min_bytes) || !base::is.finite(min_bytes) || min_bytes < 1) {
    stop("`min_bytes` must be a positive finite number.")
  }
  payload_valid <- function(path) {
    isTRUE(tryCatch(
      validator(path = path, min_bytes = min_bytes),
      error = function(e) FALSE
    ))
  }

  dir_path <- base::dirname(final_path)
  if (!base::dir.exists(dir_path)) {
    base::dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }

  # Keep the original extension on the temp file: graphics devices (Cairo),
  # ggplot2::ggsave and openxlsx all infer the output format from the file
  # extension, so the temp name must end in the same `.pdf`/`.png`/`.xlsx`.
  ext <- tools::file_ext(final_path)
  file_ext <- if (base::nzchar(ext)) base::paste0(".", ext) else ""
  tmp_path <- base::tempfile(
    pattern = base::paste0(
      base::basename(tools::file_path_sans_ext(final_path)),
      ".part-"
    ),
    tmpdir = dir_path,
    fileext = file_ext
  )
  backup_path <- NULL
  backup_created <- FALSE
  replacement_started <- FALSE
  committed <- FALSE
  base::on.exit(
    {
      if (base::file.exists(tmp_path)) {
        base::try(base::file.remove(tmp_path), silent = TRUE)
      }
      if (isTRUE(backup_created) && !base::is.null(backup_path) &&
        base::file.exists(backup_path)) {
        if (isTRUE(committed) || !isTRUE(replacement_started)) {
          base::try(base::file.remove(backup_path), silent = TRUE)
        } else {
          if (base::file.exists(final_path)) {
            base::try(base::file.remove(final_path), silent = TRUE)
          }
          restored <- base::suppressWarnings(base::file.rename(backup_path, final_path))
          if (!isTRUE(restored)) {
            restored <- base::file.copy(backup_path, final_path, overwrite = TRUE)
            if (isTRUE(restored)) {
              base::try(base::file.remove(backup_path), silent = TRUE)
            }
          }
          if (!isTRUE(restored)) {
            base::warning(
              "Could not restore the previous output after a failed write. ",
              "The backup remains at '", backup_path, "'.",
              call. = FALSE
            )
          }
        }
      } else if (!isTRUE(committed) && isTRUE(replacement_started) &&
        base::file.exists(final_path)) {
        base::try(base::file.remove(final_path), silent = TRUE)
      }
    },
    add = TRUE
  )

  producer(tmp_path)

  if (!payload_valid(tmp_path)) {
    stop(base::sprintf(
      "Temporary output '%s' was not written or has an invalid payload.", tmp_path
    ))
  }

  if (base::file.exists(final_path)) {
    backup_path <- base::tempfile(
      pattern = base::paste0(
        base::basename(tools::file_path_sans_ext(final_path)),
        ".backup-"
      ),
      tmpdir = dir_path,
      fileext = file_ext
    )
    backup_created <- base::file.copy(final_path, backup_path, overwrite = FALSE)
    if (!isTRUE(backup_created) ||
      !base::file.exists(backup_path) ||
      base::file.info(backup_path)$size != base::file.info(final_path)$size) {
      stop("Could not create a backup of the existing output '", final_path, "'.")
    }
  } else {
    replacement_started <- TRUE
  }

  moved <- base::suppressWarnings(base::file.rename(tmp_path, final_path))
  if (!isTRUE(moved)) {
    if (base::file.exists(final_path)) {
      removed <- base::suppressWarnings(base::file.remove(final_path))
      if (!isTRUE(removed)) {
        stop("Could not replace the existing output '", final_path, "'.")
      }
      replacement_started <- TRUE
    }
    moved <- base::suppressWarnings(base::file.rename(tmp_path, final_path))
    if (!isTRUE(moved)) {
      moved <- base::file.copy(tmp_path, final_path, overwrite = FALSE)
    }
  } else {
    replacement_started <- TRUE
  }
  if (!isTRUE(moved) || !payload_valid(final_path)) {
    stop(base::sprintf(
      "Failed to finalize output '%s' (write to a synced folder such as Sciebo/OneDrive may have been interrupted).",
      final_path
    ))
  }

  committed <- TRUE
  invisible(final_path)
}

.hc_file_matches_reference <- function(path,
                                       reference,
                                       min_bytes = 1,
                                       attempts = 4L,
                                       delay_sec = 0.15) {
  path <- base::as.character(path[[1]])
  reference <- base::as.character(reference[[1]])
  attempts <- base::max(1L, base::as.integer(attempts[[1]]))
  delay_sec <- base::max(0, base::as.numeric(delay_sec[[1]]))

  if (!base::file.exists(reference)) {
    return(FALSE)
  }
  reference_size <- base::file.info(reference)$size
  if (base::is.na(reference_size) || reference_size < min_bytes) {
    return(FALSE)
  }
  reference_md5 <- base::unname(tools::md5sum(reference))
  if (base::length(reference_md5) != 1L ||
    base::is.na(reference_md5) ||
    !base::nzchar(reference_md5)) {
    return(FALSE)
  }

  for (attempt in base::seq_len(attempts)) {
    matches <- tryCatch(
      {
        size <- base::file.info(path)$size
        md5 <- base::unname(tools::md5sum(path))
        base::file.exists(path) &&
          !base::is.na(size) &&
          base::identical(base::as.numeric(size), base::as.numeric(reference_size)) &&
          base::length(md5) == 1L &&
          !base::is.na(md5) &&
          base::identical(md5, reference_md5)
      },
      error = function(e) FALSE
    )
    if (isTRUE(matches)) {
      return(TRUE)
    }
    if (attempt < attempts && delay_sec > 0) {
      base::Sys.sleep(delay_sec * attempt)
    }
  }
  FALSE
}

# Verify a file that should already exist is present and non-empty; warn (do not
# stop) so a failed export becomes visible in the log instead of silent.
.hc_verify_output_file <- function(path, label = NULL, min_bytes = 1) {
  if (base::is.null(path) || base::length(path) == 0) {
    return(invisible(FALSE))
  }
  path <- base::as.character(path[[1]])
  if (base::is.na(path) || !base::nzchar(path)) {
    return(invisible(FALSE))
  }
  ok <- .hc_output_payload_valid(path, min_bytes = min_bytes)
  if (!ok) {
    base::warning(
      "Expected output ",
      if (!base::is.null(label)) label else "",
      if (!base::is.null(label)) " " else "",
      "'", path, "' is missing, empty, or invalid after writing. ",
      "If the save folder is on a synced drive (Sciebo/OneDrive), the sync ",
      "client may have interrupted the write; try a local output folder.",
      call. = FALSE
    )
  }
  invisible(ok)
}

.hc_stage_and_publish_xlsx <- function(file, producer) {
  if (!base::is.function(producer)) {
    stop("`producer` must be a function of one argument (the local XLSX path).")
  }
  local_stage <- base::tempfile(
    pattern = "hc-xlsx-stage-",
    tmpdir = base::tempdir(),
    fileext = ".xlsx"
  )
  base::on.exit(
    {
      if (base::file.exists(local_stage)) {
        base::file.remove(local_stage)
      }
    },
    add = TRUE
  )

  producer(local_stage)
  .hc_xlsx_repair_dangling_parts(local_stage)
  if (!.hc_output_payload_valid(local_stage)) {
    stop(
      "The locally staged XLSX workbook is invalid before publishing to '",
      file,
      "'."
    )
  }

  .hc_write_atomic(
    final_path = file,
    producer = function(tmp) {
      copied <- base::file.copy(local_stage, tmp, overwrite = FALSE)
      if (!isTRUE(copied)) {
        stop("Could not copy the locally staged XLSX workbook to '", tmp, "'.")
      }
    },
    validator = function(path, min_bytes) {
      .hc_file_matches_reference(
        path = path,
        reference = local_stage,
        min_bytes = min_bytes
      )
    }
  )
}

# Atomic drop-in for `openxlsx::write.xlsx(x = ..., file = ..., ...)`: keep the
# call site identical except for the function name. Build and validate the ZIP
# archive on R's local temporary filesystem first. Only then copy the complete,
# byte-verified archive to a temporary sibling of `file` and atomically rename
# it, so a synced or bind-mounted output folder cannot interfere while openxlsx
# is still assembling the workbook.
.hc_write_xlsx_atomic <- function(x, file, ...) {
  dots <- base::list(...)
  x <- .hc_xlsx_prepare_payload(x)
  .hc_stage_and_publish_xlsx(
    file = file,
    producer = function(local_stage) {
      base::do.call(
        openxlsx::write.xlsx,
        base::c(base::list(x = x, file = local_stage), dots)
      )
    }
  )
}

.hc_save_workbook_atomic <- function(wb, file, overwrite = TRUE) {
  if (!isTRUE(overwrite) && base::file.exists(file)) {
    stop("File already exists and `overwrite = FALSE`: '", file, "'.")
  }
  .hc_stage_and_publish_xlsx(
    file = file,
    producer = function(local_stage) {
      openxlsx::saveWorkbook(wb, local_stage, overwrite = TRUE)
    }
  )
}

.hc_export_sanitize_stem <- function(x, default = "page") {
  x <- base::as.character(x[[1]])
  if (base::is.na(x) || !base::nzchar(x)) {
    return(default)
  }
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if (!base::nzchar(x)) {
    x <- default
  }
  x
}

.hc_output_root <- function() {
  out_dir <- tryCatch(hcobject[["working_directory"]][["dir_output"]], error = function(e) NULL)
  save_folder <- tryCatch(hcobject[["global_settings"]][["save_folder"]], error = function(e) NULL)

  if (base::is.null(out_dir) || base::length(out_dir) == 0) {
    stop("No output directory configured in `hcobject`.")
  }
  out_dir <- base::as.character(out_dir[[1]])
  if (base::is.na(out_dir) || !base::nzchar(out_dir)) {
    stop("No output directory configured in `hcobject`.")
  }

  save_folder_value <- ""
  if (!base::is.null(save_folder) && base::length(save_folder) > 0) {
    save_folder_value <- base::as.character(save_folder[[1]])
    if (base::is.na(save_folder_value) || save_folder_value %in% c("", "FALSE", "false")) {
      save_folder_value <- ""
    }
  }

  if (!base::nzchar(save_folder_value)) {
    return(out_dir)
  }
  base::file.path(out_dir, save_folder_value)
}

.hc_output_dir <- function(...) {
  parts <- base::list(...)
  path <- .hc_output_root()
  for (part in parts) {
    if (base::is.null(part) || !base::nzchar(base::as.character(part[[1]]))) {
      next
    }
    path <- base::file.path(path, base::as.character(part[[1]]))
  }
  if (!base::dir.exists(path)) {
    base::dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

.hc_output_file <- function(filename, ..., create_dirs = TRUE) {
  dir_path <- .hc_output_root()
  parts <- base::list(...)
  for (part in parts) {
    if (base::is.null(part) || !base::nzchar(base::as.character(part[[1]]))) {
      next
    }
    dir_path <- base::file.path(dir_path, base::as.character(part[[1]]))
  }
  if (isTRUE(create_dirs) && !base::dir.exists(dir_path)) {
    base::dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  base::file.path(dir_path, filename)
}

.hc_open_pdf_device <- function(file,
                                width,
                                height,
                                pointsize = 11,
                                dpi = NULL,
                                onefile = TRUE,
                                bg = "white") {
  if (requireNamespace("Cairo", quietly = TRUE)) {
    if (!base::is.null(dpi)) {
      Cairo::Cairo(
        file = file,
        width = width,
        height = height,
        pointsize = pointsize,
        dpi = dpi,
        type = "pdf",
        units = "in",
        bg = bg,
        canvas = bg
      )
    } else {
      Cairo::CairoPDF(
        file = file,
        width = width,
        height = height,
        pointsize = pointsize,
        onefile = onefile,
        bg = bg
      )
    }
  } else {
    grDevices::pdf(
      file = file,
      width = width,
      height = height,
      pointsize = pointsize,
      onefile = onefile,
      bg = bg
    )
  }
}

.hc_open_png_device <- function(file,
                                width,
                                height,
                                res = 300,
                                pointsize = 11,
                                bg = "white") {
  grDevices::png(
    filename = file,
    width = width,
    height = height,
    units = "in",
    res = res,
    pointsize = pointsize,
    bg = bg
  )
}

.hc_draw_white_page_background <- function() {
  grid::grid.rect(
    x = 0.5,
    y = 0.5,
    width = 1,
    height = 1,
    gp = grid::gpar(fill = "white", col = NA)
  )
  invisible(NULL)
}

.hc_export_single_page_plot <- function(file,
                                        width,
                                        height,
                                        png_width = NULL,
                                        png_height = NULL,
                                        pointsize = 11,
                                        res = 300,
                                        pdf_dpi = NULL,
                                        draw_fun) {
  if (!base::is.character(file) || base::length(file) != 1 || !base::nzchar(file)) {
    stop("`file` must be a non-empty file path.")
  }
  if (!base::is.function(draw_fun)) {
    stop("`draw_fun` must be a function.")
  }

  pdf_file <- file
  png_file <- .hc_export_path_with_ext(file, "png")
  if (base::is.null(png_width)) {
    png_width <- width
  }
  if (base::is.null(png_height)) {
    png_height <- height
  }

  # Render to a temporary sibling file, then atomically move it onto the final
  # path, so a synced output folder (Sciebo/OneDrive) cannot leave a truncated
  # PDF/PNG behind.
  render_page <- function(target, open_device) {
    .hc_write_atomic(target, function(tmp) {
      open_device(tmp)
      on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
      .hc_draw_white_page_background()
      draw_fun()
      invisible(NULL)
    })
  }

  render_page(pdf_file, function(f) {
    .hc_open_pdf_device(
      file = f,
      width = width,
      height = height,
      pointsize = pointsize,
      dpi = pdf_dpi
    )
  })
  render_page(png_file, function(f) {
    .hc_open_png_device(
      file = f,
      width = png_width,
      height = png_height,
      res = res,
      pointsize = pointsize
    )
  })

  list(
    pdf = pdf_file,
    png = png_file
  )
}

.hc_export_multi_page_plot <- function(file,
                                       page_labels,
                                       width,
                                       height,
                                       pointsize = 11,
                                       res = 300,
                                       draw_page_fun,
                                       display = FALSE) {
  if (!base::is.character(file) || base::length(file) != 1 || !base::nzchar(file)) {
    stop("`file` must be a non-empty file path.")
  }
  if (!base::is.function(draw_page_fun)) {
    stop("`draw_page_fun` must be a function.")
  }

  page_labels <- base::as.character(page_labels)
  if (base::length(page_labels) == 0) {
    return(list(pdf = NULL, png = base::character(0)))
  }

  pdf_file <- file
  png_files <- stats::setNames(
    base::character(base::length(page_labels)),
    page_labels
  )

  # Multi-page PDF: render every page to a temporary sibling, then atomically
  # move it onto the final path (protects synced output folders from truncation).
  .hc_write_atomic(pdf_file, function(tmp) {
    .hc_open_pdf_device(
      file = tmp,
      width = width,
      height = height,
      pointsize = pointsize,
      onefile = TRUE
    )
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    for (idx in base::seq_along(page_labels)) {
      .hc_draw_white_page_background()
      draw_page_fun(idx, page_labels[[idx]])
    }
    invisible(NULL)
  })

  stem <- tools::file_path_sans_ext(pdf_file)
  for (idx in base::seq_along(page_labels)) {
    png_file <- base::sprintf(
      "%s_%03d_%s.png",
      stem,
      idx,
      .hc_export_sanitize_stem(page_labels[[idx]], default = base::sprintf("page_%03d", idx))
    )
    .hc_write_atomic(png_file, function(tmp) {
      .hc_open_png_device(
        file = tmp,
        width = width,
        height = height,
        res = res,
        pointsize = pointsize
      )
      on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
      .hc_draw_white_page_background()
      draw_page_fun(idx, page_labels[[idx]])
      invisible(NULL)
    })
    png_files[[idx]] <- png_file
  }

  # Optionally replay each page on the active graphics device so the figures
  # also appear inline (e.g. under an R Markdown chunk / in a notebook), in
  # addition to the exported files. Gated to contexts where a display target
  # exists, so batch runs do not spawn a stray Rplots.pdf.
  if (isTRUE(display) &&
    (base::interactive() || isTRUE(base::getOption("knitr.in.progress", FALSE)))) {
    for (idx in base::seq_along(page_labels)) {
      draw_page_fun(idx, page_labels[[idx]])
    }
  }

  list(
    pdf = pdf_file,
    png = png_files
  )
}

.hc_export_ggplot_file <- function(file,
                                   plot,
                                   width,
                                   height,
                                   pointsize = 11,
                                   res = 300) {
  if (!inherits(plot, "ggplot")) {
    stop("`plot` must be a ggplot object.")
  }
  .hc_export_single_page_plot(
    file = file,
    width = width,
    height = height,
    pointsize = pointsize,
    res = res,
    draw_fun = function() {
      .hc_display_object(plot)
    }
  )
}

# Save a ggplot to both a (cairo) PDF and a PNG companion in one call.
#
# Drop-in replacement for `ggplot2::ggsave(filename = "...pdf", ...)`: the PDF
# is written with `cairo_pdf` (matching the rest of the package) and a PNG of
# the same dimensions is written next to it at `res` dpi. Extra `...` arguments
# are forwarded to both saves; any `device` is ignored (PDF forces cairo_pdf,
# PNG infers from the `.png` extension). The PNG failing only warns, so a PDF is
# still produced.
.hc_ggsave_pdf_png <- function(filename,
                               plot,
                               width,
                               height,
                               units = "in",
                               res = 300,
                               ...) {
  if (!base::is.character(filename) || base::length(filename) != 1 || !base::nzchar(filename)) {
    stop("`filename` must be a non-empty file path.")
  }
  dots <- base::list(...)
  dots[["device"]] <- NULL
  if (!base::is.null(dots[["dpi"]])) {
    res <- dots[["dpi"]]
    dots[["dpi"]] <- NULL
  }

  pdf_file <- .hc_export_path_with_ext(filename, "pdf")
  png_file <- .hc_export_path_with_ext(filename, "png")

  .hc_write_atomic(pdf_file, function(tmp) {
    pdf_args <- base::c(
      base::list(
        filename = tmp,
        plot = plot,
        width = width,
        height = height,
        units = units,
        device = grDevices::cairo_pdf
      ),
      dots
    )
    if (base::is.null(pdf_args[["bg"]])) {
      pdf_args[["bg"]] <- "white"
    }
    base::do.call(ggplot2::ggsave, pdf_args)
  })

  tryCatch(
    .hc_write_atomic(png_file, function(tmp) {
      png_args <- base::c(
        base::list(
          filename = tmp,
          plot = plot,
          width = width,
          height = height,
          units = units,
          dpi = res
        ),
        dots
      )
      if (base::is.null(png_args[["bg"]])) {
        png_args[["bg"]] <- "white"
      }
      base::do.call(ggplot2::ggsave, png_args)
    }),
    error = function(e) {
      base::warning(
        "Could not write PNG companion ", png_file, ": ",
        base::conditionMessage(e),
        call. = FALSE
      )
    }
  )

  base::invisible(base::list(pdf = pdf_file, png = png_file))
}
