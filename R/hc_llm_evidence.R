# Evidence for hc_llm_enrichment(): statistical results on each set
# (functional enrichment, upstream regulators) and literature passages from a
# DoRAG retrieval server, rendered as prompt blocks.

.hc_llm_evidence_sources <- function() c("enrichment", "upstream", "rag")

# "none" -> character(0); "all" -> every source; otherwise validated subset,
# returned in a fixed order so that the same choice always gives the same
# slot name.
.hc_llm_evidence_levels <- function(evidence) {
  if (base::is.null(evidence)) {
    return(base::character(0))
  }
  evidence <- base::unique(base::tolower(base::trimws(base::as.character(evidence))))
  evidence <- evidence[!base::is.na(evidence) & base::nzchar(evidence)]
  if (base::length(evidence) == 0 || base::identical(evidence, "none")) {
    return(base::character(0))
  }
  if ("all" %in% evidence) {
    return(.hc_llm_evidence_sources())
  }
  bad <- base::setdiff(evidence, .hc_llm_evidence_sources())
  if (base::length(bad) > 0) {
    stop("Unknown `evidence`: ", base::paste(bad, collapse = ", "),
         ". Use \"none\", \"all\" or any of \"enrichment\", \"upstream\", \"rag\".",
         call. = FALSE)
  }
  .hc_llm_evidence_sources()[.hc_llm_evidence_sources() %in% evidence]
}

# Explains the evidence and how to weigh it; placed before the sets.
.hc_llm_evidence_preamble <- function(levels) {
  if (base::length(levels) == 0) {
    return(NULL)
  }
  base::c(
    "Evidence is given for each set. Weigh it as follows:",
    if (any(c("enrichment", "upstream") %in% levels)) {
      "- Enriched terms and upstream regulators are statistical tests on exactly this set (best first). They are the strongest evidence, but terms can be redundant or generic."
    },
    if ("rag" %in% levels) {
      "- Literature passages were retrieved by text similarity and need not concern this set. Use them only where they fit the features."
    },
    "- The features remain the primary evidence. Do not invent support that is not given; if the evidence does not fit, say so in supporting_evidence."
  )
}

# Build the evidence block for every set. Returns a list with `text` and
# `citations` (named character vectors, "" where nothing was found) and
# `problems` (messages for sets where a source failed).
.hc_llm_build_evidence <- function(hc, sets, levels, context_text,
                                   evidence_text = NULL,
                                   top = 10, qval = 0.05,
                                   rag_url = NULL, rag_top_k = 10,
                                   verbose = TRUE) {
  labels <- base::names(sets)
  text <- stats::setNames(base::rep("", base::length(sets)), labels)
  citations <- text
  problems <- base::character(0)
  if (base::length(levels) == 0 && base::is.null(evidence_text)) {
    return(list(text = text, citations = citations, problems = problems))
  }

  blocks <- stats::setNames(lapply(labels, function(x) base::character(0)), labels)
  if (!base::is.null(evidence_text)) {
    if (!base::is.list(evidence_text) && !base::is.character(evidence_text)) {
      stop("`evidence_text` must be a named list or named character vector.", call. = FALSE)
    }
    for (nm in base::intersect(base::names(evidence_text), labels)) {
      val <- base::paste(base::as.character(base::unlist(evidence_text[[nm]])), collapse = "\n")
      if (base::nzchar(val)) blocks[[nm]] <- c(blocks[[nm]], base::paste0("Supplied evidence:\n", val))
    }
  }

  if ("enrichment" %in% levels) {
    terms <- .hc_llm_collect_enrichment_terms(hc, top = top, qval = qval)
    for (nm in labels) {
      txt <- .hc_llm_format_enrichment_context(terms[[nm]])
      blocks[[nm]] <- c(blocks[[nm]], if (base::is.null(txt)) {
        "Enriched terms: none significant for this set."
      } else {
        txt
      })
    }
  }

  if ("upstream" %in% levels) {
    regs <- .hc_llm_collect_upstream(hc, top = top)
    for (nm in labels) {
      txt <- .hc_llm_format_upstream_context(regs[[nm]])
      blocks[[nm]] <- c(blocks[[nm]], if (base::is.null(txt)) {
        "Upstream regulators: none significant for this set."
      } else {
        txt
      })
    }
  }

  if ("rag" %in% levels) {
    url <- .hc_llm_rag_url(rag_url)
    for (nm in labels) {
      query <- .hc_llm_auto_rag_query(label = nm, genes = sets[[nm]], context_text = context_text)
      if (isTRUE(verbose)) message("Literature retrieval for ", nm, " ...")
      rag <- tryCatch(
        .hc_llm_request_rag(query = query, url = url, top_k = rag_top_k),
        error = function(e) e
      )
      if (inherits(rag, "error")) {
        problems <- c(problems, base::paste0(nm, ": ", base::conditionMessage(rag)))
        blocks[[nm]] <- c(blocks[[nm]], "Literature passages: retrieval failed.")
        next
      }
      rag <- .hc_llm_prepare_rag_result(rag)
      txt <- .hc_llm_format_rag_context(rag)
      citations[[nm]] <- .hc_llm_rag_citations_text(rag)
      blocks[[nm]] <- c(blocks[[nm]], if (base::nzchar(txt)) txt else "Literature passages: none retrieved.")
    }
    if (base::length(problems) > 0) {
      warning("Literature retrieval failed for ", base::length(problems), " set(s): ",
              problems[[1]], call. = FALSE)
    }
  }

  for (nm in labels) {
    text[[nm]] <- base::paste(blocks[[nm]], collapse = "\n")
  }
  list(text = text, citations = citations, problems = problems)
}

# ---- functional enrichment ------------------------------------------------

# Significant enrichment terms per module, keyed by module label and colour.
.hc_llm_collect_enrichment_terms <- function(hc, top = 10, qval = 0.05) {
  if (base::is.null(hc)) {
    stop("`evidence = \"enrichment\"` needs `hc` with results of ",
         "hc_functional_enrichment(); for `genes =`, supply `evidence_text`.",
         call. = FALSE)
  }
  enr <- tryCatch(hc_satellite(hc, "enrichments"), error = function(e) NULL)
  tbl <- if (base::is.list(enr)) enr[["significant_enrichments_all_dbs"]] else NULL
  if (base::is.null(tbl) || base::nrow(base::as.data.frame(tbl)) == 0) {
    stop("No significant functional enrichment stored in `hc`. Run ",
         "hc_functional_enrichment() first.", call. = FALSE)
  }
  tbl <- base::as.data.frame(tbl, stringsAsFactors = FALSE)
  if ("qvalue" %in% base::colnames(tbl)) {
    q <- .hc_as_numeric_safely(tbl$qvalue)
    tbl <- tbl[!base::is.na(q) & q <= qval, , drop = FALSE]
    tbl <- tbl[base::order(.hc_as_numeric_safely(tbl$qvalue)), , drop = FALSE]
  }
  out <- list()
  for (col in base::intersect(c("module_label", "cluster"), base::colnames(tbl))) {
    for (m in base::unique(base::as.character(tbl[[col]]))) {
      if (base::is.na(m) || !base::nzchar(m)) next
      df <- tbl[base::as.character(tbl[[col]]) %in% m, , drop = FALSE]
      out[[m]] <- utils::head(df, top)
    }
  }
  out
}

.hc_llm_format_enrichment_context <- function(terms, max_chars = 4000) {
  if (base::is.null(terms) || base::nrow(base::as.data.frame(terms)) == 0) {
    return(NULL)
  }
  df <- base::as.data.frame(terms, stringsAsFactors = FALSE)
  get <- function(nm) {
    if (nm %in% base::colnames(df)) base::as.character(df[[nm]]) else base::rep("", base::nrow(df))
  }
  lines <- base::vapply(base::seq_len(base::nrow(df)), function(i) {
    q <- .hc_as_numeric_safely(get("qvalue")[[i]])
    base::paste0(
      "- ", get("term")[[i]],
      if (base::nzchar(get("database")[[i]])) base::paste0(" [", get("database")[[i]], "]") else "",
      if (base::is.finite(q)) base::paste0("  q=", base::format(q, digits = 2, scientific = TRUE)) else "",
      if (base::nzchar(get("GeneRatio")[[i]])) base::paste0("  genes=", get("GeneRatio")[[i]]) else ""
    )
  }, base::character(1))
  txt <- base::paste(c("Enriched terms (hypergeometric test on this set, best first):", lines),
                     collapse = "\n")
  if (base::nchar(txt) > max_chars) {
    txt <- base::paste0(base::substr(txt, 1, max_chars), "\n[truncated]")
  }
  txt
}

# ---- upstream regulators --------------------------------------------------

# Significant module x regulator links of hc_upstream_inference(), keyed by
# module label and colour; representatives (not redundant) first.
.hc_llm_collect_upstream <- function(hc, top = 10) {
  if (base::is.null(hc)) {
    stop("`evidence = \"upstream\"` needs `hc` with results of ",
         "hc_upstream_inference().", call. = FALSE)
  }
  up <- tryCatch(hc_satellite(hc, "upstream_inference"), error = function(e) NULL)
  tbl <- if (base::is.list(up)) up[["significant_upstream_all"]] else NULL
  if (base::is.null(tbl)) {
    stop("No upstream results stored in `hc`. Run hc_upstream_inference() first.",
         call. = FALSE)
  }
  tbl <- base::as.data.frame(tbl, stringsAsFactors = FALSE)
  if (base::nrow(tbl) == 0) {
    return(list())
  }
  red <- if ("redundant_with" %in% base::colnames(tbl)) base::nzchar(base::as.character(tbl$redundant_with)) else FALSE
  tbl <- tbl[base::order(red, .hc_as_numeric_safely(tbl$qvalue)), , drop = FALSE]
  out <- list()
  for (col in base::intersect(c("module_label", "cluster"), base::colnames(tbl))) {
    for (m in base::unique(base::as.character(tbl[[col]]))) {
      if (base::is.na(m) || !base::nzchar(m)) next
      out[[m]] <- utils::head(tbl[base::as.character(tbl[[col]]) %in% m, , drop = FALSE], top)
    }
  }
  out
}

.hc_llm_format_upstream_context <- function(regs, max_chars = 3000) {
  if (base::is.null(regs) || base::nrow(regs) == 0) {
    return(NULL)
  }
  get <- function(nm) {
    if (nm %in% base::colnames(regs)) base::as.character(regs[[nm]]) else base::rep("", base::nrow(regs))
  }
  lines <- base::vapply(base::seq_len(base::nrow(regs)), function(i) {
    q <- .hc_as_numeric_safely(get("qvalue")[[i]])
    in_mod <- base::identical(base::toupper(get("regulator_in_module")[[i]]), "TRUE")
    base::paste0(
      "- ", get("term")[[i]], " [", get("resource")[[i]],
      if (base::nzchar(get("database")[[i]])) base::paste0(", ", get("database")[[i]]) else "", "]",
      if (base::is.finite(q)) base::paste0("  q=", base::format(q, digits = 2, scientific = TRUE)) else "",
      if (base::nzchar(get("n_overlap")[[i]])) base::paste0("  targets in set=", get("n_overlap")[[i]]) else "",
      if (base::nzchar(get("direction")[[i]])) base::paste0("  ", get("direction")[[i]]) else "",
      if (base::nzchar(get("peak_condition")[[i]])) base::paste0(" (strongest in ", get("peak_condition")[[i]], ")") else "",
      if (in_mod) "  regulator itself is in this set" else ""
    )
  }, base::character(1))
  txt <- base::paste(c(
    "Upstream regulators (over-representation of their targets in this set, best first):",
    lines
  ), collapse = "\n")
  if (base::nchar(txt) > max_chars) {
    txt <- base::paste0(base::substr(txt, 1, max_chars), "\n[truncated]")
  }
  txt
}

# ---- literature (DoRAG) ---------------------------------------------------

.hc_llm_rag_url <- function(rag_url) {
  if (base::is.null(rag_url)) {
    rag_url <- base::Sys.getenv("HCOCENA_RAG_URL", unset = "")
  }
  if (!base::is.character(rag_url) || base::length(rag_url) != 1 || !base::nzchar(rag_url)) {
    stop("`evidence = \"rag\"` needs the retrieval endpoint: pass `rag_url` or set ",
         "the environment variable HCOCENA_RAG_URL.", call. = FALSE)
  }
  rag_url
}

.hc_llm_auto_rag_query <- function(label, genes, context_text, max_genes = 60) {
  genes <- .hc_gemini_normalize_genes(genes)
  parts <- c(
    "Biological function of a co-regulated feature set",
    if (base::nzchar(context_text)) base::paste0("Biological context: ", context_text) else NULL,
    if (base::length(genes) > 0) base::paste0("Features: ", base::paste(utils::head(genes, max_genes), collapse = ", ")) else NULL
  )
  stringr::str_squish(base::paste(parts, collapse = "\n"))
}

.hc_llm_request_rag <- function(query, url, top_k = 10, category = "DoRAG",
                                timeout_sec = 120, connect_timeout_sec = 30) {
  payload <- list(query = query, category = category, top_k = base::as.integer(top_k))
  resp <- tryCatch(
    httr::POST(
      url = url,
      body = payload,
      encode = "json",
      httr::add_headers("Content-Type" = "application/json"),
      httr::timeout(timeout_sec),
      httr::config(connecttimeout = connect_timeout_sec)
    ),
    error = function(e) stop("Literature retrieval failed: ", base::conditionMessage(e), call. = FALSE)
  )
  txt <- base::paste(httr::content(resp, as = "text", encoding = "UTF-8"), collapse = "")
  if (httr::http_error(resp)) {
    stop("Literature retrieval failed with HTTP ", httr::status_code(resp), ".", call. = FALSE)
  }
  out <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE), error = function(e) NULL)
  if (!base::is.list(out)) {
    stop("Literature retrieval returned no JSON object.", call. = FALSE)
  }
  out
}

.hc_llm_prepare_rag_result <- function(rag_result) {
  if (base::is.null(rag_result) || !base::is.list(rag_result)) {
    return(NULL)
  }
  entries <- lapply(.hc_llm_rag_context_list(rag_result[["context"]]), .hc_llm_normalize_rag_context_entry)
  entries <- entries[base::vapply(entries, function(x) base::nzchar(x$chunk), base::logical(1))]
  list(
    context = entries,
    cited_papers = .hc_llm_unique_rag_papers(lapply(entries, function(x) x$paper))
  )
}

.hc_llm_rag_context_list <- function(x) {
  if (base::is.null(x)) {
    return(list())
  }
  if (base::is.data.frame(x)) {
    return(lapply(base::seq_len(base::nrow(x)), function(i) base::as.list(x[i, , drop = FALSE])))
  }
  if (base::is.list(x) && any(base::names(x) %in% c("chunk", "section", "relevance", "paper"))) {
    return(list(x))
  }
  if (base::is.list(x)) x else list(x)
}

.hc_llm_normalize_rag_context_entry <- function(entry) {
  if (!base::is.list(entry)) {
    entry <- list(chunk = entry)
  }
  list(
    chunk = .hc_llm_rag_scalar(entry[["chunk"]]),
    section = .hc_llm_rag_scalar(entry[["section"]]),
    relevance = suppressWarnings(base::as.numeric(.hc_llm_rag_scalar(entry[["relevance"]], NA))),
    paper = .hc_llm_normalize_rag_paper(entry[["paper"]])
  )
}

.hc_llm_normalize_rag_paper <- function(paper) {
  if (base::is.data.frame(paper)) {
    paper <- if (base::nrow(paper) > 0) base::as.list(paper[1, , drop = FALSE]) else list()
  }
  if (base::is.null(paper)) paper <- list()
  if (!base::is.list(paper)) paper <- list(apa_citation = paper)
  list(
    title = .hc_llm_rag_scalar(paper[["title"]]),
    apa_citation = .hc_llm_rag_scalar(paper[["apa_citation"]]),
    doi = .hc_llm_rag_scalar(paper[["doi"]])
  )
}

.hc_llm_rag_scalar <- function(x, default = "") {
  if (base::is.null(x) || base::length(x) == 0) {
    return(default)
  }
  if (base::is.list(x)) x <- x[[1]]
  val <- base::as.character(x[[1]])
  if (base::length(val) == 0 || base::is.na(val)) default else base::trimws(val)
}

.hc_llm_unique_rag_papers <- function(papers) {
  if (base::length(papers) == 0) {
    return(list())
  }
  keys <- base::vapply(papers, function(p) base::paste(p$apa_citation, p$title, p$doi, sep = "\r"), base::character(1))
  has <- base::vapply(papers, function(p) any(base::nzchar(c(p$apa_citation, p$title, p$doi))), base::logical(1))
  papers[!base::duplicated(keys) & has]
}

.hc_llm_format_rag_context <- function(rag_result, max_chars = 8000) {
  entries <- if (base::is.null(rag_result)) list() else rag_result$context
  if (base::length(entries) == 0) {
    return("")
  }
  lines <- base::vapply(base::seq_along(entries), function(i) {
    e <- entries[[i]]
    title <- if (base::nzchar(e$paper$title)) e$paper$title else e$paper$apa_citation
    base::paste0("[", i, "] ", if (base::nzchar(title)) base::paste0(title, ": ") else "", e$chunk)
  }, base::character(1))
  txt <- base::paste(c("Literature passages (retrieved by similarity, may not concern this set):", lines),
                     collapse = "\n")
  if (base::nchar(txt) > max_chars) {
    txt <- base::paste0(base::substr(txt, 1, max_chars), "\n[truncated]")
  }
  txt
}

.hc_llm_rag_citations_text <- function(rag_result) {
  papers <- if (base::is.null(rag_result)) list() else rag_result$cited_papers
  if (base::length(papers) == 0) {
    return("")
  }
  cit <- base::vapply(papers, function(p) {
    out <- if (base::nzchar(p$apa_citation)) p$apa_citation else p$title
    if (base::nzchar(p$doi)) base::paste0(out, " DOI: ", p$doi) else out
  }, base::character(1))
  base::paste(base::unique(cit[base::nzchar(cit)]), collapse = " | ")
}
