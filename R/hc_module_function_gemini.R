#' AI-assisted module function summary via Gemini, Claude, ChatGPT, or vLLM
#'
#' Uses a large language model API to summarize the likely overarching
#' biological function of a module or free gene set. This is an AI-assisted
#' interpretation step, not a statistical enrichment test.
#'
#' The function can use genes from an `HCoCenaExperiment` module or a manually
#' supplied gene vector. Optional free-text context can be added to steer the
#' interpretation towards a disease, tissue, perturbation, timepoint, or other
#' experimental background.
#'
#' By default, results are stored in `hc@satellite[[slot_name]][[label]]` when
#' an `hc` object is provided. A compact table is also written to
#' `hc@satellite[[paste0(slot_name, "_summary")]]`.
#'
#' The function name is retained for backward compatibility, but it can now use
#' Gemini, OpenAI-compatible APIs, and local vLLM servers via `llm`.
#'
#' @param hc Optional `HCoCenaExperiment`. Required when `module` is used or
#'   when `save_to_hc = TRUE`.
#' @param module Optional character vector naming one or more modules. Use
#'   `"all"` to summarize all modules present in
#'   `hc@satellite$module_gene_list`.
#' @param genes Optional character vector of gene symbols. Use this instead of
#'   `module`.
#' @param context Optional character string or vector with biological context
#'   for backward compatibility.
#' @param biological_context Optional character string or vector with biological
#'   context for the interpretation. Preferred over `context`.
#' @param label Optional label used for storing the result. Defaults to the
#'   module name or `"custom_geneset"`. Can only be used for a single module or
#'   a free gene set.
#' @param provider Provider-neutral alias for `llm`. When supplied it takes
#'   precedence over `llm`. Accepts the same values as `llm`.
#' @param llm LLM provider. Supported values are `"gemini"`, `"claude"`,
#'   `"anthropic"` (alias for `"claude"`), `"openai"`, `"chatgpt"` (alias for
#'   `"openai"`), or `"vllm"` for a local OpenAI-compatible vLLM server.
#' @param api_key API key for the selected provider. If omitted, the function
#'   uses `GEMINI_API_KEY` for Gemini, `ANTHROPIC_API_KEY` for Claude,
#'   `OPENAI_API_KEY` for OpenAI, and `VLLM_API_KEY` for vLLM. For local vLLM,
#'   the fallback is `"EMPTY"`.
#' @param gemini_model Optional Gemini model code. Only used when
#'   `llm = "gemini"`. Defaults to `"gemini-2.5-pro"`.
#' @param claude_model Optional Claude model code. Only used when
#'   `llm = "claude"`. Defaults to `"claude-sonnet-4-6"`.
#' @param vllm_model Optional vLLM model code. Only used when `llm = "vllm"`.
#'   Defaults to `"Qwen/Qwen2.5-VL-32B-Instruct"`.
#' @param vllm_base_url Base URL for the local OpenAI-compatible vLLM server.
#'   Only used when `llm = "vllm"`. Defaults to `"http://localhost:8000/v1"`.
#' @param base_url Provider-neutral alias for the local server endpoint. When
#'   supplied with `llm = "vllm"` it overrides `vllm_base_url`; for cloud
#'   providers it is ignored (with a warning).
#' @param model Optional generic model code. Mainly useful for OpenAI, or as a
#'   provider-agnostic override.
#' @param max_genes Maximum number of genes sent to the API. Use `Inf` or
#'   `NULL` to send all genes. Default is `Inf`.
#' @param temperature Generation temperature. Default is `0.2`.
#' @param timeout_sec Request timeout in seconds. Use `0` or `NULL` to disable
#'   the timeout. Default is `60` for hosted APIs and `300` for local `vllm`
#'   requests unless overridden explicitly.
#' @param pause_sec Pause in seconds between module requests. Useful for
#'   `module = "all"`. Use `0` or `NULL` for no pause. Default is `0`.
#' @param use_rag Logical. If `TRUE`, retrieves literature passages from the
#'   DoRAG raw retrieval API and injects them into the LLM prompt as optional
#'   supporting context. Default is `FALSE`.
#' @param rag_query Optional retrieval query. If `NULL`, a query is built from
#'   the biological context, module label, and submitted genes. A single string
#'   is reused for all modules; a named character vector can map module labels
#'   to queries; a character vector with one entry per module is used in order.
#'   Advanced users may pass a function with arguments `label`, `module`,
#'   `genes`, and `context_text`.
#' @param rag_url DoRAG raw retrieval API endpoint. If `NULL`, uses
#'   `HCOCENA_RAG_URL` or
#'   `"https://limesbcnr-007901.iaas.uni-bonn.de/api/rag/query"`.
#' @param rag_category Retrieval category sent to DoRAG. Default is `"DoRAG"`.
#' @param rag_top_k Number of passages requested after reranking. Default is
#'   `10`.
#' @param rag_themes Optional biological theme filter passed to DoRAG, for
#'   example `"Microbiome Science"` or `c("Microbiome Science", "Innate Immunity")`.
#' @param rag_timeout_sec DoRAG request timeout in seconds. Use `0` or `NULL`
#'   to disable the timeout. Default is `120`.
#' @param rag_connect_timeout_sec DoRAG connection timeout in seconds. This
#'   controls how long to wait while establishing the TCP connection to the
#'   RAG server. Use `0` or `NULL` to use the system default. Default is `30`.
#' @param rag_min_relevance Optional minimum rerank score. Retrieved passages
#'   below this score are omitted from the prompt after retrieval.
#' @param rag_max_context_chars Maximum number of characters of formatted DoRAG
#'   context injected into each LLM prompt. Use `Inf` or `NULL` to disable this
#'   limit. Default is `12000`.
#' @param rag_continue_on_error Logical. If `TRUE`, a failed DoRAG request is
#'   stored in the result and the LLM interpretation continues without RAG
#'   passages. Default is `FALSE`.
#' @param continue_on_error Logical. If `TRUE`, continue with the next module
#'   when one request fails and store the error in the result summary.
#' @param save_to_hc Logical. If `TRUE`, store the result in `hc@satellite` and
#'   return updated `hc`. Default is `TRUE` when `hc` is provided.
#' @param use_enrichment Logical. If `TRUE`, a further interpretation is run in
#'   which the module's statistically significant functional-enrichment terms
#'   are put into the prompt alongside the genes, the biological context and -
#'   when `use_rag = TRUE` - the retrieved literature. This is the most
#'   grounded of the levels: the enrichment terms are a hypergeometric test on
#'   exactly the gene list being interpreted, whereas the RAG passages are
#'   background retrieved by similarity and need not concern this module. The
#'   prompt states that ranking explicitly. Requires `hc` with results from
#'   [hc_functional_enrichment()], and is stored in `enrichment_response`
#'   next to `response` (genes only) and `rag_response` (genes + literature),
#'   so the three can be compared. Default is `FALSE`.
#' @param enrichment_top Maximum number of enrichment terms per module put into
#'   the prompt, best q-value first. Default is `10`.
#' @param enrichment_qval Significance cutoff applied to the stored enrichment
#'   terms before they are used. Default is `0.05`.
#' @param enrichment_terms Optional terms supplied directly, bypassing the
#'   lookup in `hc`. Accepts a character vector, a data frame with `term` and
#'   `qvalue` columns, or a list of either named by module. This is the only
#'   way to combine enrichment grounding with `genes =`, where there is no
#'   module to look up.
#' @param enrichment_max_context_chars Maximum number of characters of
#'   formatted enrichment terms injected into the prompt. Default is `4000`.
#' @param slot_name Satellite slot name used for storage. Default is
#'   `"llm_module_function"`.
#' @param system_instruction Optional custom system instruction for the model.
#' @param verbose Logical. If `TRUE`, print short progress messages.
#' @param ... Used by the backward-compatible wrapper aliases
#'   `hc_module_function_gemini()` and `hc_module_function_vllm()`.
#'
#' @return Updated `HCoCenaExperiment` if `save_to_hc = TRUE`; otherwise either
#'   a single result list or, for multiple modules, a list with `results` and
#'   `summary`.
#' @export
#'
#' @examples
#' \donttest{
#' if (nzchar(Sys.getenv("OPENAI_API_KEY"))) {
#'   res <- hc_module_function_llm(
#'     genes = c("STAT1", "IRF7", "CXCL10", "GBP1", "IFI44L"),
#'     context = "Interferon-driven blood module in acute viral infection",
#'     llm = "openai",
#'     api_key = Sys.getenv("OPENAI_API_KEY"),
#'     save_to_hc = FALSE
#'   )
#' }
#' }
hc_module_function_llm <- function(hc = NULL,
                                   module = NULL,
                                   genes = NULL,
                                   context = NULL,
                                   biological_context = NULL,
                                   label = NULL,
                                   provider = NULL,
                                   llm = c("gemini", "claude", "openai", "chatgpt", "vllm"),
                                   api_key = NULL,
                                   gemini_model = NULL,
                                   claude_model = NULL,
                                   vllm_model = NULL,
                                   vllm_base_url = NULL,
                                   base_url = NULL,
                                   model = NULL,
                                   max_genes = Inf,
                                   temperature = 0.2,
                                   timeout_sec = 60,
                                   pause_sec = 0,
                                   use_rag = FALSE,
                                   rag_query = NULL,
                                   rag_url = NULL,
                                   rag_category = "DoRAG",
                                   rag_top_k = 10,
                                   rag_themes = NULL,
                                   rag_timeout_sec = 120,
                                   rag_connect_timeout_sec = 30,
                                   rag_min_relevance = NULL,
                                   rag_max_context_chars = 12000,
                                   rag_continue_on_error = FALSE,
                                   use_enrichment = FALSE,
                                   enrichment_top = 10,
                                   enrichment_qval = 0.05,
                                   enrichment_terms = NULL,
                                   enrichment_max_context_chars = 4000,
                                   continue_on_error = FALSE,
                                   save_to_hc = !base::is.null(hc),
                                   slot_name = "llm_module_function",
                                   system_instruction = NULL,
                                   verbose = TRUE) {
  timeout_was_missing <- missing(timeout_sec)

  if (!base::is.null(module) && !base::is.null(genes)) {
    stop("Use either `module` or `genes`, not both.")
  }
  if (base::is.null(module) && base::is.null(genes)) {
    stop("Provide either `module` or `genes`.")
  }
  if (isTRUE(save_to_hc) && base::is.null(hc)) {
    stop("`hc` must be provided when `save_to_hc = TRUE`.")
  }

  # `provider` is a sprechender alias for `llm`; when supplied it wins.
  if (!base::is.null(provider) && base::nzchar(base::as.character(provider[[1]]))) {
    llm <- provider
  }
  llm <- base::tolower(base::as.character(llm[[1]]))
  if (!(llm %in% c("gemini", "claude", "anthropic", "openai", "chatgpt", "vllm"))) {
    stop("`llm` must be one of `gemini`, `claude`, `anthropic`, `openai`, `chatgpt`, or `vllm`.")
  }
  if (llm == "chatgpt") {
    llm <- "openai"
  } else if (llm == "anthropic") {
    llm <- "claude"
  }

  if (llm == "gemini" &&
    !base::is.null(gemini_model) &&
    base::nzchar(base::as.character(gemini_model[[1]]))) {
    model <- base::as.character(gemini_model[[1]])
  } else if (llm == "claude" &&
    !base::is.null(claude_model) &&
    base::nzchar(base::as.character(claude_model[[1]]))) {
    model <- base::as.character(claude_model[[1]])
  } else if (llm == "vllm" &&
    !base::is.null(vllm_model) &&
    base::nzchar(base::as.character(vllm_model[[1]]))) {
    model <- base::as.character(vllm_model[[1]])
  } else if (base::is.null(model) || !base::nzchar(base::as.character(model[[1]]))) {
    model <- if (llm == "gemini") {
      "gemini-2.5-pro"
    } else if (llm == "claude") {
      "claude-sonnet-4-6"
    } else if (llm == "vllm") {
      "Qwen/Qwen2.5-VL-32B-Instruct"
    } else {
      "gpt-4o-mini"
    }
  } else {
    model <- base::as.character(model[[1]])
  }
  if (!base::is.character(model) || base::length(model) != 1 || !base::nzchar(model)) {
    stop("`model` must be a non-empty character scalar.")
  }
  if (!base::is.null(max_genes)) {
    if (!base::is.numeric(max_genes) || base::length(max_genes) != 1 ||
      base::is.na(max_genes[[1]]) || max_genes[[1]] <= 0 ||
      (base::is.finite(max_genes[[1]]) && max_genes[[1]] != base::floor(max_genes[[1]]))) {
      stop("`max_genes` must be NULL, Inf, or a single positive integer.")
    }
  }
  if (!base::is.numeric(temperature) || base::length(temperature) != 1 || !base::is.finite(temperature)) {
    stop("`temperature` must be a single finite numeric value.")
  }
  if (!base::is.null(timeout_sec)) {
    if (!base::is.numeric(timeout_sec) || base::length(timeout_sec) != 1 ||
      !base::is.finite(timeout_sec) || timeout_sec < 0) {
      stop("`timeout_sec` must be NULL or a single finite numeric value >= 0.")
    }
  }
  if (!base::is.null(pause_sec)) {
    if (!base::is.numeric(pause_sec) || base::length(pause_sec) != 1 || !base::is.finite(pause_sec) || pause_sec < 0) {
      stop("`pause_sec` must be NULL or a single non-negative number.")
    }
  }
  continue_flag_label <- "`continue_on_error`"
  if (!base::is.logical(continue_on_error) || base::length(continue_on_error) != 1 || base::is.na(continue_on_error)) {
    stop(continue_flag_label, " must be TRUE or FALSE.")
  }
  if (!base::is.logical(verbose) || base::length(verbose) != 1 || base::is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.")
  }
  # Fourth interpretation level: significant enrichment terms for the module.
  # `enrichment_terms` lets a caller supply them directly, which is the only
  # way this works together with `genes =` (there is no module to look up).
  enrichment_lookup <- NULL
  if (!base::is.null(enrichment_terms)) {
    enrichment_lookup <- if (base::is.character(enrichment_terms)) {
      base::paste(enrichment_terms, collapse = "
")
    } else if (base::is.data.frame(enrichment_terms)) {
      .hc_llm_format_enrichment_context(
        terms = enrichment_terms,
        max_chars = enrichment_max_context_chars
      )
    } else {
      base::lapply(enrichment_terms, function(x) {
        if (base::is.character(x)) {
          base::paste(x, collapse = "
")
        } else {
          .hc_llm_format_enrichment_context(
            terms = x, max_chars = enrichment_max_context_chars
          )
        }
      })
    }
  } else if (isTRUE(use_enrichment)) {
    per_module <- .hc_llm_collect_enrichment_terms(
      hc = hc, top = enrichment_top, qval = enrichment_qval
    )
    enrichment_lookup <- base::lapply(per_module, function(df) {
      .hc_llm_format_enrichment_context(
        terms = df, max_chars = enrichment_max_context_chars
      )
    })
    if (isTRUE(verbose)) {
      message(
        "Enrichment grounding: ", base::length(per_module),
        " module keys with significant terms (q <= ", enrichment_qval, ")."
      )
    }
  }

  rag_options <- .hc_llm_resolve_rag_options(
    use_rag = use_rag,
    rag_query = rag_query,
    rag_url = rag_url,
    rag_category = rag_category,
    rag_top_k = rag_top_k,
    rag_themes = rag_themes,
    rag_timeout_sec = rag_timeout_sec,
    rag_connect_timeout_sec = rag_connect_timeout_sec,
    rag_min_relevance = rag_min_relevance,
    rag_max_context_chars = rag_max_context_chars,
    rag_continue_on_error = rag_continue_on_error
  )

  api_key <- .hc_llm_resolve_api_key(api_key = api_key, llm = llm)
  # `base_url` is a provider-neutral alias for the local server endpoint; it
  # currently maps onto the OpenAI-compatible `vllm` provider.
  if (!base::is.null(base_url) && base::nzchar(base::as.character(base_url[[1]]))) {
    if (llm == "vllm") {
      vllm_base_url <- base_url
    } else {
      warning(
        "`base_url` currently applies only to local `vllm` models; ignoring it for provider `",
        llm, "`.",
        call. = FALSE
      )
    }
  }
  vllm_base_url <- .hc_llm_resolve_vllm_base_url(vllm_base_url = vllm_base_url, llm = llm)
  timeout_sec <- .hc_llm_resolve_timeout(
    timeout_sec = timeout_sec,
    llm = llm,
    timeout_was_missing = timeout_was_missing
  )
  pause_sec <- .hc_llm_normalize_pause(pause_sec)

  if (base::is.null(system_instruction) || !base::nzchar(base::as.character(system_instruction[[1]]))) {
    system_instruction <- if (llm == "vllm") {
      paste(
        "You are a careful transcriptomics analyst.",
        "Infer the main biological program of the provided gene module.",
        "Use the biological context as background information, but do not simply repeat it.",
        "Be concise and module-specific.",
        "For contextual_state, describe the immediate transcriptional or cellular state, not the study design or cohort.",
        "Avoid repeating generic context words such as immature, maturation, preterm, infant, or monocyte unless they are essential and specifically supported by the genes.",
        "Do not write sentence starters like 'this module' or 'the genes'.",
        "Use the supplied genes and optional context only as evidence.",
        "Do not claim statistical enrichment was performed.",
        "Return valid JSON only."
      )
    } else {
      paste(
        "You are a careful transcriptomics analyst.",
        "Infer the main biological program of the provided gene module.",
        "Use the biological context only as an interpretation frame, not as text to repeat.",
        "Prioritize module-specific biology over generic cohort-level wording.",
        "Avoid generic answers such as monocyte maturation, immune differentiation, or immune activation unless the genes strongly support nothing more specific.",
        "Prefer concrete programs such as interferon signaling, antigen presentation, cell cycle, erythroid/megakaryocytic bias, platelet program, mitochondrial metabolism, phagolysosome, inflammatory signaling, ribosome biogenesis, stress response, chemotaxis, or tissue contamination if clearly supported.",
        "Make different modules distinguishable from each other.",
        "Return compact labels and compact descriptions only.",
        "Do not write sentence starters like 'this module', 'this program', or 'is best described as'.",
        "Use the supplied genes and optional context only as evidence.",
        "Do not claim statistical enrichment was performed.",
        "State uncertainty when the signal is weak or mixed.",
        "Return valid JSON only."
      )
    }
  } else {
    system_instruction <- base::as.character(system_instruction[[1]])
  }

  if (base::is.null(module)) {
    gene_infos <- list(
      list(
        label = "custom_geneset",
        module = NULL,
        genes = .hc_gemini_normalize_genes(genes)
      )
    )
  } else {
    if (base::is.null(hc)) {
      stop("`hc` must be provided when `module` is used.")
    }
    modules_to_run <- .hc_llm_resolve_modules(hc = hc, module = module)
    gene_infos <- lapply(modules_to_run, function(x) .hc_gemini_get_module_genes(hc = hc, module = x))
  }

  gene_infos <- gene_infos[!base::vapply(gene_infos, base::is.null, FUN.VALUE = base::logical(1))]
  if (base::length(gene_infos) == 0) {
    stop("No valid module or gene inputs were resolved.")
  }
  if (base::length(gene_infos) > 1 && !base::is.null(label) && base::nzchar(base::as.character(label[[1]]))) {
    stop("`label` can only be used for a single module or gene set.")
  }

  biological_context <- if (base::is.null(biological_context) || !base::length(biological_context)) {
    context
  } else {
    biological_context
  }
  context_text <- .hc_gemini_normalize_context(biological_context)
  results <- .hc_llm_run_sequential(
    gene_infos = gene_infos,
    label = label,
    context_text = context_text,
    llm = llm,
    api_key = api_key,
    model = model,
    vllm_base_url = vllm_base_url,
    max_genes = max_genes,
    temperature = temperature,
    timeout_sec = timeout_sec,
    pause_sec = pause_sec,
    continue_on_error = continue_on_error,
    system_instruction = system_instruction,
    response_schema = .hc_llm_response_schema(),
    rag_options = rag_options,
    enrichment_lookup = enrichment_lookup,
    verbose = verbose
  )

  if (!isTRUE(save_to_hc)) {
    if (base::length(results) == 1) {
      return(results[[1]])
    }
    return(list(
      results = stats::setNames(results, base::vapply(results, function(x) x$label, FUN.VALUE = base::character(1))),
      summary = .hc_llm_summary_from_results(results = results, hc = hc)
    ))
  }

  sat <- tryCatch(base::as.list(hc@satellite), error = function(e) list())
  slot_obj <- sat[[slot_name]]
  if (base::is.null(slot_obj) || !base::is.list(slot_obj)) {
    slot_obj <- list()
  }
  for (res in results) {
    slot_obj[[res$label]] <- res
  }
  sat[[slot_name]] <- slot_obj
  summary_tbl <- .hc_llm_summary_from_results(results = slot_obj, hc = hc)
  sat[[base::paste0(slot_name, "_summary")]] <- summary_tbl
  hc@satellite <- S4Vectors::SimpleList(sat)
  .hc_llm_export_results_excel(
    hc = hc,
    results = slot_obj,
    summary_tbl = summary_tbl,
    slot_name = slot_name
  )
  hc
}

#' @rdname hc_module_function_llm
#' @export
hc_module_function_gemini <- function(...) {
  hc_module_function_llm(...)
}

#' @rdname hc_module_function_llm
#' @export
hc_module_function_vllm <- function(...) {
  hc_module_function_llm(llm = "vllm", ...)
}

#' Pick the enrichment prompt block belonging to one module
#' @noRd
.hc_llm_enrichment_text_for <- function(lookup, gene_info, label) {
  if (base::is.null(lookup) || base::length(lookup) == 0) {
    return(NULL)
  }
  if (base::is.character(lookup) && base::length(lookup) == 1) {
    return(lookup)
  }
  keys <- base::c(gene_info$module, gene_info$label, label)
  keys <- base::as.character(keys)
  keys <- keys[!base::is.na(keys) & base::nzchar(keys)]
  for (k in keys) {
    if (!base::is.null(lookup[[k]])) {
      return(lookup[[k]])
    }
  }
  NULL
}

.hc_llm_run_sequential <- function(gene_infos,
                                   label,
                                   context_text,
                                   llm,
                                   api_key,
                                   model,
                                   vllm_base_url,
                                   max_genes,
                                   temperature,
                                   timeout_sec,
                                   pause_sec,
                                   continue_on_error,
                                   system_instruction,
                                   response_schema,
                                   rag_options,
                                   enrichment_lookup = NULL,
                                   verbose) {
  results <- vector("list", length = base::length(gene_infos))
  for (i in base::seq_along(gene_infos)) {
    gene_info <- gene_infos[[i]]
    this_label <- if (base::is.null(label) || !base::nzchar(base::as.character(label[[1]]))) {
      gene_info$label
    } else {
      base::as.character(label[[1]])
    }
    results[[i]] <- tryCatch(
      .hc_llm_summarize_one(
        gene_info = gene_info,
        label = this_label,
        context_text = context_text,
        llm = llm,
        api_key = api_key,
        model = model,
        vllm_base_url = vllm_base_url,
        max_genes = max_genes,
        temperature = temperature,
        timeout_sec = timeout_sec,
        system_instruction = system_instruction,
        response_schema = response_schema,
        rag_options = rag_options,
        enrichment_text = .hc_llm_enrichment_text_for(
          lookup = enrichment_lookup,
          gene_info = gene_info,
          label = this_label
        ),
        index = i,
        n_inputs = base::length(gene_infos),
        verbose = verbose
      ),
      error = function(e) {
        if (!isTRUE(continue_on_error)) {
          stop(e)
        }
        if (isTRUE(verbose)) {
          message("LLM module summary failed for `", this_label, "`: ", base::conditionMessage(e))
        }
        .hc_llm_error_result(
          gene_info = gene_info,
          label = this_label,
          context_text = context_text,
          llm = llm,
          model = model,
          error_message = base::conditionMessage(e)
        )
      }
    )
    if (pause_sec > 0 && i < base::length(gene_infos)) {
      if (isTRUE(verbose)) {
        message("Waiting ", pause_sec, " seconds before next module request.")
      }
      Sys.sleep(pause_sec)
    }
  }
  results
}

.hc_llm_summarize_many <- function(gene_infos,
                                   context_text,
                                   llm,
                                   api_key,
                                   model,
                                   vllm_base_url,
                                   max_genes,
                                   temperature,
                                   timeout_sec,
                                   system_instruction,
                                   verbose) {
  batch_inputs <- lapply(gene_infos, function(gene_info) {
    genes_all <- .hc_gemini_normalize_genes(gene_info$genes)
    genes_use <- .hc_llm_limit_genes(genes = genes_all, max_genes = max_genes)
    list(
      label = gene_info$label,
      module = gene_info$module,
      genes_all = genes_all,
      genes_use = genes_use,
      truncated = base::length(genes_use) < base::length(genes_all)
    )
  })

  if (isTRUE(verbose)) {
    message(
      "LLM module summary: sending ",
      base::length(batch_inputs),
      " modules in one combined request using provider `",
      llm,
      "` and model `",
      model,
      "`."
    )
  }

  prompt <- .hc_llm_build_batch_prompt(batch_inputs = batch_inputs, context_text = context_text)
  response_schema <- .hc_llm_batch_response_schema(n_modules = base::length(batch_inputs))

  req_result <- if (llm == "gemini") {
    .hc_llm_request_gemini(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      response_schema = response_schema,
      temperature = temperature,
      timeout_sec = timeout_sec
    )
  } else if (llm == "claude") {
    .hc_llm_request_claude(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      temperature = temperature,
      timeout_sec = timeout_sec
    )
  } else if (llm == "vllm") {
    .hc_llm_request_vllm(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      temperature = temperature,
      timeout_sec = timeout_sec,
      base_url = vllm_base_url
    )
  } else {
    .hc_llm_request_openai(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      response_schema = response_schema,
      temperature = temperature,
      timeout_sec = timeout_sec
    )
  }

  parsed_result <- tryCatch(
    jsonlite::fromJSON(req_result$result_text, simplifyVector = FALSE),
    error = function(e) {
      stop("Could not parse combined LLM response JSON: ", base::conditionMessage(e), call. = FALSE)
    }
  )

  module_entries <- parsed_result[["modules"]]
  if (base::is.null(module_entries) || !base::is.list(module_entries) || base::length(module_entries) == 0) {
    stop("Combined LLM response did not contain a `modules` list.", call. = FALSE)
  }

  entry_labels <- base::vapply(
    module_entries,
    function(x) {
      lbl <- x[["label"]]
      if (base::is.null(lbl)) "" else base::as.character(lbl[[1]])
    },
    FUN.VALUE = base::character(1)
  )

  out <- lapply(batch_inputs, function(inp) {
    idx <- base::match(inp$label, entry_labels)
    if (base::is.na(idx)) {
      return(.hc_llm_error_result(
        gene_info = list(module = inp$module, genes = inp$genes_all),
        label = inp$label,
        context_text = context_text,
        llm = llm,
        model = model,
        error_message = "Combined LLM response did not return an entry for this module."
      ))
    }
    entry <- module_entries[[idx]]
    list(
      label = inp$label,
      module = inp$module,
      genes_input = inp$genes_all,
      genes_sent = inp$genes_use,
      gene_count_input = base::length(inp$genes_all),
      gene_count_sent = base::length(inp$genes_use),
      truncated = inp$truncated,
      context = if (base::nzchar(context_text)) context_text else NULL,
      llm = llm,
      model = model,
      status = "ok",
      error_message = NA_character_,
      prompt = prompt,
      response = entry,
      response_text = jsonlite::toJSON(entry, auto_unbox = TRUE),
      raw_response_text = req_result$raw_response_text,
      timestamp = base::as.character(Sys.time())
    )
  })

  out
}

.hc_llm_summarize_one <- function(gene_info,
                                  label,
                                  context_text,
                                  llm,
                                  api_key,
                                  model,
                                  vllm_base_url,
                                  max_genes,
                                  temperature,
                                  timeout_sec,
                                  system_instruction,
                                  response_schema,
                                  rag_options = NULL,
                                  enrichment_text = NULL,
                                  index = 1L,
                                  n_inputs = 1L,
                                  verbose) {
  genes_all <- .hc_gemini_normalize_genes(gene_info$genes)
  if (base::length(genes_all) == 0) {
    stop("No valid genes available for LLM interpretation.")
  }

  genes_use <- .hc_llm_limit_genes(genes = genes_all, max_genes = max_genes)
  truncated <- base::length(genes_use) < base::length(genes_all)
  rag_result <- NULL
  rag_query_used <- NULL
  rag_context_text <- NULL
  rag_error_message <- NA_character_

  if (isTRUE(verbose)) {
    if (truncated) {
      message(
        "LLM module summary: sending ",
        base::length(genes_use),
        " of ", base::length(genes_all),
        " genes for `", label, "` using provider `", llm,
        "` and model `", model, "`."
      )
    } else {
      message(
        "LLM module summary: sending ",
        base::length(genes_use),
        " genes for `", label, "` using provider `", llm,
        "` and model `", model, "`."
      )
    }
  }

  if (!base::is.null(rag_options)) {
    rag_query_used <- .hc_llm_resolve_rag_query(
      rag_query = rag_options$query,
      label = label,
      module = gene_info$module,
      genes = genes_use,
      context_text = context_text,
      index = index,
      n_inputs = n_inputs
    )
    if (isTRUE(verbose)) {
      message(
        "DoRAG retrieval: requesting ",
        rag_options$top_k,
        " passages for `", label, "`."
      )
    }
    rag_raw <- tryCatch(
      .hc_llm_request_rag(
        query = rag_query_used,
        url = rag_options$url,
        category = rag_options$category,
        top_k = rag_options$top_k,
        themes = rag_options$themes,
        timeout_sec = rag_options$timeout_sec,
        connect_timeout_sec = rag_options$connect_timeout_sec
      ),
      error = function(e) {
        if (!isTRUE(rag_options$continue_on_error)) {
          stop(e)
        }
        e
      }
    )
    if (inherits(rag_raw, "error")) {
      rag_error_message <- base::conditionMessage(rag_raw)
      if (isTRUE(verbose)) {
        message(
          "DoRAG retrieval failed for `",
          label,
          "`; continuing without RAG passages: ",
          rag_error_message
        )
      }
      rag_raw <- NULL
    }
    if (base::is.null(rag_raw)) {
      rag_result <- list(
        query = rag_query_used,
        answer = "",
        context = list(),
        cited_papers = list(),
        status = "error",
        error_message = rag_error_message
      )
    } else {
      rag_result <- .hc_llm_prepare_rag_result(
        rag_result = rag_raw,
        min_relevance = rag_options$min_relevance
      )
      rag_context_text <- .hc_llm_format_rag_context(
        rag_result = rag_result,
        max_chars = rag_options$max_context_chars
      )
      rag_result$status <- "ok"
      rag_result$error_message <- NA_character_
    }
    if (isTRUE(verbose)) {
      message(
        "DoRAG retrieval: using ",
        .hc_llm_rag_context_count(rag_result),
        " retrieved passages for `", label, "`."
      )
    }
  }

  run_interpretation <- function(level, interpretation_rag_context = NULL,
                                 interpretation_enrichment = NULL) {
    .hc_llm_interpretation_level(
      level = level,
      label = label,
      genes = genes_use,
      total_gene_count = base::length(genes_all),
      context_text = context_text,
      truncated = truncated,
      rag_context_text = interpretation_rag_context,
      enrichment_text = interpretation_enrichment,
      llm = llm,
      api_key = api_key,
      model = model,
      vllm_base_url = vllm_base_url,
      temperature = temperature,
      timeout_sec = timeout_sec,
      system_instruction = system_instruction,
      response_schema = response_schema
    )
  }

  has_rag_context <- !base::is.null(rag_context_text) &&
    base::nzchar(base::as.character(rag_context_text[[1]]))

  baseline <- run_interpretation(level = "baseline")
  rag_interpretation <- NULL
  if (has_rag_context) {
    if (isTRUE(verbose)) {
      message("LLM module summary: running the additional RAG interpretation for `", label, "`.")
    }
    rag_interpretation <- run_interpretation(
      level = "rag",
      interpretation_rag_context = rag_context_text
    )
  }

  # Fourth level: genes + context + retrieved literature + the statistically
  # significant enrichment terms for this module. RAG is optional here - when
  # it was not requested this level is genes + context + enrichment.
  has_enrichment_text <- !base::is.null(enrichment_text) &&
    base::nzchar(base::as.character(enrichment_text[[1]]))
  enrichment_interpretation <- NULL
  if (has_enrichment_text) {
    if (isTRUE(verbose)) {
      message(
        "LLM module summary: running the enrichment-grounded interpretation for `",
        label, "`", if (has_rag_context) " (with RAG passages)." else ".", ""
      )
    }
    enrichment_interpretation <- run_interpretation(
      level = "enrichment",
      interpretation_rag_context = if (has_rag_context) rag_context_text else NULL,
      interpretation_enrichment = enrichment_text
    )
  }

  out <- list(
    label = label,
    module = gene_info$module,
    genes_input = genes_all,
    genes_sent = genes_use,
    gene_count_input = base::length(genes_all),
    gene_count_sent = base::length(genes_use),
    truncated = truncated,
    context = if (base::nzchar(context_text)) context_text else NULL,
    llm = llm,
    model = model,
    status = "ok",
    error_message = NA_character_,
    prompt = baseline$prompt,
    response = baseline$response,
    response_text = baseline$response_text,
    raw_response_text = baseline$raw_response_text,
    rag_prompt = if (!base::is.null(rag_interpretation)) rag_interpretation$prompt else NULL,
    rag_response = if (!base::is.null(rag_interpretation)) rag_interpretation$response else NULL,
    rag_response_text = if (!base::is.null(rag_interpretation)) rag_interpretation$response_text else NULL,
    rag_raw_response_text = if (!base::is.null(rag_interpretation)) rag_interpretation$raw_response_text else NULL,
    rag = rag_result,
    rag_query = rag_query_used,
    rag_error_message = rag_error_message,
    rag_context_text = if (!base::is.null(rag_context_text) && base::nzchar(rag_context_text)) rag_context_text else NULL,
    enrichment_prompt = if (!base::is.null(enrichment_interpretation)) enrichment_interpretation$prompt else NULL,
    enrichment_response = if (!base::is.null(enrichment_interpretation)) enrichment_interpretation$response else NULL,
    enrichment_response_text = if (!base::is.null(enrichment_interpretation)) enrichment_interpretation$response_text else NULL,
    enrichment_raw_response_text = if (!base::is.null(enrichment_interpretation)) enrichment_interpretation$raw_response_text else NULL,
    enrichment_context_text = if (has_enrichment_text) base::as.character(enrichment_text[[1]]) else NULL,
    enrichment_used_rag = has_enrichment_text && has_rag_context,
    timestamp = base::as.character(Sys.time())
  )
  out
}

.hc_llm_interpretation_level <- function(level,
                                         label,
                                         genes,
                                         total_gene_count,
                                         context_text,
                                         truncated,
                                         rag_context_text,
                                         enrichment_text = NULL,
                                         llm,
                                         api_key,
                                         model,
                                         vllm_base_url,
                                         temperature,
                                         timeout_sec,
                                         system_instruction,
                                         response_schema) {
  prompt <- .hc_gemini_build_prompt(
    label = label,
    genes = genes,
    total_gene_count = total_gene_count,
    context_text = context_text,
    truncated = truncated,
    rag_context_text = rag_context_text,
    enrichment_text = enrichment_text,
    llm = llm
  )

  req_result <- .hc_llm_request_by_provider(
    llm = llm,
    api_key = api_key,
    model = model,
    prompt = prompt,
    system_instruction = system_instruction,
    response_schema = response_schema,
    temperature = temperature,
    timeout_sec = timeout_sec,
    vllm_base_url = vllm_base_url
  )
  result_text <- req_result$result_text

  list(
    level = level,
    context = if (base::nzchar(context_text)) context_text else NULL,
    rag_used = !base::is.null(rag_context_text) &&
      base::nzchar(base::as.character(rag_context_text[[1]])),
    rag_context_text = if (!base::is.null(rag_context_text) &&
      base::nzchar(base::as.character(rag_context_text[[1]]))) {
      base::as.character(rag_context_text[[1]])
    } else {
      NULL
    },
    enrichment_used = !base::is.null(enrichment_text) &&
      base::nzchar(base::as.character(enrichment_text[[1]])),
    prompt = prompt,
    response = .hc_llm_parse_module_response(result_text),
    response_text = result_text,
    raw_response_text = req_result$raw_response_text
  )
}

.hc_llm_request_by_provider <- function(llm,
                                        api_key,
                                        model,
                                        prompt,
                                        system_instruction,
                                        response_schema,
                                        temperature,
                                        timeout_sec,
                                        vllm_base_url) {
  if (llm == "gemini") {
    .hc_llm_request_gemini(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      response_schema = response_schema,
      temperature = temperature,
      timeout_sec = timeout_sec
    )
  } else if (llm == "claude") {
    .hc_llm_request_claude(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      temperature = temperature,
      timeout_sec = timeout_sec
    )
  } else if (llm == "vllm") {
    .hc_llm_request_vllm(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      temperature = temperature,
      timeout_sec = timeout_sec,
      base_url = vllm_base_url
    )
  } else {
    .hc_llm_request_openai(
      api_key = api_key,
      model = model,
      prompt = prompt,
      system_instruction = system_instruction,
      response_schema = response_schema,
      temperature = temperature,
      timeout_sec = timeout_sec
    )
  }
}

.hc_llm_parse_module_response <- function(result_text) {
  tryCatch(
    jsonlite::fromJSON(result_text, simplifyVector = TRUE),
    error = function(e) {
      list(
        general_processes = NA_character_,
        contextual_state = NA_character_,
        key_regulators = result_text,
        parse_error = base::conditionMessage(e)
      )
    }
  )
}

.hc_llm_limit_genes <- function(genes, max_genes) {
  if (base::is.null(max_genes) || !base::is.finite(max_genes[[1]])) {
    return(genes)
  }
  utils::head(genes, base::as.integer(max_genes[[1]]))
}

.hc_llm_default_rag_url <- function() {
  Sys.getenv("HCOCENA_RAG_URL", unset = "https://limesbcnr-007901.iaas.uni-bonn.de/api/rag/query")
}

.hc_llm_resolve_rag_options <- function(use_rag,
                                        rag_query,
                                        rag_url,
                                        rag_category,
                                        rag_top_k,
                                        rag_themes,
                                        rag_timeout_sec,
                                        rag_connect_timeout_sec,
                                        rag_min_relevance,
                                        rag_max_context_chars,
                                        rag_continue_on_error) {
  if (!base::is.logical(use_rag) || base::length(use_rag) != 1 || base::is.na(use_rag)) {
    stop("`use_rag` must be TRUE or FALSE.")
  }
  if (!isTRUE(use_rag)) {
    return(NULL)
  }

  if (base::is.null(rag_url) || !base::nzchar(base::as.character(rag_url[[1]]))) {
    rag_url <- .hc_llm_default_rag_url()
  } else {
    rag_url <- base::as.character(rag_url[[1]])
  }
  rag_url <- base::trimws(rag_url)
  if (!base::nzchar(rag_url)) {
    stop("`rag_url` must be a non-empty URL when `use_rag = TRUE`.")
  }

  if (base::is.null(rag_category) || !base::nzchar(base::as.character(rag_category[[1]]))) {
    stop("`rag_category` must be a non-empty character scalar.")
  }
  rag_category <- base::trimws(base::as.character(rag_category[[1]]))
  if (!base::nzchar(rag_category)) {
    stop("`rag_category` must be a non-empty character scalar.")
  }

  if (!base::is.numeric(rag_top_k) || base::length(rag_top_k) != 1 ||
    base::is.na(rag_top_k[[1]]) || !base::is.finite(rag_top_k[[1]]) ||
    rag_top_k[[1]] <= 0 || rag_top_k[[1]] != base::floor(rag_top_k[[1]])) {
    stop("`rag_top_k` must be a single positive integer.")
  }
  rag_top_k <- base::as.integer(rag_top_k[[1]])

  if (!base::is.null(rag_timeout_sec)) {
    if (!base::is.numeric(rag_timeout_sec) || base::length(rag_timeout_sec) != 1 ||
      base::is.na(rag_timeout_sec[[1]]) || !base::is.finite(rag_timeout_sec[[1]]) ||
      rag_timeout_sec[[1]] < 0) {
      stop("`rag_timeout_sec` must be NULL or a single finite numeric value >= 0.")
    }
  }
  rag_timeout_sec <- .hc_llm_normalize_timeout(rag_timeout_sec)

  if (!base::is.null(rag_connect_timeout_sec)) {
    if (!base::is.numeric(rag_connect_timeout_sec) || base::length(rag_connect_timeout_sec) != 1 ||
      base::is.na(rag_connect_timeout_sec[[1]]) || !base::is.finite(rag_connect_timeout_sec[[1]]) ||
      rag_connect_timeout_sec[[1]] < 0) {
      stop("`rag_connect_timeout_sec` must be NULL or a single finite numeric value >= 0.")
    }
  }
  rag_connect_timeout_sec <- .hc_llm_normalize_timeout(rag_connect_timeout_sec)

  if (!base::is.null(rag_min_relevance)) {
    if (!base::is.numeric(rag_min_relevance) || base::length(rag_min_relevance) != 1 ||
      base::is.na(rag_min_relevance[[1]]) || !base::is.finite(rag_min_relevance[[1]])) {
      stop("`rag_min_relevance` must be NULL or a single finite numeric value.")
    }
    rag_min_relevance <- base::as.numeric(rag_min_relevance[[1]])
  }

  if (!base::is.null(rag_max_context_chars)) {
    if (!base::is.numeric(rag_max_context_chars) || base::length(rag_max_context_chars) != 1 ||
      base::is.na(rag_max_context_chars[[1]]) || rag_max_context_chars[[1]] <= 0) {
      stop("`rag_max_context_chars` must be NULL, Inf, or a single positive number.")
    }
    if (!base::is.finite(rag_max_context_chars[[1]])) {
      rag_max_context_chars <- NULL
    } else {
      rag_max_context_chars <- base::as.integer(rag_max_context_chars[[1]])
    }
  }

  if (!base::is.null(rag_themes)) {
    rag_themes <- base::as.character(rag_themes)
    rag_themes <- base::trimws(rag_themes)
    rag_themes <- base::unique(rag_themes[!base::is.na(rag_themes) & base::nzchar(rag_themes)])
    if (base::length(rag_themes) == 0) {
      rag_themes <- NULL
    }
  }

  if (!base::is.logical(rag_continue_on_error) ||
    base::length(rag_continue_on_error) != 1 ||
    base::is.na(rag_continue_on_error)) {
    stop("`rag_continue_on_error` must be TRUE or FALSE.")
  }

  if (!base::is.null(rag_query) && !base::is.function(rag_query)) {
    query_names <- base::names(rag_query)
    rag_query <- base::as.character(rag_query)
    rag_query <- base::trimws(rag_query)
    keep <- !base::is.na(rag_query) & base::nzchar(rag_query)
    if (base::length(query_names) == base::length(rag_query)) {
      query_names <- base::as.character(query_names)
      query_names[base::is.na(query_names)] <- ""
      query_names <- query_names[keep]
    } else {
      query_names <- NULL
    }
    rag_query <- rag_query[keep]
    if (base::length(rag_query) == 0) {
      stop("`rag_query` must contain at least one non-empty query when provided.")
    }
    if (!base::is.null(query_names)) {
      base::names(rag_query) <- query_names
    }
  }

  list(
    query = rag_query,
    url = rag_url,
    category = rag_category,
    top_k = rag_top_k,
    themes = rag_themes,
    timeout_sec = rag_timeout_sec,
    connect_timeout_sec = rag_connect_timeout_sec,
    min_relevance = rag_min_relevance,
    max_context_chars = rag_max_context_chars,
    continue_on_error = rag_continue_on_error
  )
}

.hc_llm_resolve_rag_query <- function(rag_query,
                                      label,
                                      module,
                                      genes,
                                      context_text,
                                      index,
                                      n_inputs) {
  if (base::is.null(rag_query)) {
    return(.hc_llm_auto_rag_query(label = label, genes = genes, context_text = context_text))
  }

  if (base::is.function(rag_query)) {
    query <- rag_query(
      label = label,
      module = module,
      genes = genes,
      context_text = context_text
    )
    return(.hc_llm_normalize_rag_query_scalar(query, context = "`rag_query` function result"))
  }

  query_names <- base::names(rag_query)
  if (!base::is.null(query_names) && base::length(query_names) == base::length(rag_query)) {
    query_names <- base::as.character(query_names)
    idx <- base::match(base::as.character(label[[1]]), query_names)
    if (base::is.na(idx) && !base::is.null(module)) {
      idx <- base::match(base::as.character(module[[1]]), query_names)
    }
    if (!base::is.na(idx)) {
      return(.hc_llm_normalize_rag_query_scalar(rag_query[[idx]], context = "`rag_query`"))
    }
  }

  if (base::length(rag_query) == 1) {
    return(.hc_llm_normalize_rag_query_scalar(rag_query[[1]], context = "`rag_query`"))
  }

  if (base::length(rag_query) == n_inputs) {
    return(.hc_llm_normalize_rag_query_scalar(rag_query[[index]], context = "`rag_query`"))
  }

  stop(
    "`rag_query` must be NULL, a single string, a named vector matching module labels, ",
    "or a vector with one query per module.",
    call. = FALSE
  )
}

.hc_llm_auto_rag_query <- function(label, genes, context_text, max_genes = 60) {
  genes <- .hc_gemini_normalize_genes(genes)
  gene_text <- if (base::length(genes) > 0) {
    base::paste(utils::head(genes, max_genes), collapse = ", ")
  } else {
    ""
  }
  parts <- base::c(
    "Transcriptomic module biological function",
    if (base::nzchar(context_text)) base::paste0("Biological context: ", context_text) else NULL,
    if (!base::is.null(label) && base::nzchar(base::as.character(label[[1]]))) base::paste0("Module label: ", base::as.character(label[[1]])) else NULL,
    if (base::nzchar(gene_text)) base::paste0("Genes: ", gene_text) else NULL
  )
  stringr::str_squish(base::paste(parts, collapse = "\n"))
}

.hc_llm_normalize_rag_query_scalar <- function(query, context = "`rag_query`") {
  if (base::is.null(query) || base::length(query) == 0) {
    stop(context, " must resolve to a non-empty character scalar.", call. = FALSE)
  }
  query <- base::as.character(query[[1]])
  query <- base::trimws(query)
  if (base::is.na(query) || !base::nzchar(query)) {
    stop(context, " must resolve to a non-empty character scalar.", call. = FALSE)
  }
  query
}

#' Pull the significant enrichment terms per module out of an hc object
#'
#' Returns a named list keyed by both module label ("M3") and module colour, so
#' the caller can look a module up either way. Each entry is a data frame of the
#' terms that passed `qval`, best first.
#' @noRd
.hc_llm_collect_enrichment_terms <- function(hc, top = 10, qval = 0.05) {
  if (base::is.null(hc)) {
    stop(
      "`use_enrichment = TRUE` needs an `hc` object. Pass `hc` together with ",
      "`module`, or supply the terms yourself via `enrichment_terms`.",
      call. = FALSE
    )
  }
  sat <- tryCatch(base::as.list(hc@satellite), error = function(e) list())
  enr <- sat[["enrichments"]]
  tbl <- if (base::is.list(enr)) enr[["significant_enrichments_all_dbs"]] else NULL
  if (base::is.null(tbl) || base::nrow(base::as.data.frame(tbl)) == 0) {
    stop(
      "No significant functional enrichment found in `hc@satellite$enrichments`. ",
      "Run `hc_functional_enrichment()` first, or lower `qval`.",
      call. = FALSE
    )
  }
  tbl <- base::as.data.frame(tbl, stringsAsFactors = FALSE)

  if ("qvalue" %in% base::colnames(tbl)) {
    keep <- .hc_as_numeric_safely(tbl$qvalue)
    tbl <- tbl[!base::is.na(keep) & keep <= qval, , drop = FALSE]
    tbl <- tbl[base::order(.hc_as_numeric_safely(tbl$qvalue)), , drop = FALSE]
  }
  if (base::nrow(tbl) == 0) {
    stop("No enrichment term passed `enrichment_qval = ", qval, "`.", call. = FALSE)
  }

  out <- list()
  for (col in base::intersect(base::c("module_label", "cluster"), base::colnames(tbl))) {
    for (m in base::unique(base::as.character(tbl[[col]]))) {
      df <- tbl[base::as.character(tbl[[col]]) %in% m, , drop = FALSE]
      if (base::nrow(df) > top) df <- df[base::seq_len(top), , drop = FALSE]
      if (!base::is.na(m) && base::nzchar(m)) out[[m]] <- df
    }
  }
  out
}

#' Render enrichment terms as a prompt block
#' @noRd
.hc_llm_format_enrichment_context <- function(terms, max_chars = 4000) {
  if (base::is.null(terms) || base::nrow(base::as.data.frame(terms)) == 0) {
    return(NULL)
  }
  df <- base::as.data.frame(terms, stringsAsFactors = FALSE)
  get <- function(nm) if (nm %in% base::colnames(df)) base::as.character(df[[nm]]) else base::rep("", base::nrow(df))
  lines <- base::vapply(base::seq_len(base::nrow(df)), function(i) {
    q <- .hc_as_numeric_safely(get("qvalue")[[i]])
    base::paste0(
      "- ", get("term")[[i]],
      if (base::nzchar(get("database")[[i]])) base::paste0(" [", get("database")[[i]], "]") else "",
      if (base::is.finite(q)) base::paste0("  q=", base::format(q, digits = 2, scientific = TRUE)) else "",
      if (base::nzchar(get("GeneRatio")[[i]])) base::paste0("  genes=", get("GeneRatio")[[i]]) else ""
    )
  }, character(1))

  txt <- base::paste(
    base::paste0(
      "Statistically significant over-representation for exactly this gene list ",
      "(hypergeometric test, q-value cutoff applied), best first:"
    ),
    base::paste(lines, collapse = "\n"),
    sep = "\n"
  )
  if (base::nchar(txt) > max_chars) {
    txt <- base::paste0(base::substr(txt, 1, max_chars), "\n[truncated]")
  }
  txt
}

.hc_llm_request_rag <- function(query,
                                url,
                                category,
                                top_k,
                                themes,
                                timeout_sec,
                                connect_timeout_sec = NULL) {
  if (!base::requireNamespace("httr", quietly = TRUE)) {
    stop("Package `httr` is required for DoRAG retrieval.", call. = FALSE)
  }

  payload <- list(
    query = query,
    category = category,
    top_k = base::as.integer(top_k)
  )
  if (!base::is.null(themes) && base::length(themes) > 0) {
    payload$themes <- themes
  }

  configs <- list(httr::add_headers("Content-Type" = "application/json"))
  if (!base::is.null(timeout_sec)) {
    configs <- base::c(configs, list(httr::timeout(base::as.numeric(timeout_sec[[1]]))))
  }
  if (!base::is.null(connect_timeout_sec)) {
    configs <- base::c(configs, list(httr::config(connecttimeout = base::as.numeric(connect_timeout_sec[[1]]))))
  }

  resp <- tryCatch(
    do.call(
      httr::POST,
      base::c(
        list(
          url = url,
          body = payload,
          encode = "json"
        ),
        configs
      )
    ),
    error = function(e) {
      stop("DoRAG retrieval request failed: ", base::conditionMessage(e), call. = FALSE)
    }
  )

  txt <- tryCatch(
    httr::content(resp, as = "text", encoding = "UTF-8"),
    error = function(e) {
      stop("Could not read DoRAG retrieval response: ", base::conditionMessage(e), call. = FALSE)
    }
  )
  txt <- base::paste(base::as.character(txt), collapse = "")

  if (httr::http_error(resp)) {
    stop(
      "DoRAG retrieval failed with HTTP ",
      httr::status_code(resp),
      if (base::nzchar(txt)) ": " else "",
      if (base::nzchar(txt)) stringr::str_trunc(stringr::str_squish(txt), 300) else "",
      call. = FALSE
    )
  }
  if (!base::nzchar(txt)) {
    stop("DoRAG retrieval response was empty.", call. = FALSE)
  }

  out <- tryCatch(
    jsonlite::fromJSON(txt, simplifyVector = FALSE),
    error = function(e) {
      stop("Could not parse DoRAG retrieval JSON: ", base::conditionMessage(e), call. = FALSE)
    }
  )
  if (!base::is.list(out)) {
    stop("DoRAG retrieval response did not contain a JSON object.", call. = FALSE)
  }
  out$raw_response_text <- txt
  out
}

.hc_llm_prepare_rag_result <- function(rag_result, min_relevance = NULL) {
  if (base::is.null(rag_result) || !base::is.list(rag_result)) {
    return(NULL)
  }

  context_entries <- .hc_llm_rag_context_list(rag_result[["context"]])
  context_entries <- lapply(context_entries, .hc_llm_normalize_rag_context_entry)
  context_entries <- context_entries[base::vapply(
    context_entries,
    function(x) base::nzchar(x$chunk),
    FUN.VALUE = base::logical(1)
  )]

  if (!base::is.null(min_relevance)) {
    context_entries <- context_entries[base::vapply(
      context_entries,
      function(x) !base::is.na(x$relevance) && x$relevance >= min_relevance,
      FUN.VALUE = base::logical(1)
    )]
  }

  cited_papers <- .hc_llm_unique_rag_papers(lapply(context_entries, function(x) x$paper))
  if (base::length(cited_papers) == 0) {
    cited_papers <- .hc_llm_unique_rag_papers(
      lapply(.hc_llm_rag_context_list(rag_result[["cited_papers"]]), .hc_llm_normalize_rag_paper)
    )
  }

  out <- list(
    query = .hc_llm_rag_scalar(rag_result[["query"]]),
    answer = .hc_llm_rag_scalar(rag_result[["answer"]]),
    context = context_entries,
    cited_papers = cited_papers
  )
  if (!base::is.null(rag_result[["raw_response_text"]])) {
    out$raw_response_text <- .hc_llm_rag_scalar(rag_result[["raw_response_text"]])
  }
  out
}

.hc_llm_rag_context_list <- function(x) {
  if (base::is.null(x)) {
    return(list())
  }
  if (base::is.data.frame(x)) {
    return(lapply(base::seq_len(base::nrow(x)), function(i) {
      row <- base::as.list(x[i, , drop = FALSE])
      lapply(row, function(v) {
        if (base::length(v) == 1) v[[1]] else v
      })
    }))
  }
  if (base::is.list(x) && base::length(x) > 0 &&
    any(base::names(x) %in% c("chunk", "section", "relevance", "paper", "apa_citation", "title", "doi"))) {
    return(list(x))
  }
  if (base::is.list(x)) {
    return(x)
  }
  list(x)
}

.hc_llm_normalize_rag_context_entry <- function(entry) {
  if (!base::is.list(entry)) {
    entry <- list(chunk = entry)
  }
  paper <- .hc_llm_normalize_rag_paper(entry[["paper"]])
  list(
    chunk = .hc_llm_rag_scalar(entry[["chunk"]]),
    section = .hc_llm_rag_scalar(entry[["section"]]),
    relevance = .hc_llm_rag_numeric(entry[["relevance"]]),
    paper = paper
  )
}

.hc_llm_normalize_rag_paper <- function(paper) {
  if (base::is.null(paper)) {
    paper <- list()
  }
  if (base::is.data.frame(paper)) {
    paper <- if (base::nrow(paper) > 0) {
      base::as.list(paper[1, , drop = FALSE])
    } else {
      list()
    }
  }
  if (!base::is.list(paper)) {
    paper <- list(apa_citation = paper)
  }
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
  if (base::is.data.frame(x)) {
    if (base::nrow(x) == 0 || base::ncol(x) == 0) {
      return(default)
    }
    x <- x[[1]][[1]]
  } else if (base::is.list(x) && base::length(x) == 1 && !base::is.list(x[[1]])) {
    x <- x[[1]]
  }
  val <- base::as.character(x[[1]])
  if (base::length(val) == 0 || base::is.na(val)) {
    return(default)
  }
  base::trimws(val)
}

.hc_llm_rag_numeric <- function(x) {
  if (base::is.null(x) || base::length(x) == 0) {
    return(NA_real_)
  }
  val <- suppressWarnings(base::as.numeric(x[[1]]))
  if (base::length(val) == 0 || base::is.na(val)) {
    return(NA_real_)
  }
  val[[1]]
}

.hc_llm_unique_rag_papers <- function(papers) {
  if (base::is.null(papers) || base::length(papers) == 0) {
    return(list())
  }
  papers <- lapply(papers, .hc_llm_normalize_rag_paper)
  keys <- base::vapply(
    papers,
    function(p) base::paste(p$apa_citation, p$title, p$doi, sep = "\r"),
    FUN.VALUE = base::character(1)
  )
  has_value <- base::vapply(
    papers,
    function(p) base::any(base::nzchar(base::c(p$apa_citation, p$title, p$doi))),
    FUN.VALUE = base::logical(1)
  )
  keep <- !base::duplicated(keys) & has_value
  papers[keep]
}

.hc_llm_format_rag_context <- function(rag_result, max_chars = 12000) {
  entries <- .hc_llm_rag_context_list(rag_result[["context"]])
  if (base::length(entries) == 0) {
    return("")
  }

  formatted_entries <- base::vapply(
    base::seq_along(entries),
    function(i) {
      entry <- entries[[i]]
      rel <- if (base::is.na(entry$relevance)) {
        "NA"
      } else {
        base::format(base::round(entry$relevance, 3), nsmall = 3, trim = TRUE)
      }
      paper <- entry$paper
      meta <- base::c(
        base::paste0("[", i, "] relevance=", rel),
        if (base::nzchar(entry$section)) base::paste0("section=", entry$section) else NULL,
        if (base::nzchar(paper$title)) base::paste0("title=", paper$title) else NULL
      )
      base::paste(
        base::paste(meta, collapse = " | "),
        if (base::nzchar(paper$apa_citation)) base::paste0("Citation: ", paper$apa_citation) else NULL,
        if (base::nzchar(paper$doi)) base::paste0("DOI: ", paper$doi) else NULL,
        base::paste0("Passage: ", entry$chunk),
        sep = "\n"
      )
    },
    FUN.VALUE = base::character(1)
  )

  citations <- .hc_llm_rag_citations_text(rag_result)
  out <- base::paste(
    "Use these passages only as supporting literature evidence. The gene list remains the primary evidence.",
    base::paste(formatted_entries, collapse = "\n\n"),
    if (base::nzchar(citations)) base::paste0("Unique cited papers: ", citations) else NULL,
    sep = "\n\n"
  )

  if (!base::is.null(max_chars) && base::nchar(out, type = "chars") > max_chars) {
    out <- base::paste0(
      base::substr(out, 1L, max_chars),
      "\n[DoRAG context truncated to ", max_chars, " characters.]"
    )
  }
  out
}

.hc_llm_rag_context_count <- function(rag_result) {
  base::length(.hc_llm_rag_context_list(rag_result[["context"]]))
}

.hc_llm_rag_citations_text <- function(rag_result) {
  papers <- .hc_llm_rag_context_list(rag_result[["cited_papers"]])
  if (base::length(papers) == 0) {
    return("")
  }
  citations <- base::vapply(
    papers,
    function(p) {
      p <- .hc_llm_normalize_rag_paper(p)
      out <- if (base::nzchar(p$apa_citation)) p$apa_citation else p$title
      if (base::nzchar(p$doi)) {
        out <- base::paste0(out, " DOI: ", p$doi)
      }
      out
    },
    FUN.VALUE = base::character(1)
  )
  citations <- base::unique(citations[base::nzchar(citations)])
  base::paste(citations, collapse = " | ")
}

.hc_llm_error_result <- function(gene_info,
                                 label,
                                 context_text,
                                 llm,
                                 model,
                                 error_message) {
  genes_all <- .hc_gemini_normalize_genes(gene_info$genes)
  list(
    label = label,
    module = gene_info$module,
    genes_input = genes_all,
    genes_sent = base::character(0),
    gene_count_input = base::length(genes_all),
    gene_count_sent = 0L,
    truncated = FALSE,
    context = if (base::nzchar(context_text)) context_text else NULL,
    llm = llm,
    model = model,
    status = "error",
    error_message = error_message,
    prompt = NULL,
    response = list(
      general_processes = NA_character_,
      contextual_state = NA_character_,
      key_regulators = error_message
    ),
    response_text = NA_character_,
    raw_response_text = NA_character_,
    timestamp = base::as.character(Sys.time())
  )
}

.hc_llm_natural_module_order <- function(modules) {
  modules <- base::unique(base::as.character(modules))
  modules <- modules[!base::is.na(modules) & base::nzchar(modules)]
  if (base::length(modules) <= 1) {
    return(modules)
  }

  parsed <- lapply(modules, function(x) {
    hit <- regexec("^([^0-9]*?)([0-9]+)([^0-9]*)$", x, perl = TRUE)
    parts <- regmatches(x, hit)[[1]]
    if (base::length(parts) == 4) {
      return(list(
        prefix = base::tolower(parts[[2]]),
        number = .hc_first_numeric_value(.hc_as_integer_safely(parts[[3]])),
        suffix = base::tolower(parts[[4]])
      ))
    }
    list(
      prefix = base::tolower(x),
      number = Inf,
      suffix = ""
    )
  })

  ord <- base::order(
    base::vapply(parsed, function(x) x$prefix, FUN.VALUE = base::character(1)),
    base::vapply(parsed, function(x) x$number, FUN.VALUE = base::numeric(1)),
    base::vapply(parsed, function(x) x$suffix, FUN.VALUE = base::character(1)),
    base::tolower(modules),
    modules
  )
  modules[ord]
}

.hc_llm_is_quota_error <- function(msg) {
  msg <- base::tolower(base::as.character(msg[[1]]))
  any(base::grepl(
    pattern = c("quota exceeded", "resource_exhausted", "rate limit", "too many requests", "retry in"),
    x = msg,
    fixed = TRUE
  ))
}

.hc_llm_resolve_modules <- function(hc, module) {
  module <- base::as.character(module)
  module <- base::trimws(module)
  module <- module[!base::is.na(module) & module != ""]
  if (base::length(module) == 0) {
    stop("No valid `module` values provided.")
  }

  sat <- tryCatch(base::as.list(hc@satellite), error = function(e) list())
  tbl <- sat[["module_gene_list"]]
  if (base::is.null(tbl) || !base::is.data.frame(tbl) || !"module" %in% base::colnames(tbl)) {
    stop(
      "No `hc@satellite$module_gene_list` found. Run `hc_plot_cluster_heatmap()` once first ",
      "to create the module-to-gene table."
    )
  }

  available_modules <- unique(base::as.character(tbl$module))
  available_modules <- available_modules[!base::is.na(available_modules) & available_modules != ""]
  if (base::length(available_modules) == 0) {
    stop("No modules found in `hc@satellite$module_gene_list`.")
  }

  if (base::length(module) == 1 && base::tolower(module[[1]]) == "all") {
    return(.hc_llm_module_order(hc = hc))
  }

  module
}

.hc_llm_resolve_api_key <- function(api_key, llm) {
  if (base::is.null(api_key) || !base::nzchar(base::as.character(api_key[[1]]))) {
    api_key <- if (llm == "gemini") {
      Sys.getenv("GEMINI_API_KEY", unset = "")
    } else if (llm == "claude") {
      Sys.getenv("ANTHROPIC_API_KEY", unset = "")
    } else if (llm == "vllm") {
      Sys.getenv("VLLM_API_KEY", unset = "EMPTY")
    } else {
      Sys.getenv("OPENAI_API_KEY", unset = "")
    }
  } else {
    api_key <- base::as.character(api_key[[1]])
  }

  if (!base::nzchar(api_key)) {
    if (llm == "gemini") {
      stop("No Gemini API key found. Pass `api_key` or set `GEMINI_API_KEY`.")
    } else if (llm == "claude") {
      stop("No Claude API key found. Pass `api_key` or set `ANTHROPIC_API_KEY`.")
    } else if (llm == "vllm") {
      return("EMPTY")
    }
    stop("No OpenAI API key found. Pass `api_key` or set `OPENAI_API_KEY`.")
  }

  api_key
}

.hc_llm_resolve_vllm_base_url <- function(vllm_base_url, llm) {
  if (llm != "vllm") {
    return(NULL)
  }
  if (base::is.null(vllm_base_url) || !base::nzchar(base::as.character(vllm_base_url[[1]]))) {
    vllm_base_url <- Sys.getenv("VLLM_BASE_URL", unset = "http://localhost:8000/v1")
  } else {
    vllm_base_url <- base::as.character(vllm_base_url[[1]])
  }
  vllm_base_url <- sub("/+$", "", vllm_base_url)
  if (!base::nzchar(vllm_base_url)) {
    stop("No vLLM base URL found. Pass `vllm_base_url` or set `VLLM_BASE_URL`.")
  }
  vllm_base_url
}

.hc_llm_normalize_timeout <- function(timeout_sec) {
  if (base::is.null(timeout_sec)) {
    return(NULL)
  }
  timeout_sec <- base::as.numeric(timeout_sec[[1]])
  if (!base::is.finite(timeout_sec) || timeout_sec <= 0) {
    return(NULL)
  }
  timeout_sec
}

.hc_llm_resolve_timeout <- function(timeout_sec, llm, timeout_was_missing = FALSE) {
  if (isTRUE(timeout_was_missing) && identical(llm, "vllm")) {
    timeout_sec <- Sys.getenv("VLLM_TIMEOUT_SEC", unset = "300")
  }
  .hc_llm_normalize_timeout(timeout_sec)
}

.hc_llm_normalize_pause <- function(pause_sec) {
  if (base::is.null(pause_sec)) {
    return(0)
  }
  pause_sec <- base::as.numeric(pause_sec[[1]])
  if (!base::is.finite(pause_sec) || pause_sec <= 0) {
    return(0)
  }
  pause_sec
}

.hc_llm_require_ellmer <- function() {
  if (!requireNamespace("ellmer", quietly = TRUE)) {
    stop(
      "Package `ellmer` is required for LLM requests. Install it first.",
      call. = FALSE
    )
  }
}

.hc_llm_api_key_credentials <- function(api_key) {
  force(api_key)
  function() api_key
}

.hc_llm_with_ellmer_timeout <- function(timeout_sec, expr) {
  if (base::is.null(timeout_sec)) {
    return(force(expr))
  }

  old_timeout <- getOption("ellmer_timeout_s")
  options(ellmer_timeout_s = base::as.numeric(timeout_sec))
  on.exit(options(ellmer_timeout_s = old_timeout), add = TRUE)
  force(expr)
}

.hc_llm_gemini_temperature_unsupported_error <- function(msg) {
  msg <- base::tolower(base::as.character(msg[[1]]))
  base::grepl("unknown name \"temperature\"", msg, fixed = TRUE) ||
    base::grepl("cannot find field", msg, fixed = TRUE)
}

.hc_llm_request_gemini <- function(api_key,
                                   model,
                                   prompt,
                                   system_instruction,
                                   response_schema,
                                   temperature,
                                   timeout_sec) {
  .hc_llm_require_ellmer()

  run_request <- function(include_temperature = TRUE) {
    api_args <- list()
    if (isTRUE(include_temperature) &&
      !base::is.null(temperature) &&
      base::is.finite(temperature)) {
      api_args$temperature <- temperature
    }

    chat <- tryCatch(
      do.call(
        ellmer::chat_google_gemini,
        list(
          system_prompt = system_instruction,
          base_url = "https://generativelanguage.googleapis.com/v1beta/",
          credentials = .hc_llm_api_key_credentials(api_key),
          model = model,
          api_args = api_args,
          api_headers = c("x-goog-api-client" = "hcocena/1.9"),
          echo = "none"
        )
      ),
      error = function(e) {
        stop("Could not initialize Gemini chat via ellmer: ", base::conditionMessage(e), call. = FALSE)
      }
    )

    .hc_llm_with_ellmer_timeout(timeout_sec, chat$chat(prompt))
  }

  result_text <- tryCatch(
    run_request(include_temperature = TRUE),
    error = function(e) {
      msg <- base::conditionMessage(e)
      if (.hc_llm_gemini_temperature_unsupported_error(msg)) {
        return(
          tryCatch(
            run_request(include_temperature = FALSE),
            error = function(e2) {
              stop("Gemini request failed via ellmer: ", base::conditionMessage(e2), call. = FALSE)
            }
          )
        )
      }
      stop("Gemini request failed via ellmer: ", msg, call. = FALSE)
    }
  )

  result_text <- .hc_llm_strip_json_fences(result_text)
  if (!base::nzchar(result_text)) {
    stop("Gemini response contained no text payload.", call. = FALSE)
  }

  list(
    result_text = result_text,
    raw_response_text = result_text
  )
}

.hc_llm_request_openai <- function(api_key,
                                   model,
                                   prompt,
                                   system_instruction,
                                   response_schema,
                                   temperature,
                                   timeout_sec) {
  .hc_llm_require_ellmer()

  chat <- tryCatch(
    ellmer::chat_openai(
      system_prompt = system_instruction,
      base_url = "https://api.openai.com/v1",
      credentials = .hc_llm_api_key_credentials(api_key),
      model = model,
      api_args = list(
        temperature = temperature
      ),
      echo = "none"
    ),
    error = function(e) {
      stop("Could not initialize OpenAI chat via ellmer: ", base::conditionMessage(e), call. = FALSE)
    }
  )

  result_text <- tryCatch(
    .hc_llm_with_ellmer_timeout(timeout_sec, chat$chat(prompt)),
    error = function(e) {
      stop("OpenAI request failed via ellmer: ", base::conditionMessage(e), call. = FALSE)
    }
  )

  result_text <- .hc_llm_strip_json_fences(result_text)
  if (!base::nzchar(result_text)) {
    stop("OpenAI response contained no text payload.", call. = FALSE)
  }

  list(
    result_text = result_text,
    raw_response_text = result_text
  )
}

.hc_llm_request_claude <- function(api_key,
                                   model,
                                   prompt,
                                   system_instruction,
                                   temperature,
                                   timeout_sec) {
  .hc_llm_require_ellmer()

  chat <- tryCatch(
    ellmer::chat_anthropic(
      system_prompt = system_instruction,
      base_url = "https://api.anthropic.com/v1",
      credentials = .hc_llm_api_key_credentials(api_key),
      model = model,
      api_args = list(
        temperature = temperature
      ),
      echo = "none"
    ),
    error = function(e) {
      stop("Could not initialize Claude chat via ellmer: ", base::conditionMessage(e), call. = FALSE)
    }
  )

  result_text <- tryCatch(
    .hc_llm_with_ellmer_timeout(timeout_sec, chat$chat(prompt)),
    error = function(e) {
      stop("Claude request failed via ellmer: ", base::conditionMessage(e), call. = FALSE)
    }
  )

  result_text <- .hc_llm_strip_json_fences(result_text)
  if (!base::nzchar(result_text)) {
    stop("Claude response contained no text payload.", call. = FALSE)
  }

  list(
    result_text = result_text,
    raw_response_text = result_text
  )
}

.hc_llm_request_vllm <- function(api_key,
                                 model,
                                 prompt,
                                 system_instruction,
                                 temperature,
                                 timeout_sec,
                                 base_url) {
  .hc_llm_require_ellmer()

  chat <- tryCatch(
    ellmer::chat_vllm(
      base_url = base_url,
      model = model,
      credentials = .hc_llm_api_key_credentials(api_key),
      system_prompt = system_instruction,
      api_args = list(
        temperature = temperature,
        max_tokens = 8000,
        chat_template_kwargs = list(enable_thinking = FALSE)
      )
    ),
    error = function(e) {
      stop("Could not initialize vLLM chat via ellmer: ", base::conditionMessage(e), call. = FALSE)
    }
  )

  result_text <- tryCatch(
    .hc_llm_with_ellmer_timeout(timeout_sec, chat$chat(prompt)),
    error = function(e) {
      stop("vLLM request failed via ellmer: ", base::conditionMessage(e), call. = FALSE)
    }
  )

  result_text <- .hc_llm_strip_think_blocks(result_text)
  result_text <- .hc_llm_strip_json_fences(result_text)
  if (!base::nzchar(result_text)) {
    stop("vLLM response contained no text payload.", call. = FALSE)
  }

  list(
    result_text = result_text,
    raw_response_text = result_text
  )
}

.hc_gemini_normalize_genes <- function(genes) {
  genes <- base::as.character(genes)
  genes <- base::trimws(genes)
  genes <- genes[!base::is.na(genes) & genes != ""]
  genes <- genes[!duplicated(genes)]
  genes
}

.hc_gemini_normalize_context <- function(context) {
  if (base::is.null(context)) {
    return("")
  }
  context <- base::as.character(context)
  context <- base::trimws(context)
  context <- context[!base::is.na(context) & context != ""]
  if (base::length(context) == 0) {
    return("")
  }
  base::paste(context, collapse = "\n")
}

.hc_gemini_get_module_genes <- function(hc, module) {
  module <- base::as.character(module[[1]])
  sat <- tryCatch(base::as.list(hc@satellite), error = function(e) list())
  tbl <- sat[["module_gene_list"]]

  if (base::is.null(tbl) || !base::is.data.frame(tbl) || !"genes" %in% base::colnames(tbl) || !"module" %in% base::colnames(tbl)) {
    stop(
      "No `hc@satellite$module_gene_list` found. Run `hc_plot_cluster_heatmap()` once first ",
      "to create the module-to-gene table."
    )
  }

  tbl$genes <- base::as.character(tbl$genes)
  tbl$module <- base::as.character(tbl$module)

  module_lookup <- module
  if (!(module_lookup %in% tbl$module)) {
    label_map <- tryCatch(hc@integration@cluster[["module_label_map"]], error = function(e) NULL)
    if (!base::is.null(label_map) && base::length(label_map) > 0) {
      label_map <- base::as.character(label_map)
      map_names <- base::names(label_map)
      if (!base::is.null(map_names) && base::length(map_names) == base::length(label_map)) {
        base::names(label_map) <- base::as.character(map_names)
        if (module %in% base::names(label_map)) {
          module_lookup <- base::as.character(label_map[[module]])
        }
      }
    }
  }

  genes <- tbl$genes[tbl$module == module_lookup]
  genes <- .hc_gemini_normalize_genes(genes)
  if (base::length(genes) == 0) {
    stop("Could not resolve genes for module `", module, "`.")
  }

  list(
    label = module_lookup,
    module = module_lookup,
    genes = genes
  )
}

.hc_gemini_build_prompt <- function(label,
                                    genes,
                                    total_gene_count,
                                    context_text,
                                    truncated,
                                    rag_context_text = NULL,
                                    enrichment_text = NULL,
                                    llm = "gemini") {
  trunc_note <- if (isTRUE(truncated)) {
    base::paste0(
      "Only the first ", base::length(genes),
      " genes are sent here out of ", total_gene_count,
      " total input genes due to `max_genes`."
    )
  } else {
    base::paste0("All ", total_gene_count, " input genes are included.")
  }

  biological_context <- if (base::nzchar(context_text)) context_text else "none provided"
  has_rag_context <- !base::is.null(rag_context_text) &&
    base::nzchar(base::as.character(rag_context_text[[1]]))
  has_enrichment <- !base::is.null(enrichment_text) &&
    base::nzchar(base::as.character(enrichment_text[[1]]))
  # The three evidence sources carry very different weight, so say so: the
  # genes are the measurement, the enrichment terms are a statistical test on
  # exactly those genes, and the retrieved passages are background found by
  # similarity that may not concern this module at all.
  enrichment_instruction <- if (has_enrichment) {
    base::paste(
      "Statistical over-representation results for this exact gene list are supplied.",
      "Rank the evidence: the gene list first, the enrichment terms second, the retrieved literature last.",
      "Name enriched terms only where they sharpen the answer; do not list them back.",
      "Where the enrichment terms and the retrieved passages disagree, follow the enrichment terms."
    )
  } else {
    NULL
  }

  prompt_instructions <- if (llm == "vllm") {
    base::paste(
      "Infer the main biological program of this transcriptomic module.",
      base::paste0("The biological context is: ", biological_context, "."),
      "Return JSON matching the provided schema.",
      "Provide exactly three compact but informative text fields.",
      "Do not start with phrases like this module or the genes.",
      "Use the biological context as framing information, but do not simply repeat it.",
      "Make the answer module-specific and distinguishable from other modules from the same study.",
      "Avoid one-word or overly generic labels if the genes support something more specific.",
      "For `general_processes`, list 2 to 4 specific biological themes separated by ' / '.",
      "For `contextual_state`, give one compact but informative phrase of about 4 to 10 words describing the immediate module state.",
      "Do not use `contextual_state` to restate the cohort or timeline.",
      "Avoid generic phrases such as immature monocyte maturation, preterm infant maturation, or developing monocytes unless the module is truly nonspecific.",
      "Prefer specific states such as interferon-high inflammatory state, ribosome-high proliferative state, phagolysosomal activated state, antigen-presenting inflammatory state, platelet-like metabolic state, erythroid-skewed progenitor-like state, or macrophage-like transition state when supported.",
      "For `key_regulators`, list 2 to 5 likely transcription factors or signaling regulators separated by ' / '.",
      "If regulator evidence is weak, provide the most plausible regulators briefly rather than repeating the process.",
      if (has_rag_context) "Use the retrieved DoRAG passages as optional literature support, but prioritize the supplied genes and biological context. Do not claim that RAG is statistical enrichment." else NULL,
      enrichment_instruction,
      "Example style only:",
      '{"general_processes":"interferon signaling / antiviral innate immunity / antigen presentation","contextual_state":"activated interferon-high inflammatory monocyte state","key_regulators":"STAT1 / IRF7 / IRF9 / NFKB1"}'
    )
  } else {
    base::paste(
      "Infer the main biological program of this transcriptomic module.",
      base::paste0("The biological context is: ", biological_context, "."),
      "Return JSON matching the provided schema.",
      "Be extremely concise and provide exactly three short text elements.",
      "Do not start with phrases like this module or the genes.",
      "Treat the biological context as framing information, not as wording to repeat.",
      "Do not simply restate broad phrases from the context such as monocyte maturation unless the module is truly nonspecific.",
      "Make the answer module-specific and distinguishable from other modules from the same study.",
      "If the genes support a more specific program such as interferon response, phagolysosome, antigen presentation, cell cycle, platelet-like program, erythroid bias, mitochondrial metabolism, ribosome biogenesis, glycolysis, chemotaxis, inflammatory signaling, or tissue contamination, prefer that over generic immune wording.",
      "Provide `general_processes` as a short noun phrase listing the main biological program.",
      "Provide `contextual_state` as a short phrase describing the specific monocyte or transcriptional state in this study context.",
      "Provide `key_regulators` as a short phrase naming likely driving transcription factors or signaling regulators.",
      "If regulator evidence is weak, state the most plausible regulators briefly rather than repeating the biological process.",
      if (has_rag_context) "Use the retrieved DoRAG passages as optional literature support, but prioritize the supplied genes and biological context. Do not claim that RAG is statistical enrichment." else NULL,
      enrichment_instruction
    )
  }

  base::paste(
    prompt_instructions,
    "",
    base::paste0("Label: ", label),
    base::paste0("Gene-count note: ", trunc_note),
    base::paste0("Genes:\n", base::paste(genes, collapse = ", ")),
    if (has_enrichment) base::as.character(enrichment_text[[1]]) else NULL,
    if (has_rag_context) {
      base::paste0(
        "Retrieved DoRAG literature context:\n",
        base::as.character(rag_context_text[[1]])
      )
    } else {
      NULL
    },
    sep = "\n"
  )
}

.hc_llm_build_batch_prompt <- function(batch_inputs, context_text) {
  context_block <- if (base::nzchar(context_text)) {
    base::paste0("Context:\n", context_text, "\n")
  } else {
    "Context:\nNone provided.\n"
  }

  module_blocks <- base::vapply(
    batch_inputs,
    function(inp) {
      trunc_note <- if (isTRUE(inp$truncated)) {
        base::paste0(
          "Only the first ", base::length(inp$genes_use),
          " genes are listed here out of ", base::length(inp$genes_all),
          " total input genes due to `max_genes`."
        )
      } else {
        base::paste0("All ", base::length(inp$genes_all), " input genes are included.")
      }
      base::paste(
        base::paste0("Module label: ", inp$label),
        base::paste0("Gene-count note: ", trunc_note),
        base::paste0("Genes:\n", base::paste(inp$genes_use, collapse = ", ")),
        sep = "\n"
      )
    },
    FUN.VALUE = base::character(1)
  )

  base::paste(
    "Summarize the likely overarching biological function of each transcriptomic module separately.",
    "Return JSON matching the provided schema with one entry per module.",
    "For each module, provide `short_title` as a short plot-ready label with about 3 to 8 words.",
    "The provided context is the primary interpretation frame and must strongly constrain the answer.",
    "Prefer explanations that are compatible with the given biological context, cell type, cohort, and tissue.",
    "Do not assign unrelated tissue programs such as neuronal, epithelial, ciliary, muscular, or organ-specific identities unless the evidence is overwhelming and no context-compatible explanation fits.",
    "If genes look context-mismatched, keep the interpretation context-aware and mention uncertainty in `caveats` instead of drifting to an unrelated lineage.",
    "This is an interpretation task, not a statistical enrichment test.",
    "",
    context_block,
    base::paste(module_blocks, collapse = "\n\n---\n\n"),
    sep = "\n"
  )
}

.hc_llm_response_schema <- function() {
  json_schema <- list(
    type = "object",
    properties = list(
      general_processes = list(
        type = "string",
        description = "Two to four specific biological themes separated by ' / '."
      ),
      contextual_state = list(
        type = "string",
        description = "A compact but informative phrase of about four to ten words describing the specific transcriptional or cellular state."
      ),
      key_regulators = list(
        type = "string",
        description = "Two to five likely transcription factors or signaling regulators separated by ' / '."
      )
    ),
    required = c("general_processes", "contextual_state", "key_regulators")
  )
  json_schema
}

.hc_llm_batch_response_schema <- function(n_modules) {
  list(
    type = "object",
    properties = list(
      modules = list(
        type = "array",
        minItems = base::as.integer(n_modules),
        maxItems = base::as.integer(n_modules),
        items = .hc_llm_response_schema()
      )
    ),
    required = base::as.list("modules"),
    additionalProperties = FALSE
  )
}

.hc_gemini_extract_response_text <- function(resp_obj) {
  candidates <- resp_obj[["candidates"]]
  if (base::is.null(candidates) || base::length(candidates) == 0) {
    return("")
  }
  first_candidate <- candidates[[1]]
  content <- first_candidate[["content"]]
  if (base::is.null(content)) {
    return("")
  }
  parts <- content[["parts"]]
  if (base::is.null(parts) || base::length(parts) == 0) {
    return("")
  }
  text_parts <- base::vapply(
    parts,
    function(x) {
      txt <- x[["text"]]
      if (base::is.null(txt)) "" else base::as.character(txt)
    },
    FUN.VALUE = base::character(1)
  )
  text_parts <- text_parts[text_parts != ""]
  if (base::length(text_parts) == 0) {
    return("")
  }
  base::paste(text_parts, collapse = "\n")
}

.hc_openai_extract_response_text <- function(resp_obj) {
  message_obj <- tryCatch(resp_obj$choices[[1]]$message, error = function(e) NULL)
  if (base::is.null(message_obj)) {
    return("")
  }

  content <- message_obj$content
  if (base::is.null(content)) {
    return("")
  }

  if (base::is.character(content) && base::length(content) >= 1) {
    return(base::as.character(content[[1]]))
  }

  if (base::is.list(content) && base::length(content) > 0) {
    text_parts <- base::character(0)
    for (part in content) {
      part_text <- part$text
      if (base::is.character(part_text) && base::length(part_text) >= 1) {
        text_parts <- c(text_parts, base::as.character(part_text[[1]]))
      } else if (base::is.list(part_text) && !base::is.null(part_text$value)) {
        text_parts <- c(text_parts, base::as.character(part_text$value[[1]]))
      }
    }
    text_parts <- text_parts[!base::is.na(text_parts) & text_parts != ""]
    if (base::length(text_parts) > 0) {
      return(base::paste(text_parts, collapse = "\n"))
    }
  }

  ""
}

.hc_llm_strip_json_fences <- function(x) {
  x <- base::as.character(x[[1]])
  if (!base::nzchar(x)) {
    return("")
  }
  x <- sub("^\\s*```(?:json)?\\s*", "", x, perl = TRUE)
  x <- sub("\\s*```\\s*$", "", x, perl = TRUE)
  stringr::str_squish(x)
}

.hc_llm_strip_think_blocks <- function(x) {
  x <- base::as.character(x[[1]])
  if (!base::nzchar(x)) {
    return("")
  }
  # Local Qwen/vLLM deployments may emit reasoning in <think>...</think>.
  x <- gsub("(?is)<think>.*?</think>", " ", x, perl = TRUE)
  x <- sub("(?is)<think>.*$", " ", x, perl = TRUE)
  stringr::str_squish(x)
}

.hc_llm_summary_from_results <- function(results, hc = NULL) {
  if (base::is.null(results) || base::length(results) == 0) {
    return(base::data.frame(
      module = base::character(0),
      module_color = base::character(0),
      llm = base::character(0),
      model = base::character(0),
      general_processes = base::character(0),
      contextual_state = base::character(0),
      key_regulators = base::character(0),
      llm_long_output = base::character(0),
      response_json = base::character(0),
      short_title = base::character(0),
      overarching_function = base::character(0),
      confidence = base::character(0),
      gene_count_input = base::integer(0),
      gene_count_sent = base::integer(0),
      truncated = base::logical(0),
      rag_used = base::logical(0),
      rag_query = base::character(0),
      rag_context_count = base::integer(0),
      rag_citations = base::character(0),
      rag_error_message = base::character(0),
      rag_general_processes = base::character(0),
      rag_contextual_state = base::character(0),
      rag_key_regulators = base::character(0),
      enrichment_used = base::logical(0),
      enrichment_general_processes = base::character(0),
      enrichment_contextual_state = base::character(0),
      enrichment_key_regulators = base::character(0),
      status = base::character(0),
      error_message = base::character(0),
      timestamp = base::character(0),
      stringsAsFactors = FALSE
    ))
  }

  results <- results[base::vapply(results, base::is.list, FUN.VALUE = base::logical(1))]
  if (base::length(results) == 0) {
    return(base::data.frame(
      module = base::character(0),
      module_color = base::character(0),
      llm = base::character(0),
      model = base::character(0),
      general_processes = base::character(0),
      contextual_state = base::character(0),
      key_regulators = base::character(0),
      llm_long_output = base::character(0),
      response_json = base::character(0),
      short_title = base::character(0),
      overarching_function = base::character(0),
      confidence = base::character(0),
      gene_count_input = base::integer(0),
      gene_count_sent = base::integer(0),
      truncated = base::logical(0),
      rag_used = base::logical(0),
      rag_query = base::character(0),
      rag_context_count = base::integer(0),
      rag_citations = base::character(0),
      rag_error_message = base::character(0),
      rag_general_processes = base::character(0),
      rag_contextual_state = base::character(0),
      rag_key_regulators = base::character(0),
      enrichment_used = base::logical(0),
      enrichment_general_processes = base::character(0),
      enrichment_contextual_state = base::character(0),
      enrichment_key_regulators = base::character(0),
      status = base::character(0),
      error_message = base::character(0),
      timestamp = base::character(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- lapply(results, function(res) {
    base::data.frame(
      module = .hc_llm_result_scalar(res, c("label"), default = .hc_llm_result_scalar(res, c("module"), default = NA_character_)),
      module_color = .hc_llm_module_color(label = .hc_llm_result_label(res), hc = hc),
      llm = .hc_llm_result_scalar(res, c("llm")),
      model = .hc_llm_result_scalar(res, c("model")),
      general_processes = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "general_processes"))),
      contextual_state = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "contextual_state"))),
      key_regulators = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "key_regulators"))),
      llm_long_output = .hc_llm_result_long_output(res),
      response_json = .hc_llm_result_json(res),
      short_title = .hc_llm_result_short_title(res),
      overarching_function = .hc_llm_result_overarching(res),
      confidence = NA_character_,
      gene_count_input = .hc_llm_result_scalar(res, c("gene_count_input"), default = NA_integer_, mode = "integer"),
      gene_count_sent = .hc_llm_result_scalar(res, c("gene_count_sent"), default = NA_integer_, mode = "integer"),
      truncated = .hc_llm_result_scalar(res, c("truncated"), default = FALSE, mode = "logical"),
      rag_used = .hc_llm_result_rag_used(res),
      rag_query = .hc_llm_result_scalar(res, c("rag_query")),
      rag_context_count = .hc_llm_result_rag_context_count(res),
      rag_citations = .hc_llm_result_rag_citations(res),
      rag_error_message = .hc_llm_result_scalar(res, c("rag_error_message")),
      rag_general_processes = .hc_llm_result_rag_scalar(res, "general_processes"),
      rag_contextual_state = .hc_llm_result_rag_scalar(res, "contextual_state"),
      rag_key_regulators = .hc_llm_result_rag_scalar(res, "key_regulators"),
      enrichment_used = !base::is.null(.hc_llm_result_field(res, c("enrichment_response"))),
      enrichment_general_processes = .hc_llm_result_enrichment_scalar(res, "general_processes"),
      enrichment_contextual_state = .hc_llm_result_enrichment_scalar(res, "contextual_state"),
      enrichment_key_regulators = .hc_llm_result_enrichment_scalar(res, "key_regulators"),
      status = .hc_llm_result_scalar(res, c("status")),
      error_message = .hc_llm_result_scalar(res, c("error_message")),
      timestamp = .hc_llm_result_scalar(res, c("timestamp")),
      stringsAsFactors = FALSE
    )
  })
  out <- base::do.call(base::rbind, out)
  base::rownames(out) <- NULL

  if (!base::is.null(hc)) {
    module_order <- .hc_llm_module_order(hc = hc)
    if (base::length(module_order) > 0) {
      keep <- module_order[module_order %in% out$module]
      extras <- base::setdiff(out$module, keep)
      ord <- base::match(base::c(keep, extras), out$module)
      ord <- ord[!base::is.na(ord)]
      out <- out[ord, , drop = FALSE]
      base::rownames(out) <- NULL
    }
  }

  out
}

.hc_llm_result_label <- function(res) {
  lbl <- .hc_llm_result_scalar(res, c("label"), default = NA_character_)
  if (base::is.na(lbl) || !base::nzchar(lbl)) {
    lbl <- .hc_llm_result_scalar(res, c("module"), default = NA_character_)
  }
  lbl
}

.hc_llm_result_short_title <- function(res) {
  short_title <- .hc_llm_clean_text(.hc_llm_result_field(res, c("response", "contextual_state")))
  if (base::is.null(short_title) || base::length(short_title) == 0) {
    short_title <- NA_character_
  } else {
    short_title <- base::as.character(short_title[[1]])
  }
  if (!base::nzchar(short_title) || base::is.na(short_title)) {
    short_title <- .hc_llm_short_title_fallback(
      .hc_llm_result_field(res, c("response", "general_processes"))
    )
  }
  short_title
}

.hc_llm_result_overarching <- function(res) {
  val <- .hc_llm_clean_text(.hc_llm_result_field(res, c("response", "general_processes")))
  if (base::is.null(val) || base::length(val) == 0) {
    return(NA_character_)
  }
  base::as.character(val[[1]])
}

.hc_llm_result_long_output <- function(res) {
  gp <- .hc_llm_clean_text(.hc_llm_result_field(res, c("response", "general_processes")))
  cs <- .hc_llm_clean_text(.hc_llm_result_field(res, c("response", "contextual_state")))
  kr <- .hc_llm_clean_text(.hc_llm_result_field(res, c("response", "key_regulators")))
  parts <- c(
    if (!base::is.null(gp) && base::nzchar(base::as.character(gp[[1]]))) base::paste0("General processes: ", base::as.character(gp[[1]])) else NULL,
    if (!base::is.null(cs) && base::nzchar(base::as.character(cs[[1]]))) base::paste0("Contextual state: ", base::as.character(cs[[1]])) else NULL,
    if (!base::is.null(kr) && base::nzchar(base::as.character(kr[[1]]))) base::paste0("Key regulators: ", base::as.character(kr[[1]])) else NULL
  )
  if (base::length(parts) == 0) {
    return(NA_character_)
  }
  base::paste(parts, collapse = " | ")
}

.hc_llm_result_json <- function(res) {
  resp <- .hc_llm_result_field(res, c("response"))
  if (base::is.null(resp)) {
    return(NA_character_)
  }
  out <- tryCatch(
    jsonlite::toJSON(resp, auto_unbox = TRUE, null = "null"),
    error = function(e) NA_character_
  )
  base::as.character(out[[1]])
}

.hc_llm_result_field <- function(x, path) {
  cur <- x
  for (nm in path) {
    if (base::is.null(cur) || !(nm %in% base::names(cur))) {
      return(NULL)
    }
    cur <- cur[[nm]]
  }
  cur
}

.hc_llm_result_scalar <- function(res, path, default = NA_character_, mode = c("character", "integer", "logical")) {
  mode <- base::match.arg(mode)
  val <- .hc_llm_result_field(res, path)
  if (base::is.null(val) || base::length(val) == 0) {
    return(default)
  }

  if (identical(mode, "character")) {
    val <- base::as.character(val[[1]])
    if (base::length(val) == 0) {
      return(base::as.character(default[[1]]))
    }
    return(val)
  }
  if (identical(mode, "integer")) {
    val <- .hc_as_integer_safely(val[[1]])
    if (base::length(val) == 0 || base::is.na(val)) {
      return(.hc_as_integer_safely(default[[1]]))
    }
    return(val)
  }

  val <- base::as.logical(val[[1]])
  if (base::length(val) == 0 || base::is.na(val)) {
    return(base::as.logical(default[[1]]))
  }
  val
}

.hc_llm_result_enrichment_scalar <- function(res, field, default = NA_character_) {
  val <- .hc_llm_result_field(res, c("enrichment_response", field))
  val <- .hc_llm_clean_text(val)
  if (base::is.null(val) || base::length(val) == 0) {
    return(default)
  }
  val <- base::as.character(val[[1]])
  if (base::length(val) == 0 || base::is.na(val) || !base::nzchar(val)) {
    return(default)
  }
  val
}

.hc_llm_result_rag_scalar <- function(res, field, default = NA_character_) {
  val <- .hc_llm_result_field(res, c("rag_response", field))
  val <- .hc_llm_clean_text(val)
  if (base::is.null(val) || base::length(val) == 0) {
    return(default)
  }
  val <- base::as.character(val[[1]])
  if (base::length(val) == 0 || base::is.na(val) || !base::nzchar(val)) {
    return(default)
  }
  val
}

.hc_llm_result_rag_used <- function(res) {
  rag_status <- .hc_llm_result_field(res, c("rag", "status"))
  status_ok <- !base::is.null(rag_status) &&
    base::length(rag_status) > 0L &&
    base::identical(base::as.character(rag_status[[1]]), "ok")
  status_ok &&
    .hc_llm_result_rag_context_count(res) > 0L &&
    !base::is.null(.hc_llm_result_field(res, c("rag_response")))
}

.hc_llm_result_rag_context_count <- function(res) {
  rag <- .hc_llm_result_field(res, c("rag"))
  if (base::is.null(rag)) {
    return(0L)
  }
  base::as.integer(.hc_llm_rag_context_count(rag))
}

.hc_llm_result_rag_citations <- function(res) {
  rag <- .hc_llm_result_field(res, c("rag"))
  if (base::is.null(rag)) {
    return(NA_character_)
  }
  citations <- .hc_llm_rag_citations_text(rag)
  if (!base::nzchar(citations)) {
    return(NA_character_)
  }
  citations
}

.hc_llm_short_title_fallback <- function(term, max_chars = 64) {
  if (base::is.null(term) || base::length(term) == 0) {
    return("No interpretation available")
  }
  term <- .hc_llm_clean_text(base::as.character(term[[1]]))
  if (!base::nzchar(term) || base::is.na(term)) {
    return("No interpretation available")
  }

  term <- stringr::str_squish(term)
  term <- sub("^This module is best described as\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module can be summarized as\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module primarily reflects\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module reflects\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module primarily captures\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module captures\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module is characterized by\\s+", "", term, ignore.case = TRUE)
  term <- sub("^This module represents\\s+", "", term, ignore.case = TRUE)
  term <- sub("^processes related to\\s+", "", term, ignore.case = TRUE)
  term <- sub("^processes involving\\s+", "", term, ignore.case = TRUE)
  term <- sub("^the coordinated regulation of\\s+", "", term, ignore.case = TRUE)
  term <- sub("^the orchestration of\\s+", "", term, ignore.case = TRUE)
  term <- sub("^the activation of\\s+", "", term, ignore.case = TRUE)
  term <- sub("^the diverse functions of\\s+", "", term, ignore.case = TRUE)
  term <- sub("\\.$", "", term)

  pieces <- unlist(strsplit(
    term,
    "\\s*,\\s*|\\s*;\\s*|\\s+and\\s+|\\s+with\\s+|\\s+linked to\\s+|\\s+coupled to\\s+",
    perl = TRUE
  ))
  pieces <- trimws(pieces)
  pieces <- pieces[nzchar(pieces)]
  if (length(pieces) == 0) {
    pieces <- term
  }

  chosen <- character(0)
  for (piece in pieces) {
    candidate <- paste(c(chosen, piece), collapse = " / ")
    if (nchar(candidate) > max_chars && length(chosen) > 0) {
      break
    }
    chosen <- c(chosen, piece)
    if (nchar(candidate) >= (max_chars - 8)) {
      break
    }
  }

  if (length(chosen) == 0) {
    chosen <- pieces[[1]]
  }

  short_title <- paste(chosen, collapse = " / ")
  short_title <- .hc_llm_clean_text(stringr::str_squish(short_title))
  stringr::str_trunc(short_title, width = max_chars)
}

.hc_llm_clean_text <- function(x) {
  if (base::is.null(x) || base::length(x) == 0) {
    return(x)
  }
  x <- base::as.character(x)
  x <- stringr::str_squish(x)
  x <- gsub("\\s*/\\s*/+\\s*", " / ", x, perl = TRUE)
  x <- gsub("\\s*\\|\\s*\\|+\\s*", " | ", x, perl = TRUE)
  x <- gsub("\\s{2,}", " ", x, perl = TRUE)
  x <- base::trimws(x)
  x
}

.hc_llm_export_results_excel <- function(hc,
                                         results,
                                         summary_tbl,
                                         slot_name) {
  if (base::is.null(hc) || !requireNamespace("openxlsx", quietly = TRUE)) {
    return(invisible(NULL))
  }

  details_tbl <- base::do.call(
    base::rbind,
    lapply(results, function(res) {
      base::data.frame(
        module = .hc_llm_result_label(res),
        llm = .hc_llm_result_scalar(res, c("llm")),
        model = .hc_llm_result_scalar(res, c("model")),
        general_processes = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "general_processes"))),
        contextual_state = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "contextual_state"))),
        key_regulators = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "key_regulators"))),
        llm_long_output = .hc_llm_result_long_output(res),
        response_json = .hc_llm_result_json(res),
        gene_count_input = .hc_llm_result_scalar(res, c("gene_count_input"), default = NA_integer_, mode = "integer"),
        gene_count_sent = .hc_llm_result_scalar(res, c("gene_count_sent"), default = NA_integer_, mode = "integer"),
        truncated = .hc_llm_result_scalar(res, c("truncated"), default = FALSE, mode = "logical"),
        rag_used = .hc_llm_result_rag_used(res),
        rag_query = .hc_llm_result_scalar(res, c("rag_query")),
        rag_context_count = .hc_llm_result_rag_context_count(res),
        rag_citations = .hc_llm_result_rag_citations(res),
        rag_error_message = .hc_llm_result_scalar(res, c("rag_error_message")),
        rag_general_processes = .hc_llm_result_rag_scalar(res, "general_processes"),
        rag_contextual_state = .hc_llm_result_rag_scalar(res, "contextual_state"),
        rag_key_regulators = .hc_llm_result_rag_scalar(res, "key_regulators"),
        enrichment_general_processes = .hc_llm_result_enrichment_scalar(res, "general_processes"),
        enrichment_contextual_state = .hc_llm_result_enrichment_scalar(res, "contextual_state"),
        enrichment_key_regulators = .hc_llm_result_enrichment_scalar(res, "key_regulators"),
        status = .hc_llm_result_scalar(res, c("status")),
        error_message = .hc_llm_result_scalar(res, c("error_message")),
        prompt = .hc_llm_result_scalar(res, c("prompt")),
        rag_prompt = .hc_llm_result_scalar(res, c("rag_prompt")),
        enrichment_prompt = .hc_llm_result_scalar(res, c("enrichment_prompt")),
        enrichment_context_text = .hc_llm_result_scalar(res, c("enrichment_context_text")),
        timestamp = .hc_llm_result_scalar(res, c("timestamp")),
        stringsAsFactors = FALSE
      )
    })
  )
  base::rownames(details_tbl) <- NULL

  out_dir <- .hc_resolve_output_dir(hc)
  file <- base::file.path(out_dir, base::paste0(slot_name, "_summary.xlsx"))
  tryCatch(
    .hc_write_xlsx_atomic(
      x = list(
        summary = summary_tbl,
        details = details_tbl
      ),
      file = file,
      overwrite = TRUE
    ),
    error = function(e) warning("Could not write LLM Excel summary: ", base::conditionMessage(e))
  )
  invisible(file)
}

.hc_llm_module_order <- function(hc) {
  heatmap_info <- tryCatch(.hc_heatmap_cache_info(hc@integration@cluster), error = function(e) NULL)
  label_map <- tryCatch(hc@integration@cluster[["module_label_map"]], error = function(e) NULL)
  if (!base::is.null(heatmap_info) &&
    !base::is.null(heatmap_info$row_order) &&
    base::length(heatmap_info$row_order) > 0 &&
    !base::is.null(label_map) &&
    base::length(label_map) > 0) {
    label_map <- base::as.character(label_map)
    map_names <- base::names(hc@integration@cluster[["module_label_map"]])
    if (!base::is.null(map_names) && base::length(map_names) == base::length(label_map)) {
      base::names(label_map) <- base::as.character(map_names)
    }
    mapped <- base::as.character(label_map[base::as.character(heatmap_info$row_order)])
    mapped <- mapped[!base::is.na(mapped) & base::nzchar(mapped)]
    if (base::length(mapped) > 0) {
      return(base::unique(mapped))
    }
  }

  sat <- tryCatch(base::as.list(hc@satellite), error = function(e) list())
  tbl <- sat[["module_gene_list"]]
  if (base::is.null(tbl) || !base::is.data.frame(tbl) || !"module" %in% base::colnames(tbl)) {
    return(base::character(0))
  }
  .hc_llm_natural_module_order(tbl$module)
}

.hc_llm_module_color <- function(label, hc = NULL) {
  if (base::is.null(hc)) {
    return("grey70")
  }
  label <- base::as.character(label[[1]])
  label_map <- tryCatch(hc@integration@cluster[["module_label_map"]], error = function(e) NULL)
  if (base::is.null(label_map) || base::length(label_map) == 0) {
    return("grey70")
  }
  label_map <- base::as.character(label_map)
  map_names <- base::names(hc@integration@cluster[["module_label_map"]])
  if (!base::is.null(map_names) && base::length(map_names) == base::length(label_map)) {
    base::names(label_map) <- base::as.character(map_names)
  }
  inv_map <- stats::setNames(base::names(label_map), label_map)
  col <- inv_map[[label]]
  if (base::is.null(col) || !base::nzchar(base::as.character(col))) {
    return("grey70")
  }
  base::as.character(col)
}
