#' LLM-based interpretation of modules or feature sets
#'
#' Asks a large language model what biological program each module (or any
#' other set of genes, proteins or other features) most likely represents.
#' This is an AI-assisted interpretation, not a statistical test.
#'
#' Every set is interpreted in up to three separate requests, so that each
#' answer sees only the input it is meant to see:
#' \itemize{
#'   \item `general_processes`: from the features alone, without context;
#'   \item `contextual_state` and `key_regulators`: from the features and the
#'     `context`;
#'   \item `grounded_processes` and `supporting_evidence`: from the features,
#'     the `context` and the `evidence` - by default the significant terms of
#'     [hc_functional_enrichment()]; optionally also the regulators of
#'     [hc_upstream_inference()] and literature passages from a DoRAG server.
#' }
#' [hc_plot_llm_enrichment()] draws one plot per field, so the plots can be
#' compared directly: what the model makes of the genes alone, with context,
#' and with statistical evidence. Without stored enrichment results the third
#' request is left out (with a message).
#'
#' The requests can be sent in two ways:
#' * through an API (`provider = "claude"`, `"openai"`, `"gemini"`, or
#'   `"local"` for an OpenAI-compatible vLLM server), which needs an API key
#'   for the hosted providers;
#' * without an API key (`provider = "manual"`): one prompt per request type is
#'   prepared (the first is copied to the clipboard); paste each into a new
#'   chat of any provider and read every answer back with [hc_llm_import()].
#'   The results are stored exactly like API results.
#'
#' Results are stored in `hc_satellite(hc, slot_name)` (one entry per set) and
#' `hc_satellite(hc, paste0(slot_name, "_summary"))` (a table).
#'
#' @param hc Optional `HCoCenaExperiment`. Required to annotate modules, to use
#'   stored enrichment results and to store the results; not needed for
#'   `genes`.
#' @param genes Optional feature set(s) to interpret instead of modules: a
#'   character vector (one set) or a named list of character vectors.
#' @param modules Modules of `hc` to interpret, or `"all"` (default). Needs the
#'   module table created by [hc_plot_cluster_heatmap()]. Ignored when `genes`
#'   is given.
#' @param system_prompt Optional instruction that tells the model who it is
#'   and how to answer, e.g. `"You are an expert proteomics analyst."`. A
#'   neutral default suited to any kind of features is used when `NULL`.
#' @param context Optional free text describing the study, e.g. tissue,
#'   disease, platform or time point. Used for all fields except
#'   `general_processes`.
#' @param provider `"claude"`, `"openai"`, `"gemini"`, `"local"` (vLLM or any
#'   OpenAI-compatible server) or `"manual"` (no API key; copy and paste via
#'   the web chat of any provider).
#' @param api_key API key. If `NULL`, it is read from `ANTHROPIC_API_KEY`,
#'   `OPENAI_API_KEY`, `GEMINI_API_KEY` or `VLLM_API_KEY`. Not used for
#'   `"manual"`.
#' @param model Optional model name. Defaults to `"claude-sonnet-4-6"`,
#'   `"gpt-4o-mini"`, `"gemini-2.5-pro"` or
#'   `"Qwen/Qwen2.5-VL-32B-Instruct"` (local). For `"manual"` it is only
#'   recorded with the result.
#' @param base_url Server address for `provider = "local"`. If `NULL`, read
#'   from `VLLM_BASE_URL`, falling back to `"http://localhost:8000/v1"`.
#' @param evidence Evidence for `grounded_processes`: `"enrichment"`
#'   (default; needs [hc_functional_enrichment()]), any combination with
#'   `"upstream"` (needs [hc_upstream_inference()]) and `"rag"` (literature
#'   passages, needs `rag_url`), `"all"`, or `"none"` to skip this request.
#' @param evidence_text Optional evidence supplied directly, as a named list or
#'   named character vector (names = set names). Added to the evidence above;
#'   the only way to give evidence together with `genes`.
#' @param evidence_top Maximum number of enriched terms and of regulators per
#'   set. Default 10.
#' @param evidence_qval Adjusted p-value cutoff for enriched terms. Default
#'   0.05.
#' @param rag_url Retrieval endpoint of a DoRAG server for `"rag"`. If `NULL`,
#'   read from the environment variable `HCOCENA_RAG_URL`.
#' @param rag_top_k Number of literature passages per set. Default 10.
#' @param slot_name Name under which results are stored in `hc`. Default
#'   `"llm_enrichment"`.
#' @param verbose Logical; print progress messages. Default `TRUE`.
#'
#' @return With `hc`: the updated `HCoCenaExperiment` (for `"manual"`, with
#'   the pending request stored). Without `hc`: a data frame with one row per
#'   set (for `"manual"`, an `hc_llm_request` object to pass to
#'   [hc_llm_import()]).
#'
#' @seealso [hc_llm_import()], [hc_plot_llm_enrichment()],
#'   [hc_list_llm_models()]
#' @examples
#' # Without an API key: build the requests, answer them in a web chat, import.
#' req <- hc_llm_enrichment(
#'   genes = list(IFN = c("STAT1", "IRF7", "CXCL10", "GBP1", "IFI44L")),
#'   context = "Whole blood of patients with acute viral infection",
#'   evidence_text = list(IFN = "Enriched: type I interferon signaling (q = 1e-6)"),
#'   provider = "manual",
#'   verbose = FALSE
#' )
#' length(req$parts)  # genes only / with context / with evidence
#' cat(substr(req$parts[[1]], 1, 300))
#'
#' # The answers copied from the web chat (here written by hand):
#' answer <- paste0(
#'   '{"request_id": "', req$id, '", "results": [{"set": "IFN", ',
#'   '"general_processes": "type I interferon signaling / antiviral response",',
#'   ' "contextual_state": "interferon-high antiviral state",',
#'   ' "key_regulators": "STAT1 / IRF7 / IRF9",',
#'   ' "grounded_processes": "type I interferon response",',
#'   ' "supporting_evidence": "type I interferon signaling"}]}'
#' )
#' hc_llm_import(req, text = answer, verbose = FALSE)
#'
#' \donttest{
#' # Through an API (needs a key):
#' if (nzchar(Sys.getenv("ANTHROPIC_API_KEY"))) {
#'   hc_llm_enrichment(
#'     genes = c("STAT1", "IRF7", "CXCL10", "GBP1", "IFI44L"),
#'     provider = "claude"
#'   )
#' }
#' }
#' @export
hc_llm_enrichment <- function(hc = NULL,
                              genes = NULL,
                              modules = "all",
                              system_prompt = NULL,
                              context = NULL,
                              provider = c("claude", "openai", "gemini", "local", "manual"),
                              api_key = NULL,
                              model = NULL,
                              base_url = NULL,
                              evidence = "enrichment",
                              evidence_text = NULL,
                              evidence_top = 10,
                              evidence_qval = 0.05,
                              rag_url = NULL,
                              rag_top_k = 10,
                              slot_name = "llm_enrichment",
                              verbose = TRUE) {
  provider <- base::match.arg(provider)
  if (!base::is.null(hc) && !methods::is(hc, "HCoCenaExperiment")) {
    stop("`hc` must be an `HCoCenaExperiment` or NULL.", call. = FALSE)
  }
  if (base::is.null(hc) && base::is.null(genes)) {
    stop("Provide `hc` (to interpret its modules) or `genes`.", call. = FALSE)
  }
  slot_name <- .hc_llm_check_scalar(slot_name, "slot_name")
  evidence_default <- base::missing(evidence)
  levels <- .hc_llm_evidence_levels(evidence)

  sets <- .hc_llm_gene_sets(hc = hc, genes = genes, modules = modules)
  system_prompt <- .hc_llm_system_prompt(system_prompt)
  context_text <- .hc_gemini_normalize_context(context)
  ev <- .hc_llm_prepare_evidence(
    hc = hc, sets = sets, levels = levels, default = evidence_default,
    context_text = context_text, evidence_text = evidence_text,
    top = evidence_top, qval = evidence_qval,
    rag_url = rag_url, rag_top_k = rag_top_k, verbose = verbose
  )
  tasks <- .hc_llm_tasks(with_evidence = !base::is.null(ev))

  if (identical(provider, "manual")) {
    return(.hc_llm_manual_request(
      hc = hc,
      sets = sets,
      tasks = tasks,
      system_prompt = system_prompt,
      context_text = context_text,
      evidence = ev,
      model = model,
      slot_name = slot_name,
      verbose = verbose
    ))
  }

  llm <- if (identical(provider, "local")) "vllm" else provider
  model <- .hc_llm_default_model(model, llm)
  api_key <- .hc_llm_key_or_hint(api_key, llm)
  vllm_base_url <- .hc_llm_resolve_vllm_base_url(base_url, llm)
  timeout_sec <- if (identical(llm, "vllm")) 300 else 120

  results <- base::vector("list", base::length(sets))
  for (i in base::seq_along(sets)) {
    label <- base::names(sets)[[i]]
    features <- sets[[i]]
    if (isTRUE(verbose)) {
      message("[", i, "/", base::length(sets), "] ", label, " (",
              base::length(features), " features) -> ", provider, " / ", model)
    }
    response <- list()
    prompts <- base::character(0)
    raw <- base::character(0)
    errors <- base::character(0)
    for (task in tasks) {
      prompt <- .hc_llm_task_prompt(task, label, features, context_text, ev)
      prompts[[task]] <- prompt
      reply <- tryCatch(
        .hc_llm_request_by_provider(
          llm = llm,
          api_key = api_key,
          model = model,
          prompt = prompt,
          system_instruction = system_prompt,
          response_schema = .hc_llm_response_schema(fields = .hc_llm_task_fields(task)),
          temperature = 0.2,
          timeout_sec = timeout_sec,
          vllm_base_url = vllm_base_url
        ),
        error = function(e) e
      )
      if (inherits(reply, "error")) {
        errors <- c(errors, base::paste0(task, ": ", base::conditionMessage(reply)))
        next
      }
      parsed <- .hc_llm_parse_module_response(reply$result_text)
      if (!base::is.null(parsed$parse_error)) {
        errors <- c(errors, base::paste0(task, ": ", parsed$parse_error))
      }
      for (f in .hc_llm_task_fields(task)) response[[f]] <- parsed[[f]]
      raw[[task]] <- reply$raw_response_text
    }
    if (base::length(errors) > 0) {
      warning("LLM request for `", label, "` failed: ", errors[[1]], call. = FALSE)
    }
    results[[i]] <- .hc_llm_result(
      label = label,
      features = features,
      context_text = context_text,
      llm = provider,
      model = model,
      system_prompt = system_prompt,
      prompt = prompts,
      response = response,
      raw_text = raw,
      error = if (base::length(errors) > 0) base::paste(errors, collapse = "; ") else NULL,
      evidence = ev,
      tasks = tasks
    )
  }
  base::names(results) <- base::names(sets)

  answered <- base::vapply(results, function(r) {
    any(!base::is.na(base::unlist(r$response)))
  }, base::logical(1))
  if (!any(answered)) {
    stop(
      "All LLM requests failed. First problem: ", results[[1]]$error_message,
      call. = FALSE
    )
  }
  .hc_llm_store_results(hc = hc, results = results, slot_name = slot_name)
}

#' Import an answer from a provider's web chat
#'
#' Reads an answer to a request created with
#' `hc_llm_enrichment(provider = "manual")`, checks that it belongs to that
#' request, and stores it exactly like an API result. Markdown code fences,
#' text around the JSON and typographic quotes are tolerated.
#'
#' A request has one part per request type (features only, with context,
#' with evidence), and long requests are split further. Import every answer
#' with this function; after each import the next open part is copied to the
#' clipboard. An answer may also contain the fields of several parts.
#'
#' @param x The `HCoCenaExperiment` returned by
#'   `hc_llm_enrichment(provider = "manual")`, or the `hc_llm_request` object
#'   it returned when no `hc` was given.
#' @param text Optional answer text. If `NULL` and `file` is `NULL`, the
#'   clipboard is read (needs the `clipr` package and a clipboard, i.e. not
#'   inside a remote RStudio Server or Docker session).
#' @param file Optional path to a text file containing the answer.
#' @param model Optional name of the model used in the web chat, recorded with
#'   the result.
#' @param slot_name Slot used in `hc_llm_enrichment()`. Default
#'   `"llm_enrichment"`.
#' @param verbose Logical; print progress messages. Default `TRUE`.
#'
#' @return For an `HCoCenaExperiment`, the updated object. For an
#'   `hc_llm_request`, a data frame with one row per set answered so far.
#' @seealso [hc_llm_enrichment()]
#' @examples
#' req <- hc_llm_enrichment(
#'   genes = list(IFN = c("STAT1", "IRF7", "CXCL10")),
#'   provider = "manual", verbose = FALSE
#' )
#' answer <- paste0(
#'   "Here is the result:\n```json\n",
#'   '{"request_id": "', req$id, '", "results": [{"set": "IFN", ',
#'   '"general_processes": "interferon signaling", ',
#'   '"contextual_state": "antiviral state", "key_regulators": "STAT1"}]}',
#'   "\n```"
#' )
#' hc_llm_import(req, text = answer, model = "web chat", verbose = FALSE)
#' @export
hc_llm_import <- function(x,
                          text = NULL,
                          file = NULL,
                          model = NULL,
                          slot_name = "llm_enrichment",
                          verbose = TRUE) {
  is_hc <- methods::is(x, "HCoCenaExperiment")
  if (is_hc) {
    request <- hc_satellite(x, base::paste0(slot_name, "_pending"))
    if (base::is.null(request)) {
      stop(
        "No pending request found in `", slot_name, "_pending`. Run ",
        "`hc_llm_enrichment(hc, provider = \"manual\")` first.",
        call. = FALSE
      )
    }
  } else if (inherits(x, "hc_llm_request")) {
    request <- x
  } else {
    stop("`x` must be an `HCoCenaExperiment` or an `hc_llm_request`.",
         call. = FALSE)
  }

  answer <- .hc_llm_read_answer(text = text, file = file)
  parsed <- .hc_llm_extract_answer_json(answer)

  answer_id <- parsed$request_id
  if (base::is.null(answer_id) || !base::nzchar(base::as.character(answer_id))) {
    warning("The answer has no `request_id`; it could not be checked against ",
            "the request.", call. = FALSE)
  } else if (!identical(base::as.character(answer_id), request$id)) {
    stop(
      "This answer belongs to request `", answer_id, "`, not to `",
      request$id, "`. Copy the answer to the current request.",
      call. = FALSE
    )
  }

  entries <- parsed$results
  if (base::is.null(entries) || base::length(entries) == 0) {
    stop("The answer contains no `results`.", call. = FALSE)
  }
  set_names <- base::names(request$sets)
  key <- function(z) base::tolower(base::trimws(base::as.character(z)))
  model <- if (base::is.null(model)) "web chat" else .hc_llm_check_scalar(model, "model")
  all_fields <- base::unlist(lapply(request$tasks, .hc_llm_task_fields))
  answers <- request$answers
  if (base::is.null(answers)) answers <- list()

  touched <- base::character(0)
  unknown <- base::character(0)
  for (entry in entries) {
    nm <- entry$set
    if (base::is.null(nm)) nm <- entry$module
    hit <- base::match(key(nm), key(set_names))
    if (base::length(nm) != 1 || base::is.na(hit)) {
      unknown <- c(unknown, base::as.character(nm)[1])
      next
    }
    label <- set_names[[hit]]
    got <- answers[[label]]
    if (base::is.null(got)) got <- list()
    for (f in all_fields) {
      if (!base::is.null(entry[[f]])) {
        got[[f]] <- base::paste(base::as.character(base::unlist(entry[[f]])), collapse = " / ")
      }
    }
    answers[[label]] <- got
    touched <- c(touched, label)
  }
  if (base::length(unknown) > 0) {
    warning("Ignored entries for unknown sets: ",
            base::paste(unknown, collapse = ", "), call. = FALSE)
  }
  if (base::length(touched) == 0) {
    stop("None of the answered sets match the request (expected: ",
         base::paste(set_names, collapse = ", "), ").", call. = FALSE)
  }
  request$answers <- answers
  request$done <- .hc_llm_done_parts(request)

  results <- lapply(base::names(answers), function(label) {
    .hc_llm_result(
      label = label,
      features = request$sets[[label]],
      context_text = request$context,
      llm = "manual",
      model = model,
      system_prompt = request$system_prompt,
      prompt = request$parts[base::vapply(request$part_sets, function(x) label %in% x, base::logical(1))],
      response = answers[[label]],
      raw_text = answer,
      error = NULL,
      evidence = request$evidence,
      tasks = request$tasks
    )
  })
  base::names(results) <- base::names(answers)

  if (!is_hc) {
    open_parts <- base::setdiff(base::seq_along(request$parts), request$done)
    if (base::length(open_parts) > 0 && isTRUE(verbose)) {
      message("Parts not yet answered: ", base::paste(open_parts, collapse = ", "),
              ". They are in `request$parts`.")
    }
    return(.hc_llm_summary_from_results(results = results, hc = NULL))
  }

  hc <- .hc_llm_store_results(
    hc = x,
    results = results,
    slot_name = slot_name,
    request = request
  )
  if (isTRUE(verbose)) {
    message("Imported ", base::length(base::unique(touched)), " set(s); ",
            base::length(request$done), " of ", base::length(request$parts),
            " part(s) answered.")
  }
  .hc_llm_announce_next_part(hc, slot_name = slot_name, verbose = verbose)
  hc
}

# ---- internals -------------------------------------------------------------

# Request types: what each request sees and which fields it answers.
.hc_llm_tasks <- function(with_evidence = TRUE) {
  c("features", "context", if (isTRUE(with_evidence)) "evidence")
}

.hc_llm_task_fields <- function(task) {
  switch(task,
    features = "general_processes",
    context = c("contextual_state", "key_regulators"),
    evidence = c("grounded_processes", "supporting_evidence")
  )
}

.hc_llm_task_title <- function(task) {
  switch(task,
    features = "features only",
    context = "features + context",
    evidence = "features + context + evidence"
  )
}

.hc_llm_fields <- function(evidence = FALSE) {
  base::unlist(lapply(.hc_llm_tasks(with_evidence = evidence), .hc_llm_task_fields))
}

.hc_llm_field_spec <- function(fields = .hc_llm_fields(evidence = TRUE)) {
  spec <- c(
    general_processes = "- general_processes: 2 to 4 specific biological processes or pathways as short phrases, separated by \" / \", at most about 90 characters in total; do not name individual genes or proteins and do not list them in parentheses",
    contextual_state = "- contextual_state: one phrase of about 4 to 10 words describing the specific biological or cellular state the features indicate in this context",
    key_regulators = "- key_regulators: 2 to 5 likely upstream regulators (e.g. transcription factors or signalling molecules), separated by \" / \"",
    grounded_processes = "- grounded_processes: 2 to 4 specific biological processes or pathways as short phrases, separated by \" / \", at most about 90 characters in total, based on the features, the context and the evidence; do not name individual genes or proteins",
    supporting_evidence = "- supporting_evidence: the given enriched terms, regulators or numbered literature passages that support grounded_processes, separated by \" / \"; \"none\" if none fits"
  )
  base::unname(spec[fields])
}

.hc_llm_system_prompt <- function(system_prompt) {
  if (!base::is.null(system_prompt)) {
    return(.hc_llm_check_scalar(system_prompt, "system_prompt"))
  }
  base::paste(
    "You are an expert in molecular and systems biology.",
    "You interpret sets of co-regulated features such as genes, proteins or metabolites.",
    "Infer the biological program each set most likely represents.",
    "Be specific: prefer concrete pathways, processes and cell states over generic labels, and make different sets distinguishable from each other.",
    "Use only the information given in the request and do not claim results that were not supplied.",
    "If the signal is weak or mixed, say so briefly."
  )
}

.hc_llm_check_scalar <- function(x, arg) {
  if (!base::is.character(x) || base::length(x) != 1 || base::is.na(x) ||
      !base::nzchar(base::trimws(x))) {
    stop("`", arg, "` must be a single non-empty character string.", call. = FALSE)
  }
  x
}

.hc_llm_default_model <- function(model, llm) {
  if (!base::is.null(model)) {
    return(.hc_llm_check_scalar(model, "model"))
  }
  switch(llm,
    claude = "claude-sonnet-4-6",
    openai = "gpt-4o-mini",
    gemini = "gemini-2.5-pro",
    vllm = "Qwen/Qwen2.5-VL-32B-Instruct"
  )
}

.hc_llm_key_or_hint <- function(api_key, llm) {
  tryCatch(
    .hc_llm_resolve_api_key(api_key, llm),
    error = function(e) {
      stop(base::conditionMessage(e),
           " Without an API key, use `provider = \"manual\"`.", call. = FALSE)
    }
  )
}

.hc_llm_gene_sets <- function(hc, genes, modules) {
  if (!base::is.null(genes)) {
    if (base::is.atomic(genes)) {
      genes <- list(custom_geneset = genes)
    }
    if (!base::is.list(genes) || base::length(genes) == 0) {
      stop("`genes` must be a character vector or a named list of them.", call. = FALSE)
    }
    nms <- base::names(genes)
    if (base::is.null(nms)) nms <- base::rep("", base::length(genes))
    nms[!base::nzchar(nms)] <- base::paste0("set_", base::which(!base::nzchar(nms)))
    if (base::anyDuplicated(nms)) {
      stop("Names in `genes` must be unique.", call. = FALSE)
    }
    sets <- lapply(genes, .hc_gemini_normalize_genes)
    base::names(sets) <- nms
  } else {
    module_ids <- .hc_llm_resolve_modules(hc = hc, module = modules)
    infos <- lapply(module_ids, function(m) .hc_gemini_get_module_genes(hc = hc, module = m))
    sets <- lapply(infos, function(i) i$genes)
    base::names(sets) <- base::vapply(infos, function(i) i$label, base::character(1))
  }
  empty <- base::names(sets)[base::lengths(sets) == 0]
  if (base::length(empty) > 0) {
    stop("These sets contain no features: ", base::paste(empty, collapse = ", "),
         call. = FALSE)
  }
  sets
}

# Evidence for the third request. With the default `evidence` a missing
# source only drops that request (message); an explicit choice must work.
.hc_llm_prepare_evidence <- function(hc, sets, levels, default, context_text,
                                     evidence_text, top, qval, rag_url,
                                     rag_top_k, verbose) {
  if (base::length(levels) == 0 && base::is.null(evidence_text)) {
    return(NULL)
  }
  ev <- tryCatch(
    .hc_llm_build_evidence(
      hc = hc, sets = sets, levels = levels, context_text = context_text,
      evidence_text = evidence_text, top = top, qval = qval,
      rag_url = rag_url, rag_top_k = rag_top_k, verbose = verbose
    ),
    error = function(e) e
  )
  if (inherits(ev, "error")) {
    if (!isTRUE(default)) {
      stop(base::conditionMessage(ev), call. = FALSE)
    }
    if (base::is.null(evidence_text)) {
      if (isTRUE(verbose)) {
        message("No evidence request (grounded_processes): ", base::conditionMessage(ev))
      }
      return(NULL)
    }
    ev <- .hc_llm_build_evidence(
      hc = hc, sets = sets, levels = base::character(0), context_text = context_text,
      evidence_text = evidence_text, verbose = verbose
    )
    levels <- base::character(0)
  }
  ev$levels <- c(levels, if (!base::is.null(evidence_text)) "supplied")
  ev
}

.hc_llm_task_prompt <- function(task, label, features, context_text, evidence) {
  ctx <- if (base::nzchar(context_text)) context_text else "none provided"
  ev_text <- if (identical(task, "evidence")) evidence$text[[label]] else ""
  base::paste(c(
    "Interpret the following feature set.",
    if (!identical(task, "features")) base::paste0("Context: ", ctx),
    if (identical(task, "evidence")) c("", .hc_llm_evidence_preamble(evidence$levels)),
    "",
    base::paste0("Set: ", label),
    base::paste0("Features (", base::length(features), "): ", base::paste(features, collapse = ", ")),
    if (base::nzchar(ev_text)) ev_text,
    "",
    "Return a JSON object with exactly these fields:",
    .hc_llm_field_spec(.hc_llm_task_fields(task)),
    "Return only the JSON object."
  ), collapse = "\n")
}

.hc_llm_manual_prompt <- function(id, task, sets, system_prompt, context_text,
                                  evidence, part, n_parts) {
  ctx <- if (base::nzchar(context_text)) context_text else "none provided"
  fields <- .hc_llm_task_fields(task)
  blocks <- base::vapply(base::names(sets), function(nm) {
    ev_text <- if (identical(task, "evidence")) evidence$text[[nm]] else ""
    base::paste0("### Set: ", nm, "\nFeatures (", base::length(sets[[nm]]), "): ",
                 base::paste(sets[[nm]], collapse = ", "),
                 if (base::nzchar(ev_text)) base::paste0("\n", ev_text) else "")
  }, base::character(1))
  example <- base::paste0(
    '{"request_id": "', id, '", "results": [{"set": "<set name>", ',
    base::paste0('"', fields, '": "..."', collapse = ", "), "}]}"
  )
  base::paste(c(
    system_prompt,
    "",
    base::paste0("(Request ", id, ", part ", part, " of ", n_parts, ": ",
                 .hc_llm_task_title(task), ")"),
    base::paste0("Interpret each of the following ", base::length(sets),
                 " feature sets separately."),
    if (!identical(task, "features")) base::paste0("Context: ", ctx),
    if (identical(task, "evidence")) c("", .hc_llm_evidence_preamble(evidence$levels)),
    "",
    blocks,
    "",
    "Answer with exactly one JSON object and nothing else, in this form:",
    example,
    "Rules:",
    base::paste0("- request_id: copy exactly \"", id, "\""),
    "- results: one entry per set, with \"set\" spelled exactly as given above",
    .hc_llm_field_spec(fields)
  ), collapse = "\n")
}

.hc_llm_manual_request <- function(hc, sets, tasks, system_prompt, context_text,
                                   evidence, model, slot_name, verbose,
                                   max_sets = 40L, max_chars = 60000L) {
  # Time stamp plus process id; deliberately no sample(), which would advance
  # the user's random number stream.
  id <- base::paste0(
    "hcocena-", base::gsub(".", "", base::format(Sys.time(), "%Y%m%d-%H%M%OS3"), fixed = TRUE), "-",
    base::sprintf("%04d", Sys.getpid() %% 10000L)
  )
  # One series of parts per request type, each small enough for a web chat.
  part_task <- base::character(0)
  part_sets <- list()
  for (task in tasks) {
    current <- base::character(0)
    size <- 0L
    for (nm in base::names(sets)) {
      n <- base::nchar(base::paste(sets[[nm]], collapse = ", ")) + 40L +
        if (identical(task, "evidence")) base::nchar(evidence$text[[nm]]) else 0L
      if (base::length(current) > 0 &&
          (base::length(current) >= max_sets || size + n > max_chars)) {
        part_task <- c(part_task, task)
        part_sets[[base::length(part_sets) + 1L]] <- current
        current <- base::character(0)
        size <- 0L
      }
      current <- c(current, nm)
      size <- size + n
    }
    part_task <- c(part_task, task)
    part_sets[[base::length(part_sets) + 1L]] <- current
  }
  n_parts <- base::length(part_task)
  parts <- base::vapply(base::seq_len(n_parts), function(p) {
    .hc_llm_manual_prompt(id, part_task[[p]], sets[part_sets[[p]]], system_prompt,
                          context_text, evidence, part = p, n_parts = n_parts)
  }, base::character(1))
  request <- base::structure(
    list(
      id = id,
      sets = sets,
      tasks = tasks,
      parts = parts,
      part_task = part_task,
      part_sets = part_sets,
      system_prompt = system_prompt,
      context = context_text,
      evidence = evidence,
      model = model,
      answers = list(),
      done = base::integer(0),
      created = base::as.character(Sys.time())
    ),
    class = "hc_llm_request"
  )

  files <- .hc_llm_write_parts(request, hc)
  copied <- .hc_llm_to_clipboard(parts[[1]])
  where <- if (copied) {
    "Part 1 has been copied to the clipboard.\n"
  } else {
    base::paste0("Clipboard not available; part 1 is in: ", files[[1]], "\n")
  }
  if (isTRUE(verbose)) {
    message(
      "Request ", id, ": ", base::length(sets), " set(s) in ", n_parts,
      " part(s) (", base::paste(base::vapply(tasks, .hc_llm_task_title, base::character(1)), collapse = "; "), ").\n",
      where,
      "1. Paste it into a NEW chat of your provider (ChatGPT, Claude, Gemini, ...).\n",
      "2. Copy the complete answer.\n",
      if (base::is.null(hc)) "3. Run: hc_llm_import(<request>, text = <answer>)" else
        "3. Run: hc <- hc_llm_import(hc)",
      "\nUse a new chat for every part, so that one answer cannot influence the next.",
      "\nThe next part is copied automatically after each import."
    )
  }
  if (base::is.null(hc)) {
    return(request)
  }
  sat <- base::as.list(hc@satellite)
  sat[[base::paste0(slot_name, "_pending")]] <- request
  hc@satellite <- S4Vectors::SimpleList(sat)
  invisible(hc)
}

# Parts whose fields are answered for all of their sets.
.hc_llm_done_parts <- function(request) {
  base::which(base::vapply(base::seq_along(request$parts), function(p) {
    fields <- .hc_llm_task_fields(request$part_task[[p]])
    base::all(base::vapply(request$part_sets[[p]], function(nm) {
      got <- request$answers[[nm]]
      !base::is.null(got) && base::all(fields %in% base::names(got))
    }, base::logical(1)))
  }, base::logical(1)))
}

.hc_llm_write_parts <- function(request, hc) {
  dir <- if (base::is.null(hc)) {
    base::tempdir()
  } else {
    tryCatch(.hc_resolve_output_dir(hc), error = function(e) base::tempdir())
  }
  dir <- base::sub("[/\\]+$", "", dir)
  files <- base::file.path(
    dir,
    base::paste0("llm_request_", request$id, "_part", base::seq_along(request$parts), ".txt")
  )
  for (i in base::seq_along(files)) {
    tryCatch(base::writeLines(request$parts[[i]], files[[i]], useBytes = TRUE),
             error = function(e) NULL)
  }
  files
}

.hc_llm_to_clipboard <- function(text) {
  if (!requireNamespace("clipr", quietly = TRUE) || !base::interactive()) {
    return(FALSE)
  }
  ok <- tryCatch(isTRUE(clipr::clipr_available()), error = function(e) FALSE)
  if (!ok) {
    return(FALSE)
  }
  tryCatch({
    clipr::write_clip(text, object_type = "character")
    TRUE
  }, error = function(e) FALSE)
}

.hc_llm_read_answer <- function(text, file) {
  if (!base::is.null(text)) {
    return(base::paste(base::as.character(text), collapse = "\n"))
  }
  if (!base::is.null(file)) {
    if (!base::file.exists(file)) {
      stop("File not found: ", file, call. = FALSE)
    }
    return(base::paste(base::readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  }
  if (!requireNamespace("clipr", quietly = TRUE)) {
    stop("Reading the clipboard needs the `clipr` package. Install it or pass ",
         "the answer via `text =` or `file =`.", call. = FALSE)
  }
  ok <- tryCatch(isTRUE(clipr::clipr_available()), error = function(e) FALSE)
  if (!ok) {
    stop("No clipboard available (e.g. in RStudio Server or Docker). Save the ",
         "answer to a file and use `file =`, or pass it via `text =`.", call. = FALSE)
  }
  base::paste(clipr::read_clip(), collapse = "\n")
}

# Pull the JSON object out of a pasted chat answer: tolerates code fences,
# surrounding prose and typographic quotes.
.hc_llm_extract_answer_json <- function(answer) {
  txt <- base::enc2utf8(base::as.character(answer))
  txt <- base::gsub("\ufeff|\u200b", "", txt)
  txt <- base::gsub("[\u201c\u201d\u201e\u201f\u2033]", "\"", txt)
  txt <- base::gsub("[\u2018\u2019\u201a\u201b]", "'", txt)

  candidates <- base::character(0)
  fenced <- base::regmatches(
    txt,
    base::gregexpr("```(?:json|JSON)?\\s*(\\{.*?\\})\\s*```", txt, perl = TRUE)
  )[[1]]
  if (base::length(fenced) > 0) {
    candidates <- c(candidates, base::gsub("^```(?:json|JSON)?\\s*|\\s*```$", "", fenced, perl = TRUE))
  }
  starts <- base::gregexpr("\\{", txt)[[1]]
  ends <- base::gregexpr("\\}", txt)[[1]]
  if (starts[[1]] > 0 && ends[[1]] > 0) {
    last <- base::max(ends)
    for (s in starts) {
      if (s < last) candidates <- c(candidates, base::substr(txt, s, last))
    }
  }
  for (cand in candidates) {
    parsed <- tryCatch(jsonlite::fromJSON(cand, simplifyVector = FALSE),
                       error = function(e) NULL)
    if (base::is.list(parsed) && !base::is.null(parsed$results)) {
      return(parsed)
    }
  }
  stop("Could not find the JSON answer in the pasted text. Copy the complete ",
       "answer of the chat, including the part starting with ",
       "{\"request_id\".", call. = FALSE)
}

.hc_llm_result <- function(label, features, context_text, llm, model,
                           system_prompt, prompt, response, raw_text, error,
                           evidence = NULL, tasks = .hc_llm_tasks(!base::is.null(evidence))) {
  fields <- base::unlist(lapply(tasks, .hc_llm_task_fields))
  if (base::is.null(response)) response <- list()
  full <- lapply(fields, function(f) {
    v <- response[[f]]
    if (base::is.null(v) || base::length(v) == 0) NA_character_ else base::as.character(v[[1]])
  })
  base::names(full) <- fields
  ok <- base::is.null(error) && !any(base::is.na(base::unlist(full)))
  list(
    label = label,
    module = label,
    genes_input = features,
    genes_sent = features,
    gene_count_input = base::length(features),
    gene_count_sent = base::length(features),
    truncated = FALSE,
    context = if (base::nzchar(context_text)) context_text else NULL,
    llm = llm,
    model = model,
    status = if (ok) "ok" else "error",
    error_message = if (!base::is.null(error)) {
      base::as.character(error)
    } else if (!ok) {
      "not all fields answered yet"
    } else {
      NA_character_
    },
    system_prompt = system_prompt,
    prompt = base::paste(base::as.character(base::unlist(prompt)), collapse = "\n\n---\n\n"),
    response = full,
    raw_response_text = base::paste(base::as.character(base::unlist(raw_text)), collapse = "\n\n"),
    evidence = if (base::is.null(evidence)) "none" else base::paste(evidence$levels, collapse = "+"),
    evidence_text = if (base::is.null(evidence)) NA_character_ else evidence$text[[label]],
    citations = if (base::is.null(evidence)) NA_character_ else evidence$citations[[label]],
    timestamp = base::as.character(Sys.time())
  )
}

.hc_llm_store_results <- function(hc, results, slot_name, request = NULL) {
  if (base::is.null(hc)) {
    out <- .hc_llm_summary_from_results(results = results, hc = NULL)
    base::attr(out, "results") <- results
    return(out)
  }
  sat <- base::as.list(hc@satellite)
  stored <- sat[[slot_name]]
  if (!base::is.list(stored)) stored <- list()
  for (nm in base::names(results)) stored[[nm]] <- results[[nm]]
  summary_tbl <- .hc_llm_summary_from_results(results = stored, hc = hc)
  sat[[slot_name]] <- stored
  sat[[base::paste0(slot_name, "_summary")]] <- summary_tbl

  pending_name <- base::paste0(slot_name, "_pending")
  if (!base::is.null(request)) {
    if (base::length(request$done) == base::length(request$parts)) {
      sat[[pending_name]] <- NULL
    } else {
      sat[[pending_name]] <- request
    }
  }
  hc@satellite <- S4Vectors::SimpleList(sat)
  .hc_llm_export_results_excel(hc = hc, results = stored,
                               summary_tbl = summary_tbl, slot_name = slot_name)
  hc
}

.hc_llm_announce_next_part <- function(hc, slot_name, verbose) {
  request <- hc_satellite(hc, base::paste0(slot_name, "_pending"))
  if (base::is.null(request)) {
    if (isTRUE(verbose)) message("All parts of the request are imported.")
    return(invisible(NULL))
  }
  next_part <- base::min(base::setdiff(base::seq_along(request$parts), request$done))
  copied <- .hc_llm_to_clipboard(request$parts[[next_part]])
  if (isTRUE(verbose)) {
    message(
      base::length(request$parts) - base::length(request$done), " part(s) still open. Part ",
      next_part, " (", .hc_llm_task_title(request$part_task[[next_part]]), ") ",
      if (copied) "has been copied to the clipboard." else
        "is in the output folder (clipboard not available).",
      " Paste it into a NEW chat and run hc_llm_import(hc) again."
    )
  }
  invisible(NULL)
}

#' @export
print.hc_llm_request <- function(x, ...) {
  cat("<hc_llm_request>", x$id, "\n")
  cat(" sets: ", base::length(x$sets), " (", base::paste(utils::head(base::names(x$sets), 5), collapse = ", "),
      if (base::length(x$sets) > 5) ", ..." else "", ")\n", sep = "")
  cat(" parts:", base::length(x$parts), "-",
      base::paste(base::vapply(x$part_task, .hc_llm_task_title, base::character(1)), collapse = "; "), "\n")
  cat(" import each chat answer with hc_llm_import(<this object>, text = <answer>)\n")
  invisible(x)
}
