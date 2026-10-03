# Internal engine for hc_llm_enrichment(): provider requests, response
# parsing, module/gene lookup, result tables and Excel export.

# Newer models (e.g. Claude Opus/Sonnet 5, some OpenAI reasoning models) reject
# the `temperature` parameter. Retry once without it instead of failing.
.hc_llm_temperature_rejected <- function(msg) {
  msg <- base::tolower(base::as.character(msg[[1]]))
  base::grepl("temperature", msg, fixed = TRUE) &&
    base::grepl("deprecat|not supported|unsupported|does not support|only the default|not allowed", msg)
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
  call_once <- function(temp) {
    .hc_llm_request_by_provider_once(
      llm = llm, api_key = api_key, model = model, prompt = prompt,
      system_instruction = system_instruction, response_schema = response_schema,
      temperature = temp, timeout_sec = timeout_sec, vllm_base_url = vllm_base_url
    )
  }
  tryCatch(
    call_once(temperature),
    error = function(e) {
      if (!base::is.null(temperature) &&
          llm %in% c("claude", "openai") &&
          .hc_llm_temperature_rejected(base::conditionMessage(e))) {
        return(call_once(NULL))
      }
      stop(e)
    }
  )
}

.hc_llm_request_by_provider_once <- function(llm,
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
      api_args = if (base::is.null(temperature)) list() else list(temperature = temperature),
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
      api_args = if (base::is.null(temperature)) list() else list(temperature = temperature),
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

.hc_llm_response_schema <- function(fields = c("general_processes", "contextual_state", "key_regulators")) {
  processes <- "Two to four specific biological processes or pathways as short phrases separated by ' / ', at most about 90 characters in total, without naming or listing individual genes or proteins."
  properties <- list(
    general_processes = list(type = "string", description = processes),
    contextual_state = list(
      type = "string",
      description = "A compact but informative phrase of about four to ten words describing the specific biological or cellular state in this context."
    ),
    key_regulators = list(
      type = "string",
      description = "Two to five likely transcription factors or signaling regulators separated by ' / '."
    ),
    grounded_processes = list(
      type = "string",
      description = base::paste(processes, "Based on the features, the context and the supplied evidence.")
    ),
    supporting_evidence = list(
      type = "string",
      description = "The given enriched terms, regulators or numbered literature passages that support grounded_processes, separated by ' / '; 'none' if none fits."
    )
  )
  list(type = "object", properties = properties[fields], required = fields)
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
      grounded_processes = base::character(0),
      supporting_evidence = base::character(0),
      evidence = base::character(0),
      llm_long_output = base::character(0),
      response_json = base::character(0),
      short_title = base::character(0),
      overarching_function = base::character(0),
      confidence = base::character(0),
      gene_count_input = base::integer(0),
      gene_count_sent = base::integer(0),
      truncated = base::logical(0),
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
      grounded_processes = base::character(0),
      supporting_evidence = base::character(0),
      evidence = base::character(0),
      llm_long_output = base::character(0),
      response_json = base::character(0),
      short_title = base::character(0),
      overarching_function = base::character(0),
      confidence = base::character(0),
      gene_count_input = base::integer(0),
      gene_count_sent = base::integer(0),
      truncated = base::logical(0),
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
      grounded_processes = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "grounded_processes"))),
      supporting_evidence = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "supporting_evidence"))),
      evidence = .hc_llm_result_scalar(res, c("evidence"), default = "none"),
      llm_long_output = .hc_llm_result_long_output(res),
      response_json = .hc_llm_result_json(res),
      short_title = .hc_llm_result_short_title(res),
      overarching_function = .hc_llm_result_overarching(res),
      confidence = NA_character_,
      gene_count_input = .hc_llm_result_scalar(res, c("gene_count_input"), default = NA_integer_, mode = "integer"),
      gene_count_sent = .hc_llm_result_scalar(res, c("gene_count_sent"), default = NA_integer_, mode = "integer"),
      truncated = .hc_llm_result_scalar(res, c("truncated"), default = FALSE, mode = "logical"),
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
  if (is.finite(max_chars)) stringr::str_trunc(short_title, width = max_chars) else short_title
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
        grounded_processes = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "grounded_processes"))),
        supporting_evidence = .hc_llm_clean_text(.hc_llm_result_scalar(res, c("response", "supporting_evidence"))),
        evidence = .hc_llm_result_scalar(res, c("evidence"), default = "none"),
        evidence_text = .hc_llm_result_scalar(res, c("evidence_text")),
        citations = .hc_llm_result_scalar(res, c("citations")),
        llm_long_output = .hc_llm_result_long_output(res),
        response_json = .hc_llm_result_json(res),
        gene_count_input = .hc_llm_result_scalar(res, c("gene_count_input"), default = NA_integer_, mode = "integer"),
        gene_count_sent = .hc_llm_result_scalar(res, c("gene_count_sent"), default = NA_integer_, mode = "integer"),
        truncated = .hc_llm_result_scalar(res, c("truncated"), default = FALSE, mode = "logical"),
        status = .hc_llm_result_scalar(res, c("status")),
        error_message = .hc_llm_result_scalar(res, c("error_message")),
        prompt = .hc_llm_result_scalar(res, c("prompt")),
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

