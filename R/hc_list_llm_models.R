#' List Available LLM Models
#'
#' Queries the model-listing endpoint of each selected provider and returns the
#' available models as a data frame. Supports the same providers as
#' [hc_module_function_llm()]: Gemini, Claude (Anthropic), OpenAI/ChatGPT, and a
#' local OpenAI-compatible server (vLLM).
#'
#' @param provider Character vector of providers to query. One or more of
#'   `"gemini"`, `"claude"` (`"anthropic"`), `"openai"` (`"chatgpt"`), and
#'   `"vllm"` (`"local"`); or `"all"` (default), which queries the three cloud
#'   providers plus the local server when `base_url` is supplied.
#' @param api_key Optional API key. Only applied when a single provider is
#'   requested; otherwise keys are taken from the per-provider environment
#'   variables (`GEMINI_API_KEY`, `ANTHROPIC_API_KEY`, `OPENAI_API_KEY`,
#'   `VLLM_API_KEY`).
#' @param base_url Base URL of the local OpenAI-compatible server. Required for
#'   `"vllm"`; falls back to `VLLM_BASE_URL` or `"http://localhost:8000/v1"`.
#' @param timeout_sec Per-request timeout in seconds. Default `30`.
#' @return A data frame with columns `provider`, `model`, and `info`. Providers
#'   that cannot be reached (e.g. a missing API key) are skipped with a warning
#'   rather than aborting the whole call.
#' @examples
#' # Contacts the provider APIs; each provider needs its API key set in the
#' # environment (or a reachable server for "vllm").
#' \donttest{
#' if (nzchar(Sys.getenv("GEMINI_API_KEY"))) {
#'   hc_list_llm_models("gemini")
#' }
#' if (nzchar(Sys.getenv("OPENAI_API_KEY"))) {
#'   hc_list_llm_models("openai")
#' }
#' }
#' @export
hc_list_llm_models <- function(provider = "all",
                               api_key = NULL,
                               base_url = NULL,
                               timeout_sec = 30) {
  providers <- base::tolower(base::as.character(provider))
  providers[providers == "chatgpt"] <- "openai"
  providers[providers == "anthropic"] <- "claude"
  providers[providers == "local"] <- "vllm"
  if (base::any(providers == "all")) {
    providers <- base::c("gemini", "claude", "openai")
    if (!base::is.null(base_url) && base::nzchar(base::as.character(base_url[[1]]))) {
      providers <- base::c(providers, "vllm")
    }
  }
  providers <- base::unique(providers)

  valid <- base::c("gemini", "claude", "openai", "vllm")
  unknown <- base::setdiff(providers, valid)
  if (base::length(unknown) > 0) {
    stop(
      "Unknown provider(s): ", base::paste(unknown, collapse = ", "),
      ". Use gemini, claude, openai/chatgpt, or vllm/local."
    )
  }

  single <- base::length(providers) == 1L
  results <- base::list()
  for (p in providers) {
    key_arg <- if (single) api_key else NULL
    df <- tryCatch(
      .hc_llm_fetch_models(
        provider = p, api_key = key_arg, base_url = base_url, timeout_sec = timeout_sec
      ),
      error = function(e) {
        warning("Skipping `", p, "`: ", base::conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (!base::is.null(df) && base::nrow(df) > 0) {
      results[[p]] <- df
    }
  }

  if (base::length(results) == 0) {
    return(base::data.frame(
      provider = base::character(0),
      model = base::character(0),
      info = base::character(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- base::do.call(base::rbind, results)
  out <- out[base::order(out$provider, base::tolower(out$model)), , drop = FALSE]
  base::rownames(out) <- NULL
  base::message(
    "Found ", base::nrow(out), " model(s) across ", base::length(results), " provider(s)."
  )
  out
}


#' First-or-default helper for nested JSON fields.
#' @noRd
.hc_llm_models_default <- function(x, default = "") {
  if (base::is.null(x) || base::length(x) == 0) {
    return(default)
  }
  base::as.character(x[[1]])
}


#' Build the standard model-listing data frame.
#' @noRd
.hc_llm_models_df <- function(provider, model, info) {
  base::data.frame(
    provider = base::rep(provider, base::length(model)),
    model = base::as.character(model),
    info = base::as.character(info),
    stringsAsFactors = FALSE
  )
}


#' GET a URL and parse the JSON body, raising informative errors on failure.
#' @noRd
.hc_llm_http_get_json <- function(url, headers = base::character(0), timeout_sec = 30) {
  if (!base::requireNamespace("httr", quietly = TRUE)) {
    stop("Package `httr` is required to list models.")
  }
  resp <- httr::GET(
    url,
    httr::add_headers(.headers = headers),
    httr::timeout(base::max(1, base::as.numeric(timeout_sec[[1]])))
  )
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  if (httr::http_error(resp)) {
    detail <- base::substr(base::gsub("\\s+", " ", base::as.character(txt)), 1, 300)
    stop(
      "HTTP ", httr::status_code(resp),
      if (base::nzchar(detail)) base::paste0(" - ", detail) else ""
    )
  }
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}


#' Fetch the model list for a single provider.
#' @noRd
.hc_llm_fetch_models <- function(provider, api_key, base_url, timeout_sec) {
  if (provider == "gemini") {
    key <- .hc_llm_resolve_api_key(api_key = api_key, llm = "gemini")
    js <- .hc_llm_http_get_json(
      url = base::paste0(
        "https://generativelanguage.googleapis.com/v1beta/models?key=",
        utils::URLencode(key, reserved = TRUE)
      ),
      timeout_sec = timeout_sec
    )
    models <- if (base::is.null(js$models)) base::list() else js$models
    models <- base::Filter(function(m) {
      "generateContent" %in% base::unlist(m$supportedGenerationMethods)
    }, models)
    if (base::length(models) == 0) {
      return(NULL)
    }
    return(.hc_llm_models_df(
      provider = "gemini",
      model = base::vapply(models, function(m) base::sub("^models/", "", .hc_llm_models_default(m$name)), base::character(1)),
      info = base::vapply(models, function(m) .hc_llm_models_default(m$displayName), base::character(1))
    ))
  }

  if (provider == "openai") {
    key <- .hc_llm_resolve_api_key(api_key = api_key, llm = "openai")
    js <- .hc_llm_http_get_json(
      url = "https://api.openai.com/v1/models",
      headers = base::c(Authorization = base::paste("Bearer", key)),
      timeout_sec = timeout_sec
    )
    data <- if (base::is.null(js$data)) base::list() else js$data
    if (base::length(data) == 0) {
      return(NULL)
    }
    return(.hc_llm_models_df(
      provider = "openai",
      model = base::vapply(data, function(m) .hc_llm_models_default(m$id), base::character(1)),
      info = base::vapply(data, function(m) .hc_llm_models_default(m$owned_by), base::character(1))
    ))
  }

  if (provider == "claude") {
    key <- .hc_llm_resolve_api_key(api_key = api_key, llm = "claude")
    js <- .hc_llm_http_get_json(
      url = "https://api.anthropic.com/v1/models?limit=100",
      headers = base::c(`x-api-key` = key, `anthropic-version` = "2023-06-01"),
      timeout_sec = timeout_sec
    )
    data <- if (base::is.null(js$data)) base::list() else js$data
    if (base::length(data) == 0) {
      return(NULL)
    }
    return(.hc_llm_models_df(
      provider = "claude",
      model = base::vapply(data, function(m) .hc_llm_models_default(m$id), base::character(1)),
      info = base::vapply(data, function(m) .hc_llm_models_default(m$display_name), base::character(1))
    ))
  }

  if (provider == "vllm") {
    burl <- .hc_llm_resolve_vllm_base_url(vllm_base_url = base_url, llm = "vllm")
    key <- .hc_llm_resolve_api_key(api_key = api_key, llm = "vllm")
    js <- .hc_llm_http_get_json(
      url = base::paste0(burl, "/models"),
      headers = base::c(Authorization = base::paste("Bearer", key)),
      timeout_sec = timeout_sec
    )
    data <- if (base::is.null(js$data)) base::list() else js$data
    if (base::length(data) == 0) {
      return(NULL)
    }
    return(.hc_llm_models_df(
      provider = "vllm",
      model = base::vapply(data, function(m) .hc_llm_models_default(m$id), base::character(1)),
      info = base::vapply(data, function(m) .hc_llm_models_default(m$owned_by), base::character(1))
    ))
  }

  stop("Unsupported provider: ", provider)
}
