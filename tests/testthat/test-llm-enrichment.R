## hc_llm_enrichment() / hc_llm_import(): one simple interface for all
## providers, plus a copy-and-paste route through a provider's web chat.
## API calls are replaced by a fake provider, so no key or network is needed.

fake_answer <- function(id, sets, prefix = "", evidence = FALSE) {
  entries <- vapply(sets, function(s) {
    sprintf(
      '{"set": "%s", "general_processes": "%sprocess of %s", "contextual_state": "state of %s", "key_regulators": "REG1 / REG2"%s}',
      s, prefix, s, s,
      if (evidence) sprintf(', "grounded_processes": "grounded %s", "supporting_evidence": "TERM1"', s) else ""
    )
  }, character(1))
  sprintf('{"request_id": "%s", "results": [%s]}', id, paste(entries, collapse = ", "))
}

all_fields_json <- '{"general_processes": "interferon signaling", "contextual_state": "antiviral state", "key_regulators": "STAT1", "grounded_processes": "type I interferon response", "supporting_evidence": "IFN term"}'

clustered_with_modules <- function() {
  hc <- hc_example_data("clustered")
  suppressMessages(suppressWarnings(hc_plot_cluster_heatmap(hc, file_name = FALSE)))
}

# ---- manual route (no API key) ------------------------------------------------

test_that("manual mode builds one prompt per request type without an API key", {
  withr::local_envvar(c(ANTHROPIC_API_KEY = NA, OPENAI_API_KEY = NA, GEMINI_API_KEY = NA))
  req <- hc_llm_enrichment(
    genes = list(IFN = c("STAT1", "IRF7"), Bcell = c("CD19", "MS4A1")),
    system_prompt = "You are an expert proteomics analyst.",
    context = "Serum of preterm infants",
    evidence_text = list(IFN = "Enriched: interferon alpha/beta signaling (q = 1e-6)"),
    provider = "manual",
    verbose = FALSE
  )

  expect_s3_class(req, "hc_llm_request")
  expect_match(req$id, "^hcocena-")
  expect_equal(req$part_task, c("features", "context", "evidence"))

  features_only <- req$parts[[1]]
  expect_match(features_only, "You are an expert proteomics analyst.", fixed = TRUE)
  expect_match(features_only, "### Set: IFN", fixed = TRUE)
  expect_match(features_only, "STAT1, IRF7", fixed = TRUE)
  expect_match(features_only, "### Set: Bcell", fixed = TRUE)
  expect_match(features_only, req$id, fixed = TRUE)
  expect_match(features_only, "general_processes", fixed = TRUE)
  # general_processes must not see the context or the evidence
  expect_false(grepl("Serum of preterm infants", features_only, fixed = TRUE))
  expect_false(grepl("interferon alpha/beta", features_only, fixed = TRUE))
  expect_false(grepl("contextual_state", features_only, fixed = TRUE))

  with_context <- req$parts[[2]]
  expect_match(with_context, "Serum of preterm infants", fixed = TRUE)
  expect_match(with_context, "key_regulators", fixed = TRUE)
  expect_false(grepl("interferon alpha/beta", with_context, fixed = TRUE))

  with_evidence <- req$parts[[3]]
  expect_match(with_evidence, "Serum of preterm infants", fixed = TRUE)
  expect_match(with_evidence, "interferon alpha/beta signaling", fixed = TRUE)
  expect_match(with_evidence, "grounded_processes", fixed = TRUE)
  expect_match(with_evidence, "supporting_evidence", fixed = TRUE)
})

test_that("without stored enrichment the evidence request is left out", {
  msgs <- testthat::capture_messages(
    req <- hc_llm_enrichment(genes = list(A = "X"), provider = "manual")
  )
  expect_true(any(grepl("No evidence request", msgs, fixed = TRUE)))
  expect_equal(req$part_task, c("features", "context"))
  # asked for explicitly, missing evidence is an error
  expect_error(
    hc_llm_enrichment(genes = list(A = "X"), evidence = "enrichment",
                      provider = "manual", verbose = FALSE),
    "hc_functional_enrichment"
  )
  req <- hc_llm_enrichment(genes = list(A = "X"), evidence = "none",
                           provider = "manual", verbose = FALSE)
  expect_equal(req$part_task, c("features", "context"))
})

test_that("manual round trip stores results exactly like an API run", {
  hc <- clustered_with_modules()
  hc <- hc_llm_enrichment(hc, provider = "manual", verbose = FALSE)
  pending <- hc_satellite(hc, "llm_enrichment_pending")
  expect_s3_class(pending, "hc_llm_request")
  modules <- names(pending$sets)
  expect_setequal(modules, c("M1", "M2", "M3", "M4"))

  hc <- hc_llm_import(hc, text = fake_answer(pending$id, modules),
                      model = "ChatGPT web", verbose = FALSE)

  summary_tbl <- hc_satellite(hc, "llm_enrichment_summary")
  expect_equal(nrow(summary_tbl), 4)
  expect_setequal(summary_tbl$module, modules)
  expect_true(all(summary_tbl$status == "ok"))
  expect_true(all(summary_tbl$llm == "manual"))
  expect_true(all(summary_tbl$model == "ChatGPT web"))
  expect_equal(summary_tbl$general_processes[summary_tbl$module == "M1"], "process of M1")
  expect_null(hc_satellite(hc, "llm_enrichment_pending"))

  p <- hc_plot_llm_enrichment(hc, with_heatmap = FALSE, save = FALSE)
  expect_named(p, c("general_processes", "contextual_state", "key_regulators"))
})

test_that("with evidence, the parts are imported one by one and give four plots", {
  hc <- clustered_with_modules()
  ev <- list(M1 = "Enriched: B cell receptor signaling", M2 = "Enriched: T cell activation",
             M3 = "Enriched: none", M4 = "Enriched: none")
  hc <- hc_llm_enrichment(hc, context = "Blood", evidence_text = ev,
                          provider = "manual", verbose = FALSE)
  req <- hc_satellite(hc, "llm_enrichment_pending")
  modules <- names(req$sets)
  expect_length(req$parts, 3)

  only <- function(fields) {
    entries <- vapply(modules, function(m) {
      vals <- c(general_processes = paste("process of", m), contextual_state = "state",
                key_regulators = "REG1", grounded_processes = paste("grounded", m),
                supporting_evidence = "TERM1")[fields]
      paste0('{"set": "', m, '", ', paste0('"', names(vals), '": "', vals, '"', collapse = ", "), "}")
    }, character(1))
    paste0('{"request_id": "', req$id, '", "results": [', paste(entries, collapse = ", "), "]}")
  }
  hc <- hc_llm_import(hc, text = only("general_processes"), verbose = FALSE)
  expect_equal(hc_satellite(hc, "llm_enrichment_pending")$done, 1L)
  hc <- hc_llm_import(hc, text = only(c("grounded_processes", "supporting_evidence")), verbose = FALSE)
  expect_equal(hc_satellite(hc, "llm_enrichment_pending")$done, c(1L, 3L))
  hc <- hc_llm_import(hc, text = only(c("contextual_state", "key_regulators")), verbose = FALSE)
  expect_null(hc_satellite(hc, "llm_enrichment_pending"))

  summary_tbl <- hc_satellite(hc, "llm_enrichment_summary")
  expect_true(all(summary_tbl$status == "ok"))
  expect_equal(summary_tbl$grounded_processes[summary_tbl$module == "M1"], "grounded M1")
  expect_equal(unique(summary_tbl$evidence), "supplied")
  p <- hc_plot_llm_enrichment(hc, with_heatmap = FALSE, save = FALSE)
  expect_named(p, c("general_processes", "contextual_state", "key_regulators", "grounded_processes"))
})

test_that("import tolerates code fences, surrounding text and typographic quotes", {
  req <- hc_llm_enrichment(genes = list(IFN = c("STAT1", "IRF7")),
                           provider = "manual", verbose = FALSE)
  json <- fake_answer(req$id, "IFN")
  pasted <- paste0(
    "Sure! Here is the interpretation you asked for.\n\n```json\n",
    gsub('"', "\u201c", json, fixed = TRUE),
    "\n```\nLet me know if you need anything else."
  )
  # typographic quotes are opening-only here; mixing both kinds is also fine
  pasted <- gsub("\u201c(general_processes|set)", "\u201d\\1", pasted)

  out <- hc_llm_import(req, text = pasted, verbose = FALSE)
  expect_equal(out$module, "IFN")
  expect_equal(out$general_processes, "process of IFN")
  expect_equal(out$key_regulators, "REG1 / REG2")
})

test_that("import refuses an answer that belongs to another request", {
  req <- hc_llm_enrichment(genes = list(IFN = c("STAT1")), provider = "manual",
                           verbose = FALSE)
  expect_error(
    hc_llm_import(req, text = fake_answer("hcocena-someone-else", "IFN"), verbose = FALSE),
    "belongs to request"
  )
  expect_error(
    hc_llm_import(req, text = "I could not process this request.", verbose = FALSE),
    "Could not find the JSON answer"
  )
  no_id <- '{"results": [{"set": "IFN", "general_processes": "x", "contextual_state": "y", "key_regulators": "z"}]}'
  expect_warning(
    out <- hc_llm_import(req, text = no_id, verbose = FALSE),
    "no `request_id`"
  )
  expect_equal(out$module, "IFN")
})

test_that("set names are matched case-insensitively and unknown sets are reported", {
  req <- hc_llm_enrichment(genes = list(M1 = "A", M2 = "B"), provider = "manual",
                           verbose = FALSE)
  answer <- sub('"set": "M1"', '"set": " m1 "', fake_answer(req$id, c("M1", "M2", "M9")))
  expect_warning(out <- hc_llm_import(req, text = answer, verbose = FALSE), "unknown sets: M9")
  expect_setequal(out$module, c("M1", "M2"))
})

test_that("long requests are split into parts that are imported one after another", {
  sets <- stats::setNames(lapply(seq_len(45), function(i) paste0("G", i, "_", 1:5)),
                          paste0("S", seq_len(45)))
  hc <- hc_example_data("clustered")
  hc <- hc_llm_enrichment(hc, genes = sets, provider = "manual", verbose = FALSE)
  req <- hc_satellite(hc, "llm_enrichment_pending")
  expect_equal(req$part_task, c("features", "features", "context", "context"))
  first <- req$part_sets[[1]]
  second <- req$part_sets[[2]]
  expect_equal(length(first) + length(second), 45)
  expect_match(req$parts[[1]], "part 1 of 4", fixed = TRUE)

  # one answer with all fields completes both request types of these sets
  hc <- hc_llm_import(hc, text = fake_answer(req$id, first), verbose = FALSE)
  still <- hc_satellite(hc, "llm_enrichment_pending")
  expect_equal(still$done, c(1L, 3L))

  hc <- hc_llm_import(hc, text = fake_answer(req$id, second), verbose = FALSE)
  expect_null(hc_satellite(hc, "llm_enrichment_pending"))
  expect_equal(nrow(hc_satellite(hc, "llm_enrichment_summary")), 45)
})

test_that("import reads the answer from a file, e.g. inside Docker", {
  req <- hc_llm_enrichment(genes = list(IFN = "STAT1"), provider = "manual",
                           verbose = FALSE)
  f <- withr::local_tempfile(fileext = ".txt")
  writeLines(fake_answer(req$id, "IFN"), f)
  out <- hc_llm_import(req, file = f, verbose = FALSE)
  expect_equal(out$module, "IFN")
})

test_that("building a request does not touch the user's random number stream", {
  set.seed(1)
  before <- .Random.seed
  hc_llm_enrichment(genes = list(A = "X"), provider = "manual", verbose = FALSE)
  expect_identical(.Random.seed, before)
})

# ---- API route (fake provider) -------------------------------------------------

test_that("API mode sends separate requests and returns a table", {
  seen <- list()
  local_mocked_bindings(
    .hc_llm_request_by_provider = function(llm, api_key, model, prompt,
                                           system_instruction, response_schema, ...) {
      seen[[length(seen) + 1]] <<- list(llm = llm, key = api_key, model = model,
                                         prompt = prompt, system = system_instruction,
                                         fields = response_schema$required)
      list(result_text = all_fields_json, raw_response_text = "raw")
    }
  )
  out <- hc_llm_enrichment(
    genes = list(IFN = c("STAT1", "IRF7"), Bcell = c("CD19")),
    system_prompt = "You are an expert proteomics analyst.",
    context = "Serum of preterm infants",
    provider = "claude",
    api_key = "test-key",
    verbose = FALSE
  )

  expect_s3_class(out, "data.frame")
  expect_equal(out$module, c("IFN", "Bcell"))
  expect_true(all(out$status == "ok"))
  expect_equal(out$general_processes, rep("interferon signaling", 2))
  expect_true(all(is.na(out$grounded_processes)))
  # two sets x two request types (no evidence available for `genes`)
  expect_length(seen, 4)
  expect_equal(seen[[1]]$llm, "claude")
  expect_equal(seen[[1]]$key, "test-key")
  expect_equal(seen[[1]]$model, "claude-sonnet-4-6")
  expect_equal(seen[[1]]$system, "You are an expert proteomics analyst.")
  expect_equal(seen[[1]]$fields, "general_processes")
  expect_false(grepl("Serum of preterm infants", seen[[1]]$prompt, fixed = TRUE))
  expect_match(seen[[1]]$prompt, "STAT1, IRF7", fixed = TRUE)
  expect_equal(seen[[2]]$fields, c("contextual_state", "key_regulators"))
  expect_match(seen[[2]]$prompt, "Serum of preterm infants", fixed = TRUE)
})

test_that("API mode with evidence adds a third request per set", {
  seen <- list()
  local_mocked_bindings(
    .hc_llm_request_by_provider = function(llm, api_key, model, prompt, system_instruction,
                                           response_schema, ...) {
      seen[[length(seen) + 1]] <<- list(prompt = prompt, fields = response_schema$required)
      list(result_text = all_fields_json, raw_response_text = "raw")
    }
  )
  out <- hc_llm_enrichment(
    genes = list(IFN = c("STAT1", "IRF7")), context = "Blood",
    evidence_text = list(IFN = "Enriched: interferon alpha/beta signaling"),
    provider = "openai", api_key = "k", verbose = FALSE
  )
  expect_length(seen, 3)
  expect_equal(seen[[3]]$fields, c("grounded_processes", "supporting_evidence"))
  expect_match(seen[[3]]$prompt, "interferon alpha/beta signaling", fixed = TRUE)
  expect_match(seen[[3]]$prompt, "Blood", fixed = TRUE)
  expect_equal(out$grounded_processes, "type I interferon response")
  expect_equal(out$supporting_evidence, "IFN term")
  expect_equal(out$status, "ok")
})

test_that("provider = 'local' talks to a vLLM server without an API key", {
  seen <- NULL
  local_mocked_bindings(
    .hc_llm_request_by_provider = function(llm, api_key, model, prompt,
                                           system_instruction, response_schema,
                                           temperature, timeout_sec, vllm_base_url) {
      seen <<- list(llm = llm, key = api_key, url = vllm_base_url, timeout = timeout_sec)
      list(result_text = all_fields_json, raw_response_text = "raw")
    }
  )
  withr::local_envvar(c(VLLM_API_KEY = NA, VLLM_BASE_URL = NA))
  hc_llm_enrichment(genes = "STAT1", provider = "local",
                    base_url = "http://gpu-server:8000/v1/", verbose = FALSE)
  expect_equal(seen$llm, "vllm")
  expect_equal(seen$key, "EMPTY")
  expect_equal(seen$url, "http://gpu-server:8000/v1")
  expect_equal(seen$timeout, 300)
})

test_that("a missing API key points to the manual route", {
  withr::local_envvar(c(OPENAI_API_KEY = NA))
  expect_error(
    hc_llm_enrichment(genes = "STAT1", provider = "openai", verbose = FALSE),
    'provider = "manual"',
    fixed = TRUE
  )
})

test_that("a failing set is recorded and warned about; if all fail it stops", {
  local_mocked_bindings(
    .hc_llm_request_by_provider = function(llm, api_key, model, prompt, ...) {
      if (grepl("Set: bad", prompt, fixed = TRUE)) stop("HTTP 429 rate limit")
      list(result_text = all_fields_json, raw_response_text = "raw")
    }
  )
  expect_warning(
    out <- hc_llm_enrichment(genes = list(good = "A", bad = "B"), provider = "gemini",
                             api_key = "k", verbose = FALSE),
    "HTTP 429"
  )
  expect_equal(out$status, c("ok", "error"))
  expect_match(out$error_message[2], "HTTP 429")

  expect_error(
    suppressWarnings(hc_llm_enrichment(genes = list(bad = "B"), provider = "gemini",
                                       api_key = "k", verbose = FALSE)),
    "All LLM requests failed"
  )
})

test_that("API results are stored in hc and can be plotted", {
  local_mocked_bindings(
    .hc_llm_request_by_provider = function(...) {
      list(result_text = all_fields_json, raw_response_text = "raw")
    }
  )
  hc <- clustered_with_modules()
  hc <- hc_llm_enrichment(hc, modules = c("M1", "M2"), provider = "openai",
                          api_key = "k", verbose = FALSE)
  summary_tbl <- hc_satellite(hc, "llm_enrichment_summary")
  expect_setequal(summary_tbl$module, c("M1", "M2"))
  expect_true(all(summary_tbl$llm == "openai"))
  expect_true(inherits(
    hc_plot_llm_enrichment(hc, fields = "general_processes", with_heatmap = FALSE, save = FALSE),
    "ggplot"
  ))
})

# ---- inputs and defaults ----------------------------------------------------------

test_that("the default system prompt is neutral and a custom one replaces it", {
  default <- hcocena:::.hc_llm_system_prompt(NULL)
  expect_false(grepl("monocyte", default, ignore.case = TRUE))
  expect_false(grepl("transcriptomics", default, ignore.case = TRUE))
  expect_match(default, "genes, proteins or metabolites", fixed = TRUE)
  expect_identical(hcocena:::.hc_llm_system_prompt("Custom."), "Custom.")
})

test_that("feature sets are validated", {
  expect_error(hc_llm_enrichment(provider = "manual"), "Provide `hc`")
  expect_error(hc_llm_enrichment(genes = list(A = character(0)), provider = "manual"),
               "contain no features: A")
  req <- hc_llm_enrichment(genes = list("X", B = "Y"), provider = "manual", verbose = FALSE)
  expect_equal(names(req$sets), c("set_1", "B"))
  single <- hc_llm_enrichment(genes = c("X", "Y", "X", " "), provider = "manual", verbose = FALSE)
  expect_equal(single$sets$custom_geneset, c("X", "Y"))
  expect_error(hc_llm_enrichment(genes = "X", provider = "chatbot"), "should be one of")
})

test_that("models that reject `temperature` are retried without it", {
  temps <- list()
  local_mocked_bindings(
    .hc_llm_request_by_provider_once = function(llm, api_key, model, prompt,
                                                system_instruction, response_schema,
                                                temperature, timeout_sec, vllm_base_url) {
      temps[[length(temps) + 1]] <<- if (is.null(temperature)) "none" else temperature
      if (!is.null(temperature)) {
        stop("Claude request failed via ellmer: HTTP 400 - `temperature` is deprecated for this model.")
      }
      list(result_text = all_fields_json, raw_response_text = "raw")
    }
  )
  out <- hc_llm_enrichment(genes = "STAT1", provider = "claude", api_key = "k",
                           model = "claude-opus-5", verbose = FALSE)
  expect_equal(out$status, "ok")
  # retried once per request (features, context)
  expect_equal(temps, list(0.2, "none", 0.2, "none"))
})

test_that("plot text spells out Greek letters that break Windows PNG devices", {
  ascii <- hcocena:::.hc_llm_ascii_symbols
  expect_equal(ascii("NF-\u03baB / IFN-\u03b3 \u2013 TNF-\u03b1"), "NF-kappaB / IFN-gamma - TNF-alpha")
  expect_equal(ascii("plain text"), "plain text")
})

test_that("plot text leaves out gene/protein lists and wraps onto at most two lines", {
  strip <- hcocena:::.hc_llm_strip_feature_lists
  expect_equal(
    strip("Neutrophil activation (PADI4, S100A12, CEACAM8) / Type I interferon response (IFIT1 / IFIT3)"),
    "Neutrophil activation / Type I interferon response"
  )
  # parentheses that explain something in words stay
  expect_equal(strip("Lectin pathway (complement activation)"), "Lectin pathway (complement activation)")

  long <- paste(
    "Proteasomal degradation and ubiquitin-dependent protein quality control /",
    "Mitochondrial electron transport chain (SDHA, SDHB, NDUFA5) /",
    "IL-6 and LIF family JAK-STAT inflammatory signalling"
  )
  shown <- hcocena:::.hc_llm_prepare_display_title(long, max_chars = 90)
  expect_false(grepl("SDHA", shown))
  lines <- strsplit(shown, "\n", fixed = TRUE)[[1]]
  expect_lte(length(lines), 2)
  expect_true(all(nchar(lines) <= 50))
  # nothing cut mid-phrase: what is shown are whole processes
  shown_parts <- trimws(strsplit(gsub("\n", " ", shown), " / ", fixed = TRUE)[[1]])
  expect_true(all(shown_parts %in% c(
    "Proteasomal degradation and ubiquitin-dependent protein quality control",
    "Mitochondrial electron transport chain",
    "IL-6 and LIF family JAK-STAT inflammatory signalling"
  )))
  expect_false(grepl("...", shown, fixed = TRUE))

  # short text stays on one line
  expect_equal(hcocena:::.hc_llm_prepare_display_title("Complement activation / Coagulation"),
               "Complement activation / Coagulation")

  # a line break never falls inside parentheses
  paren <- paste(
    "Acute phase response / Lectin pathway (complement activation via MBL and ficolins)",
    "/ Platelet degranulation"
  )
  for (ln in strsplit(hcocena:::.hc_llm_prepare_display_title(paren), "\n", fixed = TRUE)[[1]]) {
    expect_equal(lengths(regmatches(ln, gregexpr("(", ln, fixed = TRUE))),
                 lengths(regmatches(ln, gregexpr(")", ln, fixed = TRUE))))
  }

  # a single overlong process still yields at most two lines
  one <- paste(rep("extremely long biological process description", 6), collapse = " ")
  one_lines <- strsplit(hcocena:::.hc_llm_prepare_display_title(one), "\n", fixed = TRUE)[[1]]
  expect_lte(length(one_lines), 2)
  expect_true(all(nchar(one_lines) <= 50))
})

test_that("the prompt asks for processes without individual gene or protein names", {
  req <- hc_llm_enrichment(genes = list(A = "X"), provider = "manual", verbose = FALSE)
  expect_match(req$parts[[1]], "do not name individual genes or proteins", fixed = TRUE)
  expect_match(hcocena:::.hc_llm_task_prompt("features", "A", "X", "", NULL),
               "do not name individual genes or proteins", fixed = TRUE)
  expect_match(req$parts[[1]], "at most about 90 characters", fixed = TRUE)
})

test_that("evidence choices are validated and enrichment terms are rendered", {
  levels <- hcocena:::.hc_llm_evidence_levels
  expect_equal(levels("none"), character(0))
  expect_equal(levels("all"), c("enrichment", "upstream", "rag"))
  expect_equal(levels(c("rag", "enrichment")), c("enrichment", "rag"))
  expect_error(levels("pubmed"), "Unknown `evidence`")

  terms <- data.frame(term = c("Interferon alpha response", "Inflammatory response"),
                      database = "Hallmark", qvalue = c(1e-8, 0.003), GeneRatio = c("12/40", "5/40"))
  txt <- hcocena:::.hc_llm_format_enrichment_context(terms)
  expect_match(txt, "Interferon alpha response [Hallmark]", fixed = TRUE)
  expect_match(txt, "genes=12/40", fixed = TRUE)
  expect_null(hcocena:::.hc_llm_format_enrichment_context(terms[0, ]))

  expect_error(hcocena:::.hc_llm_rag_url(NULL), "HCOCENA_RAG_URL")
})

test_that("stored enrichment and upstream results become evidence for each module", {
  hc <- clustered_with_modules()
  sat <- as.list(hc@satellite)
  sat$enrichments <- list(significant_enrichments_all_dbs = data.frame(
    module_label = c("M1", "M1", "M2"), term = c("B cell receptor signaling", "Weak term", "T cell activation"),
    database = "Kegg", qvalue = c(1e-5, 0.2, 1e-4), GeneRatio = "5/20"
  ))
  sat$upstream_inference <- list(significant_upstream_all = data.frame(
    module_label = "M1", term = "PAX5", resource = "TF", database = "CollecTRI", qvalue = 1e-3,
    n_overlap = 6, direction = "activated", peak_condition = "severe", regulator_in_module = TRUE,
    redundant_with = ""
  ))
  hc@satellite <- S4Vectors::SimpleList(sat)
  hc <- hc_llm_enrichment(hc, evidence = c("enrichment", "upstream"), provider = "manual", verbose = FALSE)
  req <- hc_satellite(hc, "llm_enrichment_pending")
  ev <- req$evidence
  expect_equal(ev$levels, c("enrichment", "upstream"))
  expect_match(ev$text[["M1"]], "B cell receptor signaling", fixed = TRUE)
  expect_false(grepl("Weak term", ev$text[["M1"]], fixed = TRUE))
  expect_match(ev$text[["M1"]], "PAX5 [TF, CollecTRI]", fixed = TRUE)
  expect_match(ev$text[["M1"]], "regulator itself is in this set", fixed = TRUE)
  expect_match(ev$text[["M3"]], "none significant", fixed = TRUE)
  expect_match(req$parts[[3]], "statistical tests on exactly this set", fixed = TRUE)
})
