## hc_upstream_inference(): module x regulator over-representation, target
## scores, regulator activity, redundancy and the regulator's own module.

quietly <- function(expr) suppressMessages(suppressWarnings(expr))

toy_universe <- function() paste0("G", 1:200)

toy_network <- function() {
  data.frame(
    source = c(rep("TF_A", 20), rep("TF_B", 20), rep("TF_C", 10)),
    target = c(paste0("G", 1:20), paste0("G", 101:120), paste0("G", c(1:6, 150:153))),
    mor = c(rep(1, 20), rep(1, 10), rep(-1, 10), rep(1, 10)),
    stringsAsFactors = FALSE
  )
}

test_that("module ORA uses the hypergeometric test against the universe", {
  mods <- list(m1 = paste0("G", 1:25), m2 = paste0("G", 51:75))
  ora <- hcocena:::.hc_ui_module_ora(mods, toy_network(), toy_universe(), minsize = 5, padj = "BH")

  a <- ora[ora$cluster == "m1" & ora$term == "TF_A", ]
  expect_equal(a$n_overlap, 20)
  expect_equal(a$n_targets, 20)
  expect_equal(a$n_genes, 25)
  expect_equal(a$pvalue, stats::phyper(19, 20, 180, 25, lower.tail = FALSE))
  expect_equal(a$fold_enrichment, (20 / 25) / (20 / 200))
  expect_lt(a$qvalue, 1e-10)

  # rows without overlap are dropped, but count for the multiple testing:
  # 2 modules x 3 regulators = 6 tests
  expect_false(any(ora$n_overlap == 0))
  all_p <- c(ora$pvalue, rep(1, 6 - nrow(ora)))
  expect_equal(sort(ora$qvalue), sort(stats::p.adjust(all_p, "BH")[seq_len(nrow(ora))]))
})

test_that("regulators with fewer than minsize targets in the universe are not tested", {
  net <- rbind(toy_network(), data.frame(source = "TF_tiny", target = c("G1", "G2"), mor = 1))
  ora <- hcocena:::.hc_ui_module_ora(list(m1 = paste0("G", 1:25)), net, toy_universe(), minsize = 5)
  expect_false("TF_tiny" %in% ora$term)
})

test_that("random modules give no significant links (calibration)", {
  set.seed(11)
  universe <- paste0("G", 1:2000)
  net <- do.call(rbind, lapply(1:60, function(i) {
    data.frame(source = paste0("TF", i), target = sample(universe, 40), mor = 1)
  }))
  hits <- vapply(1:20, function(b) {
    pool <- sample(universe)
    mods <- split(pool[1:600], rep(1:6, each = 100))
    ora <- hcocena:::.hc_ui_module_ora(mods, net, universe, minsize = 5)
    sum(ora$qvalue <= 0.05 & ora$n_overlap >= 3)
  }, numeric(1))
  expect_lte(mean(hits > 0), 0.1)

  # a planted module is found
  mods <- list(planted = c(net$target[net$source == "TF1"][1:25], sample(universe, 75)))
  ora <- hcocena:::.hc_ui_module_ora(mods, net, universe, minsize = 5)
  expect_equal(ora$term[which.min(ora$qvalue)], "TF1")
  expect_lt(min(ora$qvalue), 1e-6)
})

test_that("target score respects the mode of regulation and reports consistency", {
  ora <- data.frame(term = "TF_B", overlap_genes = paste0("G", c(101:103, 111:113), collapse = ","), stringsAsFactors = FALSE)
  values <- matrix(c(1, 1, 1, -1, -1, -1,     # all targets move as expected
                     1, 1, 1, 1, 1, 1),       # repressed targets move up
                   ncol = 2, dimnames = list(paste0("G", c(101:103, 111:113)), c("c1", "c2")))
  st <- hcocena:::.hc_ui_module_target_scores(ora, toy_network(), values)
  expect_equal(unname(st$score[1, "c1"]), 1)
  expect_equal(unname(st$score[1, "c2"]), 0)
  expect_equal(unname(st$consistency[1, "c1"]), 1)

  values[, "c1"] <- c(1, 1, -1, -1, -1, -1)
  st <- hcocena:::.hc_ui_module_target_scores(ora, toy_network(), values)
  expect_equal(unname(st$score[1, "c1"]), 4 / 6)
  expect_equal(unname(st$consistency[1, "c1"]), 5 / 6)
})

test_that("redundant regulators point to the better-ranked one", {
  ora <- data.frame(
    cluster = c("m1", "m1", "m1", "m2"),
    term = c("A", "B", "C", "B"),
    overlap_genes = c("g1,g2,g3,g4", "g1,g2,g3,g5", "g7,g8,g9", "g1,g2,g3,g5"),
    stringsAsFactors = FALSE
  )
  red <- hcocena:::.hc_ui_flag_redundant(ora, keep = rep(TRUE, 4), jaccard = 0.5)
  expect_equal(red, c("", "A", "", ""))
  red <- hcocena:::.hc_ui_flag_redundant(ora, keep = c(FALSE, TRUE, TRUE, TRUE), jaccard = 0.5)
  expect_equal(red, c("", "", "", ""))

  # a TF that is itself a gene of the module is kept as the representative
  ora$regulator_in_module <- c(FALSE, TRUE, FALSE, FALSE)
  red <- hcocena:::.hc_ui_flag_redundant(ora, keep = rep(TRUE, 4), jaccard = 0.5)
  expect_equal(red, c("B", "", "", ""))
})

test_that("regulator correlation needs the regulator and three conditions", {
  values <- matrix(c(1, 2, 3, 4, 3, 2, 1, 0), nrow = 2, byrow = TRUE,
                   dimnames = list(c("TF_A", "G9"), paste0("c", 1:4)))
  means <- matrix(c(1, 2, 3, 4), nrow = 1, dimnames = list("m1", paste0("c", 1:4)))
  expect_equal(hcocena:::.hc_ui_regulator_cor("TF_A", "m1", values, means), 1)
  expect_true(is.na(hcocena:::.hc_ui_regulator_cor("TF_X", "m1", values, means)))
  expect_true(is.na(hcocena:::.hc_ui_regulator_cor("TF_A", "m1", values[, 1:2], means)))
})

test_that("regulator activity per condition comes from ULM over all genes", {
  skip_if_not_installed("decoupleR")
  set.seed(3)
  genes <- toy_universe()
  values <- matrix(rnorm(400, sd = 0.3), ncol = 2, dimnames = list(genes, c("up", "flat")))
  values[paste0("G", 1:20), "up"] <- values[paste0("G", 1:20), "up"] + 2
  act <- quietly(hcocena:::.hc_ui_condition_activity(values, toy_network(), minsize = 5))
  a_up <- act[act$source == "TF_A" & act$condition == "up", ]
  a_flat <- act[act$source == "TF_A" & act$condition == "flat", ]
  expect_gt(a_up$activity, 0)
  expect_lt(a_up$activity_qvalue, 0.05)
  expect_gt(a_flat$activity_qvalue, 0.05)
})

test_that("PROGENy keeps the top genes per pathway", {
  skip_if_not_installed("progeny")
  net <- hcocena:::.hc_ui_load_pathway_network("human", top = 50)
  expect_true(all(table(net$source) == 50))
  expect_true(all(c("source", "target", "mor") %in% names(net)))
})

test_that("DoRothEA can be forced and is labelled", {
  skip_if_not_installed("dorothea")
  net <- hcocena:::.hc_ui_load_tf_network("human", tf_resource = "dorothea", tf_confidence = "A")
  expect_identical(attr(net, "database"), "DoRothEA")
  expect_gt(nrow(net), 100)
})

test_that("removed options explain what replaced them", {
  skip_if_not_installed("decoupleR")
  hc <- hc_example_data("clustered")
  expect_error(quietly(hc_upstream_inference(hc, method = "ulm", plot = FALSE, save_pdf = FALSE)),
               "over-representation")
  expect_error(quietly(hc_upstream_inference(hc, activity_input = "expression", plot = FALSE, save_pdf = FALSE)),
               "removed")
  expect_error(quietly(hc_upstream_inference(hc, min_overlap = 0, plot = FALSE, save_pdf = FALSE)),
               "min_overlap")
})

test_that("results carry the link, activity and regulator columns", {
  skip_if_not_installed("decoupleR")
  hc <- hc_example_data("clustered")
  gmt <- system.file("extdata", "toy_celltype_markers.gmt", package = "hcocena")
  skip_if(!nzchar(gmt), "The toy marker GMT is unavailable.")
  hc <- quietly(hc_upstream_inference(hc, resources = "Pathway", custom_pathway_gmt = c(CellTypes = gmt),
                                      minsize = 2, min_overlap = 1, plot = FALSE, save_pdf = FALSE))
  res <- hc_satellite(hc, "upstream_inference")
  by_cond <- as.data.frame(res$all_upstream_by_condition)
  expect_true(all(c("n_overlap", "fold_enrichment", "consistency", "regulator_in_module",
                    "redundant_with", "n_active_conditions") %in% names(as.data.frame(res$significant_upstream_all))))
  expect_true(all(c("activity", "activity_qvalue") %in% names(by_cond)))
  expect_identical(res$settings$method, "ora")
})
