## hc_plot_enrichment_upstream_network(): regulation sign, term redundancy,
## term -> regulator links and the page layout.

test_that("regulation follows the sign of the targets' mode of regulation", {
  net <- data.frame(
    source = c(rep("A", 3), rep("R", 3), rep("X", 4)),
    target = c("g1", "g2", "g3", "g1", "g2", "g3", "g1", "g2", "g3", "g4"),
    mor = c(1, 1, 1, -1, -1, 1, 1, 1, -1, -1)
  )
  ora <- data.frame(term = c("A", "R", "X"), overlap_genes = c("g1,g2,g3", "g1,g2,g3", "g1,g2,g3,g4"))
  expect_equal(hcocena:::.hc_ui_module_regulation(ora, net), c("activating", "repressing", "mixed"))
})

test_that("terms lying mostly inside a better term are collapsed", {
  df <- data.frame(
    cluster = "m1",
    term = c("big", "nested", "other"),
    qvalue = c(1e-8, 1e-6, 1e-4),
    geneID = c(paste0("g", 1:20, collapse = "/"), paste0("g", 1:5, collapse = "/"),
               paste0("h", 1:5, collapse = "/"))
  )
  out <- hcocena:::.hc_kn_collapse_terms(df, contained = 0.8)
  expect_equal(out$term, c("big", "other"))
  # without gene lists nothing is removed
  expect_equal(nrow(hcocena:::.hc_kn_collapse_terms(df[, c("cluster", "term", "qvalue")])), 3)
})

test_that("term -> regulator links prefer specific over large generic terms", {
  enrich <- data.frame(
    cluster = "m1", node_key = c("db||generic", "db||specific"),
    geneID = c(paste0("g", 1:90, collapse = "/"), paste0("g", 1:10, collapse = "/"))
  )
  up <- data.frame(cluster = "m1", node_key = "TF||E2F", n_genes = 100,
                   overlap_genes = paste0("g", 1:8, collapse = ","))
  links <- hcocena:::.hc_kn_term_regulator_links(enrich, up, min_share = 0.25)
  # all 8 targets are in both terms, but only the small term is enriched for them
  expect_equal(links$term_key, "db||specific")
  expect_equal(links$n_shared, 8)
  expect_gt(links$ratio, 1.5)
})

test_that("curves start and end horizontally at the given points", {
  cv <- hcocena:::.hc_kn_curves(0, 0, 10, 4, "a", n = 11)
  expect_equal(cv$x[c(1, 11)], c(0, 10))
  expect_equal(cv$y[c(1, 11)], c(0, 4))
  expect_lt(abs(cv$y[2] - cv$y[1]), abs(cv$y[6] - cv$y[5]))
  expect_equal(hcocena:::.hc_kn_pretty_term(c("HALLMARK_E2F_TARGETS", "GOBP_CELL_CYCLE", "custom")),
               c("E2F TARGETS", "CELL CYCLE", "custom"))
})

test_that("the page keeps square heatmap cells and marks regulator effects", {
  hm <- list(
    mat = matrix(c(1, -1, 0.5, -0.5), 2, dimnames = list(c("red", "blue"), c("c1", "c2"))),
    keep_clusters = c("red", "blue"), module_labels = c("M1", "M2"), gene_counts = c(40, 25),
    value_name = "GFC", gfc_colors = c("blue", "white", "red"), scale_limits = c(-2, 2),
    scale_breaks = c(-2, 0, 2), scale_labels = c("-2", "0", "2")
  )
  enrich <- data.frame(cluster = c("red", "blue"), node_key = c("Go||A", "Kegg||B"),
                       database = c("Go", "Kegg"), term = c("GOBP_A", "KEGG_B"), qvalue = c(1e-5, 1e-3))
  up <- data.frame(cluster = c("red", "blue"), node_key = c("TF||E2F1", "TF||RB1"), resource = "TF",
                   term = c("E2F1", "RB1"), qvalue = c(1e-6, 1e-2), direction = c("activated", "inhibited"),
                   regulator_module = c("M1", NA))
  layout <- list(x = c(1, 2), limits = c(0.5, 2.5))
  page <- hcocena:::.hc_kn_page_plot(hm, c("c1", "c2"), layout, enrich, up, data.frame(),
                                     title = "test")
  expect_s3_class(page$plot, "ggplot")
  expect_identical(page$plot$coordinates$ratio, 1)
  expect_gt(page$width, page$height)
  built <- ggplot2::ggplot_build(page$plot)
  fills <- unlist(lapply(built$data, function(d) if ("fill" %in% names(d)) d$fill))
  expect_true(all(c("#C0392B", "#2471A3") %in% fills))  # activating / repressing regulators
  labels <- unlist(lapply(built$data, function(d) if ("label" %in% names(d)) as.character(d$label)))
  expect_true(all(c("40", "25", "E2F1", "A") %in% labels))

  focus <- hcocena:::.hc_kn_page_plot(hm, c("c1", "c2"), layout, enrich, up, data.frame(),
                                      focus = "red", title = "focus")
  expect_s3_class(focus$plot, "ggplot")
})
