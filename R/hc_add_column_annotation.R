#' Categorical Metadata Annotation For The Module Heatmap
#'
#' Annotates the sample groups in the module heatmap with categorical metadata from the annotations.
#'  The annotations will be shown in the heatmap when calling
#'  `hc_plot_cluster_heatmap()` the next time.
#' @param variables A vector of strings. The length of the vector must equal the number of datasets.
#'  For each dataset, give the name of the column in the corresponding annotation table, that contains the metadata of interest.
#'  E.g., if you want to show sex, and in the first annotationtable the coulm name is "sex" and in the second annotation table it's by mistake labelled as "gender"",
#'  then you set variables = c("sex", "gender").
#'  If the information only exists in some but not all datasets, then set the vector slots of those that don"t have it to NA, e.g., c(NA, "gender").
#' @param variable_label A string that describes the meta information used. This string will be used to label the annotation in the heatmap.
#'  If none is provided, the first non-NA string from 'variables' is used as the label.
#' @param type A string defining whether meta variables should be converted into percentages ("percent") or absolute ("abs"; default) values
#' @noRd

.hc_col_anno_categorical_driver <- function(variables, variable_label = NULL, type = "abs") {
  if (base::is.null(variables) || base::length(variables) == 0) {
    annos <- hcobject[["satellite_outputs"]][["column_annos_categorical"]]
    if (base::is.null(variable_label) || base::length(variable_label) == 0) {
      .hc_set_bridge_hcobject_slot(c("satellite_outputs", "column_annos_categorical"), NULL)
      message("Removed all categorical heatmap column annotations.")
    } else {
      variable_label <- base::as.character(variable_label[[1]])
      if (!base::is.null(annos) && variable_label %in% base::names(annos)) {
        annos[[variable_label]] <- NULL
        if (base::length(annos) == 0) {
          annos <- NULL
        }
        .hc_set_bridge_hcobject_slot(c("satellite_outputs", "column_annos_categorical"), annos)
        message("Removed categorical heatmap column annotation '", variable_label, "'.")
      } else {
        message("No categorical heatmap column annotation named '", variable_label, "' was found.")
      }
    }
    return(invisible(NULL))
  }

  ml <- get_anno_matrix(variables = variables)
  m <- unify_mats(ml)
  if (type == "percent") {
    m <- base::apply(m, 1, function(x) {
      if (base::sum(x) == 0) {
        x
      } else {
        x / base::sum(x) * 100
      }
    }) %>%
      base::t() %>%
      base::as.data.frame()
  }
  m[base::rowSums(m) == 0, ] <- NA
  if (base::is.null(variable_label)) {
    variable_label <- variables[!base::is.null(variables)][1]
    message("No variable label was provided. Variable will automatically be labelled as '", variable_label, "' in hcobject[['satellite_outputs']][['column_annos_categorical']].")
  }
  if (base::length(variable_label) > 1) {
    variable_label <- variable_label[1]
    message("More than one label has been provided (parameter 'variable_label'). Only the first argument will be used as the label.")
  }
  attr(m, "type") <- type
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "column_annos_categorical", variable_label), m)
}


#' Numerical Metadata Annotation For The Module Heatmap
#'
#' Annotates the sample groups in the module heatmap with numerical metadata from the annotations.
#'  The annotations will be shown in the heatmap when calling
#'  `hc_plot_cluster_heatmap()` the next time.
#' @param variables A vector of strings. The length of the vector must equal the number of datasets.
#'  For each dataset, give the name of the column in the corresponding annotation table, that contains the metadata of interest.
#'  E.g., if you want to show age, and in the first annotationtable the coulm name is "age" and in the second annotation table is "Age",
#'  then you set variables = c("age", "Age").
#'  If the information only exists in some but not all datasets, then set the vector slots of those that don"t have it to NA, e.g., c(NA, "Age").
#' @param variable_label A string that describes the meta information used. This string will be used to label the annotation in the heatmap.
#'  If none is provided, the first non-NA string from 'variables' is used as the label.
#' @noRd

.hc_col_anno_numerical_driver <- function(variables, variable_label) {
  if (base::is.null(variables) || base::length(variables) == 0) {
    annos <- hcobject[["satellite_outputs"]][["column_annos_numerical"]]
    if (base::is.null(variable_label) || base::length(variable_label) == 0) {
      .hc_set_bridge_hcobject_slot(c("satellite_outputs", "column_annos_numerical"), NULL)
      message("Removed all numerical heatmap column annotations.")
    } else {
      variable_label <- base::as.character(variable_label[[1]])
      if (!base::is.null(annos) && variable_label %in% base::names(annos)) {
        annos[[variable_label]] <- NULL
        if (base::length(annos) == 0) {
          annos <- NULL
        }
        .hc_set_bridge_hcobject_slot(c("satellite_outputs", "column_annos_numerical"), annos)
        message("Removed numerical heatmap column annotation '", variable_label, "'.")
      } else {
        message("No numerical heatmap column annotation named '", variable_label, "' was found.")
      }
    }
    return(invisible(NULL))
  }

  vals <- list()
  for (i in base::seq_along(variables)) {
    anno <- hcobject[["data"]][[base::paste0("set", i, "_anno")]]
    if (!variables[i] %in% colnames(anno) & !base::is.na(variables[i])) {
      stop(variables[i], " is not a column name found in the annotation of dataset ", i, ". Please check the spelling.")
    }
    if (!hcobject[["global_settings"]][["control"]] == "none") {
      ctrl <- base::unique(dplyr::pull(anno, hcobject[["global_settings"]][["voi"]]))[base::grepl(hcobject[["global_settings"]][["control"]],
        base::unique(dplyr::pull(anno, hcobject[["global_settings"]][["voi"]])),
        ignore.case = TRUE
      ) == TRUE]
    } else {
      ctrl <- base::c()
    }
    for (j in base::unique(dplyr::pull(anno, hcobject[["global_settings"]][["voi"]]))) {
      if (j %in% ctrl) {
        next
      }
      if (base::is.na(variables[i])) {
        vals[[j]] <- NA
      } else {
        tmp <- anno[anno[hcobject[["global_settings"]][["voi"]]] == j, ]
        vals[[j]] <- base::as.numeric(base::as.character(dplyr::pull(tmp, variables[i])))
      }
    }
  }
  if (base::is.null(variable_label)) {
    variable_label <- variables[!base::is.null(variables)][1]
    message("No variable label was provided. Variable will automatically be labelled as '", variable_label, "' in hcobject[['satellite_outputs']][['column_annos_categorical']].")
  }
  if (base::length(variable_label) > 1) {
    variable_label <- variable_label[1]
    message("More than one label has been provided (parameter 'variable_label'). Only the first argument will be used as the label.")
  }
  .hc_set_bridge_hcobject_slot(c("satellite_outputs", "column_annos_numerical", variable_label), vals)
}




#' Add numeric module heatmap annotations (S4 API)
#'
#' Adds numeric annotation tracks to the columns of the module heatmap. Call
#' with `variables = NULL` to remove previously added tracks.
#'
#' @param hc A `HCoCenaExperiment`.
#' @param variables Character vector of numeric annotation column names to add.
#'   Use `NULL` to remove existing numeric annotations.
#' @param variable_label Display label for the annotation track.
#' @return Updated `HCoCenaExperiment`.
#' @examples
#' hc <- hc_example_data("clustered")
#' hc <- hc_col_anno_numerical(hc, variables = "batch", variable_label = "batch")
#' @export
hc_col_anno_numerical <- function(hc, variables, variable_label) {
  .hc_run_driver(
    hc = hc, fun = .hc_col_anno_numerical_driver,
    variables = variables, variable_label = variable_label
  )
}
