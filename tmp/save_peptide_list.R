#' Save information from a peptides.list object to files .
#'
#' Saves the full information from objects of type `peptides.list` (e.g., output
#' of functions [extract_peptides()], [make_data_splits()] or [calc_features()])
#' to a set of CSV files for easier sharing and cross-platform use.
#'
#' @param peptides.list object of class `peptides.list`.
#' @param save_folder path to folder for saving the resulting CSV files.
#'
#' @return This function returns (invisibly) the same input list,
#' `peptides.list`, to enable insertion as part of `dplyr` pipelines.
#'
#' @author Felipe Campelo (\email{f.campelo@@aston.ac.uk})
#'
#' @importFrom utils write.csv
#' @export

save_peptide_list <- function(peptides.list, save_folder){

  # ========================================================================== #
  # Sanity checks and initial definitions
  assertthat::assert_that("peptide.list" %in% class(peptides.list),
                          is.character(save_folder),
                          length(save_folder) == 1)

  if(!dir.exists(save_folder)) dir.create(save_folder, recursive = TRUE)
  # ========================================================================== #

  # save main data frame
  message("Saving main peptides data frame")
  write.csv(peptides.list$df, paste0(save_folder, "/peptides_df_main.csv"),
            row.names = FALSE)

  # save summary peptides data frame
  message("Saving summary peptides data frame")
  write.csv(peptides.list$peptides, paste0(save_folder, "/peptides_df_summary.csv"),
            row.names = FALSE)

  # save summary proteins data frame
  if("proteins" %in% names(peptides.list)){
    message("Saving proteins data frame")
    write.csv(peptides.list$proteins, paste0(save_folder, "/proteins_df.csv"),
              row.names = FALSE)
  }

  # save function call parameters
  message("Saving function call parameters")
  fcpars <- rbind(
    # filter attributes
    data.frame(Function  = "filter_epitopes",
               Parameter = names(peptides.list$filter.attrs),
               Value     = unname(unlist(peptides.list$filter.attrs))),
    # consolidate attributes
    data.frame(Function  = "consolidate_data",
               Parameter = names(peptides.list$consolidate.attrs),
               Value     = unname(unlist(peptides.list$consolidate.attrs))),
    # peptide extraction attributes
    data.frame(Function  = "extract_peptides",
               Parameter = names(peptides.list$peptide.attrs),
               Value     = unname(unlist(peptides.list$peptide.attrs))))

  # Add split parameters, if available
  if("splits.attrs" %in% names(peptides.list)){
    tmp <- peptides.list$splits.attrs[which(names(peptides.list$splits.attrs) %in% names(formals(epitopes::make_data_splits)))]
    tmp$target_props <- paste(tmp$target_props, collapse = ",")
    tmp <- lapply(tmp, function(x) ifelse(is.null(x), "NULL", x))
    fcpars <- rbind(fcpars,
                    data.frame(Function  = "make_data_splits",
                               Parameter = names(tmp),
                               Value     = unname(unlist(tmp))))
  }

  # Add feature calculation parameters, if available
  if("feature.attrs" %in% names(peptides.list)){
    tmp <- peptides.list$feature.attrs[which(names(peptides.list$feature.attrs) %in% names(formals(epitopes::calc_features)))]
    tmp$local.features <- paste(tmp$local.features, collapse = ",")
    tmp$global.features <- paste(tmp$global.features, collapse = ",")
    tmp <- lapply(tmp, function(x) ifelse(is.null(x), "NULL", x))
    fcpars <- rbind(fcpars,
                    data.frame(Function  = "calc_features",
                               Parameter = names(tmp),
                               Value     = unname(unlist(tmp))))
  }

  # Add modelling parameters, if available
  if("model.attrs" %in% names(peptides.list)){
    tmp <- peptides.list$model.attrs[which(names(peptides.list$model.attrs) %in% names(formals(epitopes::fit_model)))]
    tmp <- lapply(tmp, function(x) ifelse(is.null(x), "NULL", x))
    fcpars <- rbind(fcpars,
                    data.frame(Function  = "fit_model",
                               Parameter = names(tmp),
                               Value     = unname(unlist(tmp))))
    if(length(peptides.list$model.attrs$other.args) > 0){
      fcpars <- rbind(fcpars,
                      data.frame(Function  = "fit_model",
                                 Parameter = names(peptides.list$model.attrs$other.args),
                                 Value     = unname(unlist(peptides.list$model.attrs$other.args))))
    }
  }

  write.csv(fcpars, paste0(save_folder, "/function_call_parameters.csv"),
            row.names = FALSE)


  # Save further split outcomes, if available
  if("splits.attrs" %in% names(peptides.list)){
    message("Saving split data")
    tmp <- peptides.list$splits.attrs[which(!(names(peptides.list$splits.attrs) %in% names(formals(epitopes::make_data_splits))))]

    df1 <- data.frame(split_props    = tmp$split_props,
                      target_props   = unname(peptides.list$splits.attrs$target_props),
                      split_balance  = unname(tmp$split_balance),
                      target_balance = unname(tmp$target_balance))

    write.csv(df1, paste0(save_folder, "/split_properties.csv"),
              row.names = FALSE)

    write.csv(tmp$cluster.alloc, paste0(save_folder, "/clusters_per_split.csv"),
              row.names = FALSE)

    write.csv(tmp$diss.matrix, paste0(save_folder, "/dissimilarity_matrix.csv"))

    write.csv(tmp$SW.scores, paste0(save_folder, "/SW_scores.csv"))

    saveRDS(tmp$clusters, paste0(save_folder, "/cluster_structure.rds"))

  }


  # save model performance, if available
  if("model.perf" %in% names(peptides.list)){
    message("Saving model information")
    tmp <- peptides.list$model.perf
    tmp$roc <- NULL
    write.csv(as.data.frame(tmp), paste0(save_folder, "/model_performance.csv"),
              row.names = FALSE)
    write.csv(peptides.list$model.perf$roc, paste0(save_folder, "/model_performance_ROC.csv"),
              row.names = FALSE)
  }

  # save model performance, if available
  if("model" %in% names(peptides.list) & !is.null(peptides.list$model)){
    message("Saving model")
    saveRDS(peptides.list$model, paste0(save_folder, "/model.rds"))
  }

  invisible(peptides.list)
}
