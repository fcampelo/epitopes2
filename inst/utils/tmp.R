library(dplyr)
library(tidymodels)
library(recipes)
library(praznik)
library(doParallel)

mydata <- readRDS("./data/processed_data_example.rds")

# Register parallel back end
idx <- rowSums(mydata$splits.summary$target_id_alloc_summary) > 0
cl <- makePSOCKcluster(sum(idx))
registerDoParallel(cl = cl)

top_p <- c(10, 20, 40, 80, 160)
top_p <- top_p[length(top_p):1]

df <- mydata$df %>% mutate(Class =  factor(Class == 1,
                                           labels = c("TRUE", "FALSE"),
                                           levels = c("TRUE", "FALSE"),
                                           ordered = TRUE))

for (i in seq_along(top_p)){
  cat(sprintf("\nStarted %s\n%04d features",
              as.character(Sys.time()), top_p[i]))

  MRMR <- recipe(Class ~ ., data = df) %>%
    step_rm(!starts_with("feat")&!matches("Class|Info_split")) %>%
    update_role("Info_split", new_role = "splitting variable") %>%
    step_zv(all_numeric_predictors(),
            id = "Remove constant predictors") %>%
    step_select_mrmr(all_numeric_predictors(),
                     outcome = "Class",
                     top_p = top_p[i],
                     id = paste0("MRMR", top_p[i]),
                     threads = ncpus)

  MRMR   <- MRMR %>% prep()
  Xfs <- MRMR %>% bake(new_data = NULL)

  cat("\nTuning model parameters...")

  folds <- group_vfold_cv(Xfs, group = Info_variant_folds, v = length(unique(Xfs$Info_variant_folds)))
  new.rec <- recipe(Binding ~ ., data = Xfs) %>%
    step_upsample(matches("Binding"), over_ratio = 0.2) %>%
    step_downsample(matches("Binding"), under_ratio = 1) %>%
    step_range(all_numeric_predictors(), min = 0, max = 1, id = "Scaling")

  xb.mod <- boost_tree(tree_depth = tune(), trees = tune(), min_n = tune()) %>%
    set_engine("lightgbm") %>%
    set_mode("classification")

  my_wf <- workflow() %>%
    add_recipe(new.rec) %>%
    add_model(xb.mod)

  search_res <- my_wf %>%
    tune_bayes(
      resamples = folds,
      initial = 15, # Generate semi-random initial solutions
      iter = 30,
      # How to measure performance?
      metrics = yardstick::metric_set(mcc, f_meas, accuracy,
                                      ppv, npv, sens, spec,
                                      roc_auc, pr_auc),
      control = control_bayes(no_improve = 15,
                              parallel_over = "resamples",
                              verbose = FALSE))

  # Collect CV metrics and select best configuration according to MCC
  CV_metrics <- search_res %>% collect_metrics()
  best_params <- search_res %>% select_best(metric = "mcc")

  saveRDS(list(CV_metrics  = CV_metrics,
               best_params = best_params),
          savefile)
  cat("\nDone!")
}
