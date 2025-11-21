# Normalize per experiment

normalization_per_exp <- function(data, column_name) {
 
  column_data <- data[[column_name]]
  
  
  # Compute medians for "Functional" and "Non-functional"
  post_pre_median_syn <- median(column_data[data$expected_clinical == "Functional"], na.rm = TRUE)
  post_pre_median_ns <- median(column_data[data$expected_clinical == "Non-functional"], na.rm = TRUE)
  
  # Handle cases where Non-functional is missing
  if (is.na(post_pre_median_ns) | nrow(data[data$expected_clinical == "Non-functional", ]) < 1) {
    post_pre_median_ns <- min(data[[column_name]], na.rm = TRUE)
    Q5.low <- quantile(data[[column_name]], probs = 0.05, na.rm = TRUE)
    post_pre_median_ns <- median(data[[column_name]][data[[column_name]] < Q5.low], na.rm = TRUE)
  }
  
  normalized_column <- data[[column_name]] / (post_pre_median_syn - post_pre_median_ns)
  return(normalized_column)
}

# Extend fit for LOESS
extend_fit <- function(fit_vector) {
  if (length(which(is.na(fit_vector))) == 0) {
    out_vector <- fit_vector
  }
  else {
    min_na <- min(which(is.na(fit_vector) == TRUE))
    min_point <- min(which(is.na(fit_vector) == FALSE))
    min_point_val <- fit_vector[min_point]
    max_point <- max(which(is.na(fit_vector) == FALSE))
    max_point_val <- fit_vector[max_point]
    max_na <- max(which(is.na(fit_vector) == TRUE))
    out_vector <- fit_vector
    if (min_na == 1) {
      out_vector[min_na:min_point] <- min_point_val
    }
    if (max_na == length(fit_vector)) {
      out_vector[(max_point+1):length(fit_vector)] <- max_point_val 
    }
  }
  return(out_vector)
}

# Get positional biais
position_biais <- function(data, col){
  posrange <- range(as.numeric(as.character(data$pos)))
  posseq <- seq(from = posrange[1], to = posrange[2])
  Q1 <- quantile(data[[col]], probs = 0.25, na.rm = T)
  
  post_pre_ratio.lo <- loess(log2(data[[col]]) ~ pos, 
                             data = data, 
                             span = 0.75,
                             model = TRUE)
  post_pre_pred <- predict(post_pre_ratio.lo, newresult = data$posseq, se = TRUE)
  post_pre_pos_effects <- post_pre_pred$fit
  print(post_pre_pos_effects)
  post_pre_pos_effects.x <- extend_fit(post_pre_pos_effects)
  return(log2(data[[col]]) - sapply(data$pos, function(name) post_pre_pos_effects.x[[name]]))
}

# Normalize per exon
per_exon_norm <- function(data){
  subset_post_pre_median_syn <- median(data$post_pre_ratio_sns_loess[data$Outcome == "Benign"], na.rm = T)
  subset_post_pre_median_ns <- median(data$post_pre_ratio_sns_loess[data$Outcome == "Non-functional"], na.rm = T)
  if (is.na(subset_post_pre_median_ns) | nrow(data[data$Outcome == "Non-functional",]) < 1){
    Q5.low <- quantile(data$post_pre_ratio_sns_loess, probs = 0.05, na.rm = T)
    subset_post_pre_median_ns <- median(data$post_pre_ratio_sns_loess[data$post_pre_ratio_sns_loess < Q5.low], na.rm = T)
  }
  
  return(data$post_pre_ratio_sns_loess - subset_post_pre_median_ns) / (post_pre_median_syn - post_pre_median_ns)
}
