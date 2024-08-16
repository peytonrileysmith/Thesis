auc_EN<-readRDS("df_auc_EN_sim.rds")
auc_Lasso<-readRDS("df_auc_Lasso_sim.rds")
auc_Ridge<-readRDS("df_auc_Ridge_sim.rds")
cor_EN<-readRDS("df_cor_EN_sim.rds")
cor_Lasso<-readRDS("df_cor_Lasso_sim.rds")
cor_Ridge<-readRDS("df_cor_Ridge_sim.rds")
library(dplyr)

calculate_5_point_diff <- function(df) {
  # Create an empty vector to store the differences
  diff_vector <- numeric()
  
  # Get the AUC column (assuming it's the second column)
  auc_values <- df[[2]]
  
  # Calculate the difference between every 5 points
  for (i in seq(1, length(auc_values), by = 5)) {
    if ((i + 4) <= length(auc_values)) {
      diff <- auc_values[i + 4] - auc_values[i]
      diff_vector <- c(diff_vector, diff)
    } else {
      # Append NA if there are not enough points to calculate the difference
      diff_vector <- c(diff_vector, NA)
    }
  }
  
  # Add the difference vector as a new column to the data frame
  df$change_in_auc <- rep(diff_vector, each = 5, length.out = nrow(df))
  
  return(df)
}

auc_EN <- calculate_5_point_diff(auc_EN)
auc_Lasso <- calculate_5_point_diff(auc_Lasso)
auc_Ridge <- calculate_5_point_diff(auc_Ridge)

calculate_5_point_diff_cor <- function(df) {
  # Create an empty vector to store the differences
  diff_vector <- numeric()
  
  # Get the correlation column (assuming it's the second column)
  cor_values <- df[[2]]
  
  # Calculate the difference between every 5 points
  for (i in seq(1, length(cor_values), by = 5)) {
    if ((i + 4) <= length(cor_values)) {
      diff <- cor_values[i + 4] - cor_values[i]
      diff_vector <- c(diff_vector, diff)
    } else {
      # Append NA if there are not enough points to calculate the difference
      diff_vector <- c(diff_vector, NA)
    }
  }
  
  # Add the difference vector as a new column to the data frame
  df$change_in_cor <- rep(diff_vector, each = 5, length.out = nrow(df))
  
  return(df)
}

cor_EN <- calculate_5_point_diff_cor(cor_EN)
cor_Lasso <- calculate_5_point_diff_cor(cor_Lasso)
cor_Ridge <- calculate_5_point_diff_cor(cor_Ridge)

#FIND POINT OF DIMINISHING RETURNS TO ADDING MORE METABOLITES
find_max_slope_point <- function(df, change_column_name) {
  # Get the change values
  change_values <- df[[change_column_name]]
  
  # Find the maximum slope and its index
  max_slope <- max(change_values, na.rm = TRUE)
  max_slope_index <- which.max(change_values)
  
  # Find the point where slope stops increasing
  # Look for the first point where the change starts decreasing or remains stable
  for (i in (max_slope_index + 1):length(change_values)) {
    if (change_values[i] < change_values[i - 1]) {
      return(i - 1)  # Return the index where the change starts decreasing
    }
  }
  
  # If no decrease is found, return the last index
  return(length(change_values))
}

max_slope_point_auc_EN <- find_max_slope_point(auc_EN, "change_in_auc")
max_slope_point_auc_Lasso <- find_max_slope_point(auc_Lasso, "change_in_auc")
max_slope_point_auc_Ridge <- find_max_slope_point(auc_Ridge, "change_in_auc")

# Apply the function to each data frame for correlation
max_slope_point_cor_EN <- find_max_slope_point(cor_EN, "change_in_cor")
max_slope_point_cor_Lasso <- find_max_slope_point(cor_Lasso, "change_in_cor")
max_slope_point_cor_Ridge <- find_max_slope_point(cor_Ridge, "change_in_cor")

# Print the results
print("AUC - EN:")
print(max_slope_point_auc_EN)
print("AUC - Lasso:")
print(max_slope_point_auc_Lasso)
print("AUC - Ridge:")
print(max_slope_point_auc_Ridge)

print("Correlation - EN:")
print(max_slope_point_cor_EN)
print("Correlation - Lasso:")
print(max_slope_point_cor_Lasso)
print("Correlation - Ridge:")
print(max_slope_point_cor_Ridge)
