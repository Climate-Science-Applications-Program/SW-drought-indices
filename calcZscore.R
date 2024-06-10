# helper functions to calculate rolling mean/sum z-scores for monthly climate variables
# MAC 06/09/24


# rolling mean z-score
climZ_mean<-function(scales,dataVar,date,month){
  tempZ<-lapply(1:length(scales), function(x) zoo::rollapply(dataVar, FUN = mean, width = scales[x],
                                                             fill=NA,align="right", by.column = TRUE))
  tempZ<-do.call(cbind, tempZ)
  colnames(tempZ)<-paste0("Z_",scales)
  tempZ<-cbind.data.frame(date,month,tempZ)
  colnames(tempZ)[1:2]<-c("dates","month")
  tempZ<-tempZ %>%
    group_by(month) %>%
    mutate(Z_3 = (Z_3-mean(Z_3, na.rm=TRUE))/sd(Z_3, na.rm=TRUE),
           Z_6 = (Z_6-mean(Z_6, na.rm=TRUE))/sd(Z_6, na.rm=TRUE),
           Z_12 = (Z_12-mean(Z_12, na.rm=TRUE))/sd(Z_12, na.rm=TRUE))
}

# rolling sum z-score
climZ_sum<-function(scales,dataVar,date,month){
  tempZ<-lapply(1:length(scales), function(x) zoo::rollapply(dataVar, FUN = sum, width = scales[x],
                                                             fill=NA,align="right", by.column = TRUE))
  tempZ<-do.call(cbind, tempZ)
  colnames(tempZ)<-paste0("Z_",scales)
  tempZ<-cbind.data.frame(date,month,tempZ)
  colnames(tempZ)[1:2]<-c("dates","month")
  tempZ<-tempZ %>%
    group_by(month) %>%
    mutate(Z_3 = (Z_3-mean(Z_3, na.rm=TRUE))/sd(Z_3, na.rm=TRUE),
           Z_6 = (Z_6-mean(Z_6, na.rm=TRUE))/sd(Z_6, na.rm=TRUE),
           Z_12 = (Z_12-mean(Z_12, na.rm=TRUE))/sd(Z_12, na.rm=TRUE))
}


##### calc bias function
# Define the function
calcBias <- function(data1, data2) {
  # Check if both data sets are numeric
  if (!is.numeric(data1)) {
    stop("The first data set must be numeric")
  }
  if (!is.numeric(data2)) {
    stop("The second data set must be numeric")
  }
  
  # Calculate the means of both data sets
  mean_data1 <- mean(data1, na.rm = TRUE)
  mean_data2 <- mean(data2, na.rm = TRUE)
  
  # Calculate the bias
  bias <- mean_data1 - mean_data2
  
  # Return the bias
  return(bias)
}

# # Use the function
# data1 <- c(5.1, 5.3, 5.0, 5.2, 5.4)
# data2 <- c(4.9, 5.0, 4.8, 4.9, 5.1)
# 
# bias_result <- calculate_bias_between_sets(data1, data2)
# print(bias_result)  # Output should reflect the difference in means between data1 and data2

#### Calc RMSE
# Define the enhanced function
calcRMSE <- function(observed, predicted) {
  # Check if both datasets are numeric
  if (!is.numeric(observed)) {
    stop("The observed dataset must be numeric")
  }
  if (!is.numeric(predicted)) {
    stop("The predicted dataset must be numeric")
  }
  
  # Check if both datasets have the same length
  if (length(observed) != length(predicted)) {
    stop("The observed and predicted datasets must have the same length")
  }
  
  # Check if both datasets are non-empty
  if (length(observed) == 0 || length(predicted) == 0) {
    stop("The observed and predicted datasets must not be empty")
  }
  
  # Calculate the squared differences
  squared_diff <- (observed - predicted)^2
  
  # Calculate the mean of the squared differences
  mean_squared_diff <- mean(squared_diff, na.rm = TRUE)
  
  # Calculate the RMSE
  rmse <- sqrt(mean_squared_diff)
  
  # Return the RMSE
  return(rmse)
}

# Use the enhanced function
# observed <- c(5.1, 5.3, 5.0, 5.2, 5.4)
# predicted <- c(5.0, 5.1, 5.0, 5.1, 5.3)
# 
# rmse_result <- calculate_rmse(observed, predicted)
# print(rmse_result)  # Output should reflect the RMSE between observed and predicted


