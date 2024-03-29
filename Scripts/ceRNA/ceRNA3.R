generateDataAndCalculatePearson <- function(data_1, data_2, n = 100) {
  # n is the number of data points, default is 100
  
  # Calculate Pearson correlation
  correlation_result <- cor(data_1, data_2, method = "pearson")
  
  # Return the result
  return(correlation_result)
}
generateDataAndCalculatePearson(150)
