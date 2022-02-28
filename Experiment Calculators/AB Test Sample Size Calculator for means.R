### Calculating sample size needed for an A/B test ###
-----------------------------------------------------------------------------
  
  #Function to remove outliers from data
  #This identifies and removes any outliers above or below 1.5*IQR
  
  rmoutliers <- function (x) {
    
    #Find the outlier boundary in the vector
    outlier_lower <- median(x) - 1.5*IQR(x)
    outlier_upper <- median(x) + 1.5*IQR(x)
    
    #Remove values greater than the minimum outlier from the vector
    filtered <- x[x < outlier_upper]
    filtered <- filtered[filtered > outlier_lower]
    
    #Percentage of values remaining
    valsremaining <- length(filtered)/length(x)
    
    #Give an error message if too much data is removed. This could happen if the data had a particularly long tail.  
    if (valsremaining < 0.95){
      stop ("This function will remove more than 5% percent of your data. You need to remove outliers manually.")}
    
    else if (length(filtered)/length(x) < 0.99){
      warning("This calculation has removed between 1% and 5% of your data.") 
      filtered
    }
    
    else
    {filtered}
  }
  
  
  minimumDetectableChange <- 20
  #This is otherwise known as the practical significance level
  #What difference in the metric would you need to see in order for the change to be worthwhile for your product
  #This is an absolute change
  
  numberOfVariants <- 2
  #The number of distinct experiences including the control - e.g. for an a/b test the number would be 2.
  
  ### read in data ###
  #The data should contain 1 weeks worth of data with a row for each user and a count of the number of metrics
  #Make sure the data includes users whose count is 0
  
  data <-read.table("data_file.tsv", sep="\t", comment.char = "", stringsAsFactors = FALSE) 
  
  metric <- rmoutliers(data$V2)
  #Removes outliers from the metric
  
  standardDeviation <- sd(metric)  
  weeklyBrowsers <- length(metric)
  
  # The coefficient for significance level 0.05 and statistical power 0.8
  coefficient <- 7.9
  
  ## Number of users needed in the test
  users <- coefficient * 2 * numberOfVariants * (standardDeviation^2/minimumDetectableChange^2)
  
  ## Number of weeks for the test to run
  #This calculation assumes that you will be able to access the same number of new users each week
  #It will not be accurate for products with a large number of regular users for a test on 100% of users
  timeInWeeks <- users/weeklyBrowsers
  