### Calculating statistical significance of an A/B test ###
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
  
  ### read in data ###
  #Make sure the data includes all users who were entered in the test, even if they didn't trigger the metric
  
  control <-read.table("control_file.tsv", sep="\t", comment.char = "", stringsAsFactors = FALSE)
  variant <-read.table("variant_file.tsv", sep="\t", comment.char = "", stringsAsFactors = FALSE) 
  
  control_metric <- rmoutliers(control$V2)
  variant_metric <- rmoutliers(variant$V2)
  #Removes outliers from the metric in both control and variant
  
  t.test(control_metric, variant_metric)
  