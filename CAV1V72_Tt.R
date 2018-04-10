


###This is the code I used to calculate the intracellular transition temperature of dual transfected
#cells (CAV1-V72 + GFP-V60) and single transfected cells (GFP-V60)
setwd("/Users/davidtyrpak/Dropbox/Lab/Data/Yue")
df.yw <- read.csv(file.choose()) #/Users/davidtyrpak/Dropbox/Lab/Data/Yue/Combine CAV1-V72 ROI Analysis.csv
View(df.yw)

##clean up data to identify number of cells from each treatment, date, and well
levels(df.yw$treatment)
df.yw$cell.number <- factor(df.yw$cell.number)
df.yw$cell.number.file <- interaction(df.yw$cell.number, df.yw$filepath)
df.yw$cell.number.file <- factor(df.yw$cell.number.file)

##identify single transfected cells
single <- df.yw[which(df.yw$treatment == "single"), ]
single$treatment <- factor(single$treatment) ##get rid of unused "dual" treatment level
single$cell.number.file <- factor(single$cell.number.file)

##identify dual transfected cells
dual <- df.yw[which(df.yw$treatment == "dual"), ]
dual$treatment <- factor(dual$treatment) ##get rid of unused "single" Treatment level
dual$cell.number.file <- factor(dual$cell.number.file); dual$cell.number.file

##Sanity check
length(df.yw$cell.number.file) == length(single$cell.number.file) + length(dual$cell.number.file)

length(unique(single$cell.number.file)) ##31 single transfected cells
length(unique(dual$cell.number.file)) ## 35 dual transfected cells


##SINGLE TRANSFECTED (GFP-V60) Tt analysis
#let's get the data row indexes for each single cell from the single data frame
single_cells_st <- sapply(unique(single$cell.number.file), function(x) which(single$cell.number.file == x)); single_cells_st
names(single_cells_st) <- unique(single$cell.number.file)

stbatch <- batch_plot(single_cells_st, single, baseline_end = 10, threshold_limit = 4, temp.column.index = 4, standard.dev.column.index = 7)


quantifiable.stbatch <- na.omit(stbatch); length(quantifiable.stbatch) ##20 cells
min(na.omit(quantifiable.stbatch)) ##26.2
median(quantifiable.stbatch) ##40.1 Celsius


#DUAL Transfected (GFP-V60 + CAV1-V72) Tt analysis
dual_cells_dt <- sapply(unique(dual$cell.number.file), function(x) which(dual$cell.number.file == x)); dual_cells_dt
names(dual_cells_dt) <- unique(dual$cell.number.file)

dtbatch <- batch_plot(dual_cells_dt, dual, baseline_end = 10, threshold_limit = 4, temp.column.index = 4, standard.dev.column.index = 7)

quantifiable.dtbatch <- na.omit(dtbatch); length(quantifiable.dtbatch) ##35 cells
min(na.omit(quantifiable.dtbatch)) ##25.4
median(quantifiable.dtbatch) ##36.9 Celsius

##statistics
wilcox.test(quantifiable.dtbatch, quantifiable.stbatch)


###This is the function that calculates transition temperatures for single cells.

batch_plot <- function(list_of_lists, data_frame, threshold_limit = 5, baseline_end = 20, temp.column.index, standard.dev.column.index){
  
  
  
  ##list_of_lists: list class; this is a list containing lists of row indexes for individual wells or cells
  
  ##data_frame: data frame class; this is the data frame that the indexes in list_of_lists refer to
  
  ##threshold_limit: int class; how many standard deviations away from the mean you want your threshold to be. 
  
  ##baseline_end: int class; the number of data points you want to go into your calculation of the baseline (e.g. baseline_end = 20 means the first 20 points of the data set are used in the calculation of the mean and standard deviation)
  
  indexes <- lapply(list_of_lists, function(x) {c(start  = x[1], stop = tail(x, n = 1))}) ##retrieve first and last index for that well/cell
  
  i = 1
  Tts <- vector(mode = "numeric", length = length(list_of_lists))
  
  ##for loop to plot a scatter plot for each well/cell in list_of_lists
  for (index in indexes) {
    start = index[1]
    stop = index[2]
    
    temps <- data_frame[ ,temp.column.index][start:stop] ##Celsius data points for this well
    
    sds <- data_frame[,standard.dev.column.index][start:stop] ##corresponding standard deviation points for this well
    
    mean <- mean(sds[1:baseline_end]) ##mean of the datapoints specified by baseline_end
    
    std_dv <- sd(sds[1:baseline_end]) ##calculate the standard deviation around these points
    
    std_normals <- (sds - mean)/std_dv ## convert the stand deviation pixal values to a standard normal distribution
    
    #threshold <- threshold_limit*std_dv ##threshold is threshold_limit times the standard deviation from the baseline/mean
    
    Tt <- which(std_normals > threshold_limit)[1]  ## first data point which crosses the threshold
    
    Tts[i] <- temps[Tt]
    
    
    
    plot(temps, std_normals, xlab = "Celsius", ylab = "Z-Score", main = names(indexes)[i], ylim = c(-10, 10))
    abline(h = 0) ##mean of first 20 points. Will be zero because we converted to standard normal
    abline(v = Tts[i], col = "red", pch = 19)
    abline(h = 0 + threshold_limit, lty = 4, col = "red") ##
    # abline(h = 0 - threshold_limit, lty = 4, col = "red")
    i = i +1
    
    
  }
  names(Tts) <- names(indexes)
  Tts
}

