###########################################################################
############################ Power Functions ##############################
###########################################################################

### Please explore these functions using example data sets using the shiny app:
### https://joelarkman.shinyapps.io/PowerTools/



# Simulate log-normal data set from input
simulateLogNormal <- function(data, covType, nsamples) {
  
  message('Simulating Log Normal dataset', appendLF = F)
  
  # Create an offset to allow lof to be taken even for values close to zero.
  offset <- abs(min(data)) + 1
  
  # Add offset to data.
  offData <- data + offset
  
  # Calculate log of data.
  logData <- log(offData)
  
  # Calculate the mean of each column in the data. 
  meansLog <- colMeans(logData)
  
  if (covType=='Estimate') {
    
    # In estimate mode, compute covariance matrix for input data.
    covLog <- cov(logData)
    
  } else if (covType=='Diagonal') {
    
    # In diagonal mode produce only diagonal 1 matrix. 
    covLog <- diag(1, ncol=ncol(logData), nrow = ncol(logData))
    
  } else {
    print('ERROR')
  }
  
  # Simulate n samples from multivariate normal distribution using colmeans and covariance matrix. 
  simData <- MASS::mvrnorm(n=nsamples, mu=meansLog, Sigma = covLog)
  
  # Find exponant to return to original scale.
  simData <- exp(simData)
  
  # Remove the initial offset.
  simData <- simData - offset
  
  # Set any simulated negative numbers to 0.
  simData[simData<0] <- 0
  
  message(' ...complete!')
  
  return(simData)
}

# Calculate confusion matrix.
calcConfMatrixUniv <- function(p,corrVector,signThreshold,corrThresh){
  
  # Examine vector of p-values to identify those below threshold p-value.
  Pf <- p < signThreshold
  
  # Create storage variables for the number of..
  # True positive variables.
  TP <- 0
  # True negative variables.
  TN <- 0
  # False positive vairables.
  FP <- 0
  # False negative variables.
  FN <- 0
  
  # For each of the variables in the data set (which is equal to number of correlation 
  # values in input vector)
  for ( i in 1:length(corrVector)) {
    
    # If the currently iterated variable correlates poorly with the test variable 
    # and is not found to significantly correlate with Y..
    if (abs(corrVector[i]) < corrThresh & Pf[i]==F) {
      TN <- TN+1 # Add one to TN count.
      
      # If the currently iterated variable correlates strongly with the test variable 
      # and is also found to significantly correlate with Y.. 
    } else if (abs(corrVector[i]) > corrThresh & Pf[i]==T) {
      TP <- TP+1 # Add one to TP count.
      
      # If the currently iterated variable correlates strongly with the test variable 
      # but is not found to significantly correlate with Y..  
    } else if (abs(corrVector[i]) > corrThresh & Pf[i]==F) {
      FN <- FN+1 # Add one to FN count. 
      
      # If the currently iterated variable correlates poorly with the test variable 
      # but is found to significantly correlate with Y..  
    } else if (abs(corrVector[i]) < corrThresh & Pf[i]==T) {
      FP <- FP+1 # Add one to FP count.
    }
  }
  
  # Calculate rates of each output as a proportion of the others. 
  TNtot <- TN/(FP+TN)
  TPtot <- TP/(TP+FN)
  FPtot <- FP/(FP+TN)
  FNtot <- FN/(TP+FN)
  FDtot <- FP/(TP+FP)
  
  # Create list of rate values. 
  output <- list('TNtot'=TNtot,'TPtot'=TPtot,'FPtot'=FPtot,'FNtot'=FNtot,'FDtot'=FDtot)
  
  # Return list of confusion matrix metrics.
  return(output)
  
}

# Continous outcome power function.
PCalc_Continuous <- function(x,
                             y=NULL,
                             sampSizes=NULL,
                             grouping=T,
                             signThreshold=0.05,
                             nSimSamp=5000,
                             nRepeats=10,
                             targetVar=NULL){
  
  if (is.null(y)){stop('Y values not provided.')} # Stop if no y data provided.
  else if (is.null(sampSizes)){stop('No sample sizes provided.')} # Stop if sample sizes are not provided.
  else if (is.null(x)){stop('No x data provided.')} # Stop if no x data is provided.
  
  # Ensure no. simulated samples is larger than largest chosen sample size to evaluate.
  if (max(sampSizes) >= nSimSamp){
    print(paste('Increasing no. simulated samples to faciltiate assessment 
                of sample size: ', max(sampSizes),'',sep = ''))
    nSimSamp <- max(sampSizes) + 500
  }
  
  # Convert x data to data frame (to ensure the presense of column names).
  data <- as.data.frame(x)
  
  # Calculate the number of x variables.
  numVars <- ncol(data)
  
  # Calculate the number of sample sizes to assess each variable at.
  nSampSizes <- length(sampSizes)
  
  # Define the number of effect sizes to assess each variable at to be 1 - their actual effect size.
  nEffSizes <- 1
  
  # Estimate effect size of each variable as the absolute value of their correlation with y.
  message('Calculating true variable effect sizes', appendLF = F)
  effectSizes <- abs(round(cor(x=x, y=y), digits = 2))
  message(' ...complete!')
  
  if (grouping==T) { # If grouping is active.
    
    message('Grouping similar variables & extracting group member with largest effect size',
            appendLF = F)
    
    # Group similar variables based on hierachicial clustering of distance matrix generated from 
    # inter-variable correlations.
    corr_clustering  <- hclust(as.dist(1-abs(cor(na.omit(data)))))
    
    # Cut the tree (and thus define each group) based on height of tree (cut at height = 0.5)
    groups <- as.factor(cutree(corr_clustering, h = 0.5))
    
    # Define variables for storage.
    out <- NULL
    groupSize <- NULL
    
    # For each of the groups identified.
    for (i in 1: length(levels(groups))) {
      
      # Extract the names of the variables in the currently iterated group.
      tmp_group <- names(groups)[which(groups==levels(groups)[i])]
      
      # Extract the effect size of the variables in the currently iterated group.
      features <- effectSizes[which(rownames(effectSizes) %in% tmp_group)]
      
      # Add variable names to the new vector of effect sizes.
      names(features) <- tmp_group
      
      # Calculate the size of the current group and append this to storage vector.
      groupSize <- c(groupSize,length(features))
      
      # Order the group members accoding to their effect size and extract group 
      # member with largest effect size.
      features <- features[order(features, decreasing = T)][1]
      
      # Add the chosen variable name / effect size for currently iterated group to 
      # output vector of features to assess.
      out <- c(out,features)
    }
    
    # Update effectSizes vector to contain only the variables to be assessed.
    effectSizes <- out
    
    # Define assessVars as the column number of each of the variables to be assessed.
    assessVars <-  which(colnames(data) %in% names(out))
    
    # Define namesVars as simply the names of each of the variables to be assessed.
    namesVars <- names(out)
    
    message(' ...complete!')
    
  } else if (grouping == F) { # If grouping is inactive.
    
    # Define all variables as assessVars and store their names.
    assessVars <- 1:ncol(data)
    namesVars <- colnames(data)
    groupSize <- NA
    
  }
  
  # Use 'simulateLogNormal' function to produce a simulated dataset of samples 
  # equal to the provided 'nSimSamp' variable and based on the correlation structure 
  # of the provided xdata.
  Samples <- simulateLogNormal(data, 'Estimate', nSimSamp)
  
  # Produce correlation matrix for the simulated dataset.
  message('Calculating correlation structure', appendLF = F)
  correlationMat <- cor(Samples)
  message(' ...complete!')
  message('')
  
  if (!is.null(targetVar)) { # If a targetVar is provided.
    
    # If the provided targetVar is the name of one of the x variables.
    if (targetVar %in% colnames(data)) {
      
      # Search for the targetVar within the groups vector and return the number of the group it is within. 
      # Then, update 'effectSizes' to only contain that of the chosen variable from this group.
      effectSizes <- effectSizes[groups[names(groups) %in% targetVar]]
      
      # Similarly, update 'namesVars'  and 'assessVars' to only contain the name and column number of the 
      # chosen variable grouped with targetVar.
      namesVars <- namesVars[groups[names(groups) %in% targetVar]]
      assessVars <- assessVars[groups[names(groups) %in% targetVar]]
      
      message(paste("Targeting variable: '", targetVar,
                    "', which is grouped with '", namesVars,"'", sep = ''))
      message('')
      
      # If the provided targetVar refers to the column number of one of the x variables.  
    } else if (targetVar %in% 1:ncol(data)) {
      
      # As above, identify the group targetVar belongs to, the variable chosen
      # to represent that group, and update 'effectSizes', 'namesVars' and 'assessVars'
      # to contain this variable only.
      effectSizes <- effectSizes[groups[names(groups) %in% colnames(data)[targetVar]]]
      namesVars <- namesVars[groups[names(groups) %in% colnames(data)[targetVar]]]
      assessVars <- assessVars[groups[names(groups) %in% colnames(data)[targetVar]]]
      
      message(paste("Targeting variable: '", targetVar,
                    "', which is grouped with '", namesVars,"'", sep = ''))
      message('')
      
    } else{stop('Invalid target variable')} # If targetVar is anything else, stop.
  }
  
  # Initialise output data structure. 
  output <- NULL
  
  
  for (i in 1:length(assessVars)) { # For each of the variables to be assessed.
    
    # Store the column number of currently iterated variable as 'currVar'.
    currVar <- assessVars[i]
    
    message(paste('(',i, '/', length(assessVars),") Assessing variable '", namesVars[i],
                  "' (effect size: ", sep = ''), appendLF = F)
    
    # Create storage data structure comprised of multiple matrices with nrows equal
    # to the number of effectSizes to assess each variable at (typically 1) and ncols 
    # equal to the number of sample sizes to assess each variable at.
    
    # Each of these matrices is defined to store calculated 'true positive', 'false positive', 
    # 'true negative', etc, values for the currently iterated variable when assessed at
    # each sample size effect sizes.
    
    # In addition, create three of these large data structures: 
    
    # One to store values calculated using uncorrected p-values. 
    uncStruct <- list('TP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'FP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'TN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'FN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'FD'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'STP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'SFP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'STN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'SFN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'SFD'=matrix(0,nrow=nEffSizes,ncol=nSampSizes))
    
    # One to store values calculated using bonferroni adjusted p-values. 
    bonfStruct <- list('TP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'FP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'TN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'FN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'FD'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'STP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'SFP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'STN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'SFN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'SFD'=matrix(0,nrow=nEffSizes,ncol=nSampSizes))
    
    # One to store values calculated using benjamini-hochberg adjusted p-values. 
    bhStruct <- list('TP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'FP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'TN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'FN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'FD'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'STP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'SFP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'STN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'SFN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'SFD'=matrix(0,nrow=nEffSizes,ncol=nSampSizes))
    
    
    for (currEff in 1:nEffSizes) {  # For each of the effect sizes each variable is being assessed at.
      
      # Create storage vector equal to the number of features in data set.
      b1 <- rep(0, numVars)
      
      # Store the effect size of the curently iterated variable at the relevent position 
      # within this storage vector (e.g if variable 4 was being assessed, store its effectSize
      # in the fourth position of b1).
      b1[currVar] <- effectSizes[i]
      
      message(paste(b1[currVar]), '):')
      message('Using sample size: ', appendLF = F)
      
      for (currSampSize in 1:nSampSizes) { # For each sample size each variable is being assessed at.
        
        if(currSampSize==nSampSizes) {
          message(paste(sampSizes[currSampSize],'.', sep = ''), appendLF = F)
          message()
        } else {
          message(paste(sampSizes[currSampSize],', ', sep = ''), appendLF = F)
        }
        
        # Create tempory storage structures to store the uncorrected and corrected confusion matrix 
        # values calculated for each of the nRepeats to be conducted.
        Results <- data.frame('TP'=rep(0,nRepeats),
                              'FP'=rep(0,nRepeats),
                              'TN'=rep(0,nRepeats),
                              'FN'=rep(0,nRepeats),
                              'FD'=rep(0,nRepeats))
        
        Bonferroni <- data.frame('TP'=rep(0,nRepeats),
                                 'FP'=rep(0,nRepeats),
                                 'TN'=rep(0,nRepeats),
                                 'FN'=rep(0,nRepeats),
                                 'FD'=rep(0,nRepeats))
        
        BH <- data.frame('TP'=rep(0,nRepeats),
                         'FP'=rep(0,nRepeats),
                         'TN'=rep(0,nRepeats),
                         'FN'=rep(0,nRepeats),
                         'FD'=rep(0,nRepeats))
        
        # Combine the three lists of storage vectors (one for uncorrected, one for 
        # bonferroni corrected and one for BH corrected p-values) to create large temporary 
        # storage structure.
        multiplerepeats <- list('Results'=Results,'Bonferroni'=Bonferroni,'BH'=BH)
        
        
        # Create a progress bar to show position with respect to number of repeats.
        #pb <- txtProgressBar(min = 0, max = nRepeats, style = 3) 
        
        for (currRepeat in 1:nRepeats) { # For each of the nRepeats
          
          # Generate selection index of random integers of length equal to the currently iterated 
          # sample size (values between 1 and the number of simulated samples).
          selectIndex <- sample.int(nSimSamp, sampSizes[currSampSize])
          
          # Use the selection index to extract a random selection of samples from the simulated 
          # dataset, equal to the currently iterated sample size. 
          SelSamples <- Samples[selectIndex,]
          
          # Center and scale the new data subset.
          SelSamples <- scale(SelSamples)
          
          # Generate random noise values (mean=0,SD=1) equal to the currently iterated sample size.
          noiseLevel <- 1
          noise <- noiseLevel * rnorm(sampSizes[currSampSize])
          
          # Generate a vector of y values by multiplying the x-values of the currently iterated variable 
          # by its effect size and adding the random noise generated above. 
          Y <- SelSamples[,currVar]*b1[currVar]
          Y <- Y + noise
          
          # Create p-value storage vector. 
          p <- rep(0, numVars)
          
          # Regress each variable in the dataset against Y value vector, extract and store p-values.
          for (i in 1:numVars) {
            B <- data.frame(Y=Y,X=SelSamples[,i])
            
            if (anyNA(B$X)) {
              B <- data.frame(Y=Y,X=SelSamples[,i])
              if (anyNA(B$X)) {
                p[i] <- 1
                next()
              }
            }
            
            model <- lm(Y~.,B)
            p[i] <- summary(model)[[4]][2,4]
          }
          
          pUnc <- p # Store uncorrected p-values.
          pBonf <- p.adjust(p, method = 'bonferroni') # Apply and store bonferroni correction.
          pBH <- p.adjust(p, method = 'BH') # Apply and store BH correction.
          
          # Extract correlation vector pertaining to how each variable in the datset correlates 
          # with the currently iterated variable.
          corrVector <- correlationMat[,currVar]
          
          # Perform 'calcConfMatrixUniv' function using uncorrected, bonferoni corrected and 
          # BH corrected p-values. 
          uncorrected <- calcConfMatrixUniv(pUnc,corrVector,signThreshold,0.8)
          bonferroni <- calcConfMatrixUniv(pBonf,corrVector,signThreshold,0.8)
          bh <- calcConfMatrixUniv(pBH,corrVector,signThreshold,0.8)
          
          # Store the uncorrected confusion matrix output values for the current repeat.
          multiplerepeats$Results$TP[currRepeat] <- uncorrected$TPtot
          multiplerepeats$Results$FP[currRepeat] <- uncorrected$FPtot
          multiplerepeats$Results$TN[currRepeat] <- uncorrected$TNtot
          multiplerepeats$Results$FN[currRepeat] <- uncorrected$FNtot
          multiplerepeats$Results$FD[currRepeat] <- uncorrected$FDtot
          
          # Store the bonferroni corrected confusion matrix output values for the current repeat.
          multiplerepeats$Bonferroni$TP[currRepeat] <- bonferroni$TPtot
          multiplerepeats$Bonferroni$FP[currRepeat] <- bonferroni$FPtot
          multiplerepeats$Bonferroni$TN[currRepeat] <- bonferroni$TNtot
          multiplerepeats$Bonferroni$FN[currRepeat] <- bonferroni$FNtot
          multiplerepeats$Bonferroni$FD[currRepeat] <- bonferroni$FDtot
          
          # Store the BH corrected confusion matrix output values for the current repeat.
          multiplerepeats$BH$TP[currRepeat] <- bh$TPtot
          multiplerepeats$BH$FP[currRepeat] <- bh$FPtot
          multiplerepeats$BH$TN[currRepeat] <- bh$TNtot
          multiplerepeats$BH$FN[currRepeat] <- bh$FNtot
          multiplerepeats$BH$FD[currRepeat] <- bh$FDtot
          
          #setTxtProgressBar(pb, currRepeat) # Update progress bar.
          
        }
        # close(pb) # Terminate progress bar. 
        
        # Extract the names of each of the confusion matrix metrics (e.g. TP, FP, TN, FN and FD).
        stats <- colnames(multiplerepeats$Bonferroni)
        
        for (nstat in 1:5) { # For each of the metrics 
          
          # Define currstat as the name of the currently iterated metric.
          currstat <- stats[nstat]
          
          # Store the average uncorrected confusion matrix values (across repeats) generated 
          # using the currently iterated variable, effect size and sample size. 
          uncStruct[[which(names(uncStruct)%in%stats[nstat])]][currEff,currSampSize] <- 
            mean(multiplerepeats$Results[,which(names(multiplerepeats$Results)%in%stats[nstat])])
          
          # Store the standard deviation of the uncorrected confusion matrix values (across repeats) 
          # generated using the currently iterted variable, effect size and sample size. 
          uncStruct[[which(names(uncStruct)%in%paste('S',stats[nstat],sep = ''))]][currEff,currSampSize] <-
            sd(multiplerepeats$Results[,which(names(multiplerepeats$Results)%in%stats[nstat])])
          
          # Store the average bonferroni corrected confusion matrix values (across repeats) generated 
          # using the currently iterated variable, effect size and sample size. 
          bonfStruct[[which(names(bonfStruct)%in%stats[nstat])]][currEff,currSampSize]  <-
            mean(multiplerepeats$Bonferroni[,which(names(multiplerepeats$Bonferroni)%in%stats[nstat])])
          
          # Store the standard deviation of the bonferroni corrected confusion matrix values (across repeats) 
          # generated using the currently iterted variable, effect size and sample size. 
          bonfStruct[[which(names(bonfStruct)%in%paste('S',stats[nstat],sep = ''))]][currEff,currSampSize] <-
            sd(multiplerepeats$Bonferroni[,which(names(multiplerepeats$Bonferroni)%in%stats[nstat])])
          
          # Store the average BH corrected confusion matrix values (across repeats) generated 
          # using the currently iterated variable, effect size and sample size. 
          bhStruct[[which(names(bhStruct)%in%stats[nstat])]][currEff,currSampSize] <- 
            mean(multiplerepeats$BH[,which(names(multiplerepeats$BH)%in%stats[nstat])])
          
          # Store the standard deviation of the BH corrected confusion matrix values (across repeats) 
          # generated using the currently iterted variable, effect size and sample size. 
          bhStruct[[which(names(bhStruct)%in%paste('S',stats[nstat],sep = ''))]][currEff,currSampSize]  <-
            sd(multiplerepeats$BH[,which(names(multiplerepeats$BH)%in%stats[nstat])])
          
        }
        
      }
      
    }
    
    # Combine each of the storage data structures for the currently iterated variable 
    # into single storeVar 'list'.
    storeVar <- list('Uncorrected'=uncStruct,'Bonferroni'=bonfStruct,'BH'=bhStruct)
    
    # Append the storeVar structure for the currently iterated variable to the main output data structure.
    output <- c(output, list(storeVar))
    
    message('')
  }
  
  # Define the names of each list entry of the main output data structure as the names of each of the
  # assessed variables. 
  names(output) <- namesVars
  
  # Create a list of the input parameters provided to the function to faciliate subsequent plotting.
  input_params <- list('effectSizes'=effectSizes,
                       'sampSizes'=sampSizes,
                       'assessVars'=assessVars,
                       'namesVars'=namesVars,
                       'groupSize'=groupSize,
                       'targetVar'=!is.null(targetVar),
                       'numSamps'=nrow(data))
  
  
  # Append the list of input parameters to the main output data structure.
  output <- c(output, 'input_params'= list(input_params))
  
  # Add S3 class of 'PCalc' to main output data structure.
  class(output) <- c('PCalc','list')
  
  # Return main output data structure.
  return(output)
}

# Two-group classification outcome power function.
PCalc_2Group <- function(x,
                         y=NULL,
                         n_controls=NULL,
                         control_group =NULL,
                         sampSizes=NULL,
                         grouping=T,
                         signThreshold=0.05,
                         nSimSamp=500,
                         nRepeats=10,
                         targetVar=NULL){
  
  if (is.null(y)){stop('Y values not provided.')} # Stop if no y data provided.
  else if (is.null(sampSizes)){stop('No sample sizes provided.')} # Stop if sample sizes are not provided.
  else if (is.null(x)){stop('No x data provided.')} # Stop if no x data is provided.
  
  # Ensure no. simulated samples is larger than largest chosen sample size to evaluate.
  if (2*(max(sampSizes)) >= nSimSamp){
    print(paste('Increasing no. simulated samples to faciltiate assessment 
                of sample size: ', max(sampSizes),sep = ''))
    nSimSamp <- 2*(max(sampSizes)) + 500
  }
  
  y <- as.factor(y)
  
  if (length(levels(y))!=2){stop('More than two Y groups were provided')}
  
  if (!is.null(control_group)) {
    case_group <- as.character(y[y != control_group][1])
    y <- factor(y, levels = c(case_group, control_group))
  }
  
  x <- as.data.frame(x)
  
  data <- as.data.frame(rbind(x[y %in% levels(y)[1],],
                              x[y %in% levels(y)[2],]))
  
  y_ordered <- c(y[y %in% levels(y)[1]],
                 y[y %in% levels(y)[2]])
  
  # Calculate the number of x variables.
  numVars <- ncol(data)
  
  # Calculate the number of sample sizes to assess each variable at.
  nSampSizes <- length(sampSizes)
  
  # Define the number of effect sizes to assess each variable at to be 1 - their actual effect size.
  nEffSizes <- 1
  
  # Estimate effect size of each variable as the absolute value of their correlation with y.
  message('Calculating true variable effect sizes', appendLF = F)
  effectSizes <- NULL
  for (i in 1:numVars) {
    n1 <- length(which(y_ordered==1))
    n2 <- length(which(y_ordered==2))
    m1 <- mean(data[which(y_ordered==1),i])
    m2 <- mean(data[which(y_ordered==2),i])
    S1 <- var(data[which(y_ordered==1),i])
    S2 <- var(data[which(y_ordered==2),i])
    s <- sqrt((((n1-1)*(S1)) + ((n2-1)*(S2)))/(n1+n2-2))
    d <- abs(round((m1-m2)/s,2))
    effectSizes <- c(effectSizes,d)
  }
  names(effectSizes) <- colnames(data)
  message(' ...complete!')
  
  if (grouping==T) { # If grouping is active.
    
    message('Grouping similar variables & extracting group member with largest effect size',
            appendLF = F)
    
    # Group similar variables based on hierachicial clustering of distance matrix generated from 
    # inter-variable correlations.
    corr_clustering  <- hclust(as.dist(1-abs(cor(na.omit(data)))))
    
    # Cut the tree (and thus define each group) based on height of tree (cut at height = 0.5)
    groups <- as.factor(cutree(corr_clustering, h = 0.5))
    
    # Define variables for storage.
    out <- NULL
    groupSize <- NULL
    
    # For each of the groups identified.
    for (i in 1: length(levels(groups))) {
      
      # Extract the names of the variables in the currently iterated group.
      tmp_group <- names(groups)[which(groups==levels(groups)[i])]
      
      # Extract the effect size of the variables in the currently iterated group.
      features <- effectSizes[which(names(effectSizes) %in% tmp_group)]
      
      # Add variable names to the new vector of effect sizes.
      names(features) <- tmp_group
      
      # Calculate the size of the current group and append this to storage vector.
      groupSize <- c(groupSize,length(features))
      
      # Order the group members accoding to their effect size and extract group 
      # member with largest effect size.
      features <- features[order(features, decreasing = T)][1]
      
      # Add the chosen variable name / effect size for currently iterated group to 
      # output vector of features to assess.
      out <- c(out,features)
    }
    
    # Update effectSizes vector to contain only the variables to be assessed.
    effectSizes <- out
    
    # Define assessVars as the column number of each of the variables to be assessed.
    assessVars <-  which(colnames(data) %in% names(out))
    
    # Define namesVars as simply the names of each of the variables to be assessed.
    namesVars <- names(out)
    
    message(' ...complete!')
    
  } else if (grouping == F) { # If grouping is inactive.
    
    # Define all variables as assessVars and store their names.
    assessVars <- 1:ncol(data)
    namesVars <- colnames(data)
    groupSize <- NA
    
  }
  
  # Use 'simulateLogNormal' function to produce a simulated dataset of samples 
  # equal to the provided 'nSimSamp' variable and based on the correlation structure 
  # of the provided xdata.
  Samples <- simulateLogNormal(data, 'Estimate', nSimSamp)
  
  # Produce correlation matrix for the simulated dataset.
  message('Calculating correlation structure', appendLF = F)
  correlationMat <- cor(Samples)
  message(' ...complete!')
  message('')
  
  if (!is.null(targetVar)) { # If a targetVar is provided.
    
    # If the provided targetVar is the name of one of the x variables.
    if (targetVar %in% colnames(data)) {
      
      # Search for the targetVar within the groups vector and return the number of the group it is within. 
      # Then, update 'effectSizes' to only contain that of the chosen variable from this group.
      effectSizes <- effectSizes[groups[names(groups) %in% targetVar]]
      
      # Similarly, update 'namesVars'  and 'assessVars' to only contain the name and column number of the 
      # chosen variable grouped with targetVar.
      namesVars <- namesVars[groups[names(groups) %in% targetVar]]
      assessVars <- assessVars[groups[names(groups) %in% targetVar]]
      
      message(paste("Targeting variable: '", targetVar,
                    "', which is grouped with '", namesVars,"'", sep = ''))
      message('')
      
      # If the provided targetVar refers to the column number of one of the x variables.  
    } else if (targetVar %in% 1:ncol(data)) {
      
      # As above, identify the group targetVar belongs to, the variable chosen
      # to represent that group, and update 'effectSizes', 'namesVars' and 'assessVars'
      # to contain this variable only.
      effectSizes <- effectSizes[groups[names(groups) %in% colnames(data)[targetVar]]]
      namesVars <- namesVars[groups[names(groups) %in% colnames(data)[targetVar]]]
      assessVars <- assessVars[groups[names(groups) %in% colnames(data)[targetVar]]]
      
      message(paste("Targeting variable: '", targetVar,
                    "', which is grouped with '", namesVars,"'", sep = ''))
      message('')
      
    } else{stop('Invalid target variable')} # If targetVar is anything else, stop.
  }
  
  # Initialise output data structure. 
  output <- NULL
  
  for (i in 1:length(assessVars)) { # For each of the variables to be assessed.
    
    # Store the column number of currently iterated variable as 'currVar'.
    currVar <- assessVars[i]
    
    message(paste('(',i, '/', length(assessVars),") Assessing variable '", namesVars[i],
                  "' (effect size: ", sep = ''), appendLF = F)
    
    # Create storage data structure comprised of multiple matrices with nrows equal
    # to the number of effectSizes to assess each variable at (typically 1) and ncols 
    # equal to the number of sample sizes to assess each variable at.
    
    # Each of these matrices is defined to store calculated 'true positive', 'false positive', 
    # 'true negative', etc, values for the currently iterated variable when assessed at
    # each sample size effect sizes.
    
    # In addition, create three of these large data structures: 
    
    # One to store values calculated using uncorrected p-values. 
    uncStruct <- list('TP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'FP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'TN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'FN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'FD'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'STP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'SFP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'STN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'SFN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                      'SFD'=matrix(0,nrow=nEffSizes,ncol=nSampSizes))
    
    # One to store values calculated using bonferroni adjusted p-values. 
    bonfStruct <- list('TP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'FP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'TN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'FN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'FD'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'STP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'SFP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'STN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'SFN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                       'SFD'=matrix(0,nrow=nEffSizes,ncol=nSampSizes))
    
    # One to store values calculated using benjamini-hochberg adjusted p-values. 
    bhStruct <- list('TP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'FP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'TN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'FN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'FD'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'STP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'SFP'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'STN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'SFN'=matrix(0,nrow=nEffSizes,ncol=nSampSizes),
                     'SFD'=matrix(0,nrow=nEffSizes,ncol=nSampSizes))
    
    
    for (currEff in 1:nEffSizes) {  # For each of the effect sizes each variable is being assessed at.
      
      # Create storage vector equal to the number of features in data set.
      b1 <- rep(0, numVars)
      
      # Store the effect size of the curently iterated variable at the relevent position 
      # within this storage vector (e.g if variable 4 was being assessed, store its effectSize
      # in the fourth position of b1).
      b1[currVar] <- effectSizes[i]
      
      message(paste(b1[currVar]), '):')
      message('Using sample size: ', appendLF = F)
      
      for (currSampSize in 1:nSampSizes) { # For each sample size each variable is being assessed at.
        
        if(currSampSize==nSampSizes) {
          message(paste(sampSizes[currSampSize],'.', sep = ''), appendLF = F)
          message()
        } else {
          message(paste(sampSizes[currSampSize],', ', sep = ''), appendLF = F)
        }
        
        # Create tempory storage structures to store the uncorrected and corrected confusion matrix 
        # values calculated for each of the nRepeats to be conducted.
        Results <- data.frame('TP'=rep(0,nRepeats),
                              'FP'=rep(0,nRepeats),
                              'TN'=rep(0,nRepeats),
                              'FN'=rep(0,nRepeats),
                              'FD'=rep(0,nRepeats))
        
        Bonferroni <- data.frame('TP'=rep(0,nRepeats),
                                 'FP'=rep(0,nRepeats),
                                 'TN'=rep(0,nRepeats),
                                 'FN'=rep(0,nRepeats),
                                 'FD'=rep(0,nRepeats))
        
        BH <- data.frame('TP'=rep(0,nRepeats),
                         'FP'=rep(0,nRepeats),
                         'TN'=rep(0,nRepeats),
                         'FN'=rep(0,nRepeats),
                         'FD'=rep(0,nRepeats))
        
        # Combine the three lists of storage vectors (one for uncorrected, one for 
        # bonferroni corrected and one for BH corrected p-values) to create large temporary 
        # storage structure.
        multiplerepeats <- list('Results'=Results,'Bonferroni'=Bonferroni,'BH'=BH)
        
        
        # Create a progress bar to show position with respect to number of repeats.
        #pb <- txtProgressBar(min = 0, max = nRepeats, style = 3) 
        
        for (currRepeat in 1:nRepeats) { # For each of the nRepeats
          
          if (is.null(n_controls)) {
            m <- 2*(sampSizes[currSampSize])
          } else {
            m <- sampSizes[currSampSize] + n_controls
          }
          
          # Generate selection index of random integers of length equal to the currently iterated 
          # sample size (values between 1 and the number of simulated samples).
          #selectIndex <- sample.int(nrow(Samples), (1+controlgroupWeight)*sampSizes[currSampSize])
          selectIndex <- sample.int(nrow(Samples), m)
          
          # Use the selection index to extract a random selection of samples from the simulated 
          # dataset, equal to the currently iterated sample size. 
          SelSamples <- Samples[selectIndex,]
          
          # Assume class balanced.
          GroupId=rep(1,nrow(SelSamples))
          GroupId[(round(sampSizes[currSampSize])+1):length(GroupId)] <- 2
          
          # Extract correlation vector pertaining to how each variable in the datset correlates 
          # with the currently iterated variable.
          corrVector <- correlationMat[,currVar]
          
          # Introduce difference between classes.
          for (i in 1:numVars){
            if (abs(corrVector[i])>0.8) {
              SelSamples[GroupId==2,i] <- SelSamples[GroupId==2,i] + b1[currVar]*sd(SelSamples[,i])
            }
          }
          
          # Center and scale the new data subset.
          SelSamples <- scale(SelSamples)
          
          # Create p-value storage vector. 
          p <- rep(0, numVars)
          
          # Conduct anova for each variable in the dataset, extract and store p-values.
          for (i in 1:numVars) {
            model <- aov(x~group, data=data.frame(x=SelSamples[,i], group=GroupId))
            p[i] <- summary(model)[[1]][1,5]
          }
          
          pUnc <- p # Store uncorrected p-values.
          pBonf <- p.adjust(p, method = 'bonferroni') # Apply and store bonferroni correction.
          pBH <- p.adjust(p, method = 'BH') # Apply and store BH correction.
          
          # Perform 'calcConfMatrixUniv' function using uncorrected, bonferoni corrected and 
          # BH corrected p-values. 
          uncorrected <- calcConfMatrixUniv(pUnc,corrVector,signThreshold,0.8)
          bonferroni <- calcConfMatrixUniv(pBonf,corrVector,signThreshold,0.8)
          bh <- calcConfMatrixUniv(pBH,corrVector,signThreshold,0.8)
          
          # Store the uncorrected confusion matrix output values for the current repeat.
          multiplerepeats$Results$TP[currRepeat] <- uncorrected$TPtot
          multiplerepeats$Results$FP[currRepeat] <- uncorrected$FPtot
          multiplerepeats$Results$TN[currRepeat] <- uncorrected$TNtot
          multiplerepeats$Results$FN[currRepeat] <- uncorrected$FNtot
          multiplerepeats$Results$FD[currRepeat] <- uncorrected$FDtot
          
          # Store the bonferroni corrected confusion matrix output values for the current repeat.
          multiplerepeats$Bonferroni$TP[currRepeat] <- bonferroni$TPtot
          multiplerepeats$Bonferroni$FP[currRepeat] <- bonferroni$FPtot
          multiplerepeats$Bonferroni$TN[currRepeat] <- bonferroni$TNtot
          multiplerepeats$Bonferroni$FN[currRepeat] <- bonferroni$FNtot
          multiplerepeats$Bonferroni$FD[currRepeat] <- bonferroni$FDtot
          
          # Store the BH corrected confusion matrix output values for the current repeat.
          multiplerepeats$BH$TP[currRepeat] <- bh$TPtot
          multiplerepeats$BH$FP[currRepeat] <- bh$FPtot
          multiplerepeats$BH$TN[currRepeat] <- bh$TNtot
          multiplerepeats$BH$FN[currRepeat] <- bh$FNtot
          multiplerepeats$BH$FD[currRepeat] <- bh$FDtot
          
          #setTxtProgressBar(pb, currRepeat) # Update progress bar.
          
        }
        # close(pb) # Terminate progress bar. 
        
        # Extract the names of each of the confusion matrix metrics (e.g. TP, FP, TN, FN and FD).
        stats <- colnames(multiplerepeats$Bonferroni)
        
        for (nstat in 1:5) { # For each of the metrics 
          
          # Define currstat as the name of the currently iterated metric.
          currstat <- stats[nstat]
          
          # Store the average uncorrected confusion matrix values (across repeats) generated 
          # using the currently iterated variable, effect size and sample size. 
          uncStruct[[which(names(uncStruct)%in%stats[nstat])]][currEff,currSampSize] <- 
            mean(multiplerepeats$Results[,which(names(multiplerepeats$Results)%in%stats[nstat])])
          
          # Store the standard deviation of the uncorrected confusion matrix values (across repeats) 
          # generated using the currently iterted variable, effect size and sample size. 
          uncStruct[[which(names(uncStruct)%in%paste('S',stats[nstat],sep = ''))]][currEff,currSampSize] <-
            sd(multiplerepeats$Results[,which(names(multiplerepeats$Results)%in%stats[nstat])])
          
          # Store the average bonferroni corrected confusion matrix values (across repeats) generated 
          # using the currently iterated variable, effect size and sample size. 
          bonfStruct[[which(names(bonfStruct)%in%stats[nstat])]][currEff,currSampSize]  <-
            mean(multiplerepeats$Bonferroni[,which(names(multiplerepeats$Bonferroni)%in%stats[nstat])])
          
          # Store the standard deviation of the bonferroni corrected confusion matrix values (across repeats) 
          # generated using the currently iterted variable, effect size and sample size. 
          bonfStruct[[which(names(bonfStruct)%in%paste('S',stats[nstat],sep = ''))]][currEff,currSampSize] <-
            sd(multiplerepeats$Bonferroni[,which(names(multiplerepeats$Bonferroni)%in%stats[nstat])])
          
          # Store the average BH corrected confusion matrix values (across repeats) generated 
          # using the currently iterated variable, effect size and sample size. 
          bhStruct[[which(names(bhStruct)%in%stats[nstat])]][currEff,currSampSize] <- 
            mean(multiplerepeats$BH[,which(names(multiplerepeats$BH)%in%stats[nstat])])
          
          # Store the standard deviation of the BH corrected confusion matrix values (across repeats) 
          # generated using the currently iterted variable, effect size and sample size. 
          bhStruct[[which(names(bhStruct)%in%paste('S',stats[nstat],sep = ''))]][currEff,currSampSize]  <-
            sd(multiplerepeats$BH[,which(names(multiplerepeats$BH)%in%stats[nstat])])
          
        }
        
      }
      
    }
    
    # Combine each of the storage data structures for the currently iterated variable 
    # into single storeVar 'list'.
    storeVar <- list('Uncorrected'=uncStruct,'Bonferroni'=bonfStruct,'BH'=bhStruct)
    
    # Append the storeVar structure for the currently iterated variable to the main output data structure.
    output <- c(output, list(storeVar))
    
    message('')
  }
  
  # Define the names of each list entry of the main output data structure as the names of each of the
  # assessed variables. 
  names(output) <- namesVars
  
  # Create a list of the input parameters provided to the function to faciliate subsequent plotting.
  input_params <- list('effectSizes'=effectSizes,
                       'sampSizes'=sampSizes,
                       'assessVars'=assessVars,
                       'namesVars'=namesVars,
                       'groupSize'=groupSize,
                       'targetVar'=!is.null(targetVar),
                       'numSamps'=nrow(data))
  
  
  # Append the list of input parameters to the main output data structure.
  output <- c(output, 'input_params'= list(input_params))
  
  # Add S3 class of 'PCalc' to main output data structure.
  class(output) <- c('PCalc','list')
  
  # Return main output data structure.
  return(output)
}

# Plotting function.
power_plot <- function(output, correction, metric){
  
  # Ensure ggplot2 package is loaded.
  require('ggplot2')
  
  # Ensure reshape2 package is loaded.
  require('reshape2')
  
  # Define index for correction type and metric.
  correction_index <- c('Uncorrected', 'Bonferroni', 'BH')
  metric_index <- c('TP', 'FP','TN','FN','FD','STP','SFP','STN','SFN','SFD')
  
  # Stop if 'output' object is not of class 'PCalc'
  if (!inherits(output, 'PCalc')) {
    stop("Must provide valid output object of class 'PCalc'.")
    
    # Stop if correction value is not a valid option.  
  } else if (!correction %in% correction_index) {
    stop("Must provide valid correction value: 'Uncorrected', 'Bonferroni' or 'BH'.")
    
    # Stop if metric value is not a valid option.  
  } else if (!metric %in% metric_index) {
    stop("Must provide valid metric value: 'TP', 'FP','TN','FN' or 'FD'.")
  }
  
  # Extract input_parameters from 'PCalc' output object.
  sampSizes <- output$input_params$sampSizes  
  effectSizes <- output$input_params$effectSizes
  assessVars <- output$input_params$assessVars 
  namesVars <- output$input_params$namesVars
  groupSize <- output$input_params$groupSize
  targetVar <- output$input_params$targetVar
  numSamps <- output$input_params$numSamps
  
  if (targetVar==F) { # If specifc target variable was not provided to PCalc function.
    
    # Create storage variable.
    plotdata <- NULL
    
    # For each variable assessed.
    for (i in 1:length(assessVars)) {
      
      # Extract the final confusion matrix values for the specified metric and at the specified 
      # level of correction.  
      plotdata <- rbind(plotdata, 
                        output[[i]]
                        [[which(correction_index==correction)]]
                        [[which(metric_index==metric)]]
                        [1,])
      
      # Append this to storage object to create single matrix of metric values containing the metric values
      # for each of the assessed variables (rows) and at each of the sample sizes used (columns). 
      plotdata <- as.matrix(plotdata)
    }
    
    # If metric being explored is true positive rate then use 'power for axis title'.
    if (metric=='TP') {
      metric_title <- 'Power'
      
      # Otherwise simply append 'Rate' to the metric name.
    } else {
      metric_title <- paste(metric, ' Rate', sep = '')
    }
    
    # Convert the plot data object to a dataframe.
    plotdata <- as.data.frame(plotdata)
    
    # Label the columns with the relevent sample size information.
    colnames(plotdata) <- sampleSizes
    
    # Add extra column containing the name of each variable (this corrosponds with 
    # each row in plot data object)
    plotdata$key <- namesVars
    
    # Add extra column containing the effect size of each of the assessed variables.
    plotdata$effectSizes <- effectSizes
    
    # Convert the data to tidy format using the melt function from reshape2 package 
    # using both 'key' and 'effectSizes' columns as key variables. 
    plotdata <- reshape2::melt(plotdata, id.vars = c('key','effectSizes'))
    
    # Convert the key column to a factor and order the levels according to the effect sizes.
    # This results in the variables being plotting in order of decreasing effect size, with 
    # largest effect sizes at the top of the graph and smallest effect size at the bottom. 
    plotdata$key <- factor(plotdata$key, 
                           levels=unique(namesVars[order(effectSizes, decreasing = T)]))
    
    # Create 'facets' variable for storage.
    facets <- NULL
    
    if (is.na(groupSize[1])) {
      for (i in 1:length(assessVars)) {
        
        # Create a facet label that contains the name of the feature and the number of other 
        # variables which grouped with it. Append this label to vector of facet labels. 
        facets <- c(facets, namesVars[i])
        
        # Name each entry in vector of facet labels after the name of the relevent variable. This is
        # how ggplot2 labeller function requires facet labels to be constructed. 
        names(facets)[i] <- namesVars[i]
      }
    } else {
      for (i in 1:length(assessVars)) {
        
        # Create a facet label that contains the name of the feature and the number of other 
        # variables which grouped with it. Append this label to vector of facet labels. 
        facets <- c(facets, paste(namesVars[i], '\n(Group size: ', groupSize[i], ')', sep = ''))
        
        # Name each entry in vector of facet labels after the name of the relevent variable. This is
        # how ggplot2 labeller function requires facet labels to be constructed. 
        names(facets)[i] <- namesVars[i]
      }
    }
    
    # For each of the variables assessed. 
    for (i in 1:length(assessVars)) {
      
      # Create a facet label that contains the name of the feature and the number of other 
      # variables which grouped with it. Append this label to vector of facet labels. 
      facets <- c(facets, paste(namesVars[i], '\n(Group size: ', groupSize[i], ')', sep = ''))
      
      # Name each entry in vector of facet labels after the name of the relevent variable. This is
      # how ggplot2 labeller function requires facet labels to be constructed. 
      names(facets)[i] <- namesVars[i]
    }
    
    # Specify colour scheme to use for plot...
    # For purple, use:
    palette <- c(3,"#CBC4D5")
    
    # For green, use:
    #palette <- c(10,'#CAD5C4')
    
    # Call ggplot:  plotting sample sizes along x axis and effect sizes along the y axis.
    ggplot(plotdata, aes(x=variable, y=as.factor(effectSizes))) +
      
      # Allow custom point shapes. 
      scale_shape_identity() +
      
      # Plot metric values as the size and colour of diamond shaped points corrosponding to each
      # XY position.
      geom_point(aes(size=value, fill=value), shape=23) +
      
      # Set the colours and break points of the fill aesthetic.
      scale_fill_distiller(palette=as.numeric(palette[1]),
                           limits=c(0, 1),
                           breaks=seq(0, 1, by=0.25),
                           direction = 1) +
      
      # Set the break points of the size aesthetic. 
      scale_size_continuous(limits=c(0, 1),
                            breaks=seq(0, 1, by=0.25)) +
      
      # Matching the fill and size legend guides allow legends to be merged into single unified 
      # legend.
      guides(fill= guide_legend(), size=guide_legend()) +
      
      facet_grid(rows=vars(as.factor(key)), # Create new facet for each variable.
                 scales = 'free_y', # Allow y axis sscale to vary between facets.
                 labeller = as_labeller(facets), # Use modified facet labels.
                 switch = 'y') + # Move facet labels to left side.
      
      ylab('Effect Size') + # Modify x axis label.
      xlab('Sample Size') + # Modify y axis label.
      labs(fill=metric_title, size=metric_title) + # Modify legend label.
      
      scale_y_discrete(position = 'right') + # Move y axis ticks to right side of plot.
      
      theme_bw() + # Change plot theme.
      theme(legend.position = 'left', # Move legend to left side of plot.
            legend.justification = 'top', # Justify legend to top of plot.
            legend.title.align = 0.5, # Centre the legend title. 
            
            # Colour the legend according to chosen pallete and add black border. 
            legend.background = element_rect(colour = "black", 
                                             size=.5, 
                                             fill=palette[2]),
            legend.key = element_rect(fill=palette[2], 
                                      colour = palette[2]),
            
            # Ensure facet labels are oriented correctely.
            strip.text.y = element_text(angle = 180),
            
            # Colour the facet labels accoring to chosen pallete.
            strip.background =element_rect(fill=palette[2]))
    
  } else if (targetVar==T) { # # If specifc target variable was provided to PCalc function.
    
    # Extract relevent data based on specified metirc and correction.
    plotdata <- output[[1]][[which(correction_index==correction)]][[which(metric_index==metric)]]
    
    # Call ggplot: plotting sample sizes along x axis and metric values along the y axis.
    ggplot(data.frame(sampleSizes,plotdata), aes(x=sampleSizes,y=plotdata[1,])) +
      geom_point(alpha=0.5) + # Add points (with alpha).
      geom_smooth(method = 'loess') + # Add smooth line between points.
      
      # Change plot theme.  
      theme_bw() +
      
      # Modify y axis breaks.  
      ylim(0,1) +
      
      # Modify x axis label.
      xlab('Sample Sizes') +
      
      # Modify y axis label
      ylab('Power') 

  }
}