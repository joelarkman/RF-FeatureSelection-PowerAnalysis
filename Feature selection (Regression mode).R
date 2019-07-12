###########################################################################
################## Feature selection - Regression mode ####################
###########################################################################

# Single instance feature selection function (Regression mode).
feature_selection <- function(i, xdata, ydata) {
  
  require(Pomona)
  
  # Load Boruta package.
  require(Boruta)
  
  # Load varselRF package.
  require(varSelRF)
  
  # Modify varselRF package to work with regression RF.
  modified.varSelRF <- function (xdata, Class, c.sd = 1, mtryFactor = 1, ntree = 5000,
                                 ntreeIterat = 2000, vars.drop.num = NULL, vars.drop.frac = 0.2,
                                 whole.range = TRUE, recompute.var.imp = FALSE, verbose = FALSE,
                                 returnFirstForest = TRUE, fitted.rf = NULL, keep.forest = FALSE)
  {
    if ((is.null(vars.drop.num) & is.null(vars.drop.frac)) |
        (!is.null(vars.drop.num) & !is.null(vars.drop.frac)))
      stop("One (and only one) of vars.drop.frac and vars.drop.num must be NULL and the other set")
    max.num.steps <- dim(xdata)[2]
    num.subjects <- dim(xdata)[1]
    if (is.null(colnames(xdata)))
      colnames(xdata) <- paste("v", 1:dim(xdata)[2], sep = "")
    n.vars <- vars <- MSE.rf <- RMSE.rf <- rep(NA, max.num.steps)
    if (!is.null(fitted.rf)) {
      if (ncol(fitted.rf$importance) < 2)
        stop("The fitted rf was not fitted with importance = TRUE")
      n.ntree <- fitted.rf$ntree
      mtry <- fitted.rf$mtry
      n.mtryFactor <- mtry/sqrt(ncol(xdata))
      if ((n.ntree != ntree) | (n.mtryFactor != mtryFactor))
        warning("Using as ntree and mtry the parameters obtained from fitted.rf",
                immediate. = TRUE)
      ntree <- n.ntree
      mtryFactor <- n.mtryFactor
      rm(n.ntree, n.mtryFactor)
      rf <- fitted.rf
    }
    else {
      mtry <- floor(sqrt(ncol(xdata)) * mtryFactor)
      rf <- randomForest(x = xdata, y = Class, ntree = ntree,
                         mtry = mtry, importance = TRUE, keep.forest = keep.forest)
    }
    if (returnFirstForest)
      FirstForest <- rf
    else FirstForest <- NULL
    
    m.iterated.mse <- m.initial.mse <- rf$mse[length(rf$mse)]
    iterated.rmse <- initial.rmse <- sqrt(m.iterated.mse)
    
    
    if (verbose) {
      print(paste("Initial MSE: mean = ", round(m.initial.mse,
                                                4), "; RMSE = ", round(initial.rmse, 4), sep = ""))
    }
    
    
    importances <- randomForest::importance(rf, type = 1, scale = FALSE)
    selected.vars <- order(importances, decreasing = TRUE)
    ordered.importances <- importances[selected.vars]
    initialImportances <- importances
    initialOrderedImportances <- ordered.importances
    j <- 1
    n.vars[j] <- dim(xdata)[2]
    vars[j] <- paste(colnames(xdata), collapse = " + ")
    
    MSE.rf[j] <- m.iterated.mse
    RMSE.rf[j] <- iterated.rmse
    
    var.simplify <- TRUE
    while (var.simplify) {
      if (verbose) {
        print("gc inside loop of varSelRF")
        print(gc())
      }
      else {
        gc()
      }
      last.rf <- rf
      last.vars <- selected.vars
      
      previous.mse <- m.iterated.mse
      previous.rmse <- iterated.rmse
      
      if (length(selected.vars) <= 2) {
        var.simplify <- FALSE
        break
      }
      if (recompute.var.imp & (j > 1)) {
        importances <- randomForest::importance(rf, type = 1, scale = FALSE)
        tmp.order <- order(importances, decreasing = TRUE)
        selected.vars <- selected.vars[tmp.order]
        ordered.importances <- importances[tmp.order]
      }
      num.vars <- length(selected.vars)
      if (is.null(vars.drop.num))
        vars.drop <- round(num.vars * vars.drop.frac)
      else vars.drop <- vars.drop.num
      if (num.vars >= (vars.drop + 2)) {
        if (vars.drop == 0) {
          vars.drop <- 1
          if ((num.vars - vars.drop) < 1)
            stop("vars.drop = 0 and num.vars -vars.drop < 1!")
        }
        selected.vars <- selected.vars[1:(num.vars - vars.drop)]
        ordered.importances <- ordered.importances[1:(num.vars -
                                                        vars.drop)]
      }
      else {
        selected.vars <- selected.vars[1:2]
        ordered.importances <- ordered.importances[1:2]
      }
      if ((length(selected.vars) < 2) | (any(selected.vars <
                                             1))) {
        var.simplify <- FALSE
        break
      }
      mtry <- floor(sqrt(length(selected.vars)) * mtryFactor)
      if (mtry > length(selected.vars))
        mtry <- length(selected.vars)
      if (recompute.var.imp)
        rf <- randomForest(x = xdata[, selected.vars], y = Class,
                           importance = TRUE, ntree = ntree, mtry = mtry,
                           keep.forest = keep.forest)
      else rf <- randomForest(x = xdata[, selected.vars], y = Class,
                              importance = FALSE, ntree = ntreeIterat, mtry = mtry,
                              keep.forest = keep.forest)
      
      m.iterated.mse <- rf$mse[length(rf$mse)]
      iterated.rmse <- sqrt(m.iterated.mse)
      
      if (verbose) {
        print(paste("..... iteration ", j, "; MSE: mean = ",
                    round(m.iterated.mse, 4), "; RMSE = ", round(iterated.rmse,
                                                                 4), "; num. vars = ", length(selected.vars),
                    sep = ""))
      }
      j <- j + 1
      n.vars[j] <- length(selected.vars)
      vars[j] <- paste(colnames(xdata)[selected.vars], collapse = " + ")
      
      MSE.rf[j] <- m.iterated.mse
      RMSE.rf[j] <- iterated.rmse
      if (!whole.range & ((m.iterated.mse > (m.initial.mse +
                                             c.sd * initial.rmse)) | (m.iterated.mse >
                                                                      (previous.mse + c.sd * previous.rmse))))
        var.simplify <- FALSE
    }
    if (!whole.range) {
      if (!is.null(colnames(xdata)))
        selected.vars <- sort(colnames(xdata)[last.vars])
      else selected.vars <- last.vars
      out <- list(selec.history = data.frame(Number.Variables = n.vars,
                                             Vars.in.Forest = vars, MSE = MSE.rf, RMSE = RMSE.rf)[1:j,
                                                                                                  ], rf.model = last.rf, selected.vars = selected.vars,
                  selected.model = paste(selected.vars, collapse = " + "),
                  best.model.nvars = length(selected.vars), initialImportances = initialImportances,
                  initialOrderedImportances = initialOrderedImportances,
                  ntree = ntree, ntreeIterat = ntreeIterat, mtryFactor = mtryFactor,
                  firstForest = FirstForest)
      class(out) <- "varSelRF"
      return(out)
    }
    else {
      n.vars <- n.vars[1:j]
      vars <- vars[1:j]
      MSE.rf <- MSE.rf[1:j]
      RMSE.rf <- RMSE.rf[1:j]
      min.mse.ci <- min(MSE.rf) + c.sd * RMSE.rf[which.min(MSE.rf)]
      best.pos <- which(MSE.rf <= min.mse.ci)[which.min(n.vars[which(MSE.rf <=
                                                                       min.mse.ci)])]
      selected.vars <- sort(unlist(strsplit(vars[best.pos],
                                            " + ", fixed = TRUE)))
      out <- list(selec.history = data.frame(Number.Variables = n.vars,
                                             Vars.in.Forest = vars, MSE = MSE.rf, RMSE = RMSE.rf),
                  rf.model = NA, selected.vars = selected.vars, selected.model = paste(selected.vars,
                                                                                       collapse = " + "), best.model.nvars = n.vars[best.pos],
                  initialImportances = initialImportances, initialOrderedImportances = initialOrderedImportances,
                  ntree = ntree, ntreeIterat = ntreeIterat, mtryFactor = mtryFactor,
                  firstForest = FirstForest)
      class(out) <- "varSelRF"
      return(out)
    }
  }
  
  # Define storage DF to store which features are selected by each method.
  varselrfstore <- data.frame(corvariables=colnames(xdata),
                              num=1:ncol(xdata),
                              RFEfreq=0,
                              borutafreq=0,
                              raw.permutationfreq=0,
                              corrected.permutationfreq=0)
  
  # Print current iteration for reference.
  print(paste('Iteration: ', i, sep = ''))
  
  # Run modified varSelRF package to faciliate recursive feature elimination (backwards elimination).
  RFE <- modified.varSelRF(xdata = xdata,
                           Class= ydata,
                           verbose = T,
                           mtryFactor = 1,
                           c.sd = 0,
                           ntree = 2000)
  
  # For each of the features selected by this method this iteration add one to tally reporting its frequency of selection by this method.
  varselrfstore$RFEfreq[varselrfstore$corvariables %in% RFE$selected.vars] <- 1
  
  # Run boruta method for feature selection.
  Boruta <- Boruta(x=xdata,
                   y=ydata,
                   pValue = 0.01,
                   doTrace = 1,
                   ntree=2000)
  
  # Extract the variables selected by boruta (those attributed with 'Confirmed').
  selected.var.Boruta <- names(Boruta$finalDecision)[which(Boruta$finalDecision=='Confirmed')]
  
  # For each of the features selected by this method this iteration add one to tally reporting its frequency of selection by this method.
  varselrfstore$borutafreq[varselrfstore$corvariables %in% selected.var.Boruta] <- 1
  
  # Run permutation method for feature selection.
  Permutation <- var.sel.perm(x=xdata,
                              y=ydata,
                              ntree=2000,
                              no.threads = 1)
  
  # Extract features with a raw p-values below 0.05.
  selected.var.raw.permutation <- which(Permutation$info$pvalue<0.05)
  
  print(selected.var.raw.permutation)
  
  # Extract features with a BH corrected p-value below 0.05.
  selected.var.corrected.permutation <- which(p.adjust(Permutation$info$pvalue, method='BH')<0.05)
  
  print(selected.var.corrected.permutation)
  
  # For each of the features selected by this method this iteration (using raw or corrected p-values) add one to tally reporting its frequency of selection by this method.
  varselrfstore$raw.permutationfreq[varselrfstore$num %in% selected.var.raw.permutation] <- 1
  varselrfstore$corrected.permutationfreq[varselrfstore$num %in% selected.var.corrected.permutation] <- 1
  
  return(varselrfstore)
  
}

# Multiple iterations of feature selection (Note: Very time consuming).
iterate_feature_selection <- function(n_iterations=100, xdata, ydata) {
  
  # Create storage variable
  sims_100 <- NULL
  
  for (i in 1:n_iterations) { # For each of the specified no. iterations.
    tmp <- feature_selection(i,xdata,ydata) # Run the feature selection function.
    sims_100 <- c(sims_100,list(tmp)) # Store whether or not each variable was chosen by each method.
  }
  
  # Create dataframe for collective results.
  varselrfstore <- data.frame(corvariables=colnames(xdata),
                              num=1:ncol(xdata),
                              RFEfreq=0,
                              borutafreq=0, 
                              raw.permutationfreq=0, 
                              corrected.permutationfreq=0)
  
  # Add the tally each feature was chosen by each method across all iterations.
  for (i in 1:length(sims_100)) {
    varselrfstore$RFEfreq <- varselrfstore$RFEfreq + sims_100[[i]]$RFEfreq
    varselrfstore$borutafreq <- varselrfstore$borutafreq + sims_100[[i]]$borutafreq
    varselrfstore$raw.permutationfreq <- varselrfstore$raw.permutationfreq + sims_100[[i]]$raw.permutationfreq
    varselrfstore$corrected.permutationfreq <- varselrfstore$corrected.permutationfreq + sims_100[[i]]$corrected.permutationfreq
  }
  
  # Return the total number of times each feature was selected by each method after specified no. iterations. 
  return(varselrfstore)
}

###########################################################################
############################## Example ####################################
###########################################################################

# Iterate function 100 times for each of the inner test datasets.
# WARNING: These functions can take several hours to complete.
selected.vars.innertest1 <- iterate_feature_selection(100, innertest_1[,-1], innertest_1[,1])
selected.vars.innertest2 <- iterate_feature_selection(100, innertest_2[,-1], innertest_2[,1])
selected.vars.innertest3 <- iterate_feature_selection(100, innertest_3[,-1], innertest_3[,1])
selected.vars.innertest4 <- iterate_feature_selection(100, innertest_4[,-1], innertest_4[,1])

# Average across all four loops, observe the average number of times each feature is selected.
selected.vars.average <- data.frame(selected.vars.innertest1[,1:2],
                            'RFEfreq'=rowMeans(cbind(selected.vars.innertest1[,3],selected.vars.innertest2[,3],selected.vars.innertest3[,3],selected.vars.innertest4[,3])),
                            'borutafreq'=rowMeans(cbind(selected.vars.innertest1[,4],selected.vars.innertest2[,4],selected.vars.innertest3[,4],selected.vars.innertest4[,4])),
                            'raw.permutationfreq'=rowMeans(cbind(selected.vars.innertest1[,5],selected.vars.innertest2[,5],selected.vars.innertest3[,5],selected.vars.innertest4[,5])),
                            'corrected.permutationfreq'=rowMeans(cbind(selected.vars.innertest1[,6],selected.vars.innertest2[,6],selected.vars.innertest3[,6],selected.vars.innertest4[,6])))

###########################################################################
########################## Stacked Bar Plot ###############################
###########################################################################

# Create storage data frame. 
store <- data.frame(Type=c('Known','Novel','Known','Novel'),
                    Key=c('High','High','Low','Low'),
                    RFEfreq=0,
                    borutafreq=0, 
                    raw.permutationfreq=0, 
                    corrected.permutationfreq=0)

# Provide vector of 'known' variables of interest.
known <- c('','','')

for ( i in 1:4) { # For each of the four methods.
  # Count the number of known features selected >= 90 times (High stringency)
  store[1,i+2] <- length(which(selected.vars.average[which(selected.vars.average[i+2]>=90),1] %in% known))
  # Count the number of novel features selected >= 90 times (High stringency)
  store[2,i+2] <- length(which(!selected.vars.average[which(selected.vars.average[i+2]>=90),1] %in% known))
  # Count the number of known features selected >= 5 times (Low stringency)
  store[3,i+2] <- length(which(selected.vars.average[which(selected.vars.average[i+2]>=5),1] %in% known))
  # Count the number of novel features selected >= 5 times (Low stringency)
  store[4,i+2] <- length(which(!selected.vars.average[which(selected.vars.average[i+2]>=5),1] %in% known))
}

# Load reshape2 package.
library(reshape2)
# Load ggplot2 package.
library(ggplot2)

# Create facet labels.
labels <- as_labeller(c('RFEfreq'='RFE','borutafreq'='Boruta','raw.permutationfreq'='Permutation\n(Raw)', 'corrected.permutationfreq'='Permutation\n(Corrected)'))

# Use melt function on 'store' dataframe and plot using ggplot. 
plot3 <- ggplot(melt(store), aes(x = Key,y = value,
                                 fill=factor(Type, levels = c('Known', 'Novel')))) +
  facet_grid(~variable, labeller = labels) +
  geom_bar(stat='identity') +
  labs(fill='Type') +
  theme_bw() +
  xlab('Stringency') +
  ylab('Feature Frequency') +
  geom_text(data = subset(melt(store), value!=0), aes(label = value),
            position = position_stack(vjust = .5), size=2) +
  theme(legend.position = 'top')


# Save plot. 
pdf('plot3.pdf', width=6, height=4)
print(plot3)
dev.off()

###########################################################################
####################### Validation model & plot ###########################
###########################################################################

# Function to produce validation models and extract power values.
modelpower <- function(outertest, selected.vars) {
  
  library(randomForest)
  
  outertest <- data
  
  # Define storage variable.
  stable.variables <- NULL
  
  # For each of the feature selection methods.
  for (i in 1:4) {
    
    # Extract the features selected by each method in more than 90% of iterations (Stable features) and add them to a list of stable features.
    stable.variables <- c(stable.variables, list(
      as.character(
        varselrfstore$corvariables[which(varselrfstore[,(i+2)]>5)])))
    
    # Name each list item after the name of the feature selection method responsible.
    names(stable.variables)[i] <- colnames(varselrfstore)[(i+2)]
  }
  
  # Create storage data frame.
  power <- data.frame('MSE'=0, 'Rsquared'=0,'vars'=names(stable.variables))
  
  # For each method.
  for (i in 1:length(stable.variables)) {
    
    print(stable.variables[i])
    
    # Set temp variables to NULL
    pred.error <- NULL
    rsquared <- NULL
    
    # Create temporary dataset featuring y and the stable variables from the outer-test data subset.
    data <- cbind(y=outertest$y,outertest[,c(stable.variables[[i]])])
    
    # Run randomforest 100 times and extract MSE and Rsquared values. 
    for (j in 1:10) {
      print(paste('Iteration: ',j,sep = ''))
      RF <- randomForest(x=data[,-1], y=data[,1],num.trees = 10000,importance = T)
      pred.error <- c(pred.error,tail(RF$mse,n = 1))
      rsquared <- c(rsquared,tail(RF$rsq,n = 1))
      
    }
    
    # Store the average values for each method.
    power[i,1] <- mean(pred.error)
    power[i,2] <- mean(rsquared)
    
  }
  
  # Print results.
  power_low <- power
  
  # Define storage variable.
  stable.variables <- NULL
  
  # For each of the feature selection methods.
  for (i in 1:4) {
    
    # Extract the features selected by each method in more than 90% of iterations (Stable features) and add them to a list of stable features.
    stable.variables <- c(stable.variables, list(
      as.character(
        selected.vars$corvariables[which(selected.vars[,(i+2)]>=90)])))
    
    # Name each list item after the name of the feature selection method responsible.
    names(stable.variables)[i] <- colnames(selected.vars)[(i+2)]
  }
  
  # Create storage data frame.
  power <- data.frame('MSE'=0, 'Rsquared'=0,'vars'=names(stable.variables))
  
  # For each method.
  for (i in 1:length(stable.variables)) {
    
    print(stable.variables[i])
    
    # Set temp variables to NULL
    pred.error <- NULL
    rsquared <- NULL
    
    # Create temporary dataset featuring y and the stable variables from the outer-test data subset.
    data <- cbind(y=outertest$y,outertest[,c(stable.variables[[i]])])
    
    # Run randomforest 10 times and extract MSE and Rsquared values. 
    for (j in 1:10) {
      print(paste('Iteration: ',j,sep = ''))
      RF <- randomForest(x=data[,-1], y=data[,1],num.trees = 10000,importance = T)
      pred.error <- c(pred.error,tail(RF$mse,n = 1))
      rsquared <- c(rsquared,tail(RF$rsq,n = 1))
      
    }
    
    # Store the average values for each method.
    power[i,1] <- mean(pred.error)
    power[i,2] <- mean(rsquared)
    
  }
  
  # Print results.
  print(power)
  
  # Store HS and LS power results.
  power_high <- power
  power_high <- data.frame(power_high, key='High')
  power_low <- data.frame(power_low, key='Low')
  power <- rbind(power_high,power_low)
  
  # Return results
  return(power)
}

# Example running function four times for each of the outerloops, using outer test data.
power_1 <- modelpower(outertest$outertest_1, selected.vars.innertest1)
power_2 <- modelpower(outertest$outertest_2, selected.vars.innertest2)
power_3 <- modelpower(outertest$outertest_3, selected.vars.innertest3)
power_4 <- modelpower(outertest$outertest_4, selected.vars.innertest4)

# Combine and melt result power values. 
power <- melt(rbind(power_1,power_2,power_3,power_4))

# Create Facet labels.
labels2 <- as_labeller(c('High'='High Stringency', 'Low'='Low Stringency'))

# Produce plot using Rsquared model metric for Y value. 
plot5 <- ggplot(power[which(power$variable=='Rsquared'),], aes(x=vars, y=value)) + 
  
  # Facet wrap plots according to Variable.
  facet_wrap(~key, labeller = labels2) +
  
  # Produce boxplots (coloured according to Day)
  geom_boxplot(aes(fill=vars)) +
  
  scale_x_discrete(labels=c("borutafreq" = "Boruta", "corrected.permutationfreq" = "Permutation\n(Corrected)",
                            "raw.permutationfreq" = "Permutation\n(Raw)", 'RFEfreq'='RFE')) +
  
  scale_fill_discrete(labels=c("borutafreq" = "Boruta", "corrected.permutationfreq" = "Permutation (Corrected)",
                               "raw.permutationfreq" = "Permutation (Raw)", 'RFEfreq'='RFE')) +
  
  xlab('Method') +
  ylab('R-Squared') +
  labs(fill='Method') +
  
  theme_bw() +
  theme(legend.position = 'top')

# Save Plot
pdf('plot5.pdf', width=8, height=4.5)
print(plot5)
dev.off()


