###########################################################################
################ Feature selection - Classification mode ##################
###########################################################################

# Single instance feature selection function (Classification mode).
feature_selection <- function(i, xdata, ydata) {
  
  # Load Pomona package.
  require(Pomona)
  
  # Load Boruta package.
  require(Boruta)
  
  # Load varselRF package.
  require(varSelRF)
  
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
  RFE <- varSelRF(xdata = xdata,
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
                   doTrace = 1,
                   num.trees = 2000)
  
  # Extract the variables selected by boruta (those attributed with 'Confirmed').
  selected.var.Boruta <- names(Boruta$finalDecision)[which(Boruta$finalDecision=='Confirmed')]
  
  # For each of the features selected by this method this iteration add one to tally reporting its frequency of selection by this method.
  varselrfstore$borutafreq[varselrfstore$corvariables %in% selected.var.Boruta] <- 1
  
  # Run permutation method for feature selection.
  Permutation <- var.sel.perm(x=xdata,
                              y=ydata,
                              ntree=2000,
                              no.threads = 1,
                              type = 'classification')
  
  # Extract features with a raw p-values below 0.05.
  selected.var.raw.permutation <- which(Permutation$info$pvalue<0.05)
  
  
  # Extract features with a BH corrected p-value below 0.05.
  selected.var.corrected.permutation <- which(p.adjust(Permutation$info$pvalue, method='BH')<0.05)
  
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