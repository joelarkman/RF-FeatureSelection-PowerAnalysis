#################################################################
######################## Data Cleaning ##########################
#################################################################

# Randomely permute the index of the row numbers. 
sample <- sample.int(nrow(data))


## Generate data index for first outerloop test data.
subsample <-sample[1:(1*round(nrow(data)*.25))] # Take the first 25% index values. 
outertrain_1 <- data[-subsample, ] # Store all data besides index as first outer train.
outertest_1  <- data[subsample, ] # Store indexed data as first outer test.

# Randomely select 90% of rows from first outer train.
subsample <- sample.int(n = nrow(outertrain_1), size = floor(.9*nrow(outertrain_1)), replace = F)
innertrain_1 <- outertrain_1[-subsample, ] # Store all data besides index as first innter train.
innertest_1  <- outertrain_1[subsample, ] # Store indexed data as first inner test. 



## Generate data index for second outerloop test data.
subsample <-sample[(1*round(nrow(data)*.25)+1):(2*round(nrow(data)*.25))] # Take the next 25% index values.
outertrain_2 <- data[-subsample, ] # Store all data besides index as second outer train.
outertest_2  <- data[subsample, ] # Store indexed data as second outer test.

# Randomely select 90% of rows from second outer train.
subsample <- sample.int(n = nrow(outertrain_2), size = floor(.9*nrow(outertrain_2)), replace = F)
innertrain_2 <- outertrain_2[-subsample, ] # Store all data besides index as second innter train.
innertest_2  <- outertrain_2[subsample, ] # Store indexed data as second inner test. 



## Generate data index for second outerloop test data.
subsample <-sample[(2*round(nrow(data)*.25)+1):(3*round(nrow(data)*.25))] # Take the next 25% index values.
outertrain_3 <- data[-subsample, ] # Store all data besides index as third outer train.
outertest_3  <- data[subsample, ] # Store indexed data as third outer test.

# Randomely select 90% of rows from third outer train.
subsample <- sample.int(n = nrow(outertrain_3), size = floor(.9*nrow(outertrain_3)), replace = F)
innertrain_3 <- outertrain_3[-subsample, ] # Store all data besides index as third innter train.
innertest_3  <- outertrain_3[subsample, ] # Store indexed data as third inner test. 



## Generate data index for final outerloop test data.
subsample <-sample[(3*round(nrow(data)*.25)+1):(4*round(nrow(data)*.25))] # Take the last 25% index values.
outertrain_4 <- data[-subsample, ] # Store all data besides index as final outer train.
outertest_4  <- data[subsample, ] # Store indexed data as final outer test.

# Randomely select 90% of rows from final outer train.
subsample <- sample.int(n = nrow(outertrain_4), size = floor(.9*nrow(outertrain_4)), replace = F)
innertrain_4 <- outertrain_4[-subsample, ] # Store all data besides index as final inner train.
innertest_4  <- outertrain_4[subsample, ] # Store indexed data as final inner test. 



# Combine all outer test data into single list for easy referencing.
outertest <- list('outertest_1'=outertest_1,'outertest_2'=outertest_2,'outertest_3'=outertest_3,'outertest_4'=outertest_4)

#################################################################
####################### Regression plot #########################
#################################################################


# Provide list of stable features found by Boruta.
variables <- c('SM(39:1)','PC(34:2','SM(36:2)')


# Store relevent data for each of the selected features.
v1 <- data.frame(y=data$`Total formula`,x=data$`SM(39:1)`, v=variables[1])
v2 <- data.frame(y=data$`Total formula`,x=data$`PC(34:2)`, v=variables[2])
v3 <- data.frame(y=data$`Total formula`,x=data$`SM(36:2)`, v=variables[3])

# Combine relevent data.
data <- rbind(v1,v2,v3)

# Define y variable name.
y_variable <- 'Infant milk amount'

# Create null variables. 
plotdata <- NULL
plotlabels <- NULL
for (i in 1:length(variables)) { # For each of the selected variables. 
  x_variable <- variables[i] # Take name of current variable.
  
  lmdat <- subset(data, v==variables[i])[1:2] # Extract relevent data for current variable.
  
  fit <- lm(y~x, lmdat) # Conduct LM.
  
  tmp.data <- data.frame(fit$model,"facet"=variables[i]) # Extract results.
  
  plotdata <- rbind(plotdata,tmp.data) # Store results.
  
  # Create custom facet label.
  plotlabels <- c(plotlabels,paste(variables[i],'\n',"Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                                   " P =",signif(summary(fit)$coef[2,4], 5)))
}

# Modify names of plot label vector.
names(plotlabels) <- variables

# Load ggplot package.
library(ggplot2)

# Create plot
ggplot(plotdata, aes(x = x, y = y)) + 
  facet_wrap(~facet, labeller = as_labeller(plotlabels), scales = 'free_x') +
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  xlab(element_blank()) +
  ylab(y_variable)

# Store plot.
pdf('plot6.pdf', width=7, height=2.5)
ggplot(plotdata, aes(x = x, y = y)) + 
  facet_wrap(~facet, labeller = as_labeller(plotlabels), scales = 'free_x') +
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  xlab(element_blank()) +
  ylab(y_variable) +
  theme_bw()
dev.off()




