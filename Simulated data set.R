# Set seed to ensure simulated data is always the same.
set.seed(10)

####################### Data simulation #######################
###############################################################

# Define no. samples to simulate.
replicates <- 200

#replicates <- 740

# Create storage vector for Y variable and matrix for X variables.
y <- rep(0,replicates)
simulation1 <- matrix(0, nrow = replicates, ncol = 5000)

# For each replicate.
for (r in 1:replicates) {
  
  # Generate 6 variables from uniform distribution.
  x <- runif(6, min=0,max=1) # Six variables generated from uniform distribution.
  
  # Specify group size.
  n <- 10
  
  # Reset stimulation_data variable.
  simulation_data <- NULL
  
  # For each of the 6 uniform distribution vairbles (x), generate a group of length n of correlated variables whereby as n increases the correlation between variables decrease. 
  for (i in 1:length(x) ) {
    for (j in 1:n) {
      # Store value in the jth position within the group of length n of varibles correlated with the currently iterated x variable.
      simulation_data[j+n*(i-1)] <- x[i] + ((0.01+(0.5*(j-1))/(n-1))*rnorm(1,0,0.3))
    }
  }
  
  # Store the resultant n*6 variables as the first columns in the matrix of x variables.
  simulation1[r,1:(n*6)] <- simulation_data
  
  # Produce a Y variable that is correlated with the first three of the uniform distribution variables (x).
  # Y is most correlated with x[1], less so with x[2] and less further still with x[3].
  y[r] <- (( 0.25 * exp((4*x[1])) ) + (4 / ( 1 + exp(-20 * (x[2]-0.5)) )) + (3*x[3])) + rnorm(1,0,0.2)
  
  # Overall this code generates a value for Y and vlaues for x*n correlated features.
  # When n=10, the first 30 variables are causal for Y, with the first 10 more correlated than 11-20 and 21-30.
  # The last three values for X are just used to create correlated features and are not causal for y.
}

# Append on to the n*6 correlated variables, enough non-correlated random uniform variables to produce 5000 features overall.
for (i in ((n*6)+1):5000){
  simulation1[,i] <- runif(replicates, min=0,max=1)
}

# Convert to dataframe.
simulation1 <- as.data.frame(simulation1)
y <- as.data.frame(y)

# Create a two-group classfication outcome, splitting samples above and below the mean of Y into two groups. 
class_y <- rep(0, length(y$y))
class_y[y$y < mean(y$y)] <- 'Group 1'
class_y[y$y > mean(y$y)] <- 'Group 2'

# Convert class_y to factor.
class_y <- as.factor(class_y)

# Produce a combined dataset with continuous outcome.
data <- cbind(y,simulation1)

# Produce a combined dataset with two-group classification outcome.
# data <- cbind(y=as.factor(class_y),simulation1) 


####################### Data Partitioning #######################
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
