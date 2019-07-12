###########################################################################
####################### Hyperparameter optimisation #######################
###########################################################################

# Load caret package.
require(caret)

# Load RF package.
require(randomForest)

# Define a custom RF (Regression) method to allow optimisation of ntree parameter as well as mtry.
# customRF <- list(type = "Regression", library = "randomForest", loop = NULL)
# customRF$parameters <- data.frame(parameter = c("mtry", "ntree", "nodesize"), class = rep("numeric", 3), label = c("mtry", "ntree", "nodesize"))
# customRF$grid <- function(x, y, len = NULL, search = "grid") {}
# customRF$fit <- function(x, y, wts, param, lev, last, weights, ...) {
#   randomForest(x, y, mtry = param$mtry, ntree=param$ntree, nodesize = param$nodesize, ...)
# }
# customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
#   predict(modelFit, newdata)
# customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
#   predict(modelFit, newdata, type = "prob")
# customRF$sort <- function(x) x[order(x[,1]),]
# customRF$levels <- function(x) x$classes

# Define a custom RF (Classification) method to allow optimisation of ntree parameter as well as mtry.
customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree", "nodesize"), class = rep("numeric", 3), label = c("mtry", "ntree", "nodesize"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, nodesize = param$nodesize, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

# Define default mtry value.
mtry <- sqrt(ncol(innertrain_1)-1)

# Define values for mtry and ntree to explore. 
# For mtry assess 50%, 75%, 100%, 125% and 150% of the default value.
tunegrid <- expand.grid(.mtry=c(round(0.5*mtry),round(0.75*mtry),round(mtry),round(1.25*mtry),round(1.5*mtry)), 
                        .ntree=c(1000, 1500, 2000, 2500),
                        .nodesize=c(1,5,10))

# Apply to first inner train set.
custom1 <- train(y~.,
                 data=innertrain_1,
                 method=customRF,
                 metric='Accuracy',
                 tuneGrid=tunegrid,
                 verbose=T)

# Apply to second inner train set.
custom2 <- train(y~.,
                 data=innertrain_2,
                 method=customRF,
                 metric='Accuracy',
                 tuneGrid=tunegrid,
                 verbose=T)

# Apply to third inner train set.
custom3 <- train(y~.,
                 data=innertest_3,
                 method=customRF,
                 metric='Accuracy',
                 tuneGrid=tunegrid,
                 verbose=T)

# Apply to final inner train set. 
custom4 <- train(y~.,
                data=innertrain_4,
                method=customRF,
                metric='Accuracy',
                tuneGrid=tunegrid,
                verbose=T)


# Load data.table package.
library(data.table)

# Average the values generated for each assessed model across the four inner test datasets. 
custom <- rbindlist(list(custom1$results,
                         custom2$results,
                         custom3$results,
                         custom4$results))[,lapply(.SD,mean), list(mtry,ntree,nodesize)]

# Produce plot of results. 
labels <- c('1'='Node Size = 1', '5'='Node Size = 5', '10'='Node Size = 10')
ggplot(custom, aes(x=mtry,y=Accuracy, group=as.factor(ntree), col=as.factor(ntree))) +
  facet_grid(~nodesize, labeller = as_labeller(labels)) +
  geom_line() + geom_point() +
  theme_bw() +
  labs(col='ntree') +
  xlab('mtry') +
  ylab('Rsquared (Repeated Cross-Validation)')

# Save plot
pdf('plot1.pdf', width=8, height=4)
ggplot(custom, aes(x=mtry,y=Accuracy, group=as.factor(ntree), col=as.factor(ntree))) +
  facet_grid(~nodesize, labeller = as_labeller(labels)) +
  geom_line() + geom_point() +
  theme_bw() +
  labs(col='ntree') +
  xlab('mtry') +
  ylab('Accuracy (Repeated Cross-Validation)')
dev.off()
