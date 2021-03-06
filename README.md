# HeatWave_Health_Impacts

# Overview

We are trying to establish a prediction model to estimate "Estimate" and "post Estimate", which are referred to the impact of a specific heatwave on humanbeings. We use
certain methods to quantify the impacts into specific values. The prediction model is based on the data named "data-raw" in my repository. In this dataset, there are 25 variables and two response variables (Estimate, post.esitmate). 
The definitions of these variables could be referred in https://github.com/gzmkobe/HyperHeatwavePredictiveModels/blob/master/README.md, Dr.Brooke Anderson's repository.

After a month investigation, three models(Tree model is actually one but we separate them) have been developed 

|Model   |File.Name   |
|---|---|
|GLM   |Heat_Wave_glm.Rmd   |
|Regression Tree I   |Heat_Wave_Tree1.Rmd   |
|Regresson Tree II   |Heat_Wave_Tree2.Rmd   | 
|NeuralNetwork   |Heat_Wave_NeuralNet.Rmd   | 

Since there are two responses we want to predict. There should be 6 models in total. We actually combine two responses in glm model adn neuralnet model. However, due to saving time, we separate two resonse with respect to tree model into two models.

For glm, we use this model only for intuition. Since there are 25 predictive variables, which could be regardes as large prediction sample. It is not possible to develop a model including these 25 varibales with a high R-squared value and make sure low collinearity.
We don't focus much more on it, instead we are concentrating on tree. However, glm model did give me some intuition.

For tree models, we have developed four models based on different methods: Single regression tree, bagging, random forest, and boosting.

|Tree models   	|
|--:	|
|Single Regresion Tree   	|   	
|Bagging   	|   	
|Random Forest   	|   	
|Boosting   	|   	

For neuralnetwork, we have mse for glm and neuralnet models which are pretty colse to each other. However, since it is hard to interpret neuralnet model, we currently are not consideirng it too much.

|Model   |GLM   |NeuralNet   |
|---|---|---|
|MSE   |0.349   |0.360   |

# Formatting Dataset


The original dataset has two columes with a same name, so we change one of them to the other name to let R recognize. 

```{r}
load("data-raw/ListOfAllHeatwaves.Rdata")
fix_colnames <- which(colnames(all.hws) == "mean.temp")
colnames(all.hws)[fix_colnames[2]] <- "mean.temp.1"
```

Also, we want to check if start year is an important predictive variables, so we change the datatype by cutting off the month and date, for example, we change "1987-03-15" to "1987"
```{r}
hw <- all.hws %>% 
  mutate(start.year = year(start.date))
colnames(hw) <- gsub(" ", ".", colnames(hw))
```
# GLM

## package loaded
```{r}
library(splines)
library(ggplot2)
library(xtable)
```

Define y1 as Estimate, y2 as post estimate. Using a loop to make ggplots iteratally with each variable vs. response. After reading all the plots, we find that post estimate is better than estimate in most cases.

```{r}
y1=hw$Estimate
y2=hw$post.estimate

library(splines)
library(ggplot2)
library(xtable)
xvar1 <- c("mean.temp","max.temp","min.temp","length","start.doy","days.above.80","days.above.85","days.above.90","days.above.95","days.above.99th", "days.above.99.5th","mean.temp.quantile","max.temp.quantile","min.temp.quantile","pop100","Ppoverty","Ppoverty75p","Purban","P75p","pop.density", "mean.temp.1","mean.summer.temp")
xvar2 <- c("start.month","first.in.season")

for (x in xvar1){
  print(x)
  my.x <-hw[,x]
  my.mod1 <- glm(y1~my.x)
  print(summary(my.mod1))
  #plot(my.x, y1, main = x, xlab = x)
  to_plot <- data.frame(my.x, y1, y2)
  a <- ggplot(to_plot, aes(x = my.x, y = y1)) + geom_point(alpha = 0.3) + geom_smooth() + xlab(x) + ggtitle(x)
  print(a)
  my.mod2 <- glm(y2~my.x)
  print(summary(my.mod2))
  #plot(my.x, y2, main = x, xlab = x)
  b <- ggplot(to_plot, aes(x = my.x, y = y2)) + geom_point() + geom_smooth(span = 0.3) + xlab(x) + ggtitle(x)
  print(b)
}
```

# Tree
## Package loaded

```{r}
library(lubridate)
library(dplyr)
# For plotting
library(ggplot2)
# For regression tree
# library(MASS)
library(tree)
# For bagging and random forest models
library(randomForest)
### For the boosting model
library(gbm)
### For tuning models
library(caret)

```
Define to_fit_y1 and to_fit_y2 as responses names.

```{r}
to_fit_y1 <- hw %>%
  select(-post.estimate)

to_fit_y2 <- hw %>%
  select(-Estimate)
```

Split data into trainning and validation dataset

```{r}
hw_indices <- 1:nrow(hw)
train <- sort(sample(hw_indices, length(hw_indices) / 2))
test <- hw_indices[-train]
```

Create Regression tree
```{r}
tree.hw <- tree(Estimate ~ ., data = to_fit_y1, subset = train)
summary(tree.hw)
```

Cross validation and we find the best for y1 are 4 nodes with lowest RMSE
```{r}
cv.hw <- cv.tree(tree.hw)
plot(cv.hw$size, cv.hw$dev, type='b')
```

Same process for bagging, random forest, boosting. We always use trainning data to establish models and use validation data to see the difference between true value and predicted value. Also, we alwyas want to tune parameters to find the optimal parameters with lowest RMSE. Below is the summary of MSE and importance varibale resulting from important plots.

|RMSE   |Single   |Bagging   |Random Forest   |Boosting   |
|--:|---|---|---|---|
|y1   |0.0363327   |1.093186510^{-5}   |1.093186510^{-5}   |2.694766910^{-5}   |
|y2   |9.741373510^{-4}   |9.618297810^{-7}   |9.618297810^{-7   |2.774876610^{-6}   |


This important varibales based on bagging and random forest using %incMSE standard and based on boositng using rel.inf standard

|Most three important varibales    |bagging   |Random Forest   |Boosting   |
|--:|---|---|---|
|y1   |max.temp/min.temp/start.year   |min.temp/min.temp.quantile/mean.temp.1   |start.doy/Ppoverty/min.temp.quantile   |
|y2   |mean.temp.quantile/max.temp.quantile/mean.temp   |mean.temp.quantile/mean.summer.temparature/days.above.99th   |mean.temp.quantile/pop.density/start.doy |



# Neural Netwrok
## Package required
```{r}
library(lubridate)
library(caret)
library(MASS)
library(neuralnet)
library(nnet)
```

## Model Descripions

First, we want to make sure there is no missing data, becasue package "neuralnet" can't deal with dataset with missings. After checking, we know there is no missing data and we can definetely run the package.
```{r}
apply(data,2,function(x) sum(is.na(x))
```
Then, a data preprocessing is required. It is good practice to normalize your data before training a neural network.
```{r}
maxs <- apply(data, 2, max) 
mins <- apply(data, 2, min)
scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))
train_ <- scaled[index,]
test_ <- scaled[-index,]
```

Then since we are dealing with neural network, we need to set how many layers and how many predictors in each layers. I used the ratio 25:10:7:1 (input:hidden:hidden:output)
```{r}
n <- names(train_)
f <- as.formula(paste("Estimate ~", paste(n[!n %in% "Estimate"], collapse = " + ")))  ## neuralnet doesn't recognize y~, so I define a formula 
nn <- neuralnet(f,data=train_,hidden=c(10,7),linear.output=T) 
```

Plot the prediction plot with respect glm and neural network to see if they produce a close result
```{r}
par(mfrow=c(1,3))

plot(test$Estimate,pr.nn_,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
abline(0,1,lwd=2) 
legend('bottomright',legend='NN',pch=18,col='red', bty='n')

plot(test$Estimate,pr.lm,col='blue',main='Real vs predicted lm',pch=18, cex=0.7)
abline(0,1,lwd=2)

plot(test$Estimate,pr.nn_,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
points(test$Estimate,pr.lm,col='blue',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend=c('NN','LM'),pch=18,col=c('red','blue'))
legend('bottomright',legend='LM',pch=18,col='blue', bty='n', cex=.95)
```

Finally, MSE based on GLM and NeuralNet is close to each other.

|Model   |GLM   |NeuralNet   |
|---|---|---|
|MSE   |0.349   |0.360   |
