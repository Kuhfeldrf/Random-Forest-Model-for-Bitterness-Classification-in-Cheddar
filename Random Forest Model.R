install.packages("Matrix")
install.packages("plyr")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
install.packages("caret")
BiocManager::install("biomformat")
install.packages("randomForest")
install.packages("glmnet")
library(biomformat)
library(plyr)
library(vegan)
library(ape)
library(glmnet)
library(caret)
library(phyloseq)
library(randomForest)
library(ggplot2)
library(e1071)
library(Matrix)


## Load Dataset
OTUR<- read.csv ("Normalized Peptide Data Abundance OUT1_Sample1 Labels.csv", header = FALSE,)
OTUR1 <-sapply(OTUR,as.numeric)
OTU <- (t(OTUR1))
OTU[is.na(OTU)]<-1000
Sample_names <- read.csv("Sample_names.csv", header = TRUE) #Sample_names is file with peptides names in the first column 
taxmat <- read.csv ("Cat_phyloseq-phylo.csv", header=FALSE,)
#updated dimensions of table to 333 rows with new data
otumat <- OTU[2:15,2:1131]
rownames(otumat) <- paste0("Sample", 1:nrow(otumat))
colnames(otumat) <- Sample_names[1:1130,1]
otumat[is.na(otumat)]<-1000
str(otumat)
otumat <- as.matrix(otumat)
otumat1 <-t(otumat)
class(otumat1)
View(otumat1)
#updated dimensions of table to 333 rows with new data
taxmat<-taxmat[2:1131,2:8]
rownames(taxmat) <- Sample_names[1:1130,1]
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat=as.matrix(taxmat)
class(taxmat)
#View(taxmat)


OTU2 <- otu_table(otumat1, taxa_are_rows = TRUE)
OTU3 <- t(OTU2)
str(OTU3)

TAX <- tax_table(taxmat)
str(TAX)

physeq <- phyloseq(OTU3, TAX)
plot_bar(physeq)

random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)


sampledata <- sample_data(data.frame(
  Class = sample(LETTERS[1:2], size=nsamples(physeq), replace=TRUE),
  Depth = sample(50:1000, size=nsamples(physeq), replace=TRUE),
  row.names=sample_names(physeq),
  stringsAsFactors=FALSE))

#Updated catagorical groups 7.31.22 -RK NB = non-bitter B = bitter
sampledata$Class <- c("NB","NB","NB","NB","NB","NB","NB","B","B","B","B","B","B","B")
sampledata$TL<- c("T_0.2a","T_0.2b","T_0.2c","T_0.2d","L_2.9","L_3.3","L_6.0","","","","","","","")
sampledata$ME<-c("","","","","","","","M_5.6","E_5.7","M_6.2","E_7.2","M_7.3","E_8.7a","E_8.7b")
ps <- merge_phyloseq(physeq,sampledata , random_tree)
View(ps)

### Normally, we would check: 1. Normalize counts and 2. Remove outliers, but we've already done this in previous practicals. Let's just take a look at a PCoA of all the samples colored by hive role

ord <- ordinate(ps, method = "PCoA", distance = 'bray') #We usually use bray curtis distance for PCoA of microbiome data, rather than euclidean. For machine learning though, because we will be transforming the variable counts so that our model weights are more interpretable, bray curtis distance will end up looking like nonsense. Therefore, for this application, we'll use euclidean.
plot_ordination(ps, ord, 'samples', color = 'Class')


### As mentioned in lecture, we'll be building 3 models, each one will return the likelihood that the sample in question is a (F=Forager, W=Worker, H=Nurse) bee. Let's make a binary column in our sample_data for each hive role. Lastly, plot an ordination colored by one of those columns (this should look very familiar from lecture)
#add a column to the sample_data that is true/false this bitter or non-bitter
ps@sam_data$Class_NB = ps@sam_data$Class == "NB"
ps@sam_data$Class_B = ps@sam_data$Class == "B"
### generates PCoA with-out prediction / shape
plot_ordination(ps, ord, 'samples', color = 'Class')+ 
  geom_text(
    mapping = aes(label = (sampledata$TL)), size = 10,hjust =-.1,vjust = -.525, Font="BOLD")+
  geom_text(
    mapping = aes(label = (sampledata$ME)), size = 10,hjust =1.1,vjust = 1.2, Font="BOLD")+
  labs(color = "Categorical grouping") +
  scale_color_manual(labels =c("Bitter","Non-bitter"), values=c("#FFA500","#0096FF")) +
  theme(
    legend.justification=c(1,-0.02),
    legend.position = c(.999,0.84),
    legend.text =  element_text(size=18),
    legend.title =  element_text(size=18),
    axis.title.x = element_text(size=20),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=18,),
    axis.text.y = element_text(size=18,))+
  geom_point(size=5)+
  geom_hline(yintercept = 0,linetype='dotted')+
  geom_vline(xintercept = 0,linetype='dotted')+
  guides(color = guide_legend(override.aes = list(size=5)))+
  guides(sahpe = guide_legend(override.aes = list(size=5)))

### Center each feature around its mean and divide by its standard deviation. This will not change the predictions of our model, but it will allow us to draw conclusions about the relative importance of features from the weights the model learns. If you do not do this step, the only conclusions you can draw from the weights of the model are whether a feature has positive or negative association with the outcome metric, but you can't say anything about the magnitude

otu_table(ps) <- otu_table(apply(ps@otu_table, 2, function(x) return( (x - mean(x)) / sd(x))), taxa_are_rows = FALSE)


## Random Forest

### Set aside testing data
set.seed(1)
index_train <- createDataPartition(ps@sam_data$Class, p = 0.7)[[1]]
x_train <- ps@otu_table[index_train, ]
x_test <- ps@otu_table[-index_train, ]

#split the phyloseq objects into training and testing to make our lives easier later on
ps_train <- phyloseq(otu_table(ps@otu_table[index_train, ], taxa_are_rows = F), ps@sam_data[index_train, ])
ps_test <- phyloseq(otu_table(ps@otu_table[-index_train, ], taxa_are_rows = F), ps@sam_data[-index_train, ])

### Train model to classify Forager bees and check performance on testing dataset

#Notes:
# this function is from the glmnet package
# family = binomial means we'll be using logistic regression 
# alpha = 1 means lasso (fewer variables); alpha = 0 means ridge (lower coefficients for each variable)
# we need to pass our input data (x) in as a sample by feature matrix
# we also need to pass in the actual answers (y)
# the function will return a trained model

model_ridge <- glmnet(x = as.matrix(x_train), y = ps_train@sam_data$Class_B, family = 'binomial', alpha = 0) 


###Check the accuracy on the testing set that the model hasn’t seen before with no regularization
#Notes: 
# s is the value of lambda (regularization strength)
pred_B <- predict(model_ridge, newx = data.matrix(x_test), s = 0) 
classifications <- pred_B > 0 # pred_bee_f is the likelihood that a sample is a forager. We change this to a classification by checking if the likelihood is greater than 0

print(paste("Accuracy on test set: ", sum(classifications == ps_test@sam_data$Class_B)*100 / nsamples(ps_test), "%")) #check the accuracy on the testing set
#
##[1] "Accuracy on test set:  75 %"
#

###Visualize the prediction rule on the training data with no regularization
#We'll be doing this a lot, so we're going to write a function that we can call to be more efficient in the future:
plot_classification_rule <- function(ps_use, model, s, color = 'Class_B', title="") {
  pred_BorNB <- predict(model, newx = as.matrix(ps_use@otu_table), s = s)
  classifications <- pred_BorNB > 0
  ps_use@sam_data$predicted <- classifications
  ord <- ordinate(ps_use, method = "PCoA", distance = 'euclidean')
  plot_ordination(ps_use, ord, 'samples', color = color, shape = 'predicted', title = title) +
    geom_point(size = 4) + #geom_point is added here to modify the size 
    labs(shape = "Classifier Is it Bitter", color = "Bitter")
}


#get predictions for entire dataset for visualization
pred_bitter <- predict(model_ridge, newx = x_test, s = 0) #again, s=0 so no regularization
classifications <- pred_bitter > 0

plot_classification_rule(ps_train, model_ridge, s = 0, title = "Predictions with no regularization")

###By changing the regularization strength (s), we’ll change the classification rule. Plug in a few different values for s (regularization strength) and rerun this chunk. Watch what happens to the predictions. What happens (Hint: keep your values between 0 and 14 for this dataset)
#get predictions for entire dataset for visualization

##### CHANGE THIS ####
s = 20
######################
#plug in values for s. what happens to your performance? How does that fit in with what you know about regularization? As you increase regularization, you should see the predictions migrate toward the majority category, because your model is being forced to 'learn less' from the particulars of the data

plot_classification_rule(ps_train, model_ridge, s = s, title = paste("Predictions with regularization lambda = ", s))


###to what value of s should we use? Let’s use cross-validation over the training data to pick the regularization strength without ever using the test data.
set.seed(6) #because we're using cross validation, we could get slightly different results due to splitting. Let's freeze the results for consistency
model_ridge <- cv.glmnet(as.matrix(x_train), y = ps_train@sam_data$Class_B, family = 'binomial', alpha = 0, nfolds = 4, lambda = seq(0, 10, by = .1)) # alpha = 0 means ridge regression; lambda parameter gives the function a range of possible regularization strengths to try

plot(model_ridge)
#Call:  cv.glmnet(x = as.matrix(x_train), y = ps_train@sam_data$Class_B,      lambda = seq(0, 10, by = 0.1), nfolds = 4, family = "binomial",      alpha = 0) 
#Measure: Binomial Deviance 
#Lambda Index Measure     SE Nonzero
#min    2.6    75  0.7831 0.3979    1130
#1se   10.0     1  0.8236 0.3321    1130
###Regularize with lambda.min and lambda.1se, explained in code comment above. You should chose which one to use by looking only at the training data, to ensure accurate error estimation on testing data
plot_classification_rule(ps_train, model_ridge, s = 'lambda.min', title = paste("Predictions with regularization, Lambda min value:  2.6"))
plot_classification_rule(ps_train, model_ridge, s = 'lambda.1se', title = paste("Predictions with regularization lambda, Lambda 1se value:  10"))
print(paste("Lambda min value: ", model_ridge$lambda.min))
print(paste("Lambda 1se value: ", model_ridge$lambda.1se))
#> print(paste("Lambda min value: ", model_ridge$lambda.min))
#[1] "Lambda min value:  2.6"
#> print(paste("Lambda 1se value: ", model_ridge$lambda.1se))
#[1] "Lambda 1se value:  10"

###Question: If you recall, there were two flavors of regularization in lecture. Which one did we use above? What happens when we use the other method?
set.seed(1) #freeze results for consistency
model_lasso <- cv.glmnet(as.matrix(x_train), y = ps_train@sam_data$Class_B, family = 'binomial', alpha = 1, nfolds = 3, lambda = seq(0, 10, by = 0.1), title = paste("lasso regularization, alpha = 1 ")) # alpha = 1 means lasso regularization (force fewer variables in model)
plot(model_lasso)

###Using lasso regularization, it looks like we could pick either strength lambda.min or lambda.1se to use for the final model, as they’re the same. We’ll use lambda.min for consistency with the previous step here. Check out the model performance on our test set

#Note: we don't have to retrain our model, because glmnet stores all models using all the values of lambda that it tried. All we have to do is apply the model to make predictions on the full dataset, rather than just the training set

plot_classification_rule(ps_test, model_lasso, s = "lambda.min", title = "lambda.min")

plot_classification_rule(ps_test, model_lasso, s = "lambda.1se", title = "lambda.1se")

###Now, what taxa sequences were positively or negatively correlated with a forager hive role? We can look at the weights, or coefficients, of the models to find out. Compare weights from lasso and ridge regularized models.
print(coefficients(model_lasso))
cml<-as.data.frame(as.matrix(coefficients(model_lasso)))
write.csv(cml, file="coefficients_model_lasso.csv")
### Find the optimal tree depth (mtry) using cross validation with the caret package
#trainControl is basically a controller for the cross-validation process. It will get passed to the train command. The package we used above, glmnet, does cross validation for you. Because glmnet doesn't implement random forests, we'll be using the caret package to handle our cross-validation

#The train function in caret will want the data as a dataframe where one column is singled out at the answers. Our answers will be the "hive_role" column which we're creating here
set.seed(1)
data_train = data.frame(x_train)
data_train$R = ps_train@sam_data$Class
control <- trainControl(method='repeatedcv', 
                        number=4, 
                        repeats=3,
                        allowParallel = F)

tunegrid <- expand.grid(.mtry=c(1:20)) #mtry is the depth of each decision tree. We'll be trying out models where each tree is 3 to 20 splits deep
rf <- train(R ~., 
            data= data_train, 
            method='rf', 
            metric='Accuracy', 
            tuneGrid=tunegrid, 
            trControl=control)
print(rf)

## Accuracy is measured using the test set assigned at each fold during cross validation. These small sub-test sets are called validation sets, for the sake of unique vocabulary. Similar to how we used cv.glmnet above to find the optimal strength of regularization, checking our error on these validation sets during cross validation allows us to pick a tree depth that will likely work best on outside data. Remember that the cross validation is happening on the training set, so we still have the actual test set to check performance on. 

# In general, the deeper the trees, the more the model overfits the training data, resulting in lower accuracy on the validation set. The "Accuracy" reported here for each value of mtry is the accuracy on the 'out of bag' samples, or the temporary validation sets created by the random forest algorithm

#In this case, you may see that regardless of tree depth, we predict the validation set perfectly (also said as "the out of bag error is 0"). In general, there should be a 'sweet spot' where a deeper tree depth overfits and a shallowed tree depth does not learn enough. 

### Let's try the model performance on the held out test set, using the value for mtry (tree depth) chosen during cross validation
mtry_best = as.numeric(rf$bestTune)
model = randomForest(x_train, y = as.factor(ps_train@sam_data$Class), mtry = mtry_best)
model

#Call:
#randomForest(x = x_train, y = as.factor(ps_train@sam_data$Class),      mtry = mtry_best) 
#               Type of random forest: classification
#                     Number of trees: 500
#No. of variables tried at each split: 13
#
#        OOB estimate of  error rate: 20%
#Confusion matrix:
#   B NB class.error
#B  3  2         0.4
#NB 0  5         0.0

#Performance on test set
preds = predict(model, x_test)
print(paste("Accuracy: ", sum(preds == as.factor(ps_test@sam_data$Class)) / nsamples(ps_test)))
#Visualize on test dataset
ord <- ordinate(ps_test, method = "PCoA", distance = 'euclidean')
ps_test@sam_data$rf_predictions = predict(model, ps_test@otu_table)
plot_ordination(ps_test, ord, 'samples', color = 'Class', shape = 'rf_predictions') + geom_point(size = 4)

#One of the reasons random forests are nice is because we don't have to deal with multiple models to classify multiple output types.

### Now we'd like to know which taxa were most important in training the full model (all data). Notice that every time you train a model and take a look at the importance of variables, you get a different graph for the importance of each variable. Run this command multiple times to see this.
model = randomForest(ps@otu_table, y = as.factor(ps@sam_data$Class), mtry = mtry_best)
varImpPlot(model, type = 2)
# Run this chunk multiple times to see how the variable importance changes


### Question: How can we tell which variables are really important?

### A common technique with random forests and other models that rely on randomness is to simply do the training process a number of times and average the results. Here, we'll do it 50 times
imp_list <- list()
for(i in 1:250){
  model = randomForest(ps@otu_table, y = as.factor(ps@sam_data$Class), mtry = mtry_best)
  imp_list[i] <- varImp(model)
}

#short list of bitter peptide canadiates
top_20<- read.csv("Top_20_List.csv", header = TRUE,)
grp_top_20=c(unlist(top_20$Sequence))

#converts the protein names to greek symbols
Beta<-"\U03B2"
Alphas1<-paste0("\U03B1","s1")
Alphas2<-paste0("\U03B1","s2")
Kappa<-paste("\U03BA")
top_20$Protein[top_20$Protein == "Beta"] <- Beta
top_20$Protein[top_20$Protein == "Alphas1"] <- Alphas1
top_20$Protein[top_20$Protein == "Alphas2"] <- Alphas2
top_20$Protein[top_20$Protein == "Kappa"] <- Kappa


#organizes data to save into bar plot
imp_df <- do.call(rbind.data.frame, imp_list)
colnames(imp_df) <- colnames(x_train)
Mean_var_imp<-colMeans(imp_df)
rf250 <- data.frame(Mean_var_imp)
rf250$Sequence<-colnames(x_train)
rf250_20<-data.frame(subset(rf250,rf250$Sequence%in%grp_top_20))
rf250_20m<-merge(rf250_20, top_20)
rf250_20m$comibned<-paste(rf250_20m$Protein,rf250_20m$Positions.in.Proteins)
rownames(rf250_20m) <- rf250_20m$comibned
rf250_20m$Positions.in.Proteins<-NULL
rf250_20m$Sequence<-NULL
rf250_20m$Protein<-NULL
rf250_20m$comibned<-NULL
rf250_20m<-rf250_20m[order(rf250_20m[1]),,drop=FALSE]
rf250_20m<-(t(rf250_20m))
#The Par(mar function increases the size of the margin. for the labels on the boxplot. this needs to be run in conjunction with the barplot code
par(mar = c(4,14,4,4)) 
barplot(rf250_20m, horiz =T, las=2)
title(ylab = "Top 20 bitter peptide candidates", line = 10) 
title(xlab= "Mean variable importance (% of total variable importance)")
#export list to csv 
colMeans(imp_df)
rf250x<-sort(colMeans(imp_df),decreasing = TRUE)
write.csv(rf250x, "Mean Varriable Importance Export")
#These importance scores should not change much, because they are averages.

#one weakness of random forests is that while they return variable importance, it is difficult (but still possible) to get the directionality of each variable (positively or negatively associated with output variable). This is because random forests allow for large amounts of dependence. A low value for a taxa2 might mean forager when paired with a high value for taxa 18, but worker when paired with a high value for taxa 21. It's difficult to pull apart those inconsistencies in the model. 