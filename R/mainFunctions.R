#'Cross validation results for a model
#' 
#'@param clinical_x clinical variables that will always be included in the model
#'@param genomic_x genomic variables that will be penalized if a penalized model is used
#'@param y response variables
#'@param data dataframe if clinical formula is used
#'@param clinical_formula formula for clinical variables
#'@param family gaussian, binomial or poisson
#'@param folds predefined partions for cross validation
#'@param k number of cross validation folds. A value of k=n is leave one out cross validation. 
#'@param fit_method glm or glmnet used to fit the model 
#'@param method_name tracking variable to include in return dataframe
#'@param n.cores Number of cores to be used
#'@param ... additional commmands to glm or cv.glmnet
#'@return returns a dataframe of predicted values and observed values. In addition, method_name is recorded if that variable is defined.
#'@examples
#'x <- matrix(rnorm(800),ncol=8)
#'y <- runif(100) < exp(1 + x[,1] + x[,5])/(1+exp(1 + x[,1] + x[,5]))
#'cv_results <- cv(x,y=y,method_name="without_formula")
#'combined_data <- data.frame(y=y,x1=x[,1],x5=x[,5])
#'gx <- x[,c(2,3,4,6,7,8)]
#'cvf <- cv(genomic_x=gx,clinical_formula=y~x1+x5,data=combined_data,method_name="with_form")
#'@author Ben Sherwood <ben.sherwood@@ku.edu>
#'@export
cv <- function(clinical_x=NULL, genomic_x=NULL, y=NULL, data=NULL, clinical_formula=NULL, family="binomial",folds=NULL, k=10, fit_method="glm", method_name=NULL,n.cores=1,...){
	if(is.null(y) & is.null(clinical_formula)){
		stop("y needs to be defined in clinical formula or y")
	}
	if(is.null(clinical_x) & (is.null(data) | is.null(clinical_formula)) & is.null(genomic_x)){
		stop("No predictors defined")
	}
	if(is.null(clinical_formula)==FALSE & is.null(data)){
		stop("clinical_formula defined with no data")
	}
	if(is.null(clinical_formula)==FALSE & is.null(clinical_x)==FALSE){
		stop("can only use clinical_formula or clincal_x, both variables cannot be defined")
	}
	if(is.null(clinical_formula)==FALSE & is.null(clinical_x)==FALSE){
		stop("can only use clinical_formula or y, both variables cannot be defined")
	}
	if(is.null(clinical_formula)==FALSE){
		#clinical_formula <- update.formula(clinical_formula, y ~ .)
		clinical_model <- glm(clinical_formula,family=family,data=data)
		clinical_x <- model.matrix(clinical_model)[,-1,drop=FALSE]
		y <- clinical_model$y
		non_pen_vars <- 1:dim(clinical_x)[2]
	} else{
		non_pen_vars <- NULL
	}
	n <- length(y)
	x <- cbind(clinical_x, genomic_x)
	if(is.null(folds)){
		folds <- randomly_assign(n,k)
	}
	if(n.cores == 1){
		results <- data.frame(do.call(rbind,lapply(unique(folds),fit_pred_fold,x,y,folds,fit_method,family,non_pen_vars=non_pen_vars,...)))
	} else{
		results <- data.frame(do.call(rbind,mclapply(unique(folds),fit_pred_fold,x,y,folds,fit_method,family,non_pen_vars=non_pen_vars,...)))
	}
	if(is.null(method_name)){
		colnames(results) <- c("prediction","response")
	} else{
		results <- cbind(results, method_name)
		colnames(results) <- c("prediction","response","method")
		results$method <- paste(results$method)
	}
	results
}

#'Assigns n samples into k groups
#' 
#'@param n sample size
#'@param k number of groups
#'@return returns a vector of length n with a random assignment of entries from 1 to k
#'@examples
#'n <- 100
#'folds_10 <- randomly_assign(n,10)
#'folds_5 <- randomly_assign(n,5)
#'@author Ben Sherwood <ben.sherwood@@ku.edu>
#'@export
randomly_assign <- function(n,k){
#randomly assign n samples into k groups
   small_set <- floor(n/k)
   group_assign <- NULL
  if(n %% k == 0){
     group_assign <-  rep(seq(1,k),n/k)
   } else{
     remainder <- n %% k
     for(i in 1:remainder){
        group_assign <- c(group_assign, rep(i,small_set+1))
     }
     group_assign <- c(group_assign, rep(seq((i+1),k),small_set))
   }
   sample(group_assign)
}

#'Cross validation on fold i
#' 
#'@param i target partition
#'@param x matrix of predictors
#'@param y vector of responses
#'@param folds defines how data is seperated into folds for cross validation
#'@param fit_method model being used to fit the data
#'@param family family used to fit the data
#'@param non_pen_vars index of variables that will not be penalized if glmnet is used
#'@param ... additional commmands to glm or cv.glmnet
#'@return returns predictions for partition i
#'@examples
#'folds_10 <- randomly_assign(100,10)
#'x <- matrix(rnorm(800),ncol=8)
#'y <- runif(100) < exp(1 + x[,1] + x[,5])/(1+exp(1 + x[,1] + x[,5]))
#'fold_1_results <- fit_pred_fold(1,x,y,folds_10,"glm","binomial")
#'fold_2_results <- fit_pred_fold(2,x,y,folds_10,"glm","binomial") 
#'@author Ben Sherwood <ben.sherwood@@ku.edu>
#'@export
fit_pred_fold <- function(i,x,y,folds,fit_method,family,non_pen_vars=NULL,...){
	if(is.null(dim(x))){
		train_x <- x[folds!=i]
		test_x  <- x[folds==i]
		p <- 1
	} else{
		train_x <- x[folds!=i,]
		test_x  <- x[folds==i,]
		p <- dim(x)[2]
	}
	train_y <- y[folds!=i]
	test_y  <- y[folds==i]
	if(fit_method == "glm"){
		train_x <- data.frame(train_x)
		colnames(train_x) <- paste0("x",1:p)
		train_x <- cbind(train_x,train_y)
		train_model <- glm(train_y ~ ., family, train_x, ...)
		test_x <- data.frame(test_x)
		colnames(test_x) <- colnames(train_x)[1:p]
		preds <- predict(train_model,test_x,type="response")
	}
	if(fit_method =="glmnet"){
		if(p == 1){
			stop("cannot use glmnet with p = 1")
		}
		pen_fact <- rep(1, p)
		pen_fact[non_pen_vars] <- 0
		train_model <- cv.glmnet(train_x, train_y, family=family, penalty.factor=pen_fact,...)
		preds <- predict(train_model,test_x,type="response")
	}
	cbind(preds,test_y)
}

#'Create ROC plot from cross validation results 
#' 
#'@param plot_data dataframe with columns: response, prediction and method where method should be 3rd column
#'@param ... additional commmands plot.roc such as main 
#'@return returns ROC plot
#'@examples
#'x <- matrix(rnorm(800),ncol=8)
#'y <- runif(100) < exp(1 + x[,1] + x[,5])/(1+exp(1 + x[,1] + x[,5]))
#'cv_results <- cv(x,y=y,method_name="without_formula")
#'combined_data <- data.frame(y=y,x1=x[,1],x5=x[,5])
#'gx <- x[,c(2,3,4,6,7,8)]
#'cvf <- cv(genomic_x=gx,clinical_formula=y~x1+x5,data=combined_data,method_name="with_form")
#'total_results <- rbind(cv_results,cvf)
#'rocplot(total_results,main="rocplot test")
#'@author Ben Sherwood <ben.sherwood@@ku.edu>
#'@export
rocplot <- function(plot_data,...){
	method_count <- 1
	auc_vals <- NULL
	method_vals <- unique(plot_data$method)
	for(method_name in method_vals){
		sub_data <- plot_data[plot_data[,3]==method_name,]#subset(plot_data, method == method_name)
		#print(paste("method count", method_count))
		if(method_name != ""){
			if(method_count==1){
			  plot.roc(sub_data$response,sub_data$prediction, lty=method_count, col=method_count, direction="<", ...)
			  auc_vals <- c(auc_vals, pROC::auc(sub_data$response, sub_data$prediction,direction="<"))
			} else{
			  lines.roc(sub_data$response,sub_data$prediction, lty=method_count, col=method_count, direction="<")
			  auc_vals <- c(auc_vals, pROC::auc(sub_data$response, sub_data$prediction, direction="<"))
			}
		}
		method_count <- method_count+1
	}
	auc_vals <- round(auc_vals,2)
	legend("bottomright", legend=paste( unique(plot_data$method), auc_vals, sep = " AUC: "), lty=1:length(unique(plot_data$method)), col=1:length(unique(plot_data$method)), lwd=2,bty = "n")
}
