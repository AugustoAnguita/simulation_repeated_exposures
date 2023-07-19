#######################################################################################################################################################################################
##
## 								List and description of all functions inluded in this document
##
## #######################################################################################################################################################################################
##
##
## check.ICC1 : Function for the datasets simulation (exposures data). This function checks if the combination of Sigma and ICC will produce a valid covariance matrix. This function checks that correlation between errors does not have to be greater than 1
##
## rebuild.ICC : Function for the datasets simulation (exposures data).
## 
## simul : Function for the datasets simulation (exposures data). This is the main function for the simulation of the datasets of baseline exposures.
##
## simX :  Function for the datasets simulation (exposures data). This function simulates repeated exposure datasets according to the ICC previously built.
##
## simtestY1 :  Function for the datasets simulation (outcome data). This function generates outcome data for the scenario 1 described in the manuscritp (All time points for true exposures are associated with the outcome).
## 
## simtestY2 :  Function for the dataset simulation (outcome). This function generates outcome data for the scenario 1 described in the manuscritp (only a single time point for each true exposures is associated with the outcome).
##
## EWAS : Implemmentation of ExWAS and ExWAS.MLR methods. Please, see manuscript for mor details.
##
## testReg : Implemmentation of ExWAS and ExWAS.MLR methods. Please, see manuscript for mor details.
## 
## applyEwas : Implemmentation of ExWAS and ExWAS.MLR methods. Please, see manuscript for mor details.
##
## EvaluateModel.cv.glmnet : Implemmentation of ElasticNet method. Please, see manuscript for mor details.
##
## applyEnet : Implemmentation of ElasticNet method. Please, see manuscript for mor details.
## 
## EvaluateModel.cv.spls : Implemmentation of sPLS method. Please, see manuscript for mor details.
##
## applySPLS : Implemmentation of sPLS method. Please, see manuscript for mor details.
##
## applysNPLS : Implemmentation of sNPLS method. Please, see manuscript for mor details.
## 
## DSAreg : Implemmentation of DSA method. Please, see manuscript for mor details.
##
## applyDSA : Implemmentation of DSA method. Please, see manuscript for mor details.
##
## applyDLNMselect : Implemmentation of DLNM with penalization and feature selection. Please, see manuscript for mor details.
##
##
## Last update: 09-05-2023
## Authors: Charline Warembourg <charline.warembourg@inserm.fr> Augusto Anguita <augusto.anguita@isglobal.org>  Xavier Basagaña <xavier.basagana@isglobal.org> 
################################################################################################################################################################################################################################################





########################################################
## Load required libraries			        ##
########################################################


list.of.packages <- c(
  "Matrix","MASS","parallel","glmnet","spls","MXM","DSA","dlnm","splines","mgcv","sNPLS","plyr","stringr","foreach","doParallel","ranger","palmerpenguins","tidyverse","kableExtra","future","dplyr"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}
for(package.i in list.of.packages){ #loading packages
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}

options(max.print=100000)

print("Step 1 - Required R libraries successfully loaded")



##### --------------------------------------------------------- 1. DATA SIMULATION --------------------------------------------------------- ##### 



##########################################################
## Functions to simulate datasets of baseline exposures ##
##########################################################


check.ICC1 = function(nt,index.exposures,ICC,rho) {
	
	Sigma = Sigma[index.exposures,index.exposures]
	l = length(index.exposures)
  	rho.min = matrix(NA,nrow=l,ncol=l)
  	rows = vector()
  	cols = vector()
  	ss <- icc1 <- icc2 <- vector()
  	k = 1
	  for (i in 1:(l-1)){
	    for (j in (i+1):l) {
	      rho.min[i,j] = max( abs(Sigma[i,j])-ICC[i] , abs(Sigma[i,j])-ICC[j]) / (sqrt(1-ICC[i])*sqrt(1-ICC[j]))
	      if (rho.min[i,j]>rho) {
		rows[k] = i
		cols[k] = j
		ss[k] = Sigma[i,j]
		icc1[k] = ICC[i]
		icc2[k] = ICC[j]
		k = k +1
	      }
	    }
	  }
	problem.in = (length(rows)>0)
	  if (problem.in) {
	    res = data.frame(rows,cols,ss,icc1,icc2)
	    # increase lowest ICC in all problematic pairs
	    newICC = ICC
	    for (i in 1:nrow(res)) {
	      if (res$icc1[i]<res$icc2[i]) {
	       if (res$icc1[i]==.1) newICC[rows] = 0.5 else if (res$icc1[i]==.5) newICC[rows] = 0.9  else if (res$icc1[i]==.9) newICC[rows] = 0.99 
	      } else {
	       if (res$icc2[i]==.1) newICC[cols] = 0.5 else if (res$icc2[i]==.5) newICC[cols] = 0.9  else if (res$icc2[i]==.9) newICC[cols] = 0.99
	      }
	    }
    
	# Now check if with the corrections, it is ok
    	rho.min = matrix(NA,nrow=l,ncol=l)
    	rows = vector()
	cols = vector()
    	ss <- icc1 <- icc2 <- vector()
    	k = 1
	    for (i in 1:(l-1)){
	      for (j in (i+1):l) {
		rho.min[i,j] = max( abs(Sigma[i,j])-newICC[i] , abs(Sigma[i,j])-newICC[j]) / (sqrt(1-newICC[i])*sqrt(1-newICC[j]))
		if (rho.min[i,j]>rho) {
		  rows[k] = i
		  cols[k] = j
		  ss[k] = Sigma[i,j]
		  icc1[k] = newICC[i]
		  icc2[k] = newICC[j]
		  k = k +1
		}
	      }
	    }
    	problem.out = (length(rows)>0)
    
  	} else {
    	newICC = ICC
    	problem.in = F
    	problem.out = F
  		}
	return(list(newICC=newICC, problem.in=problem.in, problem.out=problem.out))
}


rebuild.ICC = function(nt,index.exposures,ICC,rho) {
	  problem.out1 = T
	  n.attempts = 1
	  newICC = ICC
	  while (n.attempts<10 & problem.out1==T) {
	    check = check.ICC1(nt=nt,index.exposures=index.exposures,ICC=ICC,rho)
	    ICC = check$newICC
	    problem.out1 = check$problem.out
	    n.attempts = n.attempts + 1
	  }
	  return(list(newICC=ICC, problem.out= problem.out1))
}


simul =  function(N,nt,index.exposures,ICC) {
	  # N: sample size
	  # nt: number of time points
	  # index.exposures: the index for the exposures in Sigma (correlation matrix)
	  #        e.g. if we want to generate 122 exposures, index.exposures = 1:122
	  #        e.g. if we want to generate 50 exposures highly correlated, index.exposures = 1:50
	  #        e.g. if we want to generate 50 exposures with low correlation, index.exposures = 73:122
	  #        e.g. if we want to generate 100 exposures at random, index.exposures = sample(1:122,100,replace=F)
	  # ICC: vector of same size than index.exposures that indicates the intraclass correlation of the exposures.
	  #       This is ICC[var1] in my notation
	  #       High ICC1 means high repatebility of the variable over time.

	  # Generate from a multivariate normal distribution with a correlation 
	  #    matrix calculated from Sigma and ICC

	  # Without loss of generality, we assume all variables have mean 0
	  # We assume the same variance = 1 for all exposures - this could be changed but I'm
	  #    not sure it's worth it
	  
	  Sigma = Sigma[index.exposures,index.exposures]
	  
	  l = length(index.exposures)

	  rho.min = matrix(NA,nrow=l,ncol=l)
	  for (i in 1:(l-1)){
	    for (j in (i+1):l) {
	      rho.min[i,j] <- rho.min[j,i] <- (abs(Sigma[i,j]) - min(ICC[i],ICC[j])) / (sqrt(1-ICC[i])*sqrt(1-ICC[j]))
	    }
	  }
	  diag(rho.min) = 1
	  
	  rho.max = matrix(NA,nrow=l,ncol=l)
	  for (i in 1:(l-1)){
	    for (j in (i+1):l) {
	      rho.max[i,j] <- rho.max[j,i] <- (abs(Sigma[i,j]) ) / (sqrt(1-ICC[i])*sqrt(1-ICC[j]))
	    }
	  }
	  diag(rho.max) = 1
	  
	  rho.e = matrix(NA,nrow=l,ncol=l)
	  for (i in 1:(l-1)){
	    for (j in (i+1):l) {
	      rho.e[i,j] <- rho.e[j,i] <- runif(1,max(0,rho.min[i,j]),min(1,rho.max[i,j]))
	    }
	  }
	  diag(rho.e) = 1
	  rho.e = rho.e*sign(Sigma)

	  Psi = Sigma - rho.e * (sqrt(1-ICC) %*% t(sqrt(1-ICC))) 
	  #diag(Psi) = ICC
	  
	  blocks = list()
	  blocks[[1]] = Sigma
	  for (i in 2:nt) {
	    blocks[[i]] = Psi
	  }
	  
	  m.str <- toeplitz(1:nt)
	  
	  res <- lapply(1:nt,function(k) {
	    res <- matrix(0,ncol=ncol(m.str),nrow=nrow(m.str))
	    res[m.str == k] <- 1
	    res %x% blocks[[k]]
	  })
	  
	  S = Reduce("+",res)
	  
	  if (!is.positive.definite(S)) S = as.matrix(nearPD(S,corr=T)$mat)
	  
	  X = data.frame(mvrnorm(n=N,mu=rep(0,nt*l),Sigma=S))
	  names(X) = paste0(paste0("X",1:l),rep(paste0(".",1:nt),each=l))

	  return(list(X=X,S=S))
}





######################################################################################################
# Functions to simulate longitudinal datasets with 5-times repeated exposures and the single outcome #
######################################################################################################


simX <- function(Sigma=Sigma,N = 1200){
   # Here I'm taking a random sample of 100 exposures
   index.exposures = sample(1:122,100,replace=F)
   nvars = length(index.exposures)

   # Correlation between repeated measures can be high (0.9), medium (0.5), or low (0.1)
   #    Each variable will have one of these 3 correlations
   # Here I choose at random for each variable one of the 3 correlations 
   ICC = sample(x=c(0.1,0.5,0.9),size=nvars,prob=c(1/3,1/3,1/3),replace=T)  #ICC was called rho.vec in previous code

   # number of time points
   nt = 5

   check = rebuild.ICC(nt=nt,index.exposures=index.exposures,ICC=ICC,rho=1) 
   # The following needs to be FALSE, we have to check it every time
   message(check$problem.out)
   # number of changes done in ICC
   sum(ICC!=check$newICC)
   # we will have fewer cases of low ICCs, because they are not compatible with some correlations
   table(ICC)
   table(check$newICC)

   ICC = check$newICC

   simdata = simul(N,nt,index.exposures,ICC) 

   X = simdata[[1]]   # X was called Xbas in the previous code
   S = simdata[[2]]   # S was called U in the previous code

   # in case they are needed, names of original variables used to simulate
   names.vec = rownames(Sigma)[index.exposures]

  
  X.list<-list(X=X,rho.vec=ICC)
  
  return(X.list)
}


simtestY1 <- function(X=data.X$X,n.true.exposures=3,nvars=100,N = 1200,r=5){  #data.X = exposure matrix, n.true.exposures= number of true predictor of Y, nvars = number of exposures (not repeated), N=number of individuals, r=number of repeated measures
  
  # selected at random
  true.exposures = sample(1:nvars,n.true.exposures,replace=F)
  
  # suppose the coefficients of the variables decrease with time
  true.coefs = rep(seq(2,1,by=-(1/(r-1))),n.true.exposures)
  
  Xtrue = as.matrix(cbind(X[,paste0(paste0("X",rep(true.exposures,each=r),"."),1:r)]))
  Xtrue.names = paste0(paste0("X",rep(true.exposures,each=r),"."),1:r)
  
  linpred = Xtrue %*% matrix(true.coefs,ncol=1)
  # R2 of the model = 10%
  resid.var = var(linpred)*9
  
  Y = linpred + rnorm(N,0,sd=sqrt(resid.var))
  
  formul = paste("Y ~", paste(Xtrue.names,sep=" ", collapse="+"))
  
  mod = lm(formul,data=cbind(Y,X))
  summary(mod)
  
  resu.sim1<-list(resp=Y,true.pred=Xtrue.names,sigmaE=resid.var,beta=true.coefs,R2=0.1)
  return(resu.sim1)
}


simtestY2 <- function(X=data.X,n.true.exposures=3,nvars=100,N = 1200){
  # selected at random
  true.exposures = sample(1:nvars,n.true.exposures,replace=F)
  true.coefs = rep(1,n.true.exposures) 
  
  Xtrue = as.matrix(cbind(X[,paste0("X",true.exposures,".2")]))
  Xtrue.names = paste0("X",true.exposures,".2")
  
  linpred = Xtrue %*% matrix(true.coefs,ncol=1)
  # R2 of the model = 5%
  resid.var = var(linpred)*0.95/.05
  
  Y <- linpred + rnorm(N,0,sd=sqrt(resid.var))
  formul = paste("Y ~", paste(Xtrue.names,sep=" ", collapse="+"))
  mod = lm(formul,data=cbind(Y,X))
  summary(mod)
  # 
  # dataXtrue.names<-cbind(dataXtrue.names,Xtrue.names)
  # dataResid.var<-cbind(dataResid.var,resid.var)
  # dataTrue.coefs<-cbind(dataTrue.coefs,true.coefs)
  # 
  # dataY<-Y  
  resu.sim2<-list(resp=Y,true.pred=Xtrue.names,sigmaE=resid.var,beta=true.coefs,R2=0.1)
  return(resu.sim2)
}





##### --------------------------------------------------------- 2. IMPLEMENTATION OF EACH METHOD --------------------------------------------------------- ##### 



##########################################################################
## Functions for the implemmentation of each of assessed methods: ExWAS ##
##########################################################################


EWAS <- function(resp,Exp,beta,interaction=FALSE,add=FALSE){ 
	  require(parallel) 
	  # build the name for interactions of variables in name
	  buildName<- function(name){ 
	    if(length(name)>1) 
	      for(i in 2:length(name))    
		name[1] <- paste(name[1],":",name[i],sep="") 
	      return(name[1]) 
	  } 
	  if(interaction==1) interaction <- FALSE
	  beta <- cbind(var=gsub("[:]","*",names(beta)),val=beta)
	  # verifying formats and values of inputs and computing interactions 
	  if(interaction==FALSE) interaction <- 0 
	  if(!all(is.null(colnames(Exp)))) Exp <- Exp[,order(colnames(Exp))]  
	  # EWAS  
	  if(interaction!=FALSE) if(interaction>1)  Exp <- cbind(Exp,AddInteract(Exp,max=interaction))
	  p.values <- mclapply(1:ncol(Exp), function(x,Exp){ 
	    c(colnames(Exp)[x],summary(lm("resp~var", data=data.frame(cbind(var=Exp[,x],resp=Exp["resp"]))))$coefficients[2,])}, Exp=cbind(Exp,resp=resp)) 
	  p.values <- cbind(matrix(unlist(p.values), ncol = 5, byrow = TRUE)[,-4]) 
	  colnames(p.values) <- c("var","Est","Sd","pVal") 
	  if (all(beta%in%p.values[,1])) {
	    p.values <- merge(beta,p.values,by="var",all=TRUE)
	  }
	  else {
	    p.values <- merge(beta,p.values,by="var",all.y=TRUE)
	  }
	  #    p.values <- p.values[-which(p.values$var%in%c("(Intercept)","Intercept")),]
	  return(p.values) 
} 

testReg <- function(Exp,outc, p.values,intTest=1,corr,parTest=NA,Nperm=NA, lambda=NA, interaction=NA,EWASonly=FALSE){   
	  p.values <- p.values[order(as.character(p.values$var)),]  
	  pVal <- as.numeric(as.character(p.values$pVal))  
	  Exp <- Exp[,order(colnames(Exp))]   
	  if(!all(p.values$var==colnames(Exp))) stop("issue in ordering")   
	  resp <- outc$resp 
	  # fitting multiple linear model   
	  if(corr=="None") wh <- which(pVal<=0.05)   
	  if(corr=="Bon") wh <- which(pVal<=0.05/nrow(p.values))   
	  if(corr=="BH") wh <- which(p.adjust(pVal,"BH")<=0.05)   
	  if(corr=="BY") wh <- which(p.adjust(pVal,"BY")<=0.05)   
	  if(!corr%in%c("Bon","BH","BY","None")) stop("Please specify a known correction method for multiple testing") 
	  # EWAS   
	  wh <- p.values$var[wh]   
	  p.val <- data.frame(cbind(var=as.character(p.values$var),EWAS=as.numeric(as.character(p.values$Est)),EWAS.TP=as.numeric(p.values[,1]%in%wh)))   
	  if(EWASonly==TRUE){ RMSE=NA; R2=NA } 
	  # EWAS-LM   
	  if(EWASonly==FALSE){    
	    formula <- "~1"    
	    if(length(wh)>0) for(i in 1:length(wh))  formula <- paste(formula,"+",wh[i])    
	    modelLM <- lm(as.formula(paste("resp",formula)),data=cbind(Exp,resp=resp))   
	    predLM <- predict(modelLM, Exp, se.fit = FALSE) 
	    modelLM <- cbind(var=as.character(names(modelLM$coef)),EWAS_LM=modelLM$coef,EWAS_LM.TP=as.numeric(summary(modelLM)$coef[,4]<.05))   
	    p.val <- merge(p.val,matrix(modelLM[-1,],ncol=3,dimnames=list(NULL,colnames(modelLM))),by="var",all=TRUE)   
	    R2 <- c(NA,NA,var(predLM)/var(outc$resp))    
	    RMSE <- c(NA,NA,sqrt(mean((predLM-outc$resp)^2)))   
	  }   
	  if(nrow(p.val)!=ncol(Exp)) print("erreur de dim")  
	  if(!all(p.val$var==colnames(Exp))) print("erreur d'ordre de variables")  
	  if(any(is.na(p.val[,2])| p.val[,2]=="NaN"))  print("erreur pour EWAS")  
	  names <- p.val[,1]   
	  p.val <- t(p.val[,-1])   
	  colnames(p.val) <- names   
	  p.val <- cbind(RMSE=RMSE,R2=R2,p.val)   
	  return(p.val) 
}


applyEwas<-function(data.Y,data.X) {
	  beta <- data.Y$beta
	  names(beta)<-data.Y$true.pred
	  
	  # the method to correct for multiple testing (amongst Bon=Bonferroni, BH=Benjamini-Hochberg, BY=Benjamini-Yakutieli) 
	  res <- EWAS(resp=data.Y$resp, Exp=data.X, beta=beta, interaction=1)
	  resEWAS  <- testReg(Exp=data.X, outc=data.Y, p.values= res[,c("var","pVal")], intTest=1, corr="None", par=1, Nperm=NA, lambda=NA, interaction=1,EWASonly=FALSE) 
	  RES <- merge(res[,c("var","val","Est","pVal")],data.frame(cbind(var=colnames(resEWAS)[-(1:2)],t(resEWAS[c(1,2,3),-(1:2)]))),by="var",all.y=TRUE) 
	  colnames(RES)<-c("var","val","Est","pVal","EWAS.TP.none","EWAS_LM.none","EWAS_LM.TP.none")
	  # ExWAS with Bonferoni correction
	  resEWAS.bon  <- testReg(Exp=data.X, outc=data.Y, p.values= res[,c("var","pVal")], intTest=1, corr="Bon", par=1, Nperm=NA, lambda=NA, interaction=1,EWASonly=FALSE) 
	  RES <- merge(RES,data.frame(cbind(var=colnames(resEWAS.bon)[-(1:2)],t(resEWAS.bon[c(1,2,3),-(1:2)]))),by="var",all.y=TRUE) 
	  colnames(RES)[8:10]<-c("EWAS.TP.bon","EWAS_LM.bon","EWAS_LM.TP.bon")
	  # ExWAS with Benjamini-Hochberg correction
	  resEWAS.bh  <- testReg(Exp=data.X, outc=data.Y, p.values= res[,c("var","pVal")], intTest=1, corr="BH", par=1, Nperm=NA, lambda=NA, interaction=1,EWASonly=FALSE) 
	  RES <- merge(RES,data.frame(cbind(var=colnames(resEWAS.bh)[-(1:2)],t(resEWAS.bh[c(1,2,3),-(1:2)]))),by="var",all.y=TRUE) 
	  colnames(RES)[11:13]<-c("EWAS.TP.bh","EWAS_LM.bh","EWAS_LM.TP.bh")
	  # ExWAS with Benjamini-Yakutieli correction
	  resEWAS.by  <- testReg(Exp=data.X, outc=data.Y, p.values= res[,c("var","pVal")], intTest=1, corr="BY", par=1, Nperm=NA, lambda=NA, interaction=1,EWASonly=FALSE) 
	  RES <- merge(RES,data.frame(cbind(var=colnames(resEWAS.by)[-(1:2)],t(resEWAS.by[c(1,2,3),-(1:2)]))),by="var",all.y=TRUE) 
	  colnames(RES)[14:16]<-c("EWAS.TP.by","EWAS_LM.by","EWAS_LM.TP.by")
	  
	  #summary
	  RES2<-RES
	  RES2$true.pred<-as.factor(ifelse(is.na(RES2$val),"0","1"))
	  RES2$EWAS_LM.TP.none<-as.factor(ifelse(RES2$EWAS_LM.TP.none%in%"1","1","0"))
	  RES2$EWAS_LM.TP.bon<-as.factor(ifelse(RES2$EWAS_LM.TP.bon%in%"1","1","0"))
	  RES2$EWAS_LM.TP.bh<-as.factor(ifelse(RES2$EWAS_LM.TP.bh%in%"1","1","0"))
	  RES2$EWAS_LM.TP.by<-as.factor(ifelse(RES2$EWAS_LM.TP.by%in%"1","1","0"))
	  
	  return(RES2)
}





######################################################################################
## Functions for the implemmentation of each of assessed methods: ElasticNet (ENET) ##
######################################################################################


EvaluateModel.cv.glmnet <- function(x,y,nfolds=10) { 
	  n <- nrow(x) 
	  folds <- rep(1:nfolds,length.out=n)[sample(n,n)] #step 1: do all crossvalidations for each alpha in range 0.5 - 1.0 
	  alphasOfInterest <- seq(0.5,1,0.1) 
	  cvs <- lapply(alphasOfInterest, function(curAlpha){
	    cv.glmnet(x,y,alpha=curAlpha,family="gaussian",nfolds=nfolds,foldid=folds,standardize=FALSE)
	  }) 
	  #step 2: collect the optimum lambda for each alpha 
	  optimumPerAlpha <- sapply(seq_along(alphasOfInterest), function(curi){ 
	    curcvs <- cvs[[curi]] 
	    curAlpha <- alphasOfInterest[curi] 
	    indOfMin <- match(curcvs$lambda.min, curcvs$lambda)
	    return(c(lam=curcvs$lambda.min,alph=curAlpha,cvm=curcvs$cvm[indOfMin],cvup=curcvs$cvup[indOfMin])) 
	  }) 
	  #step 3: find the overall optimum 
	  posOfOptimum <- which.min(optimumPerAlpha["cvm",]) 
	  overall.alpha.min <- optimumPerAlpha["alph",posOfOptimum] 
	  overall.lambda.min <- optimumPerAlpha["lam",posOfOptimum] 
	  overall.criterionthreshold <- optimumPerAlpha["cvup",posOfOptimum] 
	  #step 4: now check for each alpha which lambda is the best within the threshold 
	  corrected1se <- sapply(seq_along(alphasOfInterest), function(curi){ 
	    curcvs <- cvs[[curi]] 
	    lams <- curcvs$lambda 
	    lams[lams<overall.lambda.min] <- NA 
	    lams[curcvs$cvm > overall.criterionthreshold] <- NA 
	    lam1se <- max(lams, na.rm=TRUE) 
	    c(lam=lam1se, alph=alphasOfInterest[curi]) 
	  }) 
	  #step 5: find the best (lowest) of these lambdas 
	  overall.alpha.1se <- max(corrected1se["alph",is.finite(corrected1se["lam",])]) 
	  overall.lambda.1se <- corrected1se["lam",match(overall.alpha.1se, alphasOfInterest)] 
	  model <- glmnet(x=x,y=y,alpha=overall.alpha.min,lambda=overall.lambda.min,family="gaussian",standardize=FALSE) 
	  coef <- as.matrix(predict(model,type="coefficients",s=as.numeric(overall.lambda.min)))[-1] 
	  names(coef) <- paste0("C",1:length(coef)) 
	  results.min <- cbind(data.frame(model="ENET",type="MIN"),t(coef)) 
	  model <- glmnet(x=x,y=y,alpha=overall.alpha.1se,lambda=overall.lambda.1se,family="gaussian",standardize=FALSE) 
	  coef <- as.matrix(predict(model,type="coefficients",s=as.numeric(overall.lambda.1se)))[-1] 
	  names(coef) <- paste0("C",1:length(coef)) 
	  results.opt <- cbind(data.frame(model="ENET",type="OPT"),t(coef)) 
	  return(rbind(results.min,results.opt)) 
}


applyEnet <-function(data.Y,data.X) {
  beta <- data.Y$beta
  names(beta)<-data.Y$true.pred
  RES<-init
  # ENET
  resENET <- EvaluateModel.cv.glmnet(x=as.matrix(data.X),y=data.Y$resp)
  RES <- merge(RES,cbind(var=colnames(data.X),t(resENET[,-(1:2)])),by='var',all.x=TRUE) 
  colnames(RES)[ncol(RES)+1-nrow(resENET):1] <- paste(resENET[,1],"_",resENET[,2],sep="") 
  
  #summary
  RES2<-RES
  RES2$ENET_MIN.TP<-as.factor(ifelse(is.na(RES2$ENET_MIN),"0",ifelse(RES2$ENET_MIN%in%"0","0","1")))
  RES2$ENET_OPT.TP<-as.factor(ifelse(is.na(RES2$ENET_OPT),"0",ifelse(RES2$ENET_OPT%in%"0","0","1")))
  
  return(RES2)
}





######################################################################################
## Functions for the implemmentation of each of assessed methods: sPLS	       ##
######################################################################################


EvaluateModel.cv.spls <- function(x,y,max.K=5,nfolds=10) {   
	  m <- spls::cv.spls(x=x,y=y,K=1:max.K,eta=c(1e-6,seq(0.1,0.9,0.1),1-1e-6),fold=nfolds,plot.it=FALSE,scale.x=FALSE,scale.y=FALSE)   
	  RMSE_H0 <- mean((y-mean(y))^2)   
	  K.min <- m$K.opt   
	  eta.min <- m$eta.opt   
	  if (RMSE_H0<=min(m$mspemat)) {     
	    K.min <- 0   
	  }   
	  if (K.min>0) {     
	    model <- spls::spls(x=x,y=y,K=K.min,eta=eta.min)     
	    coef <- as.numeric(coef(model))   
	  } 
	  else {     
	    coef <- rep(0,ncol(x))   
	  }   
	  names(coef) <- paste0("C",1:length(coef))   
	  results <- cbind(data.frame(model="SPLS",type="MIN"),t(coef))   
	  return(results) 
}


applySPLS <-function(data.Y,data.X) {
	  beta <- data.Y$beta
	  names(beta)<-data.Y$true.pred
	  RES<-init
	  # sPLS
	  resSPLS <- EvaluateModel.cv.spls(x=as.matrix(data.X),y=data.Y$resp)
	  colnames(resSPLS)[-(1:2)] <- colnames(data.X) 
	  RES <- merge(RES,cbind(var=colnames(resSPLS),t(resSPLS)),by='var',all.x=TRUE) 
	  colnames(RES)[ncol(RES)+1-nrow(resSPLS):1] <- paste(resSPLS[,1],"_", resSPLS[,2],sep="") 
	  
	  #summary
	  RES2<-RES
	  RES2$SPLS_MIN.TP<-as.factor(ifelse(is.na(RES2$SPLS_MIN),"0",ifelse(RES2$SPLS_MIN%in%"0","0","1")))
	  
	  return(RES2)
}





######################################################################################
## Functions for the implemmentation of each of assessed methods: sNPLS             ##
######################################################################################


applysNPLS <-function(data.Y,data.X) {

	  #Code designed to 5 time points. Need to be adjusted in case a different longitudinal structure is simulated.
	  X <- data.X
	  Y <- data.Y$resp
	  nvars <- ncol(resu.sim.dataX.i[[1]]$X)/5
	  names <- str_replace(colnames(X), '.+.(.+)', '\\1')
	  Ts <- c("1","2","3","4","5")
	  arr = array(NA, dim=c(nrow(X),ncol(X)/5,5))
	  for (i in Ts) {
	 	arr[,,as.numeric(i)] <- as.matrix(X[,grep(i,names)])
			}
			
	  cv_result <- repeat_cv(arr,Y,ncomp = 1:3, keepJ = 1:nvars, keepK = 1:5, parallel = TRUE, times = 5, nfold = 5)
	  select_best_hp <- ddply(cv_result,.(ncomp,keepJ,keepK),nrow)
		if ( length(which(select_best_hp$V1 %in% max(select_best_hp$V1))) > 1 ) {
	  	selected_hp <- select_best_hp[which(select_best_hp$V1 %in% max(select_best_hp$V1)),1:3]
	  	if (any(selected_hp$ncomp > 1)){
	  	  selected_hp <- subset(selected_hp, ncomp > 1)
	  	}else{}
	  	if (nrow(selected_hp) > 1){
	  	  selected_hp <- selected_hp[sample(1:nrow(selected_hp), 1),]
	  	  }else{}
			}else{
				selected_hp <- select_best_hp[which(select_best_hp$V1 %in% max(select_best_hp$V1)),1:3]  
				}

	  if (selected_hp[1,1]<1){ COMPs=1}else{ COMPs=selected_hp[1,1]}
	  
	  fit <- sNPLS(arr,Y, ncomp = COMPs, keepJ = rep(selected_hp[1,2],COMPs), keepK = rep(selected_hp[1,3],COMPs), silent = FALSE)
	  
	  
	  output_coefs <- as.data.frame(summary(fit))
	  rownames(output_coefs) <- gsub("X.","X",rownames(output_coefs))
	  names(output_coefs) <- c(1:5)
	  output_coefs$var <- rownames(output_coefs)
	  output_coefs <- output_coefs %>% pivot_longer(!var, names_to = "time", values_to = "coef")
	  output_coefs$var <- paste(output_coefs$var,output_coefs$time, sep = ".")
	  output_coefs <- output_coefs[,c(1,3)]
	  output_coefs$sNPLS.TP <- as.factor(ifelse(is.na(output_coefs$coef),"0",ifelse(output_coefs$coef%in%"0","0","1")))
	  RES <- output_coefs
	  return(RES)
  
}





#########################################################################
## Functions for the implemmentation of each of assessed methods: DSA  ##
#########################################################################


DSAreg <- function(Exp,resp, family = gaussian,maxsize = maxsize, maxsumofpow = 1, maxorderint = 1){ 
	  Exp <- data.frame(cbind(data.frame(Exp), resp=resp)) 
	  res <- DSA(resp ~ 1, data = Exp,  family = family, maxsize = maxsize, maxsumofpow = maxsumofpow, maxorderint = maxorderint ,nsplits=1,usersplits = NULL) 
	  form <- gsub("I[(]","",colnames(coefficients(res))) 
	  form <- gsub("[*]",":",gsub("[)]","",gsub("[:^:]1","",form))) 
	  if(length(grep(":",form))>0){  
	    nam <- strsplit(form[grep("[:]",form)],":")   
	    for(j in 1:length(nam)){     
	      nam[[j]] <- gsub("[[:space:]]","",nam[[j]])     
	      name <- nam[[j]][1]     
	      for(k in 2:length(nam[[j]]))  name <- paste(name,":",nam[[j]][k],sep="")     
	      Exp <- cbind(Exp,name=apply(Exp[,nam[[j]]],1,prod)) 
	    }
	  }  
	  form2 <- "resp~1"   
	  if(length(form)>1) for(i in 2:length(form)) form2 <- paste(form2,"+",form[i])  
	  res2 <- lm(form2, data=data.frame(Exp)) 
	  pred <- predict(res2,Exp)  
	  coef <- summary(res2)$coefficients 
	  coef <- cbind(var=rownames(coef),valDSAEst=coef[,1])  
	  return(list(coef=coef, pred=pred)) 
} 


applyDSA <- function(data.Y,data.X,maxsize) {
	  beta <- data.Y$beta
	  names(beta)<-data.Y$true.pred
	  RES<-init
	  # DSA 
	  resDSA <- DSAreg(Exp=data.X, resp= data.Y$resp, maxsize = maxsize, maxsumofpow =1, maxorderint = 1) 
	  RES <- merge(RES,resDSA$coef,by='var',all.x=TRUE) 
	  RES[,ncol(RES)] <- as.numeric(as.character(RES[,ncol(RES)]))
	  RES[which(is.na(RES[,ncol(RES)])), ncol(RES)] <- 0 
	  colnames(RES)[ncol(RES)] <- "DSA" 
	  
	  #summary
	  RES$DSA.TP<-as.factor(ifelse(is.na(RES$DSA),"0",ifelse(RES$DSA%in%"0","0","1")))
	  
	  RES2<-RES
	  return(RES2)
}





##########################################################################################################
## Functions for the implemmentation of each of assessed methods: Penalized DLNM with feature selection ##
##########################################################################################################


applyDLNMselect <- function(data.Y,data.X){

	  start<-seq(1,500,by = 5)
	  pvalue<-NULL
	  results<-NULL
	  resp<-data.Y$resp
	  exp<-data.X

	  # total number of periods
	  n <- 5  
	  npred<-100
	  
	  L = matrix(rep(1:5,length(resp)),byrow=T,ncol=5)
	  
	  # create 100 subset of dataframe (1 for each exposure with the 5 time points)
	  for (j in 1:npred){

	    nam <- paste("x", j, sep = "")
	    assign(nam, as.matrix(exp[,c(start[j]:(start[j]+4))]))
	    
	  }
	  
	  # build formula that includes all npred predictors
	  form = "resp ~ s(x1, L, bs='cb', xt=list(argvar=list(fun = 'lin'),arglag=list(fun='cr',knots = 1:n))) + "
	  
	  for (j in 2:(npred-1)) {
	    
	    form = paste0(form, "s(x",j,", L, bs='cb', xt=list(argvar=list(fun = 'lin'),arglag=list(fun='cr',knots = 1:n))) + ") 
	  
	  }
	  
	  form = paste0(form, "s(x",npred,", L, bs='cb', xt=list(argvar=list(fun = 'lin'),arglag=list(fun='cr',knots = 1:n)))")
	  form = as.formula(form)  
	  
	  # Run model
	  ctrl <- list(nthreads=24)
	  mod.select <- gam(form, select=T,method="REML",cluster=ctrl)
	  sum.mod = summary(mod.select)
	  
	  # P-values for each variable
	  tab.sum<-as.data.frame(sum.mod$s.table)
	  
	  res.all<-NULL
	  
	  for (i in 1:npred){

	    expo_name<-paste0("x",i)
	    pred<-crosspred(expo_name,model=mod.select,at=1,cen=0)
	    coef<-pred$coefficients
	    ci.inf<-pred$matlow
	    ci.sup<-pred$mathigh
	    res<-t(rbind(coef,ci.inf,ci.sup))
	    res.all<-rbind(res.all,res)
	    
	  }
	  
	  tab.sum$expo_name<-ifelse(nchar(as.character(rownames(tab.sum)))==7,substring(rownames(tab.sum),first=3,last=4),
		                    ifelse(nchar(as.character(rownames(tab.sum)))==8,substring(rownames(tab.sum),first=3,last=5),
		                           substring(rownames(tab.sum),first=3,last=6)))
	  
	  res.all<-as.data.frame(res.all,ncol=3)
	  res.all$var<-rownames(res.all)
	  colnames(res.all)<-c("coef","CI.inf","CI.sup","var")
	  res.all$expo_name<-ifelse(nchar(as.character(rownames(res.all)))==9,substring(rownames(res.all),first=3,last=4),
		                    ifelse(nchar(as.character(rownames(res.all)))==10,substring(rownames(res.all),first=3,last=5),
		                           substring(rownames(res.all),first=3,last=6)))
	  
	  results<-merge(tab.sum,res.all,by="expo_name")
	  
	  return(results)
	  
}  





##############################################################################################################################
## Functions for the implemmentation of each of assessed methods: Penalized DLNM with feature selection (Two-step approach) ##
##############################################################################################################################


applyDLNMselectAVG <-function(data.Y,data.X){

  start<-seq(1,ncol(data.X),by = 5)
  pvalue<-NULL
  results<-NULL
  resp<-data.Y$resp
  exp<-data.X
    
  # Total number of periods and predictors
  n <- 5  
  npred<-ncol(data.X)/5
  varS <- unique(gsub("\\..*","",colnames(data.X)))
  varS <- as.numeric(gsub("X","",varS))
  # Create matrix of lags
  L = matrix(rep(1:5,length(resp)),byrow=T,ncol=5)
  
  # Create npred subset of dataframe (1 for each exposure with the 5 time points)
  for (j in 1:npred){
    nam <- paste("x", varS[j], sep = "")
    assign(nam, as.matrix(exp[,c(start[j]:(start[j]+4))]))
  }
  
  # Build formula that includes all npred predictors
  form = paste0("resp ~ s(x",varS[1],", L, bs='cb', xt=list(argvar=list(fun = 'lin'),arglag=list(fun='cr',knots = 1:n))) + ") 
  
  for (j in 2:(npred-1)) {
    
    form = paste0(form, "s(x",varS[j],", L, bs='cb', xt=list(argvar=list(fun = 'lin'),arglag=list(fun='cr',knots = 1:n))) + ") 
    
  }
  
  form = paste0(form, "s(x",varS[npred],", L, bs='cb', xt=list(argvar=list(fun = 'lin'),arglag=list(fun='cr',knots = 1:n)))")
  
  form = as.formula(form)  
  
   # Run DLNM model with penalization (selection of smootheness parameter by RELM) and with variable selection
  ctrl <- list(nthreads=24)
  mod.select <- gam(form, select=T,method="REML",cluster=ctrl)
  
  # Retrieve statistics from the model
  sum.mod = summary(mod.select)
  
  # Get approximated P-values for each variable
  tab.sum<-as.data.frame(sum.mod$s.table)
  
  tab.sum$expo_name<-ifelse(nchar(as.character(rownames(tab.sum)))==7,substring(rownames(tab.sum),first=3,last=4),
                            ifelse(nchar(as.character(rownames(tab.sum)))==8,substring(rownames(tab.sum),first=3,last=5),
                                   substring(rownames(tab.sum),first=3,last=6)))
                                   
                                    
  res.all<-NULL
  for (z in tab.sum$expo_name){
    expo_name<-z
    pred<-crosspred(expo_name,model=mod.select,at=1,cen=0)
    coef<-pred$coefficients
    ci.inf<-pred$matlow
    ci.sup<-pred$mathigh
    res<-t(rbind(coef,ci.inf,ci.sup))
    res.all<-rbind(res.all,res)
  }
  
   
  tab.sum$expo_name<-ifelse(nchar(as.character(rownames(tab.sum)))==7,substring(rownames(tab.sum),first=3,last=4),
                            ifelse(nchar(as.character(rownames(tab.sum)))==8,substring(rownames(tab.sum),first=3,last=5),
                                   substring(rownames(tab.sum),first=3,last=6)))
  
  res.all<-as.data.frame(res.all,ncol=3)
  res.all$var<-rownames(res.all)
  colnames(res.all)<-c("coef","CI.inf","CI.sup","var")
  res.all$expo_name<-ifelse(nchar(as.character(rownames(res.all)))==9,substring(rownames(res.all),first=3,last=4),
                            ifelse(nchar(as.character(rownames(res.all)))==10,substring(rownames(res.all),first=3,last=5),
                                   substring(rownames(res.all),first=3,last=6)))
  
  results<-merge(tab.sum,res.all,by="expo_name")
  
  return(results)
  
}  





##### --------------------------------------------------------- 3. SUMMARY OF THE PERFOMANCES OF THE METHODS --------------------------------------------------------- ##### 




##################################################################################
## Functions to create 2x2 table to compare the expected and observed selection ##
## 	This function is called within the calculateMetrics() function		##
##################################################################################

squareTable <- function(x,y) {
  x <- factor(x)
  y <- factor(y)
  commonLevels <- c("0","1")
  x <- factor(x, levels = commonLevels)
  y <- factor(y, levels = commonLevels)
  table(x,y)
}


############################################################################################################
## Function to arrange the input dataset to calculate the performance independently of the true timepoint ##
## 		This function is called within the calculateMetrics() function when bytime=0		  ##
############################################################################################################


data100 <- function(data) {
  
  df <- data
  
  # arrange df to calculate the performance whatever the true timepoint
  df$expo_name<-ifelse(nchar(as.character(df$var))==4,substring(df$var,first=1,last=2),
                       ifelse(nchar(as.character(df$var))==5,substring(df$var,first=1,last=3),
                              substring(df$var,first=1,last=4)))
  nb_col<-dim(df)[2]
  df<-df[,c(1:3,nb_col,4:(nb_col<-dim(df)[2]-1))]
  
  numcol_name<-c("true.pred","EWAS.TP.none","EWAS_LM.TP.none","EWAS.TP.bon","EWAS_LM.TP.bon","EWAS.TP.bh",
                 "EWAS_LM.TP.bh", "EWAS.TP.by","EWAS_LM.TP.by","ENET_MIN.TP","ENET_OPT.TP","SPLS_MIN.TP",
                 "sNPLS.TP","MMPC.TP","DSA.TP","DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                 "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",  
                 "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
                 "DLNMselectforward.none","DLNMselectforward.bonf","DLNMselectforward.bh","DLNMselectforward.by")
  
  for (j in which(colnames(df)%in%numcol_name)){
    df[,j]<-as.numeric(as.character(df[,j]))
  } 
  

  df_100 <- df %>% group_by(numsim,expo_name) %>% 
    summarise_at(.vars = numcol_name,
                 .funs = sum)  
  
  
  df_100<-as.data.frame(df_100)
  
  # Replace all values greater than 0 with 1 in specific columns
  df_100 <- df_100 %>% mutate(across(all_of(numcol_name), ~ifelse(. > 0, 1, .)))
  
  return(df_100)
}  



#######################################################################################################
## Function to calculate the performance of each method, both considering or not the true time point ##
#######################################################################################################


calculateMetrics <- function(data,bytime) {	# bytime=1 to calculate the performances to identify the true time point and bytime=0 to calculate the performances independently of the true time point 
  if (bytime==1) {
    # Create an empty data frame to store the detailed results
    performance_detail <- data.frame(numsim = integer(),
                                       method = character(),
                                       sensitivity = numeric(),
                                       specificity = numeric(),
                                       fpr = numeric(),
                                       fnr = numeric(),
                                       false_discovery_rate = numeric(),
                                       false_omit_rate = numeric(),
                                       stringsAsFactors = FALSE)
    
    # Create an empty data frame to store the summarized results
    performance_summary <- data.frame(method = character(),
                                  mean_sensitivity = numeric(),
                                  mean_specificity = numeric(),
                                  mean_fpr = numeric(),
                                  mean_fnr = numeric(),
                                  mean_false_discovery_rate = numeric(),
                                  mean_false_omit_rate = numeric(),
                                  stringsAsFactors = FALSE)
    
    df <- data
  }
  
  else if (bytime==0) {
    # Create an empty data frame to store the detailed results
    performance_detail <- data.frame(numsim = integer(),
                                       method = character(),
                                       sensitivity = numeric(),
                                       specificity = numeric(),
                                       fpr = numeric(),
                                       fnr = numeric(),
                                       false_discovery_rate = numeric(),
                                       false_omit_rate = numeric(),
                                       stringsAsFactors = FALSE)
    
    # Create an empty data frame to store the summarized results
    performance_summary <- data.frame(method = character(),
                                  mean_sensitivity = numeric(),
                                  mean_specificity = numeric(),
                                  mean_fpr = numeric(),
                                  mean_fnr = numeric(),
                                  mean_false_discovery_rate = numeric(),
                                  mean_false_omit_rate = numeric(),
                                  stringsAsFactors = FALSE)
    
    df <- data100(data)
  }
  
  nb_col<-dim(df)[2]

  # list of methods of interest
  numcol_method<-which(colnames(df)%in%c("EWAS.TP.none","EWAS_LM.TP.none","EWAS.TP.bon","EWAS_LM.TP.bon","EWAS.TP.bh",
                                         "EWAS_LM.TP.bh", "EWAS.TP.by","EWAS_LM.TP.by","ENET_MIN.TP","ENET_OPT.TP","SPLS_MIN.TP",
                                         "sNPLS.TP","MMPC.TP","DSA.TP","DLNMpen.none","DLNMpen.bonf","DLNMpen.bh","DLNMpen.by",
                                         "DLNMselect.none","DLNMselect.bonf","DLNMselect.bh","DLNMselect.by",  
                                         "DLNMselectback.none","DLNMselectback.bonf","DLNMselectback.bh","DLNMselectback.by",
                                         "DLNMselectforward.none","DLNMselectforward.bonf","DLNMselectforward.bh","DLNMselectforward.by"))
  
  for (i in 1:100) {
    # Iterate over each method
    for (method_col in numcol_method) {
      df2<-df[df$numsim==i,]
      # Get the column name representing the method
      method_name <- names(df2)[method_col]
      
      # Extract the true predictor and method 
      true_pred <- df2$true.pred
      method_sel <- df2[, method_col]
      
      # Calculate the contingency table using the squareTable function
      table_result <- squareTable(method_sel, true_pred)
      
      # Calculate TP, TN, FP, FN
      tp <- table_result[2, 2]
      tn <- table_result[1, 1]
      fp <- table_result[2, 1]
      fn <- table_result[1, 2]
      
      # Calculate sensitivity, specificity, false_positive_rate, false_negative_rate, false discovery rate, false omit rate
      sensitivity <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
      specificity <- ifelse((tn + fp) == 0, 0, tn / (tn + fp))
      false_positive_rate <- ifelse((fp + tn) == 0, 0, fp / (fp + tn))
      false_negative_rate <- ifelse((fn + tp) == 0, 0, fn / (fn + tp))
      false_discovery_rate <- ifelse((fp + tp) == 0, 0, fp / (fp + tp))
      false_omit_rate <- ifelse((fn + tn) == 0, 0, fn / (fn + tn))
      
      # Store detailed results
      performance_detail <- rbind(performance_detail, data.frame(numsim = i,
                                                                     method = method_name,
                                                                     sensitivity = sensitivity,
                                                                     specificity = specificity,
                                                                     false_positive_rate = false_positive_rate,
                                                                     false_negative_rate = false_negative_rate,
                                                                     false_discovery_rate = false_discovery_rate,
                                                                     false_omit_rate = false_omit_rate,
                                                                     stringsAsFactors = FALSE))
    }
    
    performance_summary<-performance_detail %>% 
      group_by(method) %>% 
      summarise(
        false.disc.rate_mean = mean(false_discovery_rate),
        false.disc.rate_sd = sd(false_discovery_rate),
        sensitivity_mean = mean(sensitivity),
        sensitivity_sd = sd(sensitivity)
      )
    
  }
  
  
  # Return both the detailed and summarized results
  return(list(performance_detail = performance_detail, performance_summary = performance_summary))
}


# Call the calculateMetrics function to identify the true exposures at the true time point (bytime=1)
results_500  <- calculateMetrics(data=RES1.all.exp3,bytime=1)
# Access the performances of each method by simulated dataset 
performance_detail_500 <- results_500$performance_detail
# Access the overall performances of each method (mean and sd of all simulated datasets)
performance_summary_500 <- results_500$performance_summary

# Call the calculateMetrics function to identify the true exposures independently of the true time point (bytime=0)
results_100  <- calculateMetrics(data=RES1.all.exp3,bytime=0)
# Access the performances of each method by simulated dataset 
performance_detail_100 <- results_100$performance_detail
# Access the overall performances of each method (mean and sd of all simulated datasets)
performance_summary_100 <- results_100$performance_summary




