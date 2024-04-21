# The functions in this file belong to the inTrees package. With the following
# details:

# Package: inTrees
# Title: Interpret Tree Ensembles
# Version: 1.3
# Date: 2022-05-31
# Imports: RRF, arules, gbm, xtable, xgboost, data.table, methods
# Author: Houtao Deng, Xin Guan, Vadim Khotilovich
# Maintainer: Houtao Deng <softwaredeng@gmail.com>
# Description: For tree ensembles such as random forests, regularized random
# forests and gradient boosted trees, this package provides functions for:
# extracting, measuring and pruning rules; selecting a compact rule set;
# summarizing rules into a learner; calculating frequent variable interactions;
# formatting rules in latex code.
# BugReports: https://github.com/softwaredeng/inTrees/issues
# License: GPL (>= 3)
# NeedsCompilation: no
# Repository: CRAN
# Date/Publication: 2022-05-31

# The inTrees package has been set to be archived. Upon uploading a new
# version of the inTrees package to the CRAN, the maintainers of the CRE
# package will remove this file and import the inTrees package directly from
# CRAN.

# Loaded libraries for inTrees:
# stats (in common with CRE)
# xtable
# data.table (in common with CRE)
# RRF from RRF
# getTree from RRF
# pretty.gbm.tree from gbm
# xgb.model.dt.tree from xgboost
# arules
# as from methods

# Used functions:
# inTrees_treeVisit
# inTrees_ruleList2Exec
# inTrees_pruneRule
# inTrees_getRuleMetric

require(xgboost)
require(data.table)
require(arules)
require(methods)

voteAllRules <-
  function(ruleMetric,X,type="r",method="median"){
    xVoteList = vector("list",nrow(X))
    predY <- rep("",nrow(X))
    for(i in 1:nrow(ruleMetric)){
      ixMatch <- eval(parse(text=paste("which(",ruleMetric[i,"condition"], ")"))  )
      if(length(ixMatch)==0) next
      for(ii in ixMatch){
        xVoteList[[ii]] = c(xVoteList[[ii]], ruleMetric[i,"pred"])
      }
    }
    for(i in 1:length(xVoteList)){
      thisV <- xVoteList[[i]]
      if(length(thisV)==0) next
      if(type == "c") predY[i] <- names(table(thisV)[which.max(table(thisV))])
      if(type == "r"){
        thisV = as.numeric(thisV)
        if(method == "median"){
          predY[i] <- median(thisV)
        }else{
          predY[i] <- mean(thisV)
        }
      }

    }
    if(type=="r") predY <- as.numeric(predY)
    return(predY)
  }


inTrees_treeVisit <-
  function(tree,rowIx,count,ruleSet,rule,levelX,length,max_length,digits=NULL)
  {
    if( tree[rowIx,"status"] == -1 | length == max_length ){
      count = count + 1
      ruleSet[[count]] = rule
      return(list(ruleSet = ruleSet, count=count))
    }
    xIx <- tree[rowIx,"split var"]
    xValue <- tree[rowIx,"split point"]
    if(is.integer(digits)) xValue <- round(tree[rowIx,"split point"], digits)

    if(is.null(levelX[[xIx]])){
      lValue <- paste("X[,",xIx, "]<=",xValue,sep="")
      rValue <- paste("X[,",xIx, "]>",xValue,sep="")
    }else{
      xValue<- which(as.integer(intToBits(as.integer(xValue)))>0)
      lValue <- levelX[[xIx]][xValue]
      rValue <- setdiff(levelX[[xIx]],lValue)
      #   lValue <- paste("X[,",xIx, "]%in% '",lValue,"'",sep="")
      #   rValue <- paste("X[,",xIx, "]%in% '",rValue,"'",sep="")
    }
    xValue <- NULL
    ruleleft <- rule
    if(length(ruleleft)==0)
    {
      ruleleft[[as.character(xIx)]] <- lValue
    }else{
      if(as.character(xIx) %in% ls(ruleleft)) {
        if(!is.null(levelX[[xIx]])){
          lValue <- intersect(ruleleft[[as.character(xIx)]],lValue)
          ruleleft[[as.character(xIx)]] <- lValue
        }else{
          ruleleft[[as.character(xIx)]] <- paste(ruleleft[[as.character(xIx)]], "&", lValue)
        }
      }else{
        ruleleft[[as.character(xIx)]] <- lValue
      }
    }

    #thisItem = paste("X[,",xIx, "] %in% ", nxValue, sep="")
    ruleright <- rule
    if(length(ruleright)==0)
    {
      ruleright[[as.character(xIx)]] <- rValue
    }else{
      if(as.character(xIx) %in% ls(ruleright)) {
        if(!is.null(levelX[[xIx]])){
          rValue <- intersect(ruleright[[as.character(xIx)]],rValue)
          ruleright[[as.character(xIx)]] <- rValue
        }else{
          ruleright[[as.character(xIx)]] <- paste(ruleright[[as.character(xIx)]], "&", rValue)
        }
      }else{
        ruleright[[as.character(xIx)]] <- rValue
      }
    }

    thisList = inTrees_treeVisit(tree, tree[rowIx,"left daughter"],count,ruleSet,ruleleft,levelX,length+1,max_length,digits)
    ruleSet = thisList$ruleSet; count = thisList$count

    thisList = inTrees_treeVisit(tree, tree[rowIx,"right daughter"],count,ruleSet,ruleright,levelX,length+1,max_length,digits)
    ruleSet = thisList$ruleSet; count = thisList$count

    return(list(ruleSet = ruleSet, count=count))
  }


measureRule <-
  function(ruleExec,X,target,pred=NULL,regMethod="mean"){
    len <- length(unlist(strsplit(ruleExec, split=" & ")))
    origRule <- ruleExec
    ruleExec <- paste("which(", ruleExec, ")")
    ixMatch <- eval(parse(text=ruleExec))
    if(length(ixMatch)==0){
      v <- c("-1","-1", "-1", "", "")
      names(v) <- c("len","freq","err","condition","pred")
      return(v)
    }
    ys <- target[ixMatch]
    freq <- round(length(ys)/nrow(X),digits=3)

    if(is.numeric(target))
    {
      if(regMethod == "median"){
        ysMost = median(ys)
      }else{
        ysMost <- mean(ys)
      }
      err <- sum((ysMost - ys)^2)/length(ys)
    }else{
      if(length(pred)>0){ #if pred of the rule is provided
        ysMost = as.character(pred)
      }else{
        ysMost <- names(which.max(  table(ys))) # get back the first max
      }
      ly <- sum(as.character(ys)==ysMost)
      conf <- round(ly/length(ys),digits=3)
      err <- 1 - conf
    }
    rule <- origRule

    v <- c(len, freq, err, rule, ysMost)
    names(v) <- c("len","freq","err","condition","pred")
    return(v)
  }


rule2Table <-
  function(ruleExec,X,target){
    I <- rep(0,nrow(X))
    ruleExec <- paste("which(", ruleExec, ")")
    ixMatch <- eval(parse(text=ruleExec))
    if(length(ixMatch)>0) I[ixMatch] <- 1
    names(I) = NULL
    return(I)
  }

inTrees_RF2List <-
  function(rf){
    treeList <- NULL
    treeList$ntree <- rf$ntree
    treeList$list <- vector("list",rf$ntree)
    for(i in 1:treeList$ntree){
      treeList$list[[i]] <- getTree(rf,i,labelVar=FALSE)
    }
    return(treeList)
  }

sortRule <-
  function(M,decreasing=TRUE){
    qIx = order((1- as.numeric(M[,"err"])),
                as.numeric(M[,"freq"]),
                -as.numeric(M[,"len"]),
                decreasing=decreasing)
    return(M[qIx,])
  }


dataSimulate <-
  function(flag=1,nCol=20,nRow=1000){

    if(nCol<=2) stop("nCol must be >= 2.")
    #only the first and the second features are needed
    X <- matrix(runif(nRow*nCol, min=-2, max=2), ncol=nCol)
    target <- rep(-1,nRow)

    #linear case
    if(flag == 3) {
      target <- (X[,1]) + (X[,2])
      ix <- which(target>quantile(target, 1/2));
      target <- target*0-1;
      target[ix] <- 1
    }

    #nonlinear case
    if(flag == 2){
      target <- (X[,1])^2 + 1*(X[,2])^2
      ix <- which(target>quantile(target, 6/10));
      ix <- c(ix,which(target<quantile(target, 1/10)));
      target <- target*0-1;
      target[ix] <- 1
    }

    # team optimization
    if(flag == 1){
      X <- matrix(0,nRow,nCol)
      for(ii in 1:nrow(X)){
        ix <- sample(1:nCol,nCol/2,replace=FALSE)
        X[ii,ix] <- 1
      }
      target <- (xor(X[,1],X[,2]))
      repStr <- function(v){v[v=="1"] <- "Y";v[v=="0"] <- "N";return(v)}
      X <- data.frame(apply(X,2,repStr))
      target[target == FALSE] <- "lose"
      target[target == TRUE] <- "win"
      target <- as.factor(target)

      # X <- data.frame(X)
      # for(jj in 1:ncol(X)){
      #  X[,jj] <- as.factor(X[,jj])
      # }
    }
    return(list(X=X,target=target))
  }


dicretizeVector <-
  function(v,K=3){
    splitV <- quantile(v, probs = seq(0, 1, 1/K), na.rm = FALSE,
                       names = TRUE, type = 3)
    splitV <- splitV[-c(1,length(splitV))]

    numSplit <- length(splitV)  # split points + 1
    if(numSplit==0) return(v)
    newV <- vector("character", length(v))
    newV[which(v<=splitV[1])] = paste("L1",sep="")
    if(numSplit>=2){
      for(jj in 2:numSplit){
        newV[which(  v> splitV[jj-1] & v<=splitV[jj]) ] = paste("L",jj,sep="")
      }
    }
    newV[which( v> splitV[numSplit] ) ] =  paste("L",(numSplit+1),sep="")
    return(newV)
  }

getTypeX <-
  function(X){
    typeX = rep(0,ncol(X))
    for(i in 1:ncol(X)){ #numeric: 1; categorical: 2s
      if(is.numeric(X[,i])){ typeX[i] = 1 }else{
        typeX[i] = 2
      }
    }
    return(typeX)
  }

inTrees_pruneRule <-
  function(rules,X,target, maxDecay = 0.05, typeDecay = 2){
    newRuleMetric <- NULL
    for(i in 1:nrow(rules)){
      newRuleMetric <- rbind(newRuleMetric, pruneSingleRule(rules[i,],X,target, maxDecay, typeDecay))
    }
    return(newRuleMetric)
  }


XGB2List<-
  function(xgb, X)
  {
    feature_names <- colnames(X)
    xt <- xgb.model.dt.tree(feature_names = as.character(1:length(feature_names)), model=xgb)
    # avoid cran note: no visible binding for global variable
    Feature=Split=Yes=No=MissingNode=Missing=Weight=Cover=Prediction=Quality=Node=NULL
    xt[Feature == 'Leaf', Feature := '-1']
    xt[, 'split var' := as.integer(Feature)]
    xt[, 'split point' := Split]
    xt[, 'left daughter' := as.integer(tstrsplit(Yes, '-')[[2]]) + 1]
    xt[, 'right daughter' := as.integer(tstrsplit(No, '-')[[2]]) + 1]
    xt[, MissingNode := as.integer(tstrsplit(Missing, '-')[[2]]) + 1]
    xt[, Weight := Cover]
    xt[, Prediction := Quality]
    xt[, Node := Node + 1]
    xt[, c('ID', 'Yes', 'No', 'Split','Missing', 'Quality', 'Cover', 'Feature') := NULL]
    for (f in c('left daughter', 'right daughter', 'MissingNode'))
      set(xt, which(is.na(xt[[f]])), f, -1)
    treeList <- NULL
    treeList$ntree <- length(unique(xt$Tree))
    treeList$list <- split(xt, by="Tree")
    formatXGB <-
      function(tree){
        rownames(tree) <- 1:nrow(tree)
        tree$status <- ifelse(tree$`split var`==-1,-1,1)
        tree$`split point` <- as.numeric(tree$`split point`)
        tree <- tree[,c("left daughter","right daughter","MissingNode","split var","split point","status")]
        # ix <- tree$MissingNode[which(tree$MissingNode>0)]
        # if(length(ix)>0)  tree$status[ix] <- 10 #missing
        tree <- tree[,c("left daughter","right daughter","split var","split point","status")]
        tree <- as.data.frame(tree)
      }
    treeList$list <- lapply(treeList$list,formatXGB)
    return(treeList)
  }

singleRuleList2Exec <-
  function(ruleList,typeX){ #numeric: 1; categorical: 2s
    #ruleExec <- "which("
    ruleExec <- ""
    vars <- ls(ruleList)
    #ruleL <- length(unique(vars))
    vars <- vars[order(as.numeric(vars))]
    for(i in 1:length(vars)){
      if(typeX[as.numeric(vars[i])]==2){
        values <- paste("c(",paste(  paste("'",ruleList[[vars[i]]],"'",sep="")    ,collapse=","),")",sep="")
        tmp = paste("X[,",vars[i], "] %in% ", values, sep="")
      }else{
        tmp = ruleList[[vars[i]]]
      }
      if(i==1)ruleExec <- paste(ruleExec, tmp,sep="")
      if(i>1)ruleExec <- paste(ruleExec, " & ", tmp, sep="")
    }
    #ruleExec <- paste(ruleExec,")",sep="")
    return(c(ruleExec))
  }


Num2Level <-
  function(rfList,splitV){
    for(i in 1:rfList$ntree){
      rfList$list[[i]] <- data.frame(rfList$list[[i]])
      rfList$list[[i]][,"prediction"] <- data.frame(dicretizeVector(rfList$list[[i]][,"prediction"],splitV))
      colnames(rfList$list[[i]]) <- c("left daughter","right daughter","split var","split point","status","prediction")
    }
    return(rfList)
  }


applyLearner <-
  function(learner,X){
    leftIx <- 1:nrow(X)
    predY <- rep("",nrow(X))
    for(i in 1:nrow(learner)){
      ixMatch <- eval(parse(text=paste("which(",learner[i,"condition"], ")"))  )
      ixMatch <- intersect(leftIx,ixMatch)
      if(length(ixMatch)>0){
        predY[ixMatch] <- learner[i,"pred"]
        leftIx <- setdiff(leftIx,ixMatch)
      }
      if(length(leftIx)==0){
        break
      }
    }
    return(predY)
  }
presentRules <-
  function(rules,colN,digits=NULL){
    for(i in 1:nrow(rules[,"condition",drop=FALSE])){
      if(is.numeric(digits)){
        digits <- as.integer(abs(digits))
        rules[,"freq"] <- round(as.numeric( rules[,"freq"]),digits=digits)
        rules[,"err"] <- round(as.numeric( rules[,"err"]),digits=digits)
      }
      A <- regexpr("X\\[,1\\]==X\\[,1\\]", rules[i,"condition"])
      thisPos <- as.numeric(A[[1]])
      thisLen <- attr(A, "match.length")
      if(thisPos > 0){
        origStr <- substr(rules[i,"condition"], thisPos, thisPos+thisLen-1)
        rules[i,"condition"] <- gsub(origStr, "Else", rules[i,"condition"], fixed=TRUE)
      }
      while(TRUE){
        A <- regexpr("X\\[,[0-9]+\\]", rules[i,"condition"])
        thisPos <- as.numeric(A[[1]])
        thisLen <- attr(A, "match.length")
        if(thisPos <= 0) break
        origStr <- substr(rules[i,"condition"], thisPos, thisPos+thisLen-1)
        ix <- as.numeric(gsub("\\D", "", origStr))
        colStr <- colN[ix]
        rules[i,"condition"] <- gsub(origStr, colStr, rules[i,"condition"], fixed=TRUE)
      }
    }
    return(rules)
  }


buildLearner <-
  function(ruleMetric,X,target,minFreq=0.01){ #Recursive
    ruleMetric <- ruleMetric[,c("len","freq","err","condition","pred"),drop=FALSE]
    learner <- NULL
    listIxInst <- vector("list", nrow(ruleMetric))
    for(i in 1:nrow(ruleMetric)){
      ixMatch <- eval(parse(text=paste("which(",ruleMetric[i,"condition"], ")"))  )
      if(length(ixMatch)==0)next
      listIxInst[[i]] = ixMatch
    }
    ixInstLeft <- 1:length(target)
    while(TRUE){
      infor = NULL
      restErr <- 1 - max(table(target[ixInstLeft]))/length(target[ixInstLeft])
      for(i in 1:length(listIxInst)){
        thisInfor <- computeRuleInfor(listIxInst[[i]], ruleMetric[i,"pred"],target)
        infor <- rbind(infor,c(thisInfor,len=as.numeric(ruleMetric[i,"len"])))
      }
      topIx <- order(infor[,"err"],-infor[,"freq"],infor[,"len"],decreasing=FALSE)
      minSupIx <- which(infor[,"freq"] < minFreq)
      if(length(minSupIx)>0)topIx <- setdiff(topIx,minSupIx)
      if(length(topIx)>0) topIx <- topIx[1]
      if(length(topIx)==0){
        restCondition <- paste("X[,1]==X[,1]")
        restPred <- names(table(target[ixInstLeft]))[which.max(table(target[ixInstLeft]))]
        restSup <- length(ixInstLeft)/length(target)
        thisRuleMetric <- c(len=1,freq=restSup,err=restErr,condition=restCondition,pred=restPred)
        learner <- rbind(learner,thisRuleMetric)
        break
      }else if( infor[topIx,"err"] >= restErr ){
        restCondition <- paste("X[,1]==X[,1]")
        restPred <- names(table(target[ixInstLeft]))[which.max(table(target[ixInstLeft]))]
        restSup <- length(ixInstLeft)/length(target)
        thisRuleMetric <- c(len=1,freq=restSup,err=restErr,condition=restCondition,pred=restPred)
        learner <- rbind(learner,thisRuleMetric)
        break
      }
      #ruleActiveList <- c(ruleActiveList,topIx)
      thisRuleMetric <- ruleMetric[topIx,,drop=FALSE]
      thisRuleMetric[,c("freq","err","len")] <- infor[topIx,c("freq","err","len")]
      learner <- rbind(learner,thisRuleMetric)
      ixInstLeft <- setdiff(ixInstLeft,listIxInst[[topIx]])
      listIxInst <- sapply(listIxInst,setdiff,listIxInst[[topIx]])

      if(length(ixInstLeft)==0) { # if every is targetified perfectly, still set a main target
        restCondition <- paste("X[,1]==X[,1]")
        restPred <- names(table(target))[which.max(table(target))]
        restSup <- 0
        restErr <- 0
        thisRuleMetric <- c(len=1,freq=restSup,err=restErr,condition=restCondition,pred=restPred)
        learner <- rbind(learner,thisRuleMetric)
        break
      }
    }
    rownames(learner) <- NULL
    return(learner)
  }


inTrees_ruleList2Exec <-
  function(X,allRulesList){
    typeX = getTypeX(X)
    ruleExec <- unique(t(sapply(allRulesList,singleRuleList2Exec,typeX=typeX)))
    ruleExec <- t(ruleExec)
    colnames(ruleExec) <- "condition"
    return(ruleExec)
  }



getFreqPattern <-
  function(ruleMetric,minsup=0.01,minconf=0.5,minlen=1,maxlen=4){
    # set up
    predY <- as.character(ruleMetric[,"pred"])
    rulesV <- strsplit(ruleMetric[,"condition"], split=" & ")
    for(i in 1:length(rulesV)){
      rulesV[[i]] = c(rulesV[[i]],paste("=>",predY[i],sep=""))
    }
    yrhs= unique(paste("=>",ruleMetric[,"pred"],sep=""))
    trans1 <- as(rulesV, "transactions")
    rules1 <- arules::apriori(
      trans1,
      parameter = list(supp=minsup,conf=minconf,minlen=minlen,maxlen=maxlen),
      control = list(verbose=FALSE),
      appearance = list(none=NULL,rhs =yrhs,default="lhs")
    )
    #rules1= sort(rules1, decreasing = FALSE, by = "confidence")
    #quality = quality(rules1)
    #qIx = order(quality[,"confidence"],quality[,"support"],decreasing=TRUE)
    #rules1=rules1[qIx]
    #quality = quality[qIx,1:2]
    #inspect(rules1)

    lhs = methods::as(lhs(rules1),"list")
    rhs = methods::as(rhs(rules1),"list")
    rhs <- gsub("=>", "", rhs)
    quality <- quality(rules1)
    ix_empty <- NULL
    freqPattern <- NULL
    for(i in 1:length(lhs)){
      length_v <- length(lhs[[i]])
      lhs[[i]] <- paste(lhs[[i]],collapse= " & ")
      if(nchar(lhs[[i]])==0){
        ix_empty <- c(ix_empty,i)
      }
      freqPattern <- rbind(freqPattern, c(len=length_v, condition=lhs[[i]], pred=rhs[i],
                                          sup=quality[i,"support"],
                                          conf=quality[i,"confidence"]) )
    }
    if(length(ix_empty)>0)freqPattern <- freqPattern[-ix_empty,]
    qIx = order(as.numeric(freqPattern[,"sup"]), as.numeric(freqPattern[,"conf"]),
                -as.numeric(freqPattern[,"len"]),
                decreasing=TRUE)
    freqPattern <- freqPattern[qIx,c("len","sup","conf","condition","pred")]
    freqPattern[,c("sup","conf")] <- as.character(round(as.numeric(freqPattern[,c("sup","conf")]),digits=3))
    return(freqPattern)
  }


GBM2List <-
  function(gbm1,X){
    treeList <- NULL
    treeList$ntree <- gbm1$n.trees
    treeList$list <- vector("list",gbm1$n.trees)
    for(i in 1:treeList$ntree){
      treeList$list[[i]] <- pretty.gbm.tree(gbm1,i.tree = i)
    }

    v2int <- function(v){sum( (-v+1)/2 * 2^seq(0,(length(v)-1),1)  )}
    #as.integer(intToBits(3616)) pretty.gbm.tree(gbm1,i.tree = 1)
    splitBin = sapply(gbm1$c.splits,v2int)

    return(formatGBM(treeList,splitBin,X))
  }


formatGBM <-
  function(gbmList,splitBin,X){
    for(j in 1:length(gbmList$list)){
      a <- gbmList$list[[j]]
      rownames(a) <- 1:nrow(a)
      a$status <- a$SplitVar
      a <- a[,c("LeftNode","RightNode","MissingNode","SplitVar","SplitCodePred","status")]
      a[which(a[,"SplitVar"]>=0),c("SplitVar","LeftNode","RightNode","MissingNode")] <- a[which(a[,"SplitVar"]>=0),c("SplitVar","LeftNode","RightNode","MissingNode")] + 1
      ix <- a$MissingNode[which(a$MissingNode>0)]
      if(length(ix)>0)  a$status[ix] <- 10 #missing #a <- a[-ix,]
      a <- a[,c("LeftNode","RightNode","SplitVar","SplitCodePred","status")]
      cat <- which(sapply(X, is.factor) & !sapply(X, is.ordered))
      ix <- which(a[,"SplitVar"] %in% cat)

      for(i in ix) a[i,"SplitCodePred"] <- splitBin[ a[i,"SplitCodePred"]+1 ]
      colnames(a) <- c("left daughter","right daughter","split var","split point","status")
      gbmList$list[[j]] <- a
    }
    return(gbmList)
  }


inTrees_getRuleMetric <-
  function(ruleExec, X, target){
    #typeX = getTypeX(X)
    #ruleExec <- unique(t(sapply(allRulesList,RuleList2Exec,typeX=typeX)))
    #colnames(ruleExec) <- c("len","condition")
    ruleMetric <- t(sapply(ruleExec[,"condition",drop=FALSE],measureRule,X,target))
    rownames(ruleMetric) = NULL;
    # ruleMetric <- cbind( ruleExec[,1] ,  ruleMetric )
    colnames(ruleMetric) <- c("len","freq","err","condition","pred")
    dIx <- which(ruleMetric[,"len"]=="-1")
    if(length(dIx)>0){
      ruleMetric <- ruleMetric[-dIx,]
      print(paste( length(dIx)," paths are ignored.",sep=""))
    }
    return(ruleMetric)
    #qIx = order((1- as.numeric(ruleMetric[,"err"])),
    #            as.numeric(ruleMetric[,"freq"]),
    #            -as.numeric(ruleMetric[,"len"]),
    #            decreasing=TRUE)
    #return(ruleMetric[qIx,])
  }


selectRuleRRF <-
  function(ruleMetric,X,target){
    ruleI = sapply(ruleMetric[,"condition"],rule2Table,X,target)
    coefReg <- 0.95 - 0.01*as.numeric(ruleMetric[,"len"])/max(as.numeric(ruleMetric[,"len"]))
    rf <- RRF(ruleI,as.factor(target), flagReg = 1, coefReg=coefReg, mtry = (ncol(ruleI)*1/2) , ntree=50, maxnodes= 10,replace=FALSE)
    imp <- rf$importance/max(rf$importance)
    feaSet <- which(imp > 0.01)
    ruleSetPrunedRRF <- cbind(ruleMetric[feaSet,,drop=FALSE],impRRF=imp[feaSet])
    ix = order(as.numeric(ruleSetPrunedRRF[,"impRRF"]),
               - as.numeric(ruleSetPrunedRRF[,"err"]),
               - as.numeric(ruleSetPrunedRRF[,"len"]),
               decreasing=TRUE)
    ruleSelect <- ruleSetPrunedRRF[ix,,drop=FALSE]
    return(ruleSelect)
  }


extractRules <-
  function(treeList,X,ntree=100,maxdepth=6,random=FALSE,digits=NULL){
    if(is.numeric(digits)) digits <- as.integer(abs(digits))

    levelX = list()
    for(iX in 1:ncol(X))
      levelX <- c(levelX,list(levels(X[,iX])))
    # X <- NULL; target <- NULL
    ntree=min(treeList$ntree,ntree)
    allRulesList = list()
    for(iTree in 1:ntree){
      if(random==TRUE){max_length = sample(1:maxdepth,1,replace=FALSE)}else{
        max_length = maxdepth}
      rule = list(); count = 0; rowIx = 1;
      # tree = getTree(rf,iTree,labelVar=FALSE)
      tree <- treeList$list[[iTree]]
      if(nrow(tree)<=1) next # skip if there is no split
      ruleSet = vector("list", length(which(tree[,"status"]==-1)))
      res = inTrees_treeVisit(tree,rowIx = rowIx,count,ruleSet,rule,levelX,length=0,max_length=max_length,digits=digits)
      allRulesList = c(allRulesList, res$ruleSet)
    }
    allRulesList <- allRulesList[!unlist(lapply(allRulesList, is.null))]
    cat(paste(length(allRulesList)," rules (length<=",
              max_length, ") were extracted from the first ", ntree," trees.","\n",sep=""))

    rulesExec <- inTrees_ruleList2Exec(X,allRulesList)
    return(rulesExec)
  }


computeRuleInfor <-
  function(instIx,pred,target){
    trueCls <- as.character(target[instIx])
    err <- 1- length(which(trueCls == pred))/length(trueCls)
    return(c(err=err,freq=length(instIx)/length(target)))
  }


lookupRule <-
  function(rules,strList){
    ix <- grep(strList[1], rules[,"condition"],fixed = TRUE)
    if(length(strList)>=2){
      for(i in 2:length(strList)){
        ix2 <- grep(strList[i], rules[,"condition"],fixed = TRUE)
        ix <- intersect(ix,ix2)
      }
    }
    if(length(ix)>=1)return(rules[ix,,drop=FALSE])
    if(length(ix)==0)return(NULL)
  }


pruneSingleRule <-
  function(rule, X, target, maxDecay, typeDecay){
    # typeDecay = 1: relative error increase; otherwise: absolute error increase

    #A <- gregexpr("X\\[,[0-9]+\\]", s)
    newRuleMetric <- measureRule(rule["condition"],X,target)
    errOrig <- as.numeric(newRuleMetric["err"])
    ruleV <- unlist(strsplit(rule["condition"],split= " & "))
    pred <- rule["pred"]
    # newRule <- NULL
    if(length(ruleV)==1) return(newRuleMetric)
    for(i in length(ruleV):1){
      restRule <- ruleV[-i]
      restRule <- paste(restRule,collapse= " & ")
      metricTmp <- measureRule(restRule,X,target,pred)
      errNew <- as.numeric(metricTmp["err"])
      if(typeDecay == 1){
        decay <- (errNew-errOrig)/max(errOrig,0.000001)
      }else{
        decay <- (errNew-errOrig)
      }
      if( decay <= maxDecay){
        #if( errNew-errOrig <= maxDecay){
        ruleV <- ruleV[-i]
        # newRule saves the last changed rule and metrics
        newRuleMetric <- metricTmp
        if(length(ruleV)<=1)break
      }
    }
    return(newRuleMetric)
    #rule["condition"] <- paste(ruleV,collapse= " & ")
    #return(rule)
  }

