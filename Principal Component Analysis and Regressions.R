##########################################################
#############  Principal Components Analysis #############    
##########################################################

#### Setting response and computing principal components (PCs) ####
  XX=H$BMI # Response, dependent variable to be used in regressions
  Z=tcrossprod(pareto.IS) # Computing an nxn matrix of lipidomic data
  ZZ=eigen(Z) # Computing Eigen values and Eigen vectors of nxn matrix of lipidomic data
  EV=(ZZ$values/sum(ZZ$values))*100 # Calculating the Eigen value % of total variation


### Creating vectors to index pvalue and beta coefficient
 PCA.mypv=numeric(ncol(Z[,1:125])) #p-value for each metabolite
 PCA.myB=numeric(ncol(Z[,1:125])) #Beta coefficent for each metabolite

### Loop of BMI regressed on each PC individually ### 
for (i in 1:125){ # For each PC
  tmp=summary(lm(XX~ H$age + H$smoking + ZZ$vectors[,i])) # List of outputs from linear model
  PCA.mypv[i]=tmp$coefficients[4,4] # Indexing the p-value from linear model
  PCA.myB[i]=tmp$coefficients[4,1] # Indexing the beta coefficients from linear model
}

###  P-value adjustments  ###
  ##  Bonferroni pvalue correction  ##
    PCA.pbon=p.adjust(PCA.mypv, method='bonferroni')


  ## Benjamini-Hochberg pvalue correction ##
    PCA.pBH=p.adjust(PCA.mypv, method='BH')

    
####################################################
##### Determining lipids driving PC4 loadings  #####
####################################################
    
#### Setting model response and creating vectors for indexing  ####
  Y=ZZ$vectors[,4] # Response, dependent variable is PC4 loadings
  b=numeric() # Numeric vector for indexing beta coefficients
  pval=numeric() # Numeric vector for indexing p-values

  ###  Loop of PC4 scores regressed on each lipid individually  ###
    for (i in 1:ncol(pareto.IS)){ # For each of the 1745 lipids
      tmp=summary(lm(Y~pareto.IS[,i])) # List of summary of each linear model
      b[i]=tmp$coefficients[2,1] # Indexing beta coefficients
      pval[i]=tmp$coefficients[2,4] # Indexing p-values
    }
            
  ###  P-value adjustments  ###
    ##  Bonferroni pvalue correction  ##
      p.bon=p.adjust(pval, method='bonferroni')
  
    ## Benjamini-Hochberg pvalue correction ##
      p.BH=p.adjust(pval, method='BH')
