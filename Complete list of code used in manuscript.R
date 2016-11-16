###########################
######  Loading Data  ########
###########################

load("human.data.rda") # Object H. Identifier data of 126 study participants, includes: age, smoking, BMI, WC, 
	# serum markers
dim(H) # 126 rows 
load("IS mean imputed.rda") # Data matrix of internal standard normalized mean imputed lipidomic data
dim(IS.mean) # 126 rows and 1745 lipids
all.equal.character(rownames(H), rownames(IS.mean)) #The identifier data for each patient is aligned with their 
# lipidomic data

#############################
######  Scaling of Data  ########
############################

#### Pareto Scaling of Internal Standard Normalize Mean imputed Data  #####
  install.packages("MetabolAnalyze")
  library(MetabolAnalyze) # Package used for pareto scaling. 
      # https://cran.r-project.org/web/packages/MetabolAnalyze/index.html
  pareto.IS = scaling(as.matrix(scale(IS.mean, center=T, scale=F)), type = "pareto") # Center and pareto scale the internal standard normalized mean imputed data.
  dim(pareto.IS) # 126 rows and 1745 lipids
  all.equal.character(rownames(H), rownames(pareto.IS)) # The identifier data for each patient is aligned with their 
	# scaled lipidomic data


###################################
###### ANOVA and Post Hoc  ########
###################################

#### Kruskall Wallis (KW) one-way ANOVA was conducted across BMI categories along
  # with Dunn's test for multiple comparison  #####

### KW one-way ANOVA ###
	kwpv=numeric() # Numeric vector for indexing KW ANOVA p-values
	## Loop for KW ANOVA of age, BMI, WC, and serum adipokines, cytokines, and C-peptide
  # across BMI category

	for (i in 1:ncol(H)){
	  kwpv[i]=kruskal.test(x=H[,i], g=H$BMI.cat)$p.value # Indexing the p-value from each ANOVA
	}

	kwbh=p.adjust(kwpv, method = "BH") # Adjusting KW p-values according to Benjamini-Hochberg

### Dunn non-parametric test for multiple comparison ###
library(PMCMR)
	tmp=matrix(NA, ncol = 4, nrow = ncol(H))
	colnames(tmp)=c("var","lean-overwt","lean-obese", "overwt-obese" )
	tmp[,1]=colnames(H)

for (i in 1:ncol(H)){
  tmpindex=posthoc.kruskal.dunn.test(x=H[,i], g=H$BMI.cat, p.adjust.method = "BH")
  tmpgetp=get.pvalues(tmpindex)
  tmp[i,2]=tmpgetp[[1]]
  tmp[i,3]=tmpgetp[[2]]
  tmp[i,4]=tmpgetp[[3]]
}

######################################
######  Bayesian Analysis of Traits  ########
######################################

#### Reproducing Kernel Hilbert Spaces (RKHS) Regressions ####
 install.packages("BGLR")
 library(BGLR)
  # For more information on using BGLR refer to https://github.com/gdlc/BGLR-R
 
 ### Creating a directory to store file outputs from BGLR ###
 dir.create("BMI"); setwd("BMI")
 getwd() # Verify working directory 

 ### Setting response, lipidomic data, iterations, and burn in  ###
 X=pareto.IS # 1,745 pareto scaled lipidomic data
 y<-H$BMI # Response: body mass index values for each respective patient
 nIter=200000 # Long Markov Chain of 200,000 iterations
 burnIn=50000 # Number of iterations to discard for burn in
 
     ## Computing the metabolomic similarity matrix ##
     L<-sum(apply(X=X,FUN=var,MARGIN=2))
     G<-tcrossprod(X)/L 
      # The G matrix is an nxn matrix of distances to measure similarities
        # between participants with respect to their lipid profiles
     
     ##  RKHS model parameters  ## 
      # ETA #
      ETA.FixMet=list(Met=list(K=G,model="RKHS"), # The G matrix represents the lipidome, 
                      # and the RKHS kernel is specified
        Fix=list(~H$age+factor(H$smoking), model="FIXED")) # The fixed effects of the model are age and smoking
     
     ## RKHS regression ##
     fmGBLUP<-BGLR(y=y,ETA=ETA.FixMet, nIter=nIter, burnIn=burnIn, saveAt="GBLUP_")
        # RKHS model is: BMI = fixed effects + lipidomic data
     
 ### Variance of Best Linear Unbiased Predictor (BLUP) and variance of error ###
   # inference was done based on one of every 5 samples of the last 150,000
    # therefore, since the burn in was 50,000 we need to remove the first 10,000 samples
 list.files()
 VarU=scan("GBLUP_ETA_Met_varU.dat") #load in variance of the lipidome BLUP model (varU). Total of 30,000 	# lines.
 VarU=VarU[-c(1:10000)] # remove the burn in iterations from varU
 VarE=scan("GBLUP_varE.dat") # load in variance of the error (varE). Total of 30,000 lines.
 VarE=VarE[-c(1:10000)] #remove the burn in iterations varE
     
   ## Calculating percent of the inter-individual differences in response variables that can be 
     # attributed to lipidome profiles ##
     tmp=(VarU/(VarU+VarE))*100 # Calculate % varU as sum of varU + varE
     round(quantile(tmp, probs = c(0.05, 0.95)),0) # Determine 95% confidence intervals
     
###### The other responses were run similarly by creating individual directories for each response and modifying 
	# the y<- for each individual response  #######


################################
#####  Single Lipid Regressions #####
###############################

### Setting the model response and independent variables ###
	X=pareto.IS # Data set
	XX=H$BMI # Response, dependent variable  

### Creating numeric vectors for indexing estimated effect and model p-value ###
	BMI.mypv=numeric(ncol(X)) # Numeric variable for indexing each p-value for each lipid
	BMI.myB=numeric(ncol(X)) #Beta coefficient for each lipid

### Loop of single marker regressions of pareto scaled data ###
	for (i in 1:ncol(X)){ # For each column of data matrix
 	 tmp=summary(lm(XX~ H$age + H$smoking + X[,i])) # list of summary of each linear model 
  	BMI.mypv[i]=tmp$coefficients[5,4]  # Indexing p-value 
  	BMI.myB[i]=tmp$coefficients[5,1]  # Indexing beta coefficent
	}

###  P-value adjustments  ###
	##  Bonferroni pvalue correction  ##
		p.bon=p.adjust(BMI.mypv, method='bonferroni')

	## Benjamini-Hochberg pvalue correction ##
		p.BH=p.adjust(BMI.mypv, method='BH')

### Creating a dataframe of regression results ###
        Complete.Regression.Report=data.frame('Metabolite'=colnames(X), BMI.stats.reportBMI.myB, BMI.mypv, 'BMI.BH'= p.BH,  'BMI.bon'=p.bon) # Compile all Response beta coefficients, p-values, and corrected p-values 

###### The other responses were run similarly by modifying the response input (i.e., either BMI, WC, leptin, etc)  
# and uniquely indexing model outputs and compiling their results into the file Complete.Regression.Report  ######


#####################################
### Singular Value Decomposition  ###
#####################################
load("beta.report.rda")
# Beta.Report is a dataframe with 1745 rows and 10 columns
  # Column 1 contains the name of each of the 1745 lipids
    # Columns 2-10 contain the 1745 estimated effects for each respective response (i.e., BMI, WC, 
    	# leptin, etc)

# Singular value decompositions (SVD) is defined as
  # A= U*D*V
    # where U is an mxm matrix, D is an mxn, and V is an nxn matrix

### Singular value decomposition of traits ###
pca.beta=princomp(Beta.Report[,2:ncol(Beta.Report)])
lambda=pca.beta$sdev*sqrt(pca.beta$n.obs)
score=t(t(pca.beta$scores)/lambda)
variables=t(t(pca.beta$loadings)*lambda)

v.1=variables
v.1[8,1]=v.1[8,1]/25 #scaling MCP-1 Pc1 score by a factor of 25


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
