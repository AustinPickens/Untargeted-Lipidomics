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
