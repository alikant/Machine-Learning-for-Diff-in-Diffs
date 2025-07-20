################################################################
# File:        Masterthesis_Economics_Bonn.R
# Author:      Ali Kananitarigh
# Description: This code is part of a Master's thesis for 
#              the M.Sc. in Economics program, written under
#              the supervision of Prof. Dr. Christoph Breunig
# Submission:  July 2025
################################################################
# Clean the environment
rm(list=ls())

# Load libraries
if (!require("DRDID")) install.packages("DRDID")
if (!require("MASS")) install.packages("MASS")
if (!require("Matrix")) install.packages("Matrix")
if (!require("DoubleML")) install.packages("DoubleML")
if (!require("rpart")) install.packages("rpart")
if (!require("neuralnet")) install.packages("neuralnet")
if (!require("glmnet")) install.packages("glmnet")
if (!require("mlr3")) install.packages("mlr3")
if (!require("mlr3learners")) install.packages("mlr3learners")
if (!require("future")) install.packages("future")
if (!require("future.apply")) install.packages("future.apply")

# Parallel processing
if (parallel::detectCores() != 1) {
    # Exclude one core
    plan(multisession, workers = parallel::detectCores() - 1)
} else {
    message('Parallel processing is not available.')
}

################################################################
## Function of Data Generating Process 
################################################################

dgps <- function(sample.size = N, parameters = p, MC.sim=B, design=design) {
    
    N <- sample.size # Number of observations
    p <- parameters  # Number of parameters
    B <- MC.sim      # Number of Monte-Carlo simulations
    
    # Mean vector for covariates
    mean_vec <- (-1)^(0:(p-1))
    
    # Covariance matrix
    cov_matrix <- matrix(0, nrow = p, ncol = p) 
    for (j in 1:p) {
        for (k in 1:p) {
            cov_matrix[j, k] <- 0.5^abs(j - k)  
        }
    }
    
    # Generate covariates 
    X <- array(0, dim = c(B, N, p))
    # Set seed for reproducible random number generation
    set.seed(123)     
    for (b in 1:B) {
        X[b,,] <- MASS::mvrnorm(N, mu = mean_vec, Sigma = cov_matrix)
    }
    
    # Storage for nuisance functions
    mu <- matrix(0, B, N)
    g <- matrix(0, B, N)
    
    # Designs I, II, III, and IV
    for (i in 1:B) {
        for (j in 1:N) {
            if (design == 1) {
                g[i,j] <- 0.5 * sum(X[i,j,]/(1:p))
                mu[i,j] <- sum(X[i,j,]/(1:p))
            }
            else if (design == 2) {
                g[i,j] <- 0.5 * sum(X[i,j,]/(1:p))
                mu[i,j] <- 0.8 * sum(X[i,j,]/(1:p)) + 0.2 * sum(X[i,j,]^2/(1:p))
            }
            else if (design == 3) {
                g[i,j] <- (0.5 * sum(X[i,j,]/(1:p)) + 0.5 * sum(X[i,j,]^2/(1:p)))/4
                mu[i,j] <- sum(X[i,j,]/(1:p))
            }
            else if (design == 4) {
                g[i,j] <- (0.5 * sum(X[i,j,]/(1:p)) + 0.5 * sum(X[i,j,]^2/(1:p)))/4
                mu[i,j] <- 0.8 * sum(X[i,j,]/(1:p)) + 0.2 * sum(X[i,j,]^2/(1:p))
            }
        }
    }
    
    # Generate treatment variable
    sig_g = 1
    Nu <- matrix(0, B, N)
    
    D <- matrix(0, B, N)
    for (i in 1:B) {
        # Generate random variable from a logistic distribution
        Nu[i,] <- rlogis(N, location = 0, scale = sig_g)
        D[i,] = as.numeric(g[i,] - Nu[i,] > 0)
    }
    
    
    # Storage for the outcome variable
    Y0 <- matrix(0, B, N)    # Pre-treatment period
    Y1 <- matrix(0, B, N)    # Post-treatment period
    Y10 <- matrix(0, B, N)   # 2nd period outcome (untreated)
    Y11 <- matrix(0, B, N)   # 2nd period outcome (treated)
    trend <- matrix(0, B, N)
    true_ATT <- numeric(B)
    
    # Time trend
    trend <- 2 + mu[i,]
    
    # Generate outcome variable
    for(i in 1:B) {
        # Generate random errors and fixed effects
        random <- MASS::mvrnorm(n = N, mu = rep(0, 4), Sigma = diag(4))
        alpha <- random[,1]
        eps_0 <- random[,2]
        eps_1_0 <- random[,3]
        eps_1_1 <- random[,4]
        
        # Generate time trend
        trend <- 2 + mu[i,]
        mreg <- mu[i,] + trend + D[i,]*0
        
        Y0[i,] <- mu[i,] + D[i,]*mu[i,] + alpha + eps_0 
        Y10[i,] <- mreg + D[i,]*mu[i,] + trend + alpha + eps_1_0 
        Y11[i,] <- mreg + D[i,]*mu[i,] + trend + alpha + eps_1_1
        
        # Observable outcome in 2nd period
        Y1[i,] <- D[i,] * Y11[i,] + (1 - D[i,]) * Y10[i,]
        true_ATT[i] <- mean(D[i,]*(Y11[i,]-Y10[i,]))/mean(D[i,]) 
    }
    
    
    return(list(
        X = X,
        D = D,
        Y0 = Y0,
        Y1 = Y1,
        Y10 = Y10,
        Y11 = Y11,
        true_ATT = true_ATT,
        mu = mu,
        g = g
    ))
}

################################################################
################################################################
## Revised Outcome Regression in DRDID Package
################################################################
################################################################
# The link to the original version of the code is available at:
# https://github.com/pedrohcgs/DRDID/blob/master/R/reg_did_panel.R
# The original function in the DRDID package has has an issue. 
# The revised function, reg_did_panel_revised(), is correct.
# The revised version includes the intercept in the covariate matrix.
# All other parts of the function remain the same as in reg_did_panel().

reg_did_panel_revised <-function(y1, y0, D, covariates, i.weights = NULL,
                                 boot = FALSE, boot.type = "weighted", nboot = NULL,
                                 inffunc = FALSE){
    #-----------------------------------------------------------------------------
    # D as vector
    D <- as.vector(D)
    # Sample size
    n <- length(D)
    # generate deltaY
    deltaY <- as.vector(y1 - y0)
    # Add constant to covariate vector
    #int.cov <- as.matrix(rep(1,n))
    # if (!is.null(covariates)){
    #   if(all(as.matrix(covariates)[,1]==rep(1,n))){
    #     int.cov <- as.matrix(covariates)
    #   } else {
    #     int.cov <- as.matrix(cbind(1, covariates))
    #   }
    # }
    # Add constant to covariate vector
    ##########################################
    ######### Here is the fixed bug ##########
    ######### Add intercept column  ##########
    ######## The rest is the same   ##########
    ##########################################
    if(is.null(covariates)){
        
        int.cov <- as.matrix(rep(1,n))
    } else{
        int.cov <- as.matrix(cbind(1, covariates))
    }
    # Weights
    if(is.null(i.weights)) {
        i.weights <- as.vector(rep(1, n))
    } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
    # Normalize weights
    i.weights <- i.weights/mean(i.weights)
    #-----------------------------------------------------------------------------
    #Compute the Outcome regression for the control group using ols.
    # reg.coeff <- stats::coef(stats::lm(deltaY ~ -1 + int.cov,
    #                                    subset = D==0,
    #                                    weights = i.weights))
    control_filter <- (D == 0)
    reg.coeff <- stats::coef(fastglm::fastglm(
        x = int.cov[control_filter, , drop = FALSE],
        y = deltaY[control_filter],
        weights = i.weights[control_filter],
        family = gaussian(link = "identity")
    ))
    if(anyNA(reg.coeff)){
        stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is probably the reason for it.")
    }
    out.delta <-   as.vector(tcrossprod(reg.coeff, int.cov))
    #-----------------------------------------------------------------------------
    #Compute the OR-DiD estimator
    # First, the weights
    w.treat <- i.weights * D
    w.cont <- i.weights * D
    
    reg.att.treat <- w.treat * deltaY
    reg.att.cont <- w.cont * out.delta
    
    eta.treat <- mean(reg.att.treat) / mean(w.treat)
    eta.cont <- mean(reg.att.cont) / mean(w.cont)
    
    reg.att <-   eta.treat - eta.cont
    #-----------------------------------------------------------------------------
    #get the influence function to compute standard error
    #-----------------------------------------------------------------------------
    # First, the influence function of the nuisance functions
    # Asymptotic linear representation of OLS parameters
    weights.ols <- i.weights * (1 - D)
    wols.x <- weights.ols * int.cov
    wols.eX <- weights.ols * (deltaY - out.delta) * int.cov
    #XpX <- opt_crossprod(wols.x, int.cov, n)
    XpX <- crossprod(wols.x, int.cov)/n
    # Check if XpX is invertible
    if ( base::rcond(XpX) < .Machine$double.eps) {
        stop("The regression design matrix is singular. Consider removing some covariates.")
    }
    XpX.inv <- solve(XpX)
    asy.lin.rep.ols <-  wols.eX %*% XpX.inv
    #-----------------------------------------------------------------------------
    # Now, the influence function of the "treat" component
    # Leading term of the influence function
    inf.treat <- (reg.att.treat - w.treat * eta.treat) / mean(w.treat)
    #-----------------------------------------------------------------------------
    # Now, get the influence function of control component
    # Leading term of the influence function: no estimation effect
    inf.cont.1 <- (reg.att.cont - w.cont * eta.cont)
    # Estimation effect from beta hat (OLS using only controls)
    # Derivative matrix (k x 1 vector)
    M1 <- base::colMeans(w.cont * int.cov)
    # Now get the influence function related to the estimation effect related to beta's
    inf.cont.2 <- asy.lin.rep.ols %*% M1
    # Influence function for the control component
    inf.control <- (inf.cont.1 + inf.cont.2) / mean(w.cont)
    #-----------------------------------------------------------------------------
    #get the influence function of the DR estimator (put all pieces together)
    reg.att.inf.func <- (inf.treat - inf.control)
    #-----------------------------------------------------------------------------
    if (boot == FALSE) {
        # Estimate of standard error
        se.reg.att <- stats::sd(reg.att.inf.func)/sqrt(n)
        # Estimate of upper boudary of 95% CI
        uci <- reg.att + 1.96 * se.reg.att
        # Estimate of lower doundary of 95% CI
        lci <- reg.att - 1.96 * se.reg.att
        #Create this null vector so we can export the bootstrap draws too.
        reg.boot <- NULL
    }
    
    if (boot == TRUE) {
        if (is.null(nboot) == TRUE) nboot = 999
        if(boot.type == "multiplier"){
            # do multiplier bootstrap
            reg.boot <- mboot.did(reg.att.inf.func, nboot)
            # get bootstrap std errors based on IQR
            se.reg.att <- stats::IQR(reg.boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
            # get symmtric critival values
            cv <- stats::quantile(abs(reg.boot/se.reg.att), probs = 0.95)
            # Estimate of upper boudary of 95% CI
            uci <- reg.att + cv * se.reg.att
            # Estimate of lower doundary of 95% CI
            lci <- reg.att - cv * se.reg.att
        } else {
            # do weighted bootstrap
            reg.boot <- unlist(lapply(1:nboot, wboot.reg.panel,
                                      n = n, deltaY = deltaY, D = D, int.cov = int.cov, i.weights = i.weights))
            # get bootstrap std errors based on IQR
            se.reg.att <- stats::IQR((reg.boot - reg.att)) / (stats::qnorm(0.75) - stats::qnorm(0.25))
            # get symmtric critival values
            cv <- stats::quantile(abs((reg.boot - reg.att)/se.reg.att), probs = 0.95)
            # Estimate of upper boudary of 95% CI
            uci <- reg.att + cv * se.reg.att
            # Estimate of lower doundary of 95% CI
            lci <- reg.att - cv * se.reg.att
            
        }
    }
    
    if(inffunc == FALSE) reg.att.inf.func <- NULL
    #---------------------------------------------------------------------
    # record the call
    call.param <- match.call()
    # Record all arguments used in the function
    argu <- mget(names(formals()), sys.frame(sys.nframe()))
    boot.type <- ifelse(argu$boot.type=="multiplier", "multiplier", "weighted")
    boot <- ifelse(argu$boot == TRUE, TRUE, FALSE)
    argu <- list(
        panel = TRUE,
        boot = boot,
        boot.type = boot.type,
        nboot = nboot,
        type = "or"
    )
    ret <- (list(ATT = reg.att,
                 se = se.reg.att,
                 uci = uci,
                 lci = lci,
                 boots = reg.boot,
                 att.inf.func = reg.att.inf.func,
                 call.param = call.param,
                 argu = argu))
    # Define a new class
    class(ret) <- "drdid"
    
    # return the list
    return(ret)
}

################################################################
################################################################
## Monte-Carlo simulations for p = 5
################################################################
################################################################

# Setup parameters
N <- c(500, 1000)           # Sample size
p <- 5                      # Number of parameters
B <- 1000                   # Number of Monte Carlo simulations                
all_results_5 <- list()     # Storage for results
#X_storage <- list().       # Check the same X for each design
start_time5 <- Sys.time()   # Record start time
for (smpl in N) {
    for (design in 1:4) {
        # Storage for this specific config
        att_DR <- numeric(B)      # Improved DR (Sant'Anna & Zhao)
        att_OR <- numeric(B)      # Outcome Regression
        att_IPW <- numeric(B)     # Inverse Probability Weighting (Abadie)
        att_TWFE <- numeric(B)    # Two-way Fixed Effects
        att_DML <- numeric(B)     # Debiased/Double Machine Learning
        
        # Store CIs
        CI_DR <- matrix(0, nrow = B, ncol = 2)
        CI_OR <- matrix(0, nrow = B, ncol = 2)
        CI_IPW <- matrix(0, nrow = B, ncol = 2)
        CI_TWFE <- matrix(0, nrow = B, ncol = 2)
        CI_lower <- numeric(B)    # For DML
        CI_upper <- numeric(B)    # For DML
        se <- numeric(B)          # For DML
        coverage <- numeric(B)    # For DML
        
        # Generate data for this config
        data <- dgps(smpl, p, B, design)
        true_ATT <- data$true_ATT  # True ATT value (zero)
        
        # Monte-Carlo simulations
        for (i in 1:B) {
            # DR estimator
            dr_result <- DRDID::drdid_imp_panel(data$Y1[i, ], data$Y0[i, ], data$D[i, ], data$X[i, , ])
            att_DR[i] <- dr_result$ATT
            CI_DR[i, ] <- c(dr_result$lci, dr_result$uci)
            
            # OR estimator
            or_result <- reg_did_panel_revised(data$Y1[i, ], data$Y0[i, ], data$D[i, ], data$X[i, , ])
            att_OR[i] <- or_result$ATT
            CI_OR[i, ] <- c(or_result$lci, or_result$uci)
            
            # IPW estimator
            ipw_result <- DRDID::std_ipw_did_panel(data$Y1[i, ], data$Y0[i, ], data$D[i, ], data$X[i, , ])
            att_IPW[i] <- ipw_result$ATT
            CI_IPW[i, ] <- c(ipw_result$lci, ipw_result$uci)
            
            # TWFE estimator
            twfe_result <- DRDID::twfe_did_panel(data$Y1[i, ], data$Y0[i, ], data$D[i, ], data$X[i, , ])
            att_TWFE[i] <- twfe_result$ATT
            CI_TWFE[i, ] <- c(twfe_result$lci, twfe_result$uci)
            
            # ML methods
            # Number of covariates in the random forest
            mtry.value<- max(5,floor(sqrt(p)))
            # ML methods for g nuisance (indeed, it is m nuisance in our DGPs)
            # Neural networks
            ml_g <- lrn("regr.nnet", size = 5, decay = 0.01, maxit = 200)
            # Random forest
            #ml_g <- lrn("regr.ranger", num.trees = 300, mtry=mtry.value, min.node.size = 2, max.depth = 5)
            
            # ML methods for m nuisance (indeed, it is g nuisance in our DGPs)
            # LASSO
            ml_m <- lrn("classif.cv_glmnet", s="lambda.min")
            
            # Data preparation for DML (needs a dataframe)
            delta_y <- data$Y1[i,] - data$Y0[i,]
            dml_data <- DoubleML::double_ml_data_from_matrix(
                X = data$X[i,,], 
                y = delta_y, 
                d = data$D[i,]
            )
            
            # DML estimator
            dml_obj <- DoubleML::DoubleMLIRM$new(
                dml_data, 
                ml_g = ml_g, 
                ml_m = ml_m, 
                dml_procedure = "dml2",
                n_folds = 3,
                n_rep = 3,
                trimming_threshold = 0.01
            )
            dml_obj$fit()
            
            # Calculation of DML results
            att_DML[i] <- dml_obj$coef[1]
            se[i] <- dml_obj$se[1]
            CI_lower[i] <- dml_obj$coef[1] - 1.96 * dml_obj$se[1]
            CI_upper[i] <- dml_obj$coef[1] + 1.96 * dml_obj$se[1]
            coverage[i] <- (true_ATT[i] >= CI_lower[i] && true_ATT[i] <= CI_upper[i])
        }
        
        # Store results for this specific config
        key <- paste0("N_", smpl, "_Design_", design)
        #saved_X[[key]] <- data$X
        all_results_5[[key]] <- list(
            DR = c(
                bias = mean(att_DR - true_ATT),
                coverage = mean(CI_DR[, 1] <= true_ATT & true_ATT <= CI_DR[, 2]),
                length = mean(CI_DR[, 2] - CI_DR[, 1])
            ),
            OR = c(
                bias = mean(att_OR - true_ATT),
                coverage = mean(CI_OR[, 1] <= true_ATT & true_ATT <= CI_OR[, 2]),
                length = mean(CI_OR[, 2] - CI_OR[, 1])
            ),
            IPW = c(
                bias = mean(att_IPW - true_ATT),
                coverage = mean(CI_IPW[, 1] <= true_ATT & true_ATT <= CI_IPW[, 2]),
                length = mean(CI_IPW[, 2] - CI_IPW[, 1])
            ),
            TWFE = c(
                bias = mean(att_TWFE - true_ATT),
                coverage = mean(CI_TWFE[, 1] <= true_ATT & true_ATT <= CI_TWFE[, 2]),
                length = mean(CI_TWFE[, 2] - CI_TWFE[, 1])
            ),
            DML = c(
                bias = mean(att_DML - true_ATT),
                coverage = mean(coverage),
                length = mean(CI_upper - CI_lower)
            )
        )
    }
}

end_time5 <- Sys.time()    # Record end time
execution_time5 <- end_time5 - start_time5
# Save the results in a csv file
df_results_5 <- as.data.frame(all_results_5)
df_results_5_transposed <- as.data.frame(t(df_results_5))
write.csv(df_results_5_transposed, file = '/Users/alikananitarigh/simulation_results_5nn.csv', row.names = TRUE)
write.csv(round(df_results_5_transposed,3), file = '/Users/alikananitarigh/simulation_results_5nn_round.csv', row.names = TRUE)
################################################################
################################################################
## Monte-Carlo simulations for p = 10
################################################################
################################################################

# Setup parameters
N <- c(500, 1000)            # Sample size
p <- 10                      # Number of parameters
B <- 1000                    # Number of Monte Carlo simulations                
all_results_10 <- list()     # Storage for results
#X_storage <- list().        # Check the same X for each design
start_time10 <- Sys.time()   # Record start time
for (smpl in N) {
    for (design in 1:4) {
        # Storage for this specific config
        att_DR <- numeric(B)      # Improved DR (Sant'Anna & Zhao)
        att_OR <- numeric(B)      # Outcome Regression
        att_IPW <- numeric(B)     # Inverse Probability Weighting (Abadie)
        att_TWFE <- numeric(B)    # Two-way Fixed Effects
        att_DML <- numeric(B)     # Debiased/Double Machine Learning
        
        # Store CIs
        CI_DR <- matrix(0, nrow = B, ncol = 2)
        CI_OR <- matrix(0, nrow = B, ncol = 2)
        CI_IPW <- matrix(0, nrow = B, ncol = 2)
        CI_TWFE <- matrix(0, nrow = B, ncol = 2)
        CI_lower <- numeric(B)    # For DML
        CI_upper <- numeric(B)    # For DML
        se <- numeric(B)          # For DML
        coverage <- numeric(B)    # For DML
        
        # Generate data for this config
        data <- dgps(smpl, p, B, design)
        true_ATT <- data$true_ATT  # True ATT value (zero)
        
        # Monte-Carlo simulations
        for (i in 1:B) {
            # DR estimator
            dr_result <- DRDID::drdid_imp_panel(data$Y1[i, ], data$Y0[i, ], data$D[i, ], data$X[i, , ])
            att_DR[i] <- dr_result$ATT
            CI_DR[i, ] <- c(dr_result$lci, dr_result$uci)
            
            # OR estimator
            or_result <- reg_did_panel_revised(data$Y1[i, ], data$Y0[i, ], data$D[i, ], data$X[i, , ])
            att_OR[i] <- or_result$ATT
            CI_OR[i, ] <- c(or_result$lci, or_result$uci)
            
            # IPW estimator
            ipw_result <- DRDID::std_ipw_did_panel(data$Y1[i, ], data$Y0[i, ], data$D[i, ], data$X[i, , ])
            att_IPW[i] <- ipw_result$ATT
            CI_IPW[i, ] <- c(ipw_result$lci, ipw_result$uci)
            
            # TWFE estimator
            twfe_result <- DRDID::twfe_did_panel(data$Y1[i, ], data$Y0[i, ], data$D[i, ], data$X[i, , ])
            att_TWFE[i] <- twfe_result$ATT
            CI_TWFE[i, ] <- c(twfe_result$lci, twfe_result$uci)
            
            # ML methods
            # Number of covariates in the random forest
            mtry.value<- max(5,floor(sqrt(p)))
            # ML methods for g nuisance (indeed, it is m nuisance in our DGPs)
            # Neural networks
            ml_g <- lrn("regr.nnet", size = 5, decay = 0.01, maxit = 200)
            # Random forest
            #ml_g <- lrn("regr.ranger", num.trees = 300, mtry=mtry.value, min.node.size = 2, max.depth = 5)
            
            # ML methods for m nuisance (indeed, it is g nuisance in our DGPs)
            # LASSO
            ml_m <- lrn("classif.cv_glmnet", s="lambda.min")
            
            # Data preparation for DML (needs a dataframe)
            delta_y <- data$Y1[i,] - data$Y0[i,]
            dml_data <- DoubleML::double_ml_data_from_matrix(
                X = data$X[i,,], 
                y = delta_y, 
                d = data$D[i,]
            )
            
            # DML estimator
            dml_obj <- DoubleML::DoubleMLIRM$new(
                dml_data, 
                ml_g = ml_g, 
                ml_m = ml_m, 
                dml_procedure = "dml2",
                n_folds = 3,
                n_rep = 3,
                trimming_threshold = 0.01
            )
            dml_obj$fit()
            
            # Calculation of DML results
            att_DML[i] <- dml_obj$coef[1]
            se[i] <- dml_obj$se[1]
            CI_lower[i] <- dml_obj$coef[1] - 1.96 * dml_obj$se[1]
            CI_upper[i] <- dml_obj$coef[1] + 1.96 * dml_obj$se[1]
            coverage[i] <- (true_ATT[i] >= CI_lower[i] && true_ATT[i] <= CI_upper[i])
        }
        
        # Store results for this specific config
        key <- paste0("N_", smpl, "_Design_", design)
        #saved_X[[key]] <- data$X
        all_results_10[[key]] <- list(
            DR = c(
                bias = mean(att_DR - true_ATT),
                coverage = mean(CI_DR[, 1] <= true_ATT & true_ATT <= CI_DR[, 2]),
                length = mean(CI_DR[, 2] - CI_DR[, 1])
            ),
            OR = c(
                bias = mean(att_OR - true_ATT),
                coverage = mean(CI_OR[, 1] <= true_ATT & true_ATT <= CI_OR[, 2]),
                length = mean(CI_OR[, 2] - CI_OR[, 1])
            ),
            IPW = c(
                bias = mean(att_IPW - true_ATT),
                coverage = mean(CI_IPW[, 1] <= true_ATT & true_ATT <= CI_IPW[, 2]),
                length = mean(CI_IPW[, 2] - CI_IPW[, 1])
            ),
            TWFE = c(
                bias = mean(att_TWFE - true_ATT),
                coverage = mean(CI_TWFE[, 1] <= true_ATT & true_ATT <= CI_TWFE[, 2]),
                length = mean(CI_TWFE[, 2] - CI_TWFE[, 1])
            ),
            DML = c(
                bias = mean(att_DML - true_ATT),
                coverage = mean(coverage),
                length = mean(CI_upper - CI_lower)
            )
        )
    }
}

end_time10 <- Sys.time()    # Record end time
execution_time10 <- end_time10 - start_time10
# Save the results in a csv file
df_results_10 <- as.data.frame(all_results_10)
df_results_10_transposed <- as.data.frame(t(df_results_10))
write.csv(df_results_10_transposed, file = '/Users/alikananitarigh/simulation_results_10nn.csv', row.names = TRUE)
write.csv(round(df_results_5_transposed,3), file = '/Users/alikananitarigh/simulation_results_10nn_round.csv', row.names = TRUE)
################################################################
################################################################
## Monte-Carlo simulations for p = 20
################################################################
################################################################

# Setup parameters
N <- c(500, 1000)            # Sample size
p <- 20                      # Number of parameters
B <- 1000                    # Number of Monte Carlo simulations                
all_results_20 <- list()     # Storage for results
#X_storage <- list().        # Check the same X for each design
start_time20 <- Sys.time()   # Record start time
for (smpl in N) {
    for (design in 1:4) {
        # Storage for this specific config
        att_DR <- numeric(B)      # Improved DR (Sant'Anna & Zhao)
        att_OR <- numeric(B)      # Outcome Regression
        att_IPW <- numeric(B)     # Inverse Probability Weighting (Abadie)
        att_TWFE <- numeric(B)    # Two-way Fixed Effects
        att_DML <- numeric(B)     # Debiased/Double Machine Learning
        
        # Store CIs
        CI_DR <- matrix(0, nrow = B, ncol = 2)
        CI_OR <- matrix(0, nrow = B, ncol = 2)
        CI_IPW <- matrix(0, nrow = B, ncol = 2)
        CI_TWFE <- matrix(0, nrow = B, ncol = 2)
        CI_lower <- numeric(B)    # For DML
        CI_upper <- numeric(B)    # For DML
        se <- numeric(B)          # For DML
        coverage <- numeric(B)    # For DML
        
        # Generate data for this config
        data <- dgps(smpl, p, B, design)
        true_ATT <- data$true_ATT  # True ATT value (zero)
        
        # Monte-Carlo simulations
        for (i in 1:B) {
            # DR estimator
            dr_result <- DRDID::drdid_imp_panel(data$Y1[i, ], data$Y0[i, ], data$D[i, ], data$X[i, , ])
            att_DR[i] <- dr_result$ATT
            CI_DR[i, ] <- c(dr_result$lci, dr_result$uci)
            
            # OR estimator
            or_result <- reg_did_panel_revised(data$Y1[i, ], data$Y0[i, ], data$D[i, ], data$X[i, , ])
            att_OR[i] <- or_result$ATT
            CI_OR[i, ] <- c(or_result$lci, or_result$uci)
            
            # IPW estimator
            ipw_result <- DRDID::std_ipw_did_panel(data$Y1[i, ], data$Y0[i, ], data$D[i, ], data$X[i, , ])
            att_IPW[i] <- ipw_result$ATT
            CI_IPW[i, ] <- c(ipw_result$lci, ipw_result$uci)
            
            # TWFE estimator
            twfe_result <- DRDID::twfe_did_panel(data$Y1[i, ], data$Y0[i, ], data$D[i, ], data$X[i, , ])
            att_TWFE[i] <- twfe_result$ATT
            CI_TWFE[i, ] <- c(twfe_result$lci, twfe_result$uci)
            
            # ML methods
            # Number of covariates in the random forest
            mtry.value<- max(5,floor(sqrt(p)))
            # ML methods for g nuisance (indeed, it is m nuisance in our DGPs)
            # Neural networks
            ml_g <- lrn("regr.nnet", size = 5, decay = 0.01, maxit = 200)
            # Random forest
            #ml_g <- lrn("regr.ranger", num.trees = 300, mtry=mtry.value, min.node.size = 2, max.depth = 5)
            
            # ML methods for m nuisance (indeed, it is g nuisance in our DGPs)
            # LASSO
            ml_m <- lrn("classif.cv_glmnet", s="lambda.min")
            
            # Data preparation for DML (needs a dataframe)
            delta_y <- data$Y1[i,] - data$Y0[i,]
            dml_data <- DoubleML::double_ml_data_from_matrix(
                X = data$X[i,,], 
                y = delta_y, 
                d = data$D[i,]
            )
            
            # DML estimator
            dml_obj <- DoubleML::DoubleMLIRM$new(
                dml_data, 
                ml_g = ml_g, 
                ml_m = ml_m, 
                dml_procedure = "dml2",
                n_folds = 3,
                n_rep = 3,
                trimming_threshold = 0.01
            )
            dml_obj$fit()
            
            # Calculation of DML results
            att_DML[i] <- dml_obj$coef[1]
            se[i] <- dml_obj$se[1]
            CI_lower[i] <- dml_obj$coef[1] - 1.96 * dml_obj$se[1]
            CI_upper[i] <- dml_obj$coef[1] + 1.96 * dml_obj$se[1]
            coverage[i] <- (true_ATT[i] >= CI_lower[i] && true_ATT[i] <= CI_upper[i])
        }
        
        # Store results for this specific config
        key <- paste0("N_", smpl, "_Design_", design)
        #saved_X[[key]] <- data$X
        all_results_20[[key]] <- list(
            DR = c(
                bias = mean(att_DR - true_ATT),
                coverage = mean(CI_DR[, 1] <= true_ATT & true_ATT <= CI_DR[, 2]),
                length = mean(CI_DR[, 2] - CI_DR[, 1])
            ),
            OR = c(
                bias = mean(att_OR - true_ATT),
                coverage = mean(CI_OR[, 1] <= true_ATT & true_ATT <= CI_OR[, 2]),
                length = mean(CI_OR[, 2] - CI_OR[, 1])
            ),
            IPW = c(
                bias = mean(att_IPW - true_ATT),
                coverage = mean(CI_IPW[, 1] <= true_ATT & true_ATT <= CI_IPW[, 2]),
                length = mean(CI_IPW[, 2] - CI_IPW[, 1])
            ),
            TWFE = c(
                bias = mean(att_TWFE - true_ATT),
                coverage = mean(CI_TWFE[, 1] <= true_ATT & true_ATT <= CI_TWFE[, 2]),
                length = mean(CI_TWFE[, 2] - CI_TWFE[, 1])
            ),
            DML = c(
                bias = mean(att_DML - true_ATT),
                coverage = mean(coverage),
                length = mean(CI_upper - CI_lower)
            )
        )
    }
}

end_time20 <- Sys.time()    # Record end time
execution_time20 <- end_time20 - start_time20
# Save the results in a csv file
df_results_20 <- as.data.frame(all_results_20)
df_results_20_transposed <- as.data.frame(t(df_results_20))
write.csv(df_results_20_transposed, file = '/Users/alikananitarigh/simulation_results_20nn.csv', row.names = TRUE)
write.csv(round(df_results_5_transposed,3), file = '/Users/alikananitarigh/simulation_results_20nn_round.csv', row.names = TRUE)
# ---------------------------------------------------------------------------
# End of the code
# ---------------------------------------------------------------------------