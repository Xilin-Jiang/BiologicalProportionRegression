cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
# set the commonly used color
blue <- cbPalette[6]
red <- cbPalette[7]
green <- cbPalette[4]
orange <- cbPalette[2]
grey <- cbPalette[1]
yellow <- cbPalette[5]
purple <- cbPalette[8]
skyblue <- cbPalette[3]
utils::globalVariables(c("prob", "section", "y"))

##################################################
# Biological-Proportion Regression functions
##################################################

######################################################
# 20240711: five different methods of latent factor analysis
######################################################
padding_X_Y_cov <- function(X, Y, Z, X_padding_num, Y_padding_num){
  X <- X %>%
    scale()
  Y <- Y %>%
    scale()
  Z <- apply(Z, 2, scale)
  ZtZ_cov <- cov(Z, use = "pairwise.complete.obs")

  # first extend the covariance matrix for paddings of Y
  Y_Z_cov <- cov(x = Z, y = Y, use = "pairwise.complete.obs")
  a <- matrix(var(Y, na.rm = T), nrow = Y_padding_num, ncol = Y_padding_num)
  b <- matrix(rep(Y_Z_cov, Y_padding_num) , nrow = length(Y_Z_cov), ncol = Y_padding_num)
  c <- ZtZ_cov
  # combine covariance abc
  ab <- rbind(a, b)
  bc <- rbind(t(b), c)
  Y_Z_padded <- cbind(ab, bc)

  # second extend the paddings to X
  X_Z_cov <- cov(x = Z, y = X, use = "pairwise.complete.obs")
  a <- matrix(var(X, na.rm = T), nrow = X_padding_num, ncol = X_padding_num)
  # b is exteded with covariance between Y and X
  b <- matrix(rep(c(rep(cov(x = X, y = Y, use = "pairwise.complete.obs"), Y_padding_num),X_Z_cov), X_padding_num),
              nrow = length(X_Z_cov) + Y_padding_num, ncol = X_padding_num)
  c <- Y_Z_padded
  ab <- rbind(a, b)
  bc <- rbind(t(b), c)
  XY_Z_padded <- cbind(ab, bc)
  return(XY_Z_padded)
}
# upweight X and Y
padding_XY_Z_cov <- function(X, Y, Z, num_padding_XY){
  XY_Z_padded <- padding_X_Y_cov(X, Y, Z, X_padding_num=num_padding_XY, Y_padding_num=num_padding_XY)
  return(XY_Z_padded)
}

# bootstrapping samples from the covariance matrix
##########
####### use jackknife instead of bootstrapping (bootstrapping is likely to cause numeric issue as it would have repeated rows)
jk_cov <- function(XY_Z_padded, X_padding_num, Y_padding_num, num_jk_block = 10){
  # divide the ids into num_jk_block numbers: num_jk_block > 1
  size_block <- floor((dim(XY_Z_padded)[1] - X_padding_num - Y_padding_num)/num_jk_block)
  jk_XY_Z_COV <- list()
  for(jk_idx in 1:num_jk_block){
    if(jk_idx == num_jk_block){
      ids_exclude <-  (1+X_padding_num + Y_padding_num + (jk_idx - 1)*size_block ):(dim(XY_Z_padded)[1])
    }else{
      ids_exclude <-  (1+X_padding_num + Y_padding_num + (jk_idx - 1)*size_block ):(X_padding_num + Y_padding_num + jk_idx *size_block )
    }
    ids_include <- setdiff(1:dim(XY_Z_padded)[1], ids_exclude)
    jk_XY_Z_COV[[jk_idx]] <- XY_Z_padded[ids_include,ids_include]
  }
  return(jk_XY_Z_COV)
}
# create jackknife samples for OLS analysis
jk_ols_x2z <- function(OLS_x2z, num_jk_block = 10){
  # divide the ids into num_jk_block numbers
  size_block <- floor((length(OLS_x2z))/num_jk_block)
  jk_ols_effects <- list()
  for(jk_idx in 1:num_jk_block){
    if(jk_idx == num_jk_block){
      ids_exclude <-  (1+ (jk_idx - 1)*size_block ):(length(OLS_x2z))
    }else{
      ids_exclude <-  (1 + (jk_idx - 1)*size_block ):(jk_idx *size_block )
    }
    ids_include <- setdiff(1:length(OLS_x2z), ids_exclude)
    jk_ols_effects[[jk_idx]] <- OLS_x2z[ids_include]
  }
  return(jk_ols_effects)
}

bootstrapping_cov <- function(XY_Z_padded,X_padding_num, Y_padding_num){
  bs_idx <- sample((1+X_padding_num + Y_padding_num):dim(XY_Z_padded)[1], replace = T)
  return(XY_Z_padded[c(1:(X_padding_num + Y_padding_num), bs_idx ),
                     c(1:(X_padding_num + Y_padding_num), bs_idx )])
}


#' Title: Inverse variance weighting combines different estimates
#'
#' @param bioProp_results output of  `EWR` function.
#'
#' @return Summarised dataframe with IVW estimate.
#' @export
#'
#' @examples
#'   full_run_results_big_wrapper <- EWR(X_EXAMPLE,Z_EXAMPLE, use_as_helper = 2)
#'   IVW_summary(full_run_results_big_wrapper)
IVW_summary <- function(bioProp_results){
  iwv_estimate_per_protein <- bioProp_results %>%
    filter(sampling_se < 0.2) %>%
    mutate(ivw_bioprop = 1/sampling_se^2) %>%
    summarise(iwv_estimate = sum(Bioprop_xi1_xi2 * ivw_bioprop)/sum(ivw_bioprop), iwv_se = sqrt(1/sum(ivw_bioprop)), number_correlated_protein = n(), se_protein_biop = sd(Bioprop_xi1_xi2))
  return(iwv_estimate_per_protein)
}

############################################################
# 2024-11-25: use the partial cov to directly estimate EWR ratios
############################################################
# note the inpute matrix has to be ordered as (X,Y,Z) and we are compute with respect to Y
partial_cor_Y <- function(cor_XY_Z, Non_noise_Y_ratio = 1){ # Non_noise_Y_ratio is 1 - Proportion_measurement_noise_Y; using equation 5 from Strauss 1981
  total_col <- dim(cor_XY_Z)[1]
  cor_Zi_Zj <- cor_XY_Z[3:total_col, 3:total_col]
  rho_Zi_Y <- cor_XY_Z[3:total_col, 2]


  one_minus_rhoZiY2 <- 1-rho_Zi_Y^2/Non_noise_Y_ratio
  partial_cor_Z_Y <- (cor_Zi_Zj - rho_Zi_Y %*% t(rho_Zi_Y) / Non_noise_Y_ratio)/sqrt(one_minus_rhoZiY2 %*% t(one_minus_rhoZiY2))

  partial_cor_X_Z_Y <- (cor_XY_Z[3:total_col, 1] - cor_XY_Z[3:total_col, 2] * cor_XY_Z[2,1]/Non_noise_Y_ratio)/
    sqrt((1-cor_XY_Z[2,1]^2/Non_noise_Y_ratio) * (1-cor_XY_Z[3:total_col, 2]^2/Non_noise_Y_ratio))
  # testing
  # X_rsid <- resid(lm(X~Y,na.action="na.exclude"))
  # Z_rsid <- resid(lm(Z[,1]~Y,na.action="na.exclude"))
  # print(c(cor(X_rsid, Z_rsid), partial_cor_X_Z_Y[1]))
  partial_matrix <- cor_XY_Z
  partial_matrix[3:total_col, 3:total_col] <- partial_cor_Z_Y
  partial_matrix[3:total_col, 1] <- partial_cor_X_Z_Y
  partial_matrix[1, 3:total_col] <- partial_cor_X_Z_Y
  return(partial_matrix)
}

direct_EWR_ratio_estimate <- function(cor_XY_Z, Bioprop_Y = 1){
  total_col <- dim(cor_XY_Z)[1]
  # compute the total amount of variance in X explained by Z
  X_prop <- cor_XY_Z[1,3:total_col] %*%
    solve( cor_XY_Z[3:total_col, 3:total_col]) %*%
    cor_XY_Z[3:total_col, 1]

  # compute the proportion of variance specific for X; note here we use proportion of variance in Y to account for measurement error in Y
  # prop_Y <- cor_XY_Z[2,3:total_col] %*%
  #   solve( cor_XY_Z[3:total_col, 3:total_col]) %*%
  #   cor_XY_Z[3:total_col, 2]
  # print(paste0("initial estimation of Y has ", 1-prop_Y, " proportion of variance being noise; if this number is large, the estimate might not be accurate." ))
  #
   # start_buffer_value <- 0.05
  # partial_cor_XZ_Y <- partial_cor_Y(cor_XY_Z, Non_noise_Y_ratio = min(prop_Y+start_buffer_value,1))
  # while(any(is.na(partial_cor_XZ_Y))){
  #   start_buffer_value <- start_buffer_value + 0.01
  #   partial_cor_XZ_Y <- partial_cor_Y(cor_XY_Z, Non_noise_Y_ratio = min(prop_Y+start_buffer_value,1))
  # }
  partial_cor_XZ_Y <- partial_cor_Y(cor_XY_Z, Non_noise_Y_ratio = 1)

  X_specific_prop <- partial_cor_XZ_Y[1,3:total_col] %*%
    solve( partial_cor_XZ_Y[3:total_col, 3:total_col]) %*%
    partial_cor_XZ_Y[3:total_col, 1]

  X_specific_prop <- X_specific_prop * (1 - cor_XY_Z[2,1]^2)

  # adjust for the noise in Y
  X_share_prop <- (X_prop - X_specific_prop)/Bioprop_Y
  # recompute the X_specific_prop
  X_specific_prop <- X_prop - X_share_prop

  # cov
  # XZ_Y_cov <- partial_cor_XZ_Y[3:total_col, 1] *
  #   sqrt(1-cor_XY_Z[3:total_col,2]^2) *
  #   sqrt(1-cor_XY_Z[2,1]^2)
  # Y_conditioned_cov <- cor_XY_Z[3:total_col, 1]  - XZ_Y_cov
  # plot(Y_conditioned_cov, cor_XY_Z[3:total_col, 2])
  adjusting_factor <- 1 + X_specific_prop/X_share_prop
  results_ratios <- list()
  results_ratios$E_xi_XY <- X_share_prop
  results_ratios$E_xi_X <- X_specific_prop
  results_ratios$adjusting_factor <- adjusting_factor
  return(results_ratios)
}

# function to detect the dilusion effect and provides an alternative estimate.
BPR_dilusion_test <- function(cor_XY_Z, Bioprop_xi1, Bioprop_Y){
  total_col <- dim(cor_XY_Z)[1]
  # compute the total amount of variance in X explained by Z
  X_prop <- cor_XY_Z[1,3:total_col] %*%
    solve( cor_XY_Z[3:total_col, 3:total_col]) %*%
    cor_XY_Z[3:total_col, 1]
  partial_cor_XZ_Y <- partial_cor_Y(cor_XY_Z, Non_noise_Y_ratio = 1)
  partial_cov_XZ <- (sqrt(pmax(0, 1-cor_XY_Z[, 2]^2)) %*% t(sqrt(pmax(0, 1-cor_XY_Z[, 2]^2)))) * partial_cor_XZ_Y
  X_specific_cov <- partial_cov_XZ[1,3:total_col] %*%
    solve( partial_cov_XZ[3:total_col, 3:total_col]) %*%
    partial_cov_XZ[3:total_col, 1]
  # compute the proportion of variance
  B_term_X_share_cov <- (X_prop - X_specific_cov)
  Py_Exi1_divBterm  <- (Bioprop_xi1 * Bioprop_Y)/B_term_X_share_cov
  T_zxi1 <- 1 - Py_Exi1_divBterm
  PyExi1 <- Bioprop_xi1 * Bioprop_Y
  Py_2_Exi1 <- (2-Bioprop_Y) * Bioprop_xi1

  # only the second estimate is used as there should be only one positive solution
  if(T_zxi1 * (PyExi1^2 - Py_Exi1_divBterm * Py_2_Exi1^2) < 0){
    E_Z_xi1_second_estimate <- 0
    X_share_cov <- (X_prop-X_specific_cov)/Bioprop_Y
    X_specific_cov <- X_prop - X_share_cov
  }else{
    # E_Z_xi1_first_estimate <- -Py_2_Exi1/2 * T_zxi1 -
    #   1/2 * sqrt(T_zxi1 * (PyExi1^2 - Py_Exi1_divBterm * Py_2_Exi1^2))
    E_Z_xi1_second_estimate <- -Py_2_Exi1/2 * T_zxi1 +
      1/2 * sqrt(T_zxi1 * (PyExi1^2 - Py_Exi1_divBterm * Py_2_Exi1^2))
    X_share_cov <- Bioprop_xi1 * Bioprop_xi1/(Bioprop_xi1 + E_Z_xi1_second_estimate)
    X_specific_cov <- X_prop - X_share_cov
  }

  adjusting_factor <- 1 + X_specific_cov/X_share_cov
  results_ratios <- list()
  results_ratios$E_xi_XY <- X_share_cov
  results_ratios$E_xi_X <- X_specific_cov
  results_ratios$dilusion_testStat <- T_zxi1
  results_ratios$dilusion_variance <- E_Z_xi1_second_estimate
  results_ratios$adjusting_factor <- adjusting_factor
  return(results_ratios)
}


# clumping covariance matrix to remove highly correlated traits
clump_cov_Z <- function(cov_Z, threshold = 0.3){
  rm.idx <- c()
  keep.idx <- c()
  diag(cov_Z) <- 0
  seed.rm <- which(abs(cov_Z) == max(abs(cov_Z)), arr.ind = TRUE)
  seed.keep <- sample(seed.rm[1,], size = 1)
  keep.idx <- unique(c(keep.idx, seed.keep))
  rm.nodes <- which(cov_Z[seed.keep,] > threshold)
  rm.idx <- unique(c(rm.idx , rm.nodes))

  cov_Z_remain <- cov_Z
  cov_Z_remain[rm.idx, ] <- 0
  cov_Z_remain[,rm.idx] <- 0
  while(max(abs(cov_Z_remain)) > threshold){
    seed.rm <- which(abs(cov_Z_remain) == max(abs(cov_Z_remain)), arr.ind = TRUE)
    seed.keep <- sample(seed.rm[1,], size = 1)
    keep.idx <- unique(c(keep.idx, seed.keep))
    rm.nodes <- which(abs(cov_Z_remain[seed.keep,]) > threshold)
    # print(length(intersect(keep.idx, rm.nodes)))
    rm.idx <- unique(c(rm.idx , rm.nodes))

    cov_Z_remain <- cov_Z
    cov_Z_remain[rm.idx,] <- 0
    cov_Z_remain[,rm.idx] <- 0
  }
  return(rm.idx)

}

#  summary level data to perform fast jackknifing for se estimates
jackknife_acrossZ_partialCor <- function(cor_XY_Z, OLS_y2z, OLS_x2z, PxPy,
                                       num_jk_block = 20, covClumpThreshold = 0.5
){

  jk_cov_list <- jk_cov(cor_XY_Z,X_padding_num = 1, Y_padding_num = 1, num_jk_block = num_jk_block)
  jk_OLS_y2z <- jk_ols_x2z(OLS_y2z, num_jk_block = num_jk_block)
  jk_OLS_x2z <- jk_ols_x2z(OLS_x2z, num_jk_block = num_jk_block)

  jk_results <- data.frame(jk_idx = as.numeric(),
                           beta_EWR = as.numeric(),
                           R2_EWR_jk = as.numeric(),
                           Bioprop_xi1 = as.numeric(),
                           Partial_beta = as.numeric(),
                           Bioprop_Y = as.numeric(),
                           adjusting_factor = as.numeric(),
                           E_xi_XY = as.numeric(),
                           E_xi_X = as.numeric(),
                           Bioprop_xi1_xi2= as.numeric(),
                           dilusion_BPR_Bioprop_xi1_xi2 = as.numeric(),
                           dilusion_test_statistics = as.numeric()
                           )
  for(jk_idx in 1:num_jk_block){
    # beta_EWR
    OLS_x2z_jk <- jk_OLS_x2z[[jk_idx]]
    OLS_y2z_jk <- jk_OLS_y2z[[jk_idx]]
    cor_XY_Z_jk <- jk_cov_list[[jk_idx]]
    EWR_model_jk <- summary(lm(OLS_y2z_jk ~ OLS_x2z_jk))
    beta_EWR_jk <- EWR_model_jk$coefficients[2,1]
    R2_EWR_jk <- EWR_model_jk$adj.r.squared
    Bioprop_xi1 <- sqrt(PxPy/beta_EWR_jk^2) * R2_EWR_jk
    Partial_beta <- beta_EWR_jk/R2_EWR_jk
    Bioprop_Y <- Bioprop_xi1 * Partial_beta^2

    total_col <- dim(cor_XY_Z_jk)[1]
    cov_Z <- cor_XY_Z_jk[3:total_col, 3:total_col]
    clumping_idx <-  clump_cov_Z(cov_Z, threshold = covClumpThreshold)
    clumping_idx_pad <- clumping_idx + 2
    clumped_cor_XY_Z <- cor_XY_Z_jk[- clumping_idx_pad, -clumping_idx_pad ]
    clump_prop <- bind_rows(lapply(direct_EWR_ratio_estimate(clumped_cor_XY_Z, Bioprop_Y = Bioprop_Y), c))

    Bioprop_xi1_xi2 <- clump_prop$adjusting_factor * Bioprop_xi1

    BPR_dilusion_estimates <- bind_rows(lapply(BPR_dilusion_test(clumped_cor_XY_Z, Bioprop_xi1 = Bioprop_xi1, Bioprop_Y = Bioprop_Y), c))
    dilusion_BPR_Bioprop_xi1_xi2 <- BPR_dilusion_estimates$adjusting_factor * Bioprop_xi1

    jk_results <- jk_results %>%
      add_row(jk_idx = jk_idx,
              beta_EWR = beta_EWR_jk,
              R2_EWR_jk = R2_EWR_jk,
              Bioprop_xi1 = Bioprop_xi1,
              Partial_beta = Partial_beta,
              Bioprop_Y = Bioprop_Y,
              adjusting_factor = clump_prop$adjusting_factor,
              E_xi_XY = clump_prop$E_xi_XY,
              E_xi_X = clump_prop$E_xi_X,
              Bioprop_xi1_xi2= Bioprop_xi1_xi2,
              dilusion_BPR_Bioprop_xi1_xi2 = dilusion_BPR_Bioprop_xi1_xi2,
              dilusion_test_statistics = BPR_dilusion_estimates$dilusion_testStat)
  }
  return(jk_results)
}


#' Title: EWR using partial covariance matrix for stage 2, faster option
#'
#' @param X a vector which denote the levels of the target protein.
#' @param Y a vector (same length of X) which denote the levels of the helper protein, should be significantly correlated with X.
#' @param Z_EWR for performing EWR regression: a matrix (each row correspond to one element in X and Y) of auxiliary proteins that tags biological variance in X. We recommend at least 100 columns.
#' @param Z_proportion_ratio for estimate ratio of variance: a matrix (each row correspond to one element in X and Y) of auxiliary proteins that tags biological variance in X. We recommend at least 100 columns.
#' @param bootstrapping_number number of bootstraps for computing confidence interval.
#' @param covClumpThreshold when computing the partial correlation matrix, clumping highly correlated pairs of Zs; default threshold cor>0.5.
#'
#' @return a list containing objects:
#' \itemize{
#'    \item Bioprop_xi1_xi2: estimate of biological proportion (i.e. 1-noise) variance in protein level X.
#'    \item sampling_se: standard error of the estimate Bioprop_xi1_xi2; estimated through sampling procedure.
#'    \item sampling_mean: mean across the sample estimate (could be used as a bias detection tool).
#'    \item jk_Zlabel_se_bioprop: se estimated by jacknifing Zj j index (as opposed to jackknifing individuals).
#'    \item PxPy
#'    \item PxPy_se: standard error of PxPy estimate.
#'    \item beta_EWR: Slope of covariance y2z on x2z.
#'    \item beta_EWR_se: standard error of beta_EWR estimate.
#'    \item Bioprop_xi1: xi_1 variance
#'    \item Partial_beta: partial causal effect
#'    \item Bioprop_Y: beta^2 * E_xi12
#'    \item sample_Bioprop_xi1: sample mean of Bioprop_xi1
#'    \item sample_Partial_beta: sample mean of sample_Partial_beta
#'    \item sample_Bioprop_Y: sample mean of Bioprop_Y
#'    \item adjusting_factor: adjusting factor computed from the correlation matrix.
#'    \item sample_adjusting_factor: adjusting factor computed from the correlation matrix, mean across bootstrap samples.
#'    \item E_xi_XY: estimated variance in X that is shared with Y from the correlation matrix.
#'    \item E_xi_X: estimated variance in X that is not shared with Y from the correlation matrix.
#'    \item sample_E_xi_XY: bootstrap mean of E_xi_XY.
#'    \item sample_E_xi_X: bootstrap mean of E_xi_X.
#'    \item X_variance_corMatrix: proportion of variance estimated directly from the Z.
#'    \item se_X_variance_corMatrix: bootstrap se of X_variance_corMatrix.
#'    \item sigma_X: directly estimated residual variance in X that is not captured by Z.
#'    \item sigma_Y: directly estimated residual variance in Y that is not captured by Z.
#'}
#'
#' @export
#'
#' @examples
#' quick_results <- bioProp_EWR_partialCor_adjustment(X_EXAMPLE, Y_EXAMPLE, Z_EWR = Z_EXAMPLE[,1:200], bootstrapping_number=20)
#'
bioProp_EWR_partialCor_adjustment <- function(X, Y, Z_EWR, Z_proportion_ratio = NULL, bootstrapping_number = 100, covClumpThreshold = 0.5){ # X is N x 1, Y is N x 1, Z is N x J
  print(paste0("Number of auxiliary proteins: ", dim(Z_EWR)[2], "; Recommended at auxiliary protein > 100"))

  X_notNA_idx <- which(!is.na(X))
  print(paste0("Keep ", length(X_notNA_idx), " rows where X is not missing"))
  X <- scale(X[X_notNA_idx])
  Y <- scale(Y[X_notNA_idx])
  Z_EWR <- apply(Z_EWR[X_notNA_idx, ], 2, scale)
  if(is.null(Z_proportion_ratio)){
    Z_proportion_ratio <- Z_EWR
  }else{
    Z_proportion_ratio <- apply(Z_proportion_ratio[X_notNA_idx, ], 2, scale)
  }

  # two step procedure

  cor_XY_Z <- padding_XY_Z_cov(X, Y, Z_proportion_ratio, num_padding_XY = 1)
  total_col <- dim(cor_XY_Z)[1]
  cov_Z <- cor_XY_Z[3:total_col, 3:total_col]

  print(paste("Perform covariance regression."))
  PxPy_model <- summary(lm(scale(X) ~ scale(Y)))
  PxPy <- PxPy_model$coefficients[2,1]^2
  # beta_EWR
  OLS_x2z <- c(cor(X, Z_EWR, use = "pairwise.complete.obs"))
  OLS_y2z <- c(cor(Y, Z_EWR, use = "pairwise.complete.obs"))

  # start a 10 random clumping and test if it is robust
  num_jk_block <- 20
  print(paste("Jackknife across Z index, using ", num_jk_block, " samples."))
  jk_initial_samples <- jackknife_acrossZ_partialCor(cor_XY_Z, OLS_y2z, OLS_x2z, PxPy,
                               num_jk_block = num_jk_block, covClumpThreshold = covClumpThreshold)
  jk_se_bioprop <- sd(jk_initial_samples$Bioprop_xi1_xi2) * sqrt(num_jk_block)
  print(paste("Jackknife se: ", jk_se_bioprop, ". Z-score = ", mean(jk_initial_samples$Bioprop_xi1_xi2)/jk_se_bioprop))

  # is the clumping procedure robust?
  jk_se_E_xi_XY <- sd(jk_initial_samples$E_xi_XY) * sqrt(num_jk_block)
  clumping_zscore <- mean(jk_initial_samples$E_xi_XY)/jk_se_E_xi_XY
  if(clumping_zscore < 3){
    print(paste("Clumping procedure not robust, might have clumped too many Zs; use higher clumping threshold r = 0.5."))
    covClumpThreshold <- 0.5
  }

  # clump the covariance matrix so only non correlated Z are preserved.
  clumping_idx <-  clump_cov_Z(cov_Z, threshold = covClumpThreshold)

  Z_clumped <- Z_proportion_ratio[, -clumping_idx]
  cor_XY_Z <- padding_XY_Z_cov(X, Y, Z_clumped, num_padding_XY = 1)
  total_col <- dim(cor_XY_Z)[1]

  EWR_model <- summary(lm(OLS_y2z ~ OLS_x2z))
  beta_EWR <- EWR_model$coefficients[2,1]
  beta_EWR_se <- EWR_model$coefficients[2,2]
  R2_EWR <- EWR_model$adj.r.squared
  # print(paste("Slope of covariance y2z on x2z is: ", beta_EWR, "(", beta_EWR_se, ")", " with z-score ", beta_EWR/beta_EWR_se, ". |Z-score| > 5 required for the estimate to be approximately unbiased."))
  Bioprop_xi1 <- sqrt(PxPy/beta_EWR^2) * R2_EWR
  Partial_beta <- beta_EWR/R2_EWR
  Bioprop_Y <- Bioprop_xi1 * Partial_beta^2

  # compute the covariance matrix
  EWR_correction <- direct_EWR_ratio_estimate(cor_XY_Z, Bioprop_Y = Bioprop_Y)

  # test X Z dilution; using equations from overleaf
  dilution_correction <- BPR_dilusion_test(cor_XY_Z = cor_XY_Z, Bioprop_xi1 = Bioprop_xi1, Bioprop_Y = Bioprop_Y)

  bootstrapping_number <- bootstrapping_number
  # cor_XY_Z_bt_all <- array(dim = c( total_col, total_col, bootstrapping_number))

  bt_PxPy <- c()
  bt_beta_EWR <- c()
  bt_R2_EWR <- c()
  bt_Bioprop_xi1 <- c()
  bt_Partial_beta <- c()
  bt_Bioprop_Y <- c()
  bt_E_xi_XY <- c()
  bt_E_xi_X <- c()
  # bt_var_eta_1j <- c()
  # bt_var_eta_2j_E2_xi_X <- c()
  bt_EWR_adjusting_factor <- c()
  bt_BPR_dilusion_adjusting_factor <- c()
  bt_dilusion_test_statistics <- c()
  for(bt in 1:bootstrapping_number){
    print(paste0("bootstrapping sampling: ", bt))
    try({
      bt_labels <- sample(1:length(X), replace = T, size = length(X))

      PxPy_bt <- cor(X[bt_labels], Y[bt_labels], use = "pairwise.complete.obs")^2
      bt_PxPy <- c(bt_PxPy, PxPy_bt) # save samples for se computation
      OLS_x2z <- c(cor(X[bt_labels], Z_EWR[bt_labels,], use = "pairwise.complete.obs"))
      OLS_y2z <- c(cor(Y[bt_labels], Z_EWR[bt_labels,], use = "pairwise.complete.obs"))

      EWR_model_bt <- summary(lm(OLS_y2z ~ OLS_x2z))
      beta_EWR_bt <- EWR_model_bt$coefficients[2,1]
      bt_beta_EWR <- c(bt_beta_EWR, beta_EWR_bt) # save samples for se computation
      bt_R2_EWR <- c(bt_R2_EWR, EWR_model_bt$adj.r.squared)
      # beta_EWR_bt <- summary(lm(cor_XY_Z_bt[3:total_col, 2] ~ cor_XY_Z_bt[3:total_col, 1]))$coefficients[2,1]
      bt_Bioprop_xi1 <- c(bt_Bioprop_xi1, sqrt(PxPy_bt/beta_EWR_bt^2) * EWR_model_bt$adj.r.squared)
      bt_Partial_beta <- c(bt_Partial_beta, beta_EWR_bt/EWR_model_bt$adj.r.squared)
      Bioprop_Y_bt <- sqrt(PxPy_bt * beta_EWR_bt^2)/EWR_model_bt$adj.r.squared
      bt_Bioprop_Y <- c(bt_Bioprop_Y, Bioprop_Y_bt)

      # perform random clumping
      clumping_idx <-  clump_cov_Z(cov_Z, threshold = covClumpThreshold)
      Z_clumped <- Z_proportion_ratio[, -clumping_idx]
      cor_XY_Z_bt <- padding_XY_Z_cov(X[bt_labels], Y[bt_labels], Z_clumped[bt_labels,], num_padding_XY = 1)
      total_col <- dim(cor_XY_Z_bt)[1]
      EWR_correction_bt <- direct_EWR_ratio_estimate(cor_XY_Z_bt, Bioprop_Y = Bioprop_Y_bt)
      # cor_XY_Z_bt_all[, , bt] <- cor_XY_Z_bt

      # test for partial tagging X -> Z
      BPR_dilusion_correction_bt <- BPR_dilusion_test(cor_XY_Z_bt, Bioprop_xi1 = sqrt(PxPy_bt/beta_EWR_bt^2) * EWR_model_bt$adj.r.squared,
                                                      Bioprop_Y = Bioprop_Y_bt)
      bt_BPR_dilusion_adjusting_factor <- c(bt_BPR_dilusion_adjusting_factor, BPR_dilusion_correction_bt$adjusting_factor)
      bt_dilusion_test_statistics <- c(bt_dilusion_test_statistics, BPR_dilusion_correction_bt$dilusion_testStat)

      # compute the total amount of variance in X explained by Z
      bt_E_xi_XY <- c(bt_E_xi_XY, EWR_correction_bt$E_xi_XY)
      bt_E_xi_X <- c(bt_E_xi_X, EWR_correction_bt$E_xi_X)
      # bt_var_eta_1j <- c(bt_var_eta_1j, EWR_correction_bt$var_eta_1j)
      # bt_var_eta_2j_E2_xi_X <- c(bt_var_eta_2j_E2_xi_X, EWR_correction_bt$var_eta_2j_E2_xi_X)
      bt_EWR_adjusting_factor <- c(bt_EWR_adjusting_factor, EWR_correction_bt$adjusting_factor)
    })
  }
  beta_EWR_se <- sd(bt_beta_EWR) # have to use this se as the covariance are not independent.
  R2_EWR_se <- sd(bt_R2_EWR)
  print(paste("R2 of covariance y2z on x2z is: ", R2_EWR, "(", R2_EWR_se, ")", " with z-score ", R2_EWR/R2_EWR_se, ". |Z-score| > 5 required for the estimate to be approximately unbiased."))
  print(paste("Slope of covariance y2z on x2z is: ", beta_EWR, "(", beta_EWR_se, ")", " with z-score ", beta_EWR/beta_EWR_se, ". |Z-score| > 5 required for the estimate to be approximately unbiased."))
  # comparison below suggest the inverse of variance is the main source of bias; probably safe to proceed with simple subtraction
  # print(c(mean(bt_approx_Bioprop_xi1_xi2 * bt_EWR_adjusting_factor), approx_Bioprop_xi1_xi2 * EWR_correction$adjusting_factor))
  # print(c(mean(bt_E_xi_XY), EWR_correction$E_xi_XY))
  # print(c(mean(bt_E_xi_X), EWR_correction$E_xi_X))
  # print(c(mean(bt_var_eta_1j), EWR_correction$var_eta_1j))
  # print(c(mean(bt_var_eta_2j_E2_xi_X), EWR_correction$var_eta_2j_E2_xi_X))
  # se needs to consider the bias correction; note this is an approximation as we ignore the impact of "approx_Bioprop_xi1_xi2" on se in the bias correction.
  sampling_se <- sd(bt_Bioprop_xi1 * bt_EWR_adjusting_factor) * (1+1/sqrt(bootstrapping_number))
  sampling_mean <- mean(bt_Bioprop_xi1 * bt_EWR_adjusting_factor)
  # use bootstrap to correct for bias in the EWR adjusting procedure
  Bioprop_xi1_xi2 <- Bioprop_xi1 * EWR_correction$adjusting_factor -
    Bioprop_xi1 * EWR_correction$adjusting_factor/sampling_mean *(sampling_mean - Bioprop_xi1 * EWR_correction$adjusting_factor )

  # compute the proportion of variance in Y that are captured by other Zs
  total_col <- dim(cor_XY_Z)[1]
  prop_Y_naive <- cor_XY_Z[2,3:total_col] %*%
    solve( cor_XY_Z[3:total_col, 3:total_col]) %*%
    cor_XY_Z[3:total_col, 2]

  # compute the adjusted proportion of variance in X explained by Z and its standard error
  X_variance_corMatrix <- (EWR_correction$E_xi_XY + EWR_correction$E_xi_X) * 2 -
    (mean(bt_E_xi_XY) + mean(bt_E_xi_X))
  se_X_variance_corMatrix <- sd(bt_E_xi_XY +bt_E_xi_X)

  # test if Z only tags part of the variance in X and provide an alternative estimate that account for this
  dilusion_T_se <- sd(bt_dilusion_test_statistics)
  dilusion_test_statistics <- dilution_correction$dilusion_testStat
  print(paste("Test if Z only tags part of X; Z-score: ", dilusion_test_statistics/dilusion_T_se, "."))

  dilusion_BPR_sampling_se <- sd(bt_Bioprop_xi1 * bt_BPR_dilusion_adjusting_factor) * (1+1/sqrt(bootstrapping_number))
  dilusion_BPR_sampling_mean <- mean(bt_Bioprop_xi1 * bt_BPR_dilusion_adjusting_factor)
  # use bootstrap to correct for bias in the EWR adjusting procedure
  dilusion_BPR_Bioprop_xi1_xi2 <- Bioprop_xi1 * dilution_correction$adjusting_factor -
    Bioprop_xi1 * dilution_correction$adjusting_factor/dilusion_BPR_sampling_mean *
    (dilusion_BPR_sampling_mean - Bioprop_xi1 * dilution_correction$adjusting_factor )


  results_biprop <- list()
  results_biprop$Bioprop_xi1_xi2 <- c(Bioprop_xi1_xi2)
  results_biprop$sampling_se <- sampling_se # jackknife (or boostraping) se
  results_biprop$sampling_mean <- sampling_mean # jackknife (or boostraping) mean, use to identify bias
  results_biprop$jk_Zlabel_se_bioprop <- jk_se_bioprop # reporting the se estimated from jacknifing Zj index
  results_biprop$PxPy <- PxPy
  results_biprop$PxPy_se <- sd(bt_PxPy)
  results_biprop$PxPy_sample_mean <- mean(bt_PxPy)
  results_biprop$beta_EWR <- beta_EWR
  results_biprop$beta_EWR_se <- sd(bt_beta_EWR)
  results_biprop$beta_EWR_sample_mean <- mean(bt_beta_EWR)
  results_biprop$R2_EWR <- R2_EWR
  results_biprop$R2_EWR_se <- R2_EWR_se
  results_biprop$R2_EWR_sample_mean <- mean(bt_R2_EWR)
  results_biprop$Bioprop_xi1 <- Bioprop_xi1
  results_biprop$Partial_beta <- Partial_beta
  results_biprop$Bioprop_Y <- Bioprop_Y
  results_biprop$sample_Bioprop_xi1 <- mean(bt_Bioprop_xi1)
  results_biprop$sample_Partial_beta <- mean(bt_Partial_beta)
  results_biprop$sample_Bioprop_Y <- mean(bt_Bioprop_Y)
  results_biprop$adjusting_factor <- EWR_correction$adjusting_factor
  results_biprop$sample_adjusting_factor <- mean(bt_EWR_adjusting_factor) # this is the factor used to detect bias in the adjusting_factor
  results_biprop$E_xi_XY <- EWR_correction$E_xi_XY
  results_biprop$E_xi_X <- EWR_correction$E_xi_X
  results_biprop$sample_E_xi_XY <- mean(bt_E_xi_XY)
  results_biprop$sample_E_xi_X <- mean(bt_E_xi_X)
  results_biprop$X_variance_corMatrix <- X_variance_corMatrix
  results_biprop$se_X_variance_corMatrix <- se_X_variance_corMatrix
  results_biprop$sigma_X <- 1 - X_variance_corMatrix
  results_biprop$sigma_Y <- 1 - prop_Y_naive

  # save the alternative estimate when Z partially tags X
  results_biprop$dilusion_T_se <- dilusion_T_se
  results_biprop$dilusion_test_statistics <- dilusion_test_statistics
  results_biprop$dilusion_variance <- dilution_correction$dilusion_variance # if this number is 0, then there is no solution for the dilusion factor
  results_biprop$dilusion_BPR_sampling_se <- dilusion_BPR_sampling_se
  results_biprop$dilusion_BPR_sampling_mean <- dilusion_BPR_sampling_mean
  # use bootstrap to correct for bias in the EWR adjusting procedure
  results_biprop$dilusion_BPR_Bioprop_xi1_xi2 <- dilusion_BPR_Bioprop_xi1_xi2

  return(results_biprop)
}


#' Title: BPR wrapper function. Input is many pairs of X and Y, and all based on the same set of Z.
#'
#' @param target_protein a vector which denote the levels of the target protein.
#' @param auxiliary_proteins matrix (each row correspond to one element in target proteins) of auxiliary proteins that tags biological variance in target proteins. We recommend at least 100 columns.
#' @param top_PCs number of PCs estimated from auxiliary proteins; default 20.
#' @param number_of_Y_chosen number of PCs used as Y for the EWR estimation; deafalt 5.
#' @param bootstrapping_number number of boostrap samples for estimate uncertainty.
#' @param threshold_PxPy a threshold for cor(PC, X)^2, default 0.
#' @param Z_pca_number number of Ys chosen to compute Y.
#'
#' @returns A dataframe contains same output as the `bioProp_EWR_partialCor_adjustment` function.
#' @export
#'
#' @examples
#' full_run_results_big_wrapper <- BPR_fast_wrapper(X_EXAMPLE,Z_EXAMPLE,number_of_Y_chosen = 2, bootstrapping_number = 5)
BPR_fast_wrapper <- function(target_protein, auxiliary_proteins, top_PCs = 50, number_of_Y_chosen = 5, bootstrapping_number = 100,
                    threshold_PxPy = 0, Z_pca_number = 500){
  target_protein <- c(target_protein)
  if(length(target_protein) != dim(auxiliary_proteins)[1]){
    stop("dimension of target protein and auxiliary proteins is not matching.")
  }

  correlation_target_aux <- cor(target_protein, auxiliary_proteins, use = "pairwise.complete.obs")

  if(max(abs(correlation_target_aux)) > 0.99){
    stop("Really high correlation between target protein and on auxiliary protein. Did you accidentally put the target protein in the helper proteins?")
  }

  if(max(abs(correlation_target_aux)) < 0.05){
    stop("target protein doesn't seem to be correlated with any of auxiliary proteins. Have you checked the ordering is consistent between target and auxiliary proteins?")
  }

  # find the set of proteins for computing PCs as Y
  if(Z_pca_number > length(correlation_target_aux)/2){
    Z_pca_number <- floor(length(correlation_target_aux)/2)
  }
  # Z_pca_idx <- colnames(correlation_target_aux)[which(rank(dplyr::desc(abs(correlation_target_aux))) <=Z_pca_number)]
  # Z_nonY_idx <- colnames(correlation_target_aux)[which(rank(dplyr::desc(abs(correlation_target_aux))) >Z_pca_number)]
  Z_pca_idx <- sample(1:length(correlation_target_aux), size = Z_pca_number, replace = F)
  Z_nonY_idx <- setdiff(1:length(correlation_target_aux), Z_pca_idx)
  Z_pca <- auxiliary_proteins[,Z_pca_idx]
  Z_proportion_ratio <- auxiliary_proteins[,Z_nonY_idx]

  X <- target_protein
  X_notNA_idx <- which(!is.na(X))
  X <- scale(X[X_notNA_idx])
  Z_proportion_ratio <- apply(Z_proportion_ratio[X_notNA_idx, ], 2, scale)
  Z_pca <- apply(Z_pca[X_notNA_idx, ], 2, scale)

  print(paste0("Compute top ", top_PCs, " PCs from the auxiliary proteins"))
  Z_pca[is.na(Z_pca)] <- 0
  pca_z <- prcomp(Z_pca)
  Z_toppcs <- Z_pca %*% pca_z$rotation[,1:top_PCs]
  Z_toppcs <- apply(Z_toppcs, 2, scale)
  lm_X_zpcs <- summary(lm(X~Z_toppcs))
  Y_pc_list <- which(rank(dplyr::desc(abs(lm_X_zpcs$coefficients[2:(dim(Z_toppcs)[2] + 1),1]))) <=number_of_Y_chosen)


  EWR_results <- list()
  for(pc_idx in 1:length(Y_pc_list)){
    pc_choice <- Y_pc_list[pc_idx]
    # use a threshold of 0.05 for minial variance explained by Y
    if(lm_X_zpcs$coefficients[pc_choice + 1,1]^2 > threshold_PxPy){
      print(paste0("EWR analysis using PC", pc_choice))
      EWR_results[[pc_idx]] <- bind_rows(bioProp_EWR_partialCor_adjustment(X,Y = Z_toppcs[,pc_choice], Z_EWR = Z_pca, Z_proportion_ratio = Z_proportion_ratio, bootstrapping_number = bootstrapping_number)) %>%
        mutate(PC_id = paste0("PC_", pc_choice))
    }
  }
  EWR_results <- bind_rows(EWR_results)
  return(EWR_results)
}







