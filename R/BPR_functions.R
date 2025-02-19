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
utils::globalVariables(c("prob", "section", "y","cor", "prcomp", "lm", "sampling_se",
                         "Bioprop_xi1_xi2", "ivw_bioprop"))

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
#' @param bioProp_results output of  `BPR` function.
#'
#' @return Summarised dataframe with IVW estimate.
#' @export
#'
#' @examples
#'   full_run_results_big_wrapper <- BPR_fast_wrapper(cbind(X_EXAMPLE, X_EXAMPLE), cbind(Y_EXAMPLE,Z_EXAMPLE[,1]), Z_EXAMPLE, bootstrapping_number = 5)
#'   IVW_summary(full_run_results_big_wrapper)
IVW_summary <- function(bioProp_results){
  iwv_estimate_per_protein <- bioProp_results %>%
    filter(sampling_se < 0.2) %>%
    mutate(ivw_bioprop = 1/sampling_se^2, ivw_dilution_BPR = 1/dilution_BPR_sampling_se^2) %>%
    summarise(iwv_estimate = sum(Bioprop_xi1_xi2 * ivw_bioprop)/sum(ivw_bioprop), iwv_se = sqrt(1/sum(ivw_bioprop)),
              iwv_dilution_BPR_estimate = sum(dilution_BPR_Bioprop_xi1_xi2 * ivw_dilution_BPR)/sum(ivw_dilution_BPR), iwv_dilution_BPR_se = sqrt(1/sum(ivw_dilution_BPR)),
              number_correlated_protein = n(),
              se_protein_biop = sd(Bioprop_xi1_xi2),
              se_BPR_Bioprop = sd(dilution_BPR_Bioprop_xi1_xi2))
  return(iwv_estimate_per_protein)
}

############################################################
# 2024-11-25: use the partial cov to directly estimate BPR ratios
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

direct_BPR_ratio_estimate <- function(cor_XY_Z, Bioprop_Y = 1){
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

# function to detect the dilution effect and provides an alternative estimate.
BPR_dilution_test <- function(cor_XY_Z, Bioprop_xi1, Bioprop_Y){
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
  results_ratios$dilution_testStat <- T_zxi1
  results_ratios$dilution_variance <- E_Z_xi1_second_estimate
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
                           beta_BPR = as.numeric(),
                           R2_BPR_jk = as.numeric(),
                           Bioprop_xi1 = as.numeric(),
                           Partial_beta = as.numeric(),
                           Bioprop_Y = as.numeric(),
                           adjusting_factor = as.numeric(),
                           E_xi_XY = as.numeric(),
                           E_xi_X = as.numeric(),
                           Bioprop_xi1_xi2= as.numeric(),
                           dilution_BPR_Bioprop_xi1_xi2 = as.numeric(),
                           dilution_test_statistics = as.numeric()
                           )
  for(jk_idx in 1:num_jk_block){
    # beta_BPR
    OLS_x2z_jk <- jk_OLS_x2z[[jk_idx]]
    OLS_y2z_jk <- jk_OLS_y2z[[jk_idx]]
    cor_XY_Z_jk <- jk_cov_list[[jk_idx]]
    BPR_model_jk <- summary(lm(OLS_y2z_jk ~ OLS_x2z_jk))
    beta_BPR_jk <- BPR_model_jk$coefficients[2,1]
    R2_BPR_jk <- BPR_model_jk$adj.r.squared
    Bioprop_xi1 <- sqrt(PxPy/beta_BPR_jk^2) * R2_BPR_jk
    Partial_beta <- beta_BPR_jk/R2_BPR_jk
    Bioprop_Y <- Bioprop_xi1 * Partial_beta^2

    total_col <- dim(cor_XY_Z_jk)[1]
    cov_Z <- cor_XY_Z_jk[3:total_col, 3:total_col]
    clumping_idx <-  clump_cov_Z(cov_Z, threshold = covClumpThreshold)
    clumping_idx_pad <- clumping_idx + 2
    clumped_cor_XY_Z <- cor_XY_Z_jk[- clumping_idx_pad, -clumping_idx_pad ]
    clump_prop <- bind_rows(lapply(direct_BPR_ratio_estimate(clumped_cor_XY_Z, Bioprop_Y = Bioprop_Y), c))

    Bioprop_xi1_xi2 <- clump_prop$adjusting_factor * Bioprop_xi1

    BPR_dilution_estimates <- bind_rows(lapply(BPR_dilution_test(clumped_cor_XY_Z, Bioprop_xi1 = Bioprop_xi1, Bioprop_Y = Bioprop_Y), c))
    dilution_BPR_Bioprop_xi1_xi2 <- BPR_dilution_estimates$adjusting_factor * Bioprop_xi1

    jk_results <- jk_results %>%
      add_row(jk_idx = jk_idx,
              beta_BPR = beta_BPR_jk,
              R2_BPR_jk = R2_BPR_jk,
              Bioprop_xi1 = Bioprop_xi1,
              Partial_beta = Partial_beta,
              Bioprop_Y = Bioprop_Y,
              adjusting_factor = clump_prop$adjusting_factor,
              E_xi_XY = clump_prop$E_xi_XY,
              E_xi_X = clump_prop$E_xi_X,
              Bioprop_xi1_xi2= Bioprop_xi1_xi2,
              dilution_BPR_Bioprop_xi1_xi2 = dilution_BPR_Bioprop_xi1_xi2,
              dilution_test_statistics = BPR_dilution_estimates$dilution_testStat)
  }
  return(jk_results)
}


#' Title: Biological-proportion-regression main function.
#'
#' @param X a vector which denote the levels of the target protein.
#' @param Y a vector (same length of X) which denote the levels of the helper protein, should be significantly correlated with X.
#' @param Z_BPR for performing BPR regression: a matrix (each row correspond to one element in X and Y) of auxiliary proteins that tags biological variance in X. We recommend at least 100 columns.
#' @param Z_proportion_ratio for estimate ratio of variance (default = Z_BPR): a matrix (each row correspond to one element in X and Y) of auxiliary proteins that tags biological variance in X. We recommend at least 100 columns.
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
#'    \item beta_BPR: Slope of covariance y2z on x2z.
#'    \item beta_BPR_se: standard error of beta_BPR estimate.
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
#'    \item dilution_T_se: se for the dilution test statistics.
#'    \item dilution_test_statistics: < 0 means X only tags part of variance in X.
#'    \item dilution_variance; estimated variance that are added to xi_XY when using Z to tag X (see supplemeetary methods for details, this is the E[Z_xi1^2] term); if this number is 0, then there is no solution for the dilution factor.
#'    \item dilution_BPR_sampling_se: bootstrapping se for the biological proportion estimate (adjusting for dilution).
#'    \item dilution_BPR_sampling_mean: bootstrapping sample mean for the biological proportion estimate (adjusting for dilution).
#'    \item dilution_BPR_Bioprop_xi1_xi2: biological proportion estimate after accounting for that Z only tags part of X.
#'    }
#'
#' @export
#'
#' @examples
#' quick_results <- bioProp_BPR_partialCor_adjustment(X_EXAMPLE, Y_EXAMPLE, Z_BPR = Z_EXAMPLE[,1:200], bootstrapping_number=20)
#'
bioProp_BPR_partialCor_adjustment <- function(X, Y, Z_BPR, Z_proportion_ratio = NULL, bootstrapping_number = 100, covClumpThreshold = 0.5){ # X is N x 1, Y is N x 1, Z is N x J
  print(paste0("Number of auxiliary proteins: ", dim(Z_BPR)[2], "; Recommended at auxiliary protein > 100"))

  X_notNA_idx <- which(!is.na(X))
  print(paste0("Keep ", length(X_notNA_idx), " rows where X is not missing"))
  X <- scale(X[X_notNA_idx])
  Y <- scale(Y[X_notNA_idx])
  Z_BPR <- apply(Z_BPR[X_notNA_idx, ], 2, scale)
  if(is.null(Z_proportion_ratio)){
    Z_proportion_ratio <- Z_BPR
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
  # beta_BPR
  OLS_x2z <- c(cor(X, Z_BPR, use = "pairwise.complete.obs"))
  OLS_y2z <- c(cor(Y, Z_BPR, use = "pairwise.complete.obs"))

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

  BPR_model <- summary(lm(OLS_y2z ~ OLS_x2z))
  beta_BPR <- BPR_model$coefficients[2,1]
  beta_BPR_se <- BPR_model$coefficients[2,2]
  R2_BPR <- BPR_model$adj.r.squared
  # print(paste("Slope of covariance y2z on x2z is: ", beta_BPR, "(", beta_BPR_se, ")", " with z-score ", beta_BPR/beta_BPR_se, ". |Z-score| > 5 required for the estimate to be approximately unbiased."))
  Bioprop_xi1 <- sqrt(PxPy/beta_BPR^2) * R2_BPR
  Partial_beta <- beta_BPR/R2_BPR
  Bioprop_Y <- Bioprop_xi1 * Partial_beta^2

  # compute the covariance matrix
  BPR_correction <- direct_BPR_ratio_estimate(cor_XY_Z, Bioprop_Y = Bioprop_Y)

  # test X Z dilution; using equations from overleaf
  dilution_correction <- BPR_dilution_test(cor_XY_Z = cor_XY_Z, Bioprop_xi1 = Bioprop_xi1, Bioprop_Y = Bioprop_Y)

  bootstrapping_number <- bootstrapping_number
  # cor_XY_Z_bt_all <- array(dim = c( total_col, total_col, bootstrapping_number))

  bt_PxPy <- c()
  bt_beta_BPR <- c()
  bt_R2_BPR <- c()
  bt_Bioprop_xi1 <- c()
  bt_Partial_beta <- c()
  bt_Bioprop_Y <- c()
  bt_E_xi_XY <- c()
  bt_E_xi_X <- c()
  # bt_var_eta_1j <- c()
  # bt_var_eta_2j_E2_xi_X <- c()
  bt_BPR_adjusting_factor <- c()
  bt_BPR_dilution_adjusting_factor <- c()
  bt_dilution_test_statistics <- c()
  for(bt in 1:bootstrapping_number){
    print(paste0("bootstrapping sampling: ", bt))
    try({
      bt_labels <- sample(1:length(X), replace = T, size = length(X))

      PxPy_bt <- cor(X[bt_labels], Y[bt_labels], use = "pairwise.complete.obs")^2
      bt_PxPy <- c(bt_PxPy, PxPy_bt) # save samples for se computation
      OLS_x2z <- c(cor(X[bt_labels], Z_BPR[bt_labels,], use = "pairwise.complete.obs"))
      OLS_y2z <- c(cor(Y[bt_labels], Z_BPR[bt_labels,], use = "pairwise.complete.obs"))

      BPR_model_bt <- summary(lm(OLS_y2z ~ OLS_x2z))
      beta_BPR_bt <- BPR_model_bt$coefficients[2,1]
      bt_beta_BPR <- c(bt_beta_BPR, beta_BPR_bt) # save samples for se computation
      bt_R2_BPR <- c(bt_R2_BPR, BPR_model_bt$adj.r.squared)
      # beta_BPR_bt <- summary(lm(cor_XY_Z_bt[3:total_col, 2] ~ cor_XY_Z_bt[3:total_col, 1]))$coefficients[2,1]
      bt_Bioprop_xi1 <- c(bt_Bioprop_xi1, sqrt(PxPy_bt/beta_BPR_bt^2) * BPR_model_bt$adj.r.squared)
      bt_Partial_beta <- c(bt_Partial_beta, beta_BPR_bt/BPR_model_bt$adj.r.squared)
      Bioprop_Y_bt <- sqrt(PxPy_bt * beta_BPR_bt^2)/BPR_model_bt$adj.r.squared
      bt_Bioprop_Y <- c(bt_Bioprop_Y, Bioprop_Y_bt)

      # perform random clumping
      clumping_idx <-  clump_cov_Z(cov_Z, threshold = covClumpThreshold)
      Z_clumped <- Z_proportion_ratio[, -clumping_idx]
      cor_XY_Z_bt <- padding_XY_Z_cov(X[bt_labels], Y[bt_labels], Z_clumped[bt_labels,], num_padding_XY = 1)
      total_col <- dim(cor_XY_Z_bt)[1]
      BPR_correction_bt <- direct_BPR_ratio_estimate(cor_XY_Z_bt, Bioprop_Y = Bioprop_Y_bt)
      # cor_XY_Z_bt_all[, , bt] <- cor_XY_Z_bt

      # test for partial tagging X -> Z
      BPR_dilution_correction_bt <- BPR_dilution_test(cor_XY_Z_bt, Bioprop_xi1 = sqrt(PxPy_bt/beta_BPR_bt^2) * BPR_model_bt$adj.r.squared,
                                                      Bioprop_Y = Bioprop_Y_bt)
      bt_BPR_dilution_adjusting_factor <- c(bt_BPR_dilution_adjusting_factor, BPR_dilution_correction_bt$adjusting_factor)
      bt_dilution_test_statistics <- c(bt_dilution_test_statistics, BPR_dilution_correction_bt$dilution_testStat)

      # compute the total amount of variance in X explained by Z
      bt_E_xi_XY <- c(bt_E_xi_XY, BPR_correction_bt$E_xi_XY)
      bt_E_xi_X <- c(bt_E_xi_X, BPR_correction_bt$E_xi_X)
      # bt_var_eta_1j <- c(bt_var_eta_1j, BPR_correction_bt$var_eta_1j)
      # bt_var_eta_2j_E2_xi_X <- c(bt_var_eta_2j_E2_xi_X, BPR_correction_bt$var_eta_2j_E2_xi_X)
      bt_BPR_adjusting_factor <- c(bt_BPR_adjusting_factor, BPR_correction_bt$adjusting_factor)
    })
  }
  beta_BPR_se <- sd(bt_beta_BPR) # have to use this se as the covariance are not independent.
  R2_BPR_se <- sd(bt_R2_BPR)
  print(paste("R2 of covariance y2z on x2z is: ", R2_BPR, "(", R2_BPR_se, ")", " with z-score ", R2_BPR/R2_BPR_se, ". |Z-score| > 5 required for the estimate to be approximately unbiased."))
  print(paste("Slope of covariance y2z on x2z is: ", beta_BPR, "(", beta_BPR_se, ")", " with z-score ", beta_BPR/beta_BPR_se, ". |Z-score| > 5 required for the estimate to be approximately unbiased."))
  # comparison below suggest the inverse of variance is the main source of bias; probably safe to proceed with simple subtraction
  # print(c(mean(bt_approx_Bioprop_xi1_xi2 * bt_BPR_adjusting_factor), approx_Bioprop_xi1_xi2 * BPR_correction$adjusting_factor))
  # print(c(mean(bt_E_xi_XY), BPR_correction$E_xi_XY))
  # print(c(mean(bt_E_xi_X), BPR_correction$E_xi_X))
  # print(c(mean(bt_var_eta_1j), BPR_correction$var_eta_1j))
  # print(c(mean(bt_var_eta_2j_E2_xi_X), BPR_correction$var_eta_2j_E2_xi_X))
  # se needs to consider the bias correction; note this is an approximation as we ignore the impact of "approx_Bioprop_xi1_xi2" on se in the bias correction.
  sampling_se <- sd(bt_Bioprop_xi1 * bt_BPR_adjusting_factor) * (1+1/sqrt(bootstrapping_number))
  sampling_mean <- mean(bt_Bioprop_xi1 * bt_BPR_adjusting_factor)
  # use bootstrap to correct for bias in the BPR adjusting procedure
  Bioprop_xi1_xi2 <- Bioprop_xi1 * BPR_correction$adjusting_factor -
    Bioprop_xi1 * BPR_correction$adjusting_factor/sampling_mean *(sampling_mean - Bioprop_xi1 * BPR_correction$adjusting_factor )

  # compute the proportion of variance in Y that are captured by other Zs
  total_col <- dim(cor_XY_Z)[1]
  prop_Y_naive <- cor_XY_Z[2,3:total_col] %*%
    solve( cor_XY_Z[3:total_col, 3:total_col]) %*%
    cor_XY_Z[3:total_col, 2]

  # compute the adjusted proportion of variance in X explained by Z and its standard error
  X_variance_corMatrix <- (BPR_correction$E_xi_XY + BPR_correction$E_xi_X) * 2 -
    (mean(bt_E_xi_XY) + mean(bt_E_xi_X))
  se_X_variance_corMatrix <- sd(bt_E_xi_XY +bt_E_xi_X)

  # test if Z only tags part of the variance in X and provide an alternative estimate that account for this
  dilution_T_se <- sd(bt_dilution_test_statistics)
  dilution_test_statistics <- dilution_correction$dilution_testStat
  print(paste("Test if Z only tags part of X; Z-score: ", dilution_test_statistics/dilution_T_se, "."))

  dilution_BPR_sampling_se <- sd(bt_Bioprop_xi1 * bt_BPR_dilution_adjusting_factor) * (1+1/sqrt(bootstrapping_number))
  dilution_BPR_sampling_mean <- mean(bt_Bioprop_xi1 * bt_BPR_dilution_adjusting_factor)
  # use bootstrap to correct for bias in the BPR adjusting procedure
  dilution_BPR_Bioprop_xi1_xi2 <- Bioprop_xi1 * dilution_correction$adjusting_factor -
    Bioprop_xi1 * dilution_correction$adjusting_factor/dilution_BPR_sampling_mean *
    (dilution_BPR_sampling_mean - Bioprop_xi1 * dilution_correction$adjusting_factor )


  results_biprop <- list()
  results_biprop$Bioprop_xi1_xi2 <- c(Bioprop_xi1_xi2)
  results_biprop$sampling_se <- sampling_se # jackknife (or boostraping) se
  results_biprop$sampling_mean <- sampling_mean # jackknife (or boostraping) mean, use to identify bias
  results_biprop$jk_Zlabel_se_bioprop <- jk_se_bioprop # reporting the se estimated from jacknifing Zj index
  results_biprop$PxPy <- PxPy
  results_biprop$PxPy_se <- sd(bt_PxPy)
  results_biprop$PxPy_sample_mean <- mean(bt_PxPy)
  results_biprop$beta_BPR <- beta_BPR
  results_biprop$beta_BPR_se <- sd(bt_beta_BPR)
  results_biprop$beta_BPR_sample_mean <- mean(bt_beta_BPR)
  results_biprop$R2_BPR <- R2_BPR
  results_biprop$R2_BPR_se <- R2_BPR_se
  results_biprop$R2_BPR_sample_mean <- mean(bt_R2_BPR)
  results_biprop$Bioprop_xi1 <- Bioprop_xi1
  results_biprop$Partial_beta <- Partial_beta
  results_biprop$Bioprop_Y <- Bioprop_Y
  results_biprop$sample_Bioprop_xi1 <- mean(bt_Bioprop_xi1)
  results_biprop$sample_Partial_beta <- mean(bt_Partial_beta)
  results_biprop$sample_Bioprop_Y <- mean(bt_Bioprop_Y)
  results_biprop$adjusting_factor <- BPR_correction$adjusting_factor
  results_biprop$sample_adjusting_factor <- mean(bt_BPR_adjusting_factor) # this is the factor used to detect bias in the adjusting_factor
  results_biprop$E_xi_XY <- BPR_correction$E_xi_XY
  results_biprop$E_xi_X <- BPR_correction$E_xi_X
  results_biprop$sample_E_xi_XY <- mean(bt_E_xi_XY)
  results_biprop$sample_E_xi_X <- mean(bt_E_xi_X)
  results_biprop$X_variance_corMatrix <- X_variance_corMatrix
  results_biprop$se_X_variance_corMatrix <- se_X_variance_corMatrix
  results_biprop$sigma_X <- 1 - X_variance_corMatrix
  results_biprop$sigma_Y <- 1 - prop_Y_naive

  # save the alternative estimate when Z partially tags X
  results_biprop$dilution_T_se <- dilution_T_se
  results_biprop$dilution_test_statistics <- dilution_test_statistics
  results_biprop$dilution_variance <- dilution_correction$dilution_variance # if this number is 0, then there is no solution for the dilution factor
  results_biprop$dilution_BPR_sampling_se <- dilution_BPR_sampling_se
  results_biprop$dilution_BPR_sampling_mean <- dilution_BPR_sampling_mean
  # use bootstrap to correct for bias in the BPR adjusting procedure
  results_biprop$dilution_BPR_Bioprop_xi1_xi2 <- dilution_BPR_Bioprop_xi1_xi2

  return(results_biprop)
}


#' Title: BPR wrapper function. Input is many pairs of X and Y, and all based on the same set of Z.
#' Paralle computation of all bootstrap samples to speed up the computation; this function is much faster for many pairs of X and Y.
#' @param target_traits a matrix (nrow= number of samples; ncolumns=number of target traits), each column is one target trait whose biological proportion is to be estimated.
#' @param helper_traits a matrix of the same size as target_traits; each column is the helper trait for the corresponding target trait; column index should match.
#' @param auxiliary_traits matrix (each row correspond to one element in target proteins) of auxiliary proteins that tags biological variance in target proteins. We recommend at least 100 columns.
#' @param bootstrapping_number number of boostrap samples for estimate uncertainty.
#' @param covClumpThreshold when computing the partial correlation matrix, clumping highly correlated pairs of Zs; default threshold cor>0.5.
#'
#' @returns A dataframe contains same output as the `bioProp_BPR_partialCor_adjustment` function.
#' @export
#'
#' @examples
#' full_run_results_big_wrapper <- BPR_fast_wrapper(cbind(X_EXAMPLE,X_EXAMPLE, X_EXAMPLE), cbind(Y_EXAMPLE,Y_EXAMPLE, Z_EXAMPLE[,1]), Z_EXAMPLE, bootstrapping_number = 5)
BPR_fast_wrapper <- function(target_traits,  helper_traits, auxiliary_traits, bootstrapping_number = 20, covClumpThreshold = 0.5){

  print("This is a fast function for computing biological proportion of X for multiple X-Y pairs using the same set of Zs;")
  print("parallelisation assumes the same missingness across X-Y pairs, therefore will be slightly inaccurate if some pairs of X-Y has ")
  print("substantially more missingness compare to other X-Y pairs. In those cases, we suggest using the single X-Y pair function. ")
  target_traits <- as.matrix(target_traits)
  helper_traits <- as.matrix(helper_traits)
  if(! identical(dim(target_traits), dim(helper_traits) )){
    stop("dimension of target traits and helper traits is not matching. Make sure you have the same number of individuals (nrow) and number of X, Y pairs (ncol).")
  }

  if(dim(target_traits)[1] != dim(target_traits)[1]){
    stop("Number of rows are not matching across target_traits,  helper_traits, and auxiliary_traits.")
  }

  correlation_target_aux <- cor(target_traits, auxiliary_traits, use = "pairwise.complete.obs")
  if(max(abs(correlation_target_aux)) > 0.99){
    stop("Really high correlation between target protein and on auxiliary protein. Did you accidentally put the target protein in the auxiliary traits?")
  }
    if(max(abs(correlation_target_aux)) > 0.99){
    stop("Really high correlation between target protein and on auxiliary protein. Did you accidentally put the target protein in the auxiliary traits?")
  }

  if(max(abs(correlation_target_aux)) < 0.05){
    stop("target protein doesn't seem to be correlated with any of auxiliary proteins. Have you checked the row ordering is consistent between target and auxiliary traits?")
  }

  # also needs to get rid of Y that is the same as Zs
  # if Y = one of the Zs, remove those from Z
  correlation_helper_aux <- cor(helper_traits, auxiliary_traits, use = "pairwise.complete.obs")
  if(max(abs(correlation_helper_aux)) > 0.99 ){
    z_BPR_rm <- c()
    for(Y_idx in 1:dim(helper_traits)[2]){
      idx_replicate_Y <- which(max(abs(correlation_helper_aux[Y_idx])) > 0.99 )
      if(length(idx_replicate_Y) > 0){
        print(paste0("Z index ", idx_replicate_Y, " is perfectly correlated with Y index ", Y_idx, ", remove it from Zs."))
        z_BPR_rm <- c(z_BPR_rm, idx_replicate_Y)
      }
    }
    auxiliary_traits <- auxiliary_traits[, -z_BPR_rm]
  }

  # perform clumping and save them
  Z_BPR <- apply(auxiliary_traits, 2, scale)
  Z_proportion_ratio <- Z_BPR
  cor_Z <- cor(Z_proportion_ratio, use = "pairwise.complete.obs")

  # for each bootstrapping samples, do the clumping of cor_Z and bootstrapping idx
  bt_labels_list <- list()
  clumping_idx_list <- list()
  cor_Z_clumped_list <- list()
  for(bt in 1:bootstrapping_number){
    bt_labels <- sample(1:dim(target_traits)[1], replace = T, size = dim(target_traits)[1])
    bt_labels_list[[bt]] <- bt_labels
    clumping_idx <-  clump_cov_Z(cor_Z, threshold = covClumpThreshold)
    clumping_idx_list[[bt]] <- clumping_idx
    Z_clumped <- Z_proportion_ratio[, -clumping_idx]
    cor_Z_clumped_list[[bt]] <- cor(Z_clumped, use = "pairwise.complete.obs")
  }

  # do the analysis for each of the X-Y pairs
  number_pairs <- dim(target_traits)[2]
  results_X_Y_pairs <- list()
  for(pair_idx in 1:number_pairs){
    X <- target_traits[, pair_idx]
    Y <- helper_traits[, pair_idx]

    X_notNA_idx <- which(!is.na(X))
    print("----------------------------------------------------------------------------- ")
    print(paste0("Running X-Y pair ", pair_idx, ". Keep ", length(X_notNA_idx), " rows where X is not missing."))
    print("----------------------------------------------------------------------------- ")
    X <- scale(X[X_notNA_idx])
    Y <- scale(Y[X_notNA_idx])
    Z_BPR <- apply(auxiliary_traits[X_notNA_idx, ], 2, scale)
    Z_proportion_ratio <- Z_BPR
    # still need to compute the point estimate accounting for the missingness; only sped up the bootstrapping step
    cor_XY_Z <- padding_XY_Z_cov(X, Y, Z_proportion_ratio, num_padding_XY = 1)
    total_col <- dim(cor_XY_Z)[1]
    cov_Z <- cor_XY_Z[3:total_col, 3:total_col]

    print(paste("Perform covariance regression."))
    PxPy_model <- summary(lm(scale(X) ~ scale(Y)))
    PxPy <- PxPy_model$coefficients[2,1]^2
    # beta_BPR
    OLS_x2z <- c(cor(X, Z_BPR, use = "pairwise.complete.obs"))
    OLS_y2z <- c(cor(Y, Z_BPR, use = "pairwise.complete.obs"))

    # start a 10 random clumping and test if it is robust
    num_jk_block <- 20
    print(paste("Jackknife across Z index, using ", num_jk_block, " samples."))
    jk_initial_samples <- jackknife_acrossZ_partialCor(cor_XY_Z, OLS_y2z, OLS_x2z, PxPy,
                                                       num_jk_block = num_jk_block, covClumpThreshold = covClumpThreshold)
    jk_se_bioprop <- sd(jk_initial_samples$dilution_BPR_Bioprop_xi1_xi2) * sqrt(num_jk_block)
    print(paste("Jackknife se: ", jk_se_bioprop, ". Z-score = ", mean(jk_initial_samples$dilution_BPR_Bioprop_xi1_xi2)/jk_se_bioprop))

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

    BPR_model <- summary(lm(OLS_y2z ~ OLS_x2z))
    beta_BPR <- BPR_model$coefficients[2,1]
    beta_BPR_se <- BPR_model$coefficients[2,2]
    R2_BPR <- BPR_model$adj.r.squared
    Bioprop_xi1 <- sqrt(PxPy/beta_BPR^2) * R2_BPR
    Partial_beta <- beta_BPR/R2_BPR
    Bioprop_Y <- Bioprop_xi1 * Partial_beta^2

    # compute the covariance matrix
    BPR_correction <- direct_BPR_ratio_estimate(cor_XY_Z, Bioprop_Y = Bioprop_Y)

    # test X Z dilution; using equations from overleaf
    dilution_correction <- BPR_dilution_test(cor_XY_Z = cor_XY_Z, Bioprop_xi1 = Bioprop_xi1, Bioprop_Y = Bioprop_Y)

    # bootstrap: using the same set of bootstrapping index for each X-Y pairs
    bootstrapping_number <- bootstrapping_number

    bt_PxPy <- c()
    bt_beta_BPR <- c()
    bt_R2_BPR <- c()
    bt_Bioprop_xi1 <- c()
    bt_Partial_beta <- c()
    bt_Bioprop_Y <- c()
    bt_E_xi_XY <- c()
    bt_E_xi_X <- c()
    # bt_var_eta_1j <- c()
    # bt_var_eta_2j_E2_xi_X <- c()
    bt_BPR_adjusting_factor <- c()
    bt_BPR_dilution_adjusting_factor <- c()
    bt_dilution_test_statistics <- c()

    # has to assign X, Y, Z using original data frame as the bootstrapping index were sampled using original data
    X <- scale(target_traits[, pair_idx])
    Y <- scale(helper_traits[, pair_idx])
    Z_BPR <- apply(auxiliary_traits, 2, scale)
    Z_proportion_ratio <- Z_BPR
    for(bt in 1:bootstrapping_number){
      print(paste0("bootstrapping sampling: ", bt))
      try({
        bt_labels <- bt_labels_list[[bt]]
        PxPy_bt <- cor(X[bt_labels], Y[bt_labels], use = "pairwise.complete.obs")^2
        bt_PxPy <- c(bt_PxPy, PxPy_bt) # save samples for se computation
        OLS_x2z <- c(cor(X[bt_labels], Z_BPR[bt_labels,], use = "pairwise.complete.obs"))
        OLS_y2z <- c(cor(Y[bt_labels], Z_BPR[bt_labels,], use = "pairwise.complete.obs"))

        BPR_model_bt <- summary(lm(OLS_y2z ~ OLS_x2z))
        beta_BPR_bt <- BPR_model_bt$coefficients[2,1]
        bt_beta_BPR <- c(bt_beta_BPR, beta_BPR_bt) # save samples for se computation
        bt_R2_BPR <- c(bt_R2_BPR, BPR_model_bt$adj.r.squared)
        # beta_BPR_bt <- summary(lm(cor_XY_Z_bt[3:total_col, 2] ~ cor_XY_Z_bt[3:total_col, 1]))$coefficients[2,1]
        bt_Bioprop_xi1 <- c(bt_Bioprop_xi1, sqrt(PxPy_bt/beta_BPR_bt^2) * BPR_model_bt$adj.r.squared)
        bt_Partial_beta <- c(bt_Partial_beta, beta_BPR_bt/BPR_model_bt$adj.r.squared)
        Bioprop_Y_bt <- sqrt(PxPy_bt * beta_BPR_bt^2)/BPR_model_bt$adj.r.squared
        bt_Bioprop_Y <- c(bt_Bioprop_Y, Bioprop_Y_bt)

        # get the clumped matrix
        clumping_idx <-  clumping_idx_list[[bt]]
        cor_Z_bt <- cor_Z_clumped_list[[bt]]
        Z_clumped <- Z_proportion_ratio[, -clumping_idx]
        cor_xy2z <- cor(cbind(X[bt_labels],Y[bt_labels]), Z_clumped[bt_labels,], use = "pairwise.complete.obs")

        a <- matrix(c(1,sqrt(PxPy_bt), sqrt(PxPy_bt), 1), nrow = 2, ncol = 2)
        ab <- rbind(a, t(cor_xy2z))
        bc <- rbind(cor_xy2z, cor_Z_bt)
        cor_XY_Z_bt <- cbind(ab, bc)

        # cor_XY_Z_bt_compute <- padding_XY_Z_cov(X[bt_labels], Y[bt_labels], Z_clumped[bt_labels,], num_padding_XY = 1)
        # plot(cor_XY_Z_bt, cor_XY_Z_bt_compute)
        # cor_XY_Z_bt_compute[1:4, 1:4]
        # cor_XY_Z_bt[1:4, 1:4]

        BPR_correction_bt <- direct_BPR_ratio_estimate(cor_XY_Z_bt, Bioprop_Y = Bioprop_Y_bt)

        # test for partial tagging X -> Z
        BPR_dilution_correction_bt <- BPR_dilution_test(cor_XY_Z_bt, Bioprop_xi1 = sqrt(PxPy_bt/beta_BPR_bt^2) * BPR_model_bt$adj.r.squared,
                                                        Bioprop_Y = Bioprop_Y_bt)
        bt_BPR_dilution_adjusting_factor <- c(bt_BPR_dilution_adjusting_factor, BPR_dilution_correction_bt$adjusting_factor)
        bt_dilution_test_statistics <- c(bt_dilution_test_statistics, BPR_dilution_correction_bt$dilution_testStat)

        # compute the total amount of variance in X explained by Z
        bt_E_xi_XY <- c(bt_E_xi_XY, BPR_correction_bt$E_xi_XY)
        bt_E_xi_X <- c(bt_E_xi_X, BPR_correction_bt$E_xi_X)
        # bt_var_eta_1j <- c(bt_var_eta_1j, BPR_correction_bt$var_eta_1j)
        # bt_var_eta_2j_E2_xi_X <- c(bt_var_eta_2j_E2_xi_X, BPR_correction_bt$var_eta_2j_E2_xi_X)
        bt_BPR_adjusting_factor <- c(bt_BPR_adjusting_factor, BPR_correction_bt$adjusting_factor)
      })
    }
    beta_BPR_se <- sd(bt_beta_BPR) # have to use this se as the covariance are not independent.
    R2_BPR_se <- sd(bt_R2_BPR)
    print(paste("R2 of covariance y2z on x2z is: ", R2_BPR, "(", R2_BPR_se, ")", " with z-score ", R2_BPR/R2_BPR_se, ". |Z-score| > 5 required for the estimate to be approximately unbiased."))
    print(paste("Slope of covariance y2z on x2z is: ", beta_BPR, "(", beta_BPR_se, ")", " with z-score ", beta_BPR/beta_BPR_se, ". |Z-score| > 5 required for the estimate to be approximately unbiased."))

    sampling_se <- sd(bt_Bioprop_xi1 * bt_BPR_adjusting_factor) * (1+1/sqrt(bootstrapping_number))
    sampling_mean <- mean(bt_Bioprop_xi1 * bt_BPR_adjusting_factor)
    # use bootstrap to correct for bias in the BPR adjusting procedure
    Bioprop_xi1_xi2 <- Bioprop_xi1 * BPR_correction$adjusting_factor -
      Bioprop_xi1 * BPR_correction$adjusting_factor/sampling_mean *(sampling_mean - Bioprop_xi1 * BPR_correction$adjusting_factor )

    # compute the proportion of variance in Y that are captured by other Zs
    total_col <- dim(cor_XY_Z)[1]
    prop_Y_naive <- cor_XY_Z[2,3:total_col] %*%
      solve( cor_XY_Z[3:total_col, 3:total_col]) %*%
      cor_XY_Z[3:total_col, 2]

    # compute the adjusted proportion of variance in X explained by Z and its standard error
    X_variance_corMatrix <- (BPR_correction$E_xi_XY + BPR_correction$E_xi_X) * 2 -
      (mean(bt_E_xi_XY) + mean(bt_E_xi_X))
    se_X_variance_corMatrix <- sd(bt_E_xi_XY +bt_E_xi_X)

    # test if Z only tags part of the variance in X and provide an alternative estimate that account for this
    dilution_T_se <- sd(bt_dilution_test_statistics)
    dilution_test_statistics <- dilution_correction$dilution_testStat
    print(paste("Test if Z only tags part of X; Z-score: ", dilution_test_statistics/dilution_T_se, "."))

    dilution_BPR_sampling_se <- sd(bt_Bioprop_xi1 * bt_BPR_dilution_adjusting_factor) * (1+1/sqrt(bootstrapping_number))
    dilution_BPR_sampling_mean <- mean(bt_Bioprop_xi1 * bt_BPR_dilution_adjusting_factor)
    # use bootstrap to correct for bias in the BPR adjusting procedure
    dilution_BPR_Bioprop_xi1_xi2 <- Bioprop_xi1 * dilution_correction$adjusting_factor -
      Bioprop_xi1 * dilution_correction$adjusting_factor/dilution_BPR_sampling_mean *
      (dilution_BPR_sampling_mean - Bioprop_xi1 * dilution_correction$adjusting_factor )


    results_biprop <- list()
    results_biprop$Bioprop_xi1_xi2 <- c(Bioprop_xi1_xi2)
    results_biprop$sampling_se <- sampling_se # jackknife (or boostraping) se
    results_biprop$sampling_mean <- sampling_mean # jackknife (or boostraping) mean, use to identify bias
    results_biprop$jk_Zlabel_se_bioprop <- jk_se_bioprop # reporting the se estimated from jacknifing Zj index
    results_biprop$PxPy <- PxPy
    results_biprop$PxPy_se <- sd(bt_PxPy)
    results_biprop$PxPy_sample_mean <- mean(bt_PxPy)
    results_biprop$beta_BPR <- beta_BPR
    results_biprop$beta_BPR_se <- sd(bt_beta_BPR)
    results_biprop$beta_BPR_sample_mean <- mean(bt_beta_BPR)
    results_biprop$R2_BPR <- R2_BPR
    results_biprop$R2_BPR_se <- R2_BPR_se
    results_biprop$R2_BPR_sample_mean <- mean(bt_R2_BPR)
    results_biprop$Bioprop_xi1 <- Bioprop_xi1
    results_biprop$Partial_beta <- Partial_beta
    results_biprop$Bioprop_Y <- Bioprop_Y
    results_biprop$sample_Bioprop_xi1 <- mean(bt_Bioprop_xi1)
    results_biprop$sample_Partial_beta <- mean(bt_Partial_beta)
    results_biprop$sample_Bioprop_Y <- mean(bt_Bioprop_Y)
    results_biprop$adjusting_factor <- BPR_correction$adjusting_factor
    results_biprop$sample_adjusting_factor <- mean(bt_BPR_adjusting_factor) # this is the factor used to detect bias in the adjusting_factor
    results_biprop$E_xi_XY <- BPR_correction$E_xi_XY
    results_biprop$E_xi_X <- BPR_correction$E_xi_X
    results_biprop$sample_E_xi_XY <- mean(bt_E_xi_XY)
    results_biprop$sample_E_xi_X <- mean(bt_E_xi_X)
    results_biprop$X_variance_corMatrix <- X_variance_corMatrix
    results_biprop$se_X_variance_corMatrix <- se_X_variance_corMatrix
    results_biprop$sigma_X <- 1 - X_variance_corMatrix
    results_biprop$sigma_Y <- 1 - prop_Y_naive

    # save the alternative estimate when Z partially tags X
    results_biprop$dilution_T_se <- dilution_T_se
    results_biprop$dilution_test_statistics <- dilution_test_statistics
    results_biprop$dilution_variance <- dilution_correction$dilution_variance # if this number is 0, then there is no solution for the dilution factor
    results_biprop$dilution_BPR_sampling_se <- dilution_BPR_sampling_se
    results_biprop$dilution_BPR_sampling_mean <- dilution_BPR_sampling_mean
    # use bootstrap to correct for bias in the BPR adjusting procedure
    results_biprop$dilution_BPR_Bioprop_xi1_xi2 <- dilution_BPR_Bioprop_xi1_xi2

    results_X_Y_pairs[[pair_idx]] <- bind_rows(results_biprop) %>%
      mutate(X_Y_pair_idx = pair_idx, X_idx = colnames(target_traits)[pair_idx], Y_idx = colnames(helper_traits)[pair_idx]) %>%
      select(X_Y_pair_idx, X_idx, Y_idx, everything())
  }


  return(bind_rows(results_X_Y_pairs))
}







