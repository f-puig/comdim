SelectFeaturesRW <- function(RW = RW, results = results, ndim = NULL, blocks = NULL, threshold_cor = 1, threshold_cov = 1, mean.RW = T, plots = "NO") {

#' SelectFeaturesRW
#'
#' Finds the important variables presenting a coordinated response in all the replicate-blocks
#' @param RW The object used as input in ComDim_PCA().
#' @param results The object obtained in ComDim_PCA().
#' @param ndim The number of the component wor which the most important variables are to be calculated.
#' @param blocks A vector with the indices for the replicate blocks of the same data type.
#' @param threshold_cor The "times" parameter used to calculate the threshold in the following formula: cor(variable) > times * sd(cor(variables)). Minimal value that can be assigned to threshold_cor is 1.
#' @param threshold_cov The "times" parameter used to calculate the threshold in the following formula: cov(variable) > times * sd(cov(variables)). Minimal value that can be assigned to threshold_cor is 1.
#' @param mean.RW Logical value to indicate whether the RW data must be mean-centered (TRUE) or not (FALSE).
#' @param plots Parameter to indicate whether graphical representation must be plotted. Possible values are "NO" for no plots, "separated" for plotting each plot individually, and "together" to plot all the plots in the same grid.
#' @return An object with 2 lists. The first list contains the important variables presenting a positive relationship with the scores, while the second list contains the most important variables presenting a negative relationship.
#' @examples
#'  b1 = matrix(rnorm(500),10,50)
#'  batch_b1 = rep(1,10)
#'  b2 = matrix(rnorm(800),30,80)
#'  batch_b2 = c(rep(1,10),rep(2,10),rep(3,10))
#'  mb <- BuildMultiBlock(b1, batches = batch_b1)
#'  mb <- BuildMultiBlock(b2, growingMB = mb, batches = batch_b2, equalSampleNumber = FALSE)
#'  rw <- SplitRW(mb)
#'  results<-ComDim_PCA(rw, 2) # In this analysis, we used 2 components.
#'  features <- SelectFeaturesRW(RW = rw, results = results, ndim = 1, blocks = c(2,3,4))
#' @export

  # Check that blocks is a numerical vector
  if (!(is.numeric(blocks))){
    stop('The parameter blocks must be a numerical vector.')
  }

  # Check data

  ntable = length(RW) # number of blocks

  # Checking that RW is a suitable multi-block structure
  properties <- c('Data', 'Variables', 'Samples')

  give_error <- 0

  sample_number <- as.vector(NULL)
  variable_number <- as.vector(NULL)

  if (length(RW) < max(blocks)){
    stop("The list in 'RW' is incorrect.")
  }

  for (i in 1:ntable){

    sample_number[i] <- length(RW[[1]]$Samples)
    variable_number[i] <- length(RW[[i]]$Variables)
    if(sum(properties %in% names(RW[[i]])) != 3){ # Checking if a field is missing
      print(sprintf('Block %s has missing fields.', i))
      give_error <- 1
    }

    if (nrow(RW[[i]]$Data) != length(RW[[i]]$Samples)) {
      print(sprintf('Block %s has incorrect number of samples.', i))
      give_error <- 1
    }

    if (ncol(RW[[i]]$Data) != length(RW[[i]]$Variables)) {
      print(sprintf('Block %s has incorrect number of variables.', i))
      give_error <- 1
    }
  }

  if (length(unique(sample_number[blocks]))!= 1){
    stop('The replicate blocks to evaluate have a different sample size.')
  }

  if (length(unique(variable_number[blocks]))!= 1){
    stop('The replicate blocks to evaluate have a different variable size.')
  }

  if (give_error) {
    stop('The data is not ready for ComDim.')
  }

  # Check if result contain the Global scores
  if(!("Q" %in% names(results))){
    stop('Resuts do not contain the Global scores.')
  } else {
    if(ncol(results$Q$Data) < ndim) {
      stop("The 'ndim' value chosen for this analysis is too high.")
    }
    if(nrow(results$Q$Data) != unique(sample_number)){
      stop("The size of the local scores do not match the size of 'RW'.")
    }
  }

  # Check if result contain the blocks of local loadings
  if (!("P_Loc" %in% names(results))){
    stop('Resuts do not contain the Local loading blocks.')
  } else {
    if (length(results$P_Loc$Data) < max(blocks)){
      stop('The list of blocks to compare is incorrect. Check the local loadings.')
    } else {
      for (i in 1:length(blocks)){
        if(ncol(results$P_Loc$Data[[blocks[i]]]) < ndim) {
          stop("The 'ndim' value chosen for this analysis is too high.")
        }
        if(nrow(results$P_Loc$Data[[blocks[i]]]) != ncol(RW[[blocks[i]]]$Data)){
          stop("The size of the local loadings do not match the size of 'RW'.")
        }
      }
    }
  }

  if (!(is.logical(mean.RW))){
    stop("The parameter 'mean.RW' must be logical.")
  }

  if(!(is.character(plots))){
    if(is.logical(plots)){
      if (plots == FALSE){
        plots <- "NO"
      } else {
        plots <- "separated"
      }
    }
  }

  if(!(requireNamespace("pracma", quietly = TRUE))){
    warning("The 'pracma' package is needed but not installed." )
  }

  if(plots != "NO" && !(requireNamespace("ggplot2", quietly = TRUE))){
    warning("The 'ggplot2' package is needed but not installed." )
    plots <- "NO"
  }

  if(plots != "NO" && !(requireNamespace("gridExtra", quietly = TRUE))){
    warning("The 'ggplot2' package is needed but not installed." )
    plots <- "separated"
  }

  mean_blocks <- list(NULL)
  for (i in 1:length(blocks)){
    # Mean-center
    if (isTRUE(mean.RW) && requireNamespace("pracma", quietly = TRUE)){
      mean_blocks[[i]] <- RW[[blocks[i]]]$Data - pracma::repmat(colMeans(RW[[blocks[i]]]$Data),unique(sample_number),1)
    } else {
      mean_blocks[[i]] <- RW[[blocks[i]]]$Data
    }
  }

  if(length(ndim) != 1){
    stop("The length of 'ndim' must be of 1.")
  }

  if(!(is.numeric(ndim))){
    stop("'ndim' must be numeric.")
  }

  s1 <- matrix(, nrow = length(blocks), ncol=unique(variable_number[blocks]))
  s2 <- matrix(, nrow = length(blocks), ncol=unique(variable_number[blocks]))

  # Scale Q in case it wasnt (if the data comes from ComDim, it was).
  if (requireNamespace("pracma", quietly = TRUE)){
    qs <- results$Q$Data[,ndim]%*%pracma::pinv(t(results$Q$Data[,ndim])%*%results$Q$Data[,ndim])
  } else {
    qs <- results$Q$Data[,ndim]
    warning("If 'results' was not scaled, the output from this function may be incorrect.")
  }

  # Calculate the corr and cov values.
  for (i in 1:length(blocks)){
    # Calculate the s-plot values.
    #s1(i,:) = (T' * X_Cent) ./(N - 1); % cov
    #s2(i,:) = s1(i,:) ./ (std(T) .* std(X_Cent)); %cor
    s1[i,] <- (t(qs) %*% mean_blocks[[i]])/(unique(sample_number)-1)
    s2[i,] <- s1[i,]/(stats::sd(qs)*apply(mean_blocks[[i]], 2, stats::sd))
  }

  s1[is.nan(s1)] <- 0
  s2[is.nan(s2)] <- 0
  s1_imp <- list(NULL)
  s2_imp <- list(NULL)
  for (i in 1:length(blocks)){
    if (is.numeric(threshold_cov) && length(threshold_cov) == 1){
      if (threshold_cov < 1) {
        threshold_cov <- 1
        warning("'threshold_cov value has been modified to 1.")
      }
      std_cov <- stats::sd(s1[i,])*threshold_cov
      s1_imp[[i]] <- c(which(s1[i,] > std_cov), which(s1[i,] < -std_cov)) # Covariate and anticovariate
    } else {
      stop("'Threshold_cov' must be a numeric object of length 1.")
    }

    if (is.numeric(threshold_cor) && length(threshold_cor) == 1){
      if (threshold_cor < 1) {
        threshold_cor <- 1
        warning("'threshold_cor value has been modified to 1.")
      }
      std_cor <- stats::sd(s2[i,])*threshold_cor
      s2_imp[[i]] <- c(which(s2[i,] > std_cor), which(s2[i,] < -std_cor))
    } else {
      stop("'Threshold_cor' must be a numeric object of length 1.")
    }
  }
  imp_pos <- list(NULL)
  imp_neg <- list(NULL)
  for(i in 1:length(blocks)){
    imp <- intersect(s1_imp[[i]],s2_imp[[i]])
    imp_pos[[i]] <- intersect(imp, which(results$P_Loc$Data[[blocks[i]]] > 0))
    imp_neg[[i]] <- intersect(imp, which(results$P_Loc$Data[[blocks[i]]] < 0))
  }

  common_RW_pos <- list(NULL)
  common_RW_neg <- list(NULL)
  for(i in 1:length(blocks)){
    if (i == 1){
      common_RW_pos <- imp_pos[[i]]
      common_RW_neg <- imp_neg[[i]]
    } else {
      common_RW_pos <- intersect(common_RW_pos, imp_pos[[i]])
      common_RW_neg <- intersect(common_RW_neg, imp_neg[[i]])
    }
  }

  names(common_RW_pos) <- RW[[blocks[1]]]$Variables[common_RW_pos]
  names(common_RW_neg) <- RW[[blocks[1]]]$Variables[common_RW_neg]
  ivs <- list (positive = common_RW_pos, negative = common_RW_neg)

  if (grepl('NO', plots, ignore.case = TRUE)) {
    # Do nothing.
  }else if (grepl('separated', plots, ignore.case = TRUE) && (requireNamespace("ggplot2", quietly = TRUE))){
    for(i in 1:length(blocks)){
      xx <- data.frame(Covariance = s1[i,], Correlation = s2[i,])
      color<-rep('#000000', variable_number[blocks[i]])
      color[c(imp_pos[[i]],imp_neg[[i]])] <- rep('#FF0000', length(c(imp_pos[[i]],imp_neg[[i]])))
      labels <- rep("",variable_number[blocks[i]])
      labels[c(imp_pos[[i]],imp_neg[[i]])] <- RW[[blocks[i]]]$Variables[c(imp_pos[[i]],imp_neg[[i]])]
      useless<-readline(prompt=sprintf("S-plot for block %d shown. Press [enter] to continue", blocks[i]))
      print(ggplot2::ggplot(xx, ggplot2::aes(x= Covariance, y= Correlation)) +
        ggplot2::geom_point(color = color) +
        ggplot2::geom_text(label=labels) +
        ggplot2::ggtitle(sprintf("S-plot for Block %d and component %i", blocks[i], ndim)))
    }
  }else if (grepl('together', plots, ignore.case = TRUE)  && (requireNamespace("ggplot2", quietly = TRUE))){
    pi <- list(NULL)
    colori <- list(NULL)
    labelsi <- list(NULL)
    xx <- list(NULL)
    for(i in 1:length(blocks)){
      colori[[i]] <- rep('#000000', variable_number[blocks[i]])
      colori[[i]][c(imp_pos[[i]],imp_neg[[i]])] <- rep('#FF0000', length(c(imp_pos[[i]],imp_neg[[i]])))
      xx[[i]] <- data.frame(Covariance = s1[i,], Correlation = s2[i,])
      labelsi[[i]] <- rep("",variable_number[blocks[i]])
      labelsi[[i]][c(imp_pos[[i]],imp_neg[[i]])] <- RW[[blocks[i]]]$Variables[c(imp_pos[[i]],imp_neg[[i]])]
      pi[[i]] <- ggplot2::ggplot(xx[[i]], ggplot2::aes(x= Covariance, y= Correlation)) +
        ggplot2::geom_point(color = colori[[i]]) +
        ggplot2::geom_text(label=labelsi[[i]]) +
        ggplot2::ggtitle(sprintf("S-plot for Block %d and component %i", blocks[i], ndim))
    }
    to_eval <- "gridExtra::grid.arrange("
    for(i in 1:length(blocks)){
      to_eval <- paste0(to_eval, sprintf("pi[[%d]], ",i))
    }
    to_eval <-paste0(to_eval, " nrow = 1)")
    eval(parse(text = to_eval))
  }else{
    warning("The 'plots' parameter was not set correctly. Plots will not be displayed.")
  }

  return(ivs) # ivs = important variables
}
