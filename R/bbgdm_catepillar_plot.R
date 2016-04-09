#' bbgdm catepillar plots function
#' @param object a \code{bbgdm} object
#' @param pars name of parameters to plot.
#' @param pars_labels vector of parameter labels.
#' @param zero_line logical. Whether or not to the zero #' medians.
#' @param horizontal logical. Whether or not you would like the lines to be
#' horizontal

plot.bbgdm.caterpillar <- function(object, pars=NULL, pars_labels = NULL,
                              horizontal = TRUE, 
                              zero_line=TRUE,
                              xlims = NULL){
  # Extract all simulations
  sims <- as.data.frame(object$all.coefs.se)
  
  # Extract only desired parameters
  par_names <- names(sims)
  if(is.null(pars)){
    sims_subset <- sims
    cat('plotting all spline bases \n')
  } else {
  sims_subset <- sims[, par_names %in% pars]
  }
  if (ncol(sims_subset) == 0) {
    stop("No parameters selected. \n", call. = FALSE)
  }
  # Gather for plotting
  gathered <- tidyr::gather(sims_subset, variable, value)
  
  # Add labels
  if (!is.null(pars_labels)) {
    message("\nEnsure that your parameter labels are in the same order as the parameters.\n")
    if (length(pars_labels) !=
        length(unique(gathered$variable))) {
      stop("pars_labels must equal the number of plotted parameters.",
           call. = FALSE)
    }
    gathered$variable <- factor(gathered$variable,
                                labels = pars_labels)
  }
  
  # Find central interval
    gathered <- dplyr::group_by(gathered, variable)
    lower95 <- dplyr::summarise(gathered, quantile(value, 0.025))
    lower90 <- dplyr::summarise(gathered, quantile(value, 0.05))
    upper90 <- dplyr::summarise(gathered, quantile(value, 0.95))
    upper95 <- dplyr::summarise(gathered, quantile(value, 0.975))
  
  # Find medians
  medians <- dplyr::summarise(gathered, median(value))
  
  # Merge
  comb <- suppressMessages(dplyr::inner_join(lower95, lower90))
  comb <- suppressMessages(dplyr::inner_join(comb, medians))
  comb <- suppressMessages(dplyr::inner_join(comb, upper90))
  comb <- suppressMessages(dplyr::inner_join(comb, upper95))
  names(comb) <- c('pars', 'lower95', 'lower90', 'medians', 'upper90',
                   'upper95')
  
  # Plot
  pp <- ggplot2::ggplot(comb, ggplot2::aes(x = medians, y = pars,
                           xmin = lower95,
                           xmax = upper95)) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_segment(ggplot2::aes(x = lower95, xend = upper95, yend = pars),
                     size = 0.5) +
        ggplot2::geom_segment(ggplot2::aes(x = lower90, xend = upper90,
                         yend = pars), size = 1.5) +
        ggplot2::xlab('parameter quantiles of bayesian bootstrap samples')+
        ggplot2::ylab('') +
        if(!is.null(xlims))ggplot2::xlim(xlims)+
        ggthemes::theme_few()+
        ggplot2::theme(plot.margin=unit(c(5,5,5,5),"mm"))
  
  if (isTRUE(zero_line)) pp <- pp + ggplot2::geom_vline(xintercept = 0)
  if (!isTRUE(horizontal)) pp <- pp + ggplot2::coord_flip()
  return(pp)
}