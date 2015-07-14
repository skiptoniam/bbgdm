#' Function for plotting GDM w/ Bayesian Bootstrap.
#' 
#' Plots a Generalised dissimilarity model with bayesian bootstrap.
#' @param gdm.bb_object As derived from gdm.bb function
#' @param preddata new environmental and spatial covariate data to predict GDM.bb model to. 
#' @param type  the type of prediction required. The default is on the scale of the linear predictors; the alternative "response" is on the scale of the response variable. Thus for a default binomial model the default predictions are of log-odds (probabilities on logit scale) and type = "response" gives the predicted probabilities. The "terms" option returns a matrix giving the fitted values of each term in the model formula on the linear predictor scale.   The value of this argument can be abbreviated.
#' @param se.fit logical switch indicating if standard errors are required.
#' @param dispersion the dispersion of the GLM fit to be assumed in computing the standard errors. If omitted, that returned by summary applied to the object is used.
#' @param terms	with type = "terms" by default all terms are returned. A character vector specifies which terms are to be returned
#' @param na.action function determining what should be done with missing values in newdata. The default is to predict NA.
#' @export
#' @examples
#' x <- matrix(rbinom(1:100,1,.6),10,10)# presence absence matrix
#' form <- ~ 1 + covar_1 + covar_2
#' test.gdm.bb <- gdm.bb(form,family="binomial", dism_metric="number_shared", nboot=200, x, sim_covar=TRUE,scale_covar=FALSE)
#' predict(test.gdm.bb,preddata)


predict.gdm.bb <- function(object, newdata = NULL, type = c("link", "response","terms"),
          se.fit = FALSE, dispersion = NULL, terms = NULL, 
          na.action = na.pass, ...) {
  type <- match.arg(type)
  na.act <- object$na.action
  object$na.action <- NULL
  if (!se.fit) {
    if (missing(newdata)) {
      pred <- switch(type, link = object$starting_gdm$linear.predictors, 
                     response = object$starting_gdm$fitted.values, terms = predict.lm(object$starting_gdm, 
                                                                         se.fit = se.fit, scale = 1, type = "terms", 
                                                                         terms = terms))
      if (!is.null(na.act)) 
        pred <- napredict(na.act, pred)
    }
    else {
      pred <- predict.lm(object, newdata, se.fit, scale = 1, 
                         type = ifelse(type == "link", "response", type), 
                         terms = terms, na.action = na.action)
      switch(type, response = {
        pred <- family(object)$linkinv(pred)
      }, link = , terms = )
    }
  }
  else {
    if (inherits(object, "survreg")) 
      dispersion <- 1
    if (is.null(dispersion) || dispersion == 0) 
      dispersion <- summary(object, dispersion = dispersion)$dispersion
    residual.scale <- as.vector(sqrt(dispersion))
    pred <- predict.lm(object, newdata, se.fit, scale = residual.scale, 
                       type = ifelse(type == "link", "response", type), 
                       terms = terms, na.action = na.action)
    fit <- pred$fit
    se.fit <- pred$se.fit
    switch(type, response = {
      se.fit <- se.fit * abs(family(object)$mu.eta(fit))
      fit <- family(object)$linkinv(fit)
    }, link = , terms = )
    if (missing(newdata) && !is.null(na.act)) {
      fit <- napredict(na.act, fit)
      se.fit <- napredict(na.act, se.fit)
    }
    pred <- list(fit = fit, se.fit = se.fit, residual.scale = residual.scale)
  }
  pred
}