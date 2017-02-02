#' Simple function for testing for overdispersion.
#'
#' @description Taken from \href{http://glmm.wikidot.com/faq}{Ben Bolker at GLMM wiki}.
#'
#' @param model A merMod object from \pkg{lme4}.
#' @details How can I test for overdispersion/compute an overdispersion factor?
#' ...with the usual caveats, plus a few extras — counting degrees of freedom, etc.
#' — the usual procedure of calculating the sum of squared Pearson residuals and
#' comparing it to the residual degrees of freedom should give at least a crude
#' idea of overdispersion. The following attempt counts each variance or covariance
#' parameter as one model degree of freedom and presents the sum of squared Pearson
#' residuals, the ratio of (SSQ residuals/rdf), the residual df, and the p-value based
#' on the (approximately!!) appropriate \eqn{χ2} distribution. Do PLEASE note the usual, and
#' extra, caveats noted here: this is an APPROXIMATE estimate of an overdispersion parameter.
#'  Even in the GLM case, the expected deviance per point equaling 1 is only true as the
#'  distribution of individual deviates approaches normality, i.e. the usual \eqn{λ>5} rules of
#'  thumb for Poisson values and \eqn{min(Np,N(1−p))>5} for binomial values (e.g. see Venables
#'  and Ripley MASS p. 209). (And that's without the extra complexities due to GLMM, i.e.
#'  the "effective" residual df should be large enough to make the sums of squares converge on a \eqn{χ2} distribution …)
#' @examples
#' library(lme4)  ## 1.0-4
#' set.seed(101)
#' d <- data.frame(y=rpois(1000,lambda=3),x=runif(1000),
#'                 f=factor(sample(1:10,size=1000,replace=TRUE)))
#'                 m1 <- glmer(y~x+(1|f),data=d,family=poisson)
#'                 overdisp_fun(m1)
#' ##        chisq        ratio          rdf            p
#' ## 1026.7780815    1.0298677  997.0000000    0.2497659
#' @export

overdispFun <- function(model) {
  ## number of variance parameters in
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
