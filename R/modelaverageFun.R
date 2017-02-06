#' Uses model avaraging to return fixed-effects parameter estimates.
#' @param model A merMod object from \pkg{lme4}. Most complicated model to evaluate.
#' @param threshold Character. \emph{"C95"}, \emph{"D2"}, \emph{"D4"}, \emph{"D6"}, \emph{"D8"}
#' @param fixed optional, either a single sided formula or a character vector giving names of terms to be included in all models
#' @param r2glmm Logical. Compute R2c and R2m for (g)lmm models using \code{\link[MuMIn]{r.squaredGLMM}}
#' @details Uses the \code{\link[MuMIn]{dredge}} function which takes the full (i.e. most complex, including all terms of interest and interactions) model and automatically
#' builds all combinations of simpler models from the provided terms. For each model, AIC values are computed. An average model is built \code{\link[MuMIn]{model.avg}}
#' using the subset of best models determined by the \emph{threshold} argument. \emph{"C95"} will select the models whose cumulative AICw ≤ 0.95. \emph{"D2"}, \emph{"D4"}, \emph{"D6"}, \emph{"D8"}
#' will select models with ΔAIC ≤ 2, ΔAIC ≤ 4, ΔAIC ≤ 6 and ΔAIC ≤ 8 respectively. If only a single model meets the threshold, it will be returned
#' instead of the averaged model.
#' @return \emph{"$all.models"} Object of class \code{\link[MuMIn]{dredge}}
#' @return \emph{"$best.models"} Object of class \code{\link[MuMIn]{model.avg}}
#' @examples
#' require(lme4)
#' d1<-cbpp
#' d1$period<-as.numeric(as.character(cbpp$period))
#' d1$response<-d1$incidence/d1$size
#' gm1 <- glmer(cbind(incidence, size - incidence) ~ period+size + (1 | herd), data = cbpp, family = binomial)
#' mm1 <- modelaverageFun(model=gm1, threshold="D2")
#' mm1$all.models
#' summary(mm1$best.models)
#' @export

modelaverageFun<-function(model,threshold, fixed=NULL, r2glmm=F){
  if (class(model)[1]=="lmerMod"){
    model<-update(model,na.action=na.fail, REML=T)
  } else {
    model<-update(model,na.action=na.fail)
  }

  mm<-MuMIn::dredge(model,rank="AIC", fixed=fixed, if (r2glmm==T) {extra=MuMIn::r.squaredGLMM} else {extra=NULL})

  if (threshold=="C95"){
    if (length(which(cumsum(mm$weight)<=0.95))>1){
      mm.av<-MuMIn::model.avg(MuMIn::get.models(mm,subset=cumsum(mm$weight)<=0.95), revised.var=T)
    } else {
      mm.av<-MuMIn::get.models(mm,subset=mm$delta==0)[[1]]
    }
  }

  if (threshold=="D2"){
    if (length(which(mm$delta<=2))>1){
      mm.av<-MuMIn::model.avg(MuMIn::get.models(mm,subset=mm$delta<=2), revised.var=T)
    } else {
      mm.av<-MuMIn::get.models(mm,subset=mm$delta==0)[[1]]
    }
  }

  if (threshold=="D4"){
    if (length(which(mm$delta<=4))>1){
      mm.av<-MuMIn::model.avg(MuMIn::get.models(mm,subset=mm$delta<=4), revised.var=T)
    } else {
      mm.av<-MuMIn::get.models(mm,subset=mm$delta==0)[[1]]
    }
  }

  if (threshold=="D6"){
    if (length(which(mm$delta<=6))>1){
      mm.av<-MuMIn::model.avg(MuMIn::get.models(mm,subset=mm$delta<=6), revised.var=T)
    } else {
      mm.av<-MuMIn::get.models(mm,subset=mm$delta==0)[[1]]
    }
  }

  if (threshold=="D8"){
    if (length(which(mm$delta<=8))>1){
      mm.av<-MuMIn::model.avg(MuMIn::get.models(mm,subset=mm$delta<=8), revised.var=T)
    } else {
      mm.av<-MuMIn::get.models(mm,subset=mm$delta==0)[[1]]
    }
  }

  output<-NULL
  output$all.models<-mm
  output$best.models<-mm.av
  return(output)
}
