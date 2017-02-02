#' A wrapper around \code{\link[lattice]{xyplot}} with showing significance from \code{\link{ciFun}}.
#'
#' @description Grouping variables that do not overlap 0 are showin with *+* and *-*
#'
#' @param model A merMod object from \pkg{lme4}.
#' @param model.variable Character. The name of the random effect variable to be computed.
#' @param data An optional data frame used in \emph{model}
#' @param group.variable Character. Grouping variable for plots
#' @param return Character. Either \emph{"raw"}, \emph{"fitted"}, or \emph{"both"}
#' @param type Character. Passed to \emph{type} in \code{\link[lattice]{xyplot}}
#' @param ... Other arguments passed to \code{\link{ciFun}}
#' @examples
#' require(lme4)
#' m1 <- lmer(Reaction ~ 1 + (Days | Subject), sleepstudy)
#' xyplotFun(model=m1, model.variable= "Days", data=sleepstudy, group.variable="Subject")
#'
#' d1<-cbpp
#' d1$period<-as.numeric(as.character(cbpp$period))
#' d1$response<-d1$incidence/d1$size
#' gm1 <- glmer(response ~ 1 + (period | herd), data = d1, weights=d1$size, family = binomial)
#' xyplotFun(model=gm1, model.variable= "period", data=d1, group.variable="herd")
#' @export

xyplotFun<-function(model,model.variable,data,group.variable,return = "both", type = c("p","r"),...) {

  temp.data<-model.frame(model)

  if (any(names(model@call)=="weights")) {
    temp.data$a<-temp.data[,colnames(model@frame)[1]]*temp.data[,"(weights)"]
    temp.data$b<-temp.data[,"(weights)"]
    temp.data<-temp.data[,colnames(temp.data)!="(weights)"]

    form<-substr(as.character(model@call["formula"]), start=regexpr("~",as.character(model@call["formula"]))[1],stop=nchar(as.character(model@call["formula"])))

    temp.model<<-glmer(paste0("cbind(a,b-a) ", form),data=temp.data, family=binomial)
    temp.data<<-temp.data
  }
  if (!any(names(model@call)=="weights") & as.character(model@call["family"])=="binomial") {
    temp.data$response<-temp.data[,1][,1]/(temp.data[,1][,2]+temp.data[,1][,1])
  }

  temp.data$Site<-eval(parse(text=paste0('temp.data$',group.variable)))

  temp.data$Site<-as.factor(as.character(temp.data$Site))

  ci.mod<-ciFun(model=model,model.variable=model.variable,group.variable=group.variable,data=data,plot=F,...)

  temp.data$Site.sig<-temp.data$Site

  if (any(ci.mod$upper>0&ci.mod$lower>0)){
    levels(temp.data$Site.sig)[which(levels(temp.data$Site.sig)%in%ci.mod$groupID[(ci.mod$upper>0&ci.mod$lower>0)])]<-unique(paste0(unlist(levels(temp.data$Site.sig)[which(levels(temp.data$Site.sig)%in%ci.mod$groupID[(ci.mod$upper>0&ci.mod$lower>0)])]), "*+*"))
  }
  if (any(ci.mod$upper<0&ci.mod$lower<0)){
    levels(temp.data$Site.sig)[which(levels(temp.data$Site.sig)%in%ci.mod$groupID[(ci.mod$upper<0&ci.mod$lower<0)])]<-unique(paste0(unlist(levels(temp.data$Site.sig)[which(levels(temp.data$Site.sig)%in%ci.mod$groupID[(ci.mod$upper<0&ci.mod$lower<0)])]), "*-*"))
  }

  if (return=="both"){
    len<-nrow(temp.data)
    temp.data$Source<-"Original"
    temp.data<-rbind(temp.data,temp.data)
    temp.data[(len+1):(len*2),colnames(model@frame)[1]]<-predict(model,type="response")
    temp.data$Source[(len+1):(len*2)]<-"Fitted"
    temp.data$Source<-as.factor(temp.data$Source)

    if (!any(names(model@call)=="weights") & as.character(model@call["family"])=="binomial"){
      p<-lattice::xyplot(formula(paste0("response ~", model.variable, "| Site.sig")),
                         data=temp.data,group=Source,auto.key=T,type=type)
    } else {
      p<-lattice::xyplot(formula(paste0(colnames(model@frame)[1],"~", model.variable, "| Site.sig")),
                         data=temp.data,group=Source,auto.key=T,type=type)
    }

  }

  if (return=="raw"){
    if (!any(names(model@call)=="weights") & as.character(model@call["family"])=="binomial"){
      p<-lattice::xyplot(formula(paste0("response ~", model.variable, "| Site.sig")),
                      data=temp.data,auto.key=T,type=type)
    } else {
      p<-lattice::xyplot(formula(paste0(colnames(model@frame)[1],"~", model.variable, "| Site.sig")),
                      data=temp.data,auto.key=T,type=type)
    }
  }

  if (return=="fitted"){
    temp.data[,colnames(model@frame)[1]]<-predict(model,type="response")

    if (!any(names(model@call)=="weights") & as.character(model@call["family"])=="binomial"){
      p<-lattice::xyplot(formula(paste0("response ~", model.variable, "| Site.sig")),
                      data=temp.data,auto.key=T,type=type)
    } else {
      p<-lattice::xyplot(formula(paste0(colnames(model@frame)[1],"~", model.variable, "| Site.sig")),
                      data=temp.data,auto.key=T,type=type)
    }
  }
  return(p)
}
