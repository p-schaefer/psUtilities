#' Identify highly influential data points.
#'
#' Use \code{\link[influence.ME]{influence.mer}} from the \pkg{influence.ME} package. Influential observations
#' and/or groups are identified as being above the critical threshold of \eqn{4/N}, where \eqn{N} is the number of
#' unique observations or groups.
#'
#' @param model A merMod object from \pkg{lme4}.
#' @param model.variable Character. The name of the random effect variable to be computed.
#' @param data An optional data frame used in \emph{model}
#' @param group.variable Character. Grouping variable for plots
#' @param inf.type Character vector. Specifying which type of influential variables are returned
#' @param plot Logical. Should a plot be produced?
#' @param ... Other arguments passed to \code{\link{ciFun}}
#'
#' @examples
#' require(lme4)
#' m1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' inf1<-influenceFun(model=m1, model.variable= "Days", data=sleepstudy, group.variable="Subject")
#' inf1
#'
#' d1<-cbpp
#' d1$period<-as.numeric(as.character(cbpp$period))
#' d1$response<-d1$incidence/d1$size
#' gm1 <- glmer(cbind(incidence, size - incidence) ~ 1 + (period | herd), data = d1, family = binomial)
#'
#' inf1<-influenceFun(model=gm1, model.variable= "period", data=d1, group.variable="herd")
#' inf1
#'
#' gm2 <- glmer(response ~ 1 + (period | herd), data = d1, weights=d1$size, family = binomial)
#' inf2<-influenceFun(model=gm2, model.variable= "period", data=d1, group.variable="herd")
#' inf2
#' @export


influenceFun<-function (model, model.variable, data, group.variable, inf.type=c("observations","groups"), plot=T,...) {

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

  temp.data$Site<-eval(parse(text=paste0('data$',group.variable)))

  temp.data$Site<-as.factor(as.character(temp.data$Site))

  print("Running... function may take a long time with large models")

  if (any(inf.type=="observations")){
    if (any(names(model@call)=="weights")) {
      mod.obs<-influence.ME::influence(temp.model,obs=T)
    } else {
      mod.obs<-influence.ME::influence(model,obs=T)
    }
  }
  if (any(inf.type=="groups")) {
    if (any(names(model@call)=="weights")) {
      mod.site<-influence.ME::influence(temp.model,group.variable)
    } else {
      mod.site<-influence.ME::influence(model,group.variable)
    }
  }

  temp.data$influence<-"Not Highly Influential"
  temp.data$influence[which(data$Site%in%unique(data$Site)[which(cooks.distance(mod.site)>4/length(unique(data$Site)))])]<-"Highly Influential Site"
  temp.data$influence[which(cooks.distance(mod.obs)>4/nrow(data))]<-"Highly Influential Observation"

  output<-NULL
  output$obs<-mod.obs
  output$site<-mod.site

  ci.mod<-ciFun(model=model,model.variable=model.variable,group.variable=group.variable,data=data,plot=F,...)

  temp.data$Site.sig<-temp.data$Site
  if (any(ci.mod$upper>0&ci.mod$lower>0)){
    levels(temp.data$Site.sig)[which(levels(temp.data$Site.sig)%in%ci.mod$groupID[(ci.mod$upper>0&ci.mod$lower>0)])]<-unique(paste0(unlist(levels(temp.data$Site.sig)[which(levels(temp.data$Site.sig)%in%ci.mod$groupID[(ci.mod$upper>0&ci.mod$lower>0)])]), "*+*"))
  }
  if (any(ci.mod$upper<0&ci.mod$lower<0)){
    levels(temp.data$Site.sig)[which(levels(temp.data$Site.sig)%in%ci.mod$groupID[(ci.mod$upper<0&ci.mod$lower<0)])]<-unique(paste0(unlist(levels(temp.data$Site.sig)[which(levels(temp.data$Site.sig)%in%ci.mod$groupID[(ci.mod$upper<0&ci.mod$lower<0)])]), "*-*"))
  }

  if (plot==T){
    if (!any(names(model@call)=="weights") & as.character(model@call["family"])=="binomial"){
      #p <- recordPlot()
      p<-lattice::xyplot(formula(paste0("response ~", model.variable, "| Site.sig")),
                         data=temp.data,group=influence,auto.key=T,type="p")
      #plot.new()
    } else {
      #p <- recordPlot()
      p<-lattice::xyplot(formula(paste0(colnames(model@frame)[1],"~", model.variable, "| Site.sig")),
                         data=temp.data,group=influence,auto.key=T,type="p")
      #plot.new()
    }
    p
  }

  output<-NULL
  output$plot<-p
  #output$obs.raw<-mod.obs
  #output$site.raw<-mod.site
  if (any(temp.data$influence=="Highly Influential Site")) {
    output$inf.sites<-temp.data[which(data$Site%in%unique(data$Site)[which(cooks.distance(mod.site)>4/length(unique(data$Site)))]),
                                -c(which(colnames(temp.data) %in% c("Site","Site.sig")))]
  }
  if (any(temp.data$influence=="Highly Influential Observation")) {
    output$inf.obs<-temp.data[which(cooks.distance(mod.obs)>4/nrow(data)),
                              -c(which(colnames(temp.data) %in% c("Site","Site.sig")))]
  }

  return(output)
}
