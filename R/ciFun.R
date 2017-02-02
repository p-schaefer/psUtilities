#' Calculates and plots confidence intervals for random effects.
#' @description Calculates and plots confidence intervals for random effects of \pkg{lme4} models using \code{\link[lme4]{REsim}}
#'
#' @param model A merMod object from \pkg{lme4}.
#' @param model.variable Character. The name of the random effect variable to be computed.
#' @param data An optional data frame used in \emph{model}
#' @param col.variable An optional variable name used to colour the plotted values
#' @param plot Logical. Should a plot be produced?
#' @param stat \emph{median} or \emph{mean}
#' @param level Numeric. Confidence interval to return
#' @param oddsRatio Logical. should paramters be converted to odds ratios?
#' @param nsim Numeric.Number of simulations to use
#' @param sig.only Logical. Should only significant random effect levels be plotted?
#' @param ... other arguments passed to REsim.
#' @details Use the Gelman sim technique to build empirical Bayes estimates. Uses the sim function in the arm package
#' @return A data frame and plot of confidence intervals.
#' @examples
#' require(lme4)
#' m1 <- lmer(Reaction ~ 1 + (Days | Subject), sleepstudy)
#' ci1 <- ciFun(model=m1, model.variable= "Days")
#'
#' d1<-cbpp
#' d1$period<-as.numeric(as.character(cbpp$period))
#' d1$response<-d1$incidence/d1$size
#' gm1 <- glmer(response ~ 1 + (period | herd), data = d1, weights=d1$size, family = binomial)
#' ci2 <- ciFun(model=gm1, model.variable= "period")
#' @export


ciFun<-function(model,model.variable,data=NULL,col.variable=NULL,plot=T,stat="median",level=0.95,oddsRatio=F,nsim=1000,sig.only=FALSE,group.variable,...){
  sim<-merTools::REsim(model,n.sim=nsim,...)

  output<-sim[sim$term==model.variable,]
  output$upper<-output[,stat]+output[,"sd"]*qnorm(1-((1-level)/2))
  output$lower<-output[,stat]-output[,"sd"]*qnorm(1-((1-level)/2))
  output$sig<- output[, "lower"] > 0 | output[, "upper"] < 0
  hlineInt<-0
  if (oddsRatio == TRUE) {
    output[, "ymax"] <- exp(output[, "upper"])
    output[, stat] <- exp(output[, stat])
    output[, "ymin"] <- exp(output[, "lower"])
    hlineInt <- 1
  }

  if (plot==T) {
    if (sig.only==T){

      plot.output<-output[order(eval(parse(text=paste0('output$',stat)))),]
      plot.output<-plot.output[plot.output$sig==T,]
      sitenum<-nrow(plot.output)
      if (!is.null(data) & !is.null(col.variable)){
        colours<-data[sapply(as.character(plot.output$groupID),function(x) match(x,as.character(eval(parse(text=paste0('data$',group.variable)))))),col.variable]
      }

      plotrix::plotCI(x=1:sitenum,y=plot.output[, stat],
             li=as.numeric(plot.output$lower),
             ui=as.numeric(plot.output$upper),
             lwd=1,xlab="",ylab="site-specific slope (95% CI)",xaxt="n",las=1,main=paste0(colnames(attr(model,"frame"))[1]),
             col=if (!is.null(data)& !is.null(col.variable)) {col=colours} else {"black"})
      axis(1,at=1:sitenum,labels=plot.output$groupID,cex.axis=0.6,las=2)
      abline(h=hlineInt,lty=2,lwd=1)
      if (!is.null(data) & !is.null(col.variable)){
        legend("topleft",c("Lower","Middle","Upper"),pch=c(21),pt.bg=c(1:3),bty="n",cex=1)
      }

    }
    if (sig.only==F) {
      plot.output<-output[order(eval(parse(text=paste0('output$',stat)))),]
      sitenum<-nrow(plot.output)
      if (!is.null(data) & !is.null(col.variable)){
        colours<-data[sapply(as.character(plot.output$groupID),function(x) match(x,as.character(eval(parse(text=paste0('data$',group.variable)))))),col.variable]
      }

      plotCI(x=1:sitenum,y=plot.output[, stat],
             li=as.numeric(plot.output$lower),
             ui=as.numeric(plot.output$upper),
             lwd=1,xlab="",ylab="site-specific slope (95% CI)",xaxt="n",las=1,main=paste0(colnames(attr(model,"frame"))[1]),
             col=if (!is.null(data) & !is.null(col.variable)) {col=colours} else {"black"})
      axis(1,at=1:sitenum,labels=plot.output$groupID,cex.axis=0.6,las=2)
      abline(h=hlineInt,lty=2,lwd=1)
      if (!is.null(data) & !is.null(col.variable)){
        legend("topleft",c("Lower","Middle","Upper"),pch=c(21),pt.bg=c(1:3),bty="n",cex=1)
      }
    }
  }
  return(output)
}
