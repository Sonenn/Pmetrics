#' Non-compartmental analysis
#'
#' Performs a non-compartmental analysis on individual Bayesian posterior predicted 
#' time-observation profiles generated after an NPAG run by the \code{\link{makePost}}
#' command. 
#'
#' @param post An \emph{PMpost} object created by \code{\link{makePost}} or loaded with \code{\link{NPload}}
#' @param data An \emph{NPAG} object created by \code{\link{NPparse}} or loaded with \code{\link{NPload}}
#' @param input The number of the input (e.g. drug) to analyze; default 1.
#' @param icen Only relevant for PMpost or PMpop objects which have predictions based on median or mean of each
#' subject's Bayesian posterior parameter distribution.  Default is "median", but could be "mean".
#' @param outeq The number of the output equation to analyze; default 1
#' @param block The number of the observation block within subjects, with each block delimited by EVID=4 in the data file; default 1
#' @param start If the \code{start} time is not 0 (default), then it is assumed that steady state (multiple dose) conditions apply. 
#' @param end Set this equal to the dosing interval for steady state (multiple dose) analysis.
#' @return A dataframe of class \emph{PMnca} with columns
#'  \item{id }{Subject identification}
#'  \item{auc }{Area under the time-observation curve, using the trapezoidal approximation, from time 0 until the second dose, 
#' or if only one dose, until the last observation}
#'  \item{aumc }{Area under the first moment curve}
#'  \item{k }{Slope by least-squares linear regression of the final 6 log-transformed observations vs. time}
#'  \item{auclast }{Area under the curve from the time of the last observation to infinity, calculated as [Final obs]/k}
#'  \item{aumclast }{Area under the first moment curve from the time of the last observation to infinity}
#'  \item{aucinf }{Area under the curve from time 0 to infinity, caluculated as auc + auclast}
#'  \item{aumcinf }{Area under the first moment curve from time 0 to infinity}
#'  \item{mrt }{Mean residence time, calculated as 1/k}
#'  \item{cmax }{Maximum predicted concentration after the first dose}
#'  \item{tmax }{Time to cmax}
#'  \item{cl }{Clearance, calculated as dose/aucinf}
#'  \item{vdss }{Volume of distribution at steady state, calculated as cl*mrt}
#'  \item{thalf }{Half life of elimination, calculated as ln(2)/k}
#'  \item{dose }{First dose amount for each subject}
#' @author Michael Neely

makeNCA <- function(post,data,input=1,icen="median",outeq=1,block=1,start=0,end=Inf){

  post2 <- post[post$block==block,]
  post2 <- post2[post2$icen==icen,]
  post2 <- post2[post2$outeq==outeq,]
  post2 <- post2[!is.na(post2$pred),]

  inputCol <- input*2+1
  
  if(data$ncov>0) {
    data$dosecov <- cbind(data$dosecov[,1:2],data$dosecov[,inputCol:(inputCol+1)],data$dosecov[,(inputCol+2):ncol(data$dosecov)])
  } else {data$dosecov <- cbind(data$dosecov[,1:2],data$dosecov[,inputCol:(inputCol+1)]) }
  nsub <- data$nsub  
  
  #convert IV doses
  doseMatrix <- matrix(NA,ncol=ncol(data$dosecov))
  j <- 1
  for(i in 1:(nrow(data$dosecov))){
    j <- j+1
    doseMatrix <- rbind(doseMatrix,data$dosecov[i,])
    if(data$dosecov[i,4]==0){ #iv dose
      doseMatrix[j,4] <- data$dosecov[i,3]*(data$dosecov[i+1,2]-data$dosecov[i,2])
      i <- i+1
      if(i==nrow(data$dosecov)) break  #we've reached the end
      }  

  }
  iv <- which(doseMatrix[,4]==0)
  if(length(iv)>0) doseMatrix <- doseMatrix[-iv,]
  doseMatrix <- doseMatrix[-1,]
  
  ndoses <- table(doseMatrix[,4]!=0,doseMatrix[,1])
  
  if(end == Inf){
    #get time of second dose (if any) for each subject
    end <- vector("numeric",nsub)  
    second <- ndoses
    second[second>2] <- 2
    for (i in 1:nsub){
      index <- sum(ndoses[1:i])-ndoses[i] + second[i]
      end[i] <- doseMatrix[index,2]
    }
    end[end==0] <- Inf
  } else {
    post2 <- post2[post2$time >= start & post2$time <= end,]
    doseMatrix <- doseMatrix[doseMatrix[,2] >= start,]
    end <- rep(end,nsub)
  }
  
  #get dose amounts
  doses <- doseMatrix[!duplicated(doseMatrix[,1]),4]
  

  #make data.frame
  NCA <- matrix(NA,nrow=nsub,ncol=14)

  #cycle through each subject
  for(i in 1:nsub){
    temp <- post2[post2$id==unique(post2$id)[i] & post2$time<=end[i],]
    if(nrow(temp)==0) next
    NCA[i,10] <- max(temp$pred) #cmax
    NCA[i,11] <- temp$time[which(temp$pred==NCA[i,10])][1] - start #tmax
    NCA[i,2] <- makeAUC(temp,outeq=outeq,block=block)[,2] #auc
    temp2 <- data.frame(id=temp$id,time=temp$time,pred=temp$time*temp$pred)
    NCA[i,3] <- makeAUC(temp2,pred~time)[,2] #aumc
  
    temp <- tail(temp,6)
    temp <- temp[temp$pred>0,]
    NCA[i,4] <- -coef(lm(log(temp$pred)~temp$time))[2] #k
    NCA[i,5] <- ifelse(start==0,temp[6,"pred"]/NCA[i,4],NA) #auclast
    NCA[i,6] <- ifelse(start==0,temp[6,"pred"]*temp$time[6]/NCA[i,4] + temp[6,"pred"]/NCA[i,4]^2,NA) #aumclast

  }
  
  NCA <- data.frame(NCA)
  names(NCA) <- c("id","auc","aumc","k","auclast","aumclast","aucinf","aumcinf","mrt","cmax","tmax","cl","vdss","thalf")
  
  start <- rep(start,nsub)
  NCA$id <- unique(post$id)
  NCA$aucinf <- ifelse(start==0,NCA$auc + NCA$auclast,NA)
  NCA$aumcinf <- ifelse(start==0,NCA$aumc + NCA$aumclast,NA)
  NCA$mrt <- ifelse(start==0,NCA$aumcinf/NCA$aucinf,NCA$aumc/NCA$auc)

  NCA$dose <- doses

  NCA$cl <- ifelse(start==0,NCA$dose/NCA$aucinf,NCA$dose/NCA$auc)
  NCA$vdss <- NCA$cl*NCA$mrt
  NCA$thalf <- ifelse(start==0,log(2)/NCA$k,NCA$mrt/1.44)
  class(NCA) <- c("PMnca","data.frame")
  return(NCA)
}


