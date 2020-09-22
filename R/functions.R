#' Beta function
#'
#' Represents a tool able to accommodate various lag structures for the 
#' additional MIDAS variable observed each "low-frequency" period \eqn{t}. It can have a monotonically increasing, 
#' decreasing weighting scheme or a hump-shaped weighting scheme. The Beta function is:
#' \deqn{\delta_k(\omega)=\frac{(k/K)^{\omega_1-1} (1-k/K)^{\omega_2-1}}{\sum_{j=1}^K (j/K)^{\omega_1-1}(1-j/K)^{\omega_2-1}}.}
#' For additional details, see \insertCite{ghysels_2007;textual}{rumidas}.
#' @param k Lag of interest.
#' @param K Number of (lagged) realizations to consider.
#' @param w1,w2 Parameters governing the weights of each \eqn{k} lag.
#' @return The weights associated to each lag \eqn{k}, with \eqn{k=1,\cdots,K}.
#' @importFrom Rdpack reprompt
#' @references
#'  \insertAllCited{} 
#' @examples
#' # suppose to have four lags: 
#' # K<-5 # Note: the number of lags has to be increased by one
#' # w1<-1	# by setting w1=1, only a monotonically decreasing weighting scheme is allowed
#'		#(more recent observations weigh more)
#' # w2<-5
#' beta_function(1:5,K=5,w1=1,w2=5)
#' @export

## Beta function

beta_function<-function(k,K,w1,w2){
j<-1:K
num<-((k/K)^(w1-1))*(1-k/K)^(w2-1)
den<-sum(
((j/K)^(w1-1))*(1-j/K)^(w2-1)
)
beta_func<-num/den
return(beta_func)
}

#' MIDAS variable matrix transformation
#'
#' Implements the trasformation of the MIDAS variable into a matrix, whose dimension is 
#' \eqn{(K+1) \times N}, where \eqn{K} is the number of lagged realizations to consider and
#' \eqn{N} is the length of the variable \eqn{x}. 
#' @param x Variable according to which the MIDAS term has to be aligned. It must be an 'xts' object.
#' @param mv MIDAS variable, observed each period \eqn{t}. It must be an 'xts' object.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param type The frequency of the period of observations for the MIDAS variable. It can be 'weekly', 'monthly' or 'quarterly'.
#' @return The resulting matrix has as many rows as the number of lagged realizations (plus one) of the MIDAS variable
#' to consider, and as many columns as the length of \eqn{x}. 
#' @importFrom Rdpack reprompt
#' @importFrom xts apply.weekly
#' @importFrom xts apply.monthly
#' @importFrom xts apply.quarterly
#' @import roll
#' @import xts
#' @import tseries
#' @import lubridate
#' @import zoo
#' @examples
#' 
#' # weekly frequency
#' # obtain weekly MIDAS variable after daily aggregation
#' # RV_weekly_sum<-apply.weekly(rv5^0.5,sum) #realized volatility
#' # then allocate correctly the information
#' # RV_weekly<-as.xts(coredata(RV_weekly_sum),seq(as.Date("2000-01-10"), 
#' # by = "week", length.out = length(RV_weekly_sum)))
#' # use mv_into_mat (two cases, the second one does not work)
#' # mv_into_mat(sp500['2002/2003-12-26'],diff(RV_weekly['/2003-12']),K=4,"weekly")
#' # mv_into_mat(sp500['2002/2003-12-26'],diff(RV_weekly['/2005-12']),K=4,"weekly") #does not work
#'
#' # monthly frequency
#' # r_t<-sp500['2005/2010']
#' # mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#'
#' # quarterly frequency
#' # RV_quarterly_sum<-apply.quarterly(rv5,sum)
#' # RV_quarterly<-as.xts(coredata(RV_quarterly_sum),seq(as.Date("2000-04-01"), 
#' # by = "quarter", length.out = length(RV_quarterly_sum)))
#' # mv_into_mat(sp500['2004/2010'],diff(RV_quarterly),K=10,"quarterly")
#' 
#' @export

## mv_into_mat_f

mv_into_mat<-function(x,mv,K,type){

N<-length(x)
N_mv<-length(mv)

cond_r_t<- class(x)[1]
cond_mv<- class(mv)[1]

first_obs_r_t<- (first(index(x)))
first_obs_mv<- (first(index(mv)))
diff_time<-as.numeric(difftime(first_obs_r_t,first_obs_mv,units="day"))
num_obs_k<-ifelse(type=="weekly",5,ifelse(type=="monthly",30,120))


w_days_r_t<- wday(index(last(x)))


last_mv<-paste(last(year(mv)),last(month(mv)),"20",sep="-")
last_r_t<-paste(last(year(x)),last(month(x)),"20",sep="-")



############# checks 

if(cond_r_t != "xts") { stop("x must be an xts object. Please provide it in the correct form")}
if(cond_mv != "xts") { stop("mv must be an xts object. Please provide it in the correct form")}


if(diff_time<0|((diff_time - K*num_obs_k)<0)) { stop("mv must start at least
'K+1' periods before x. Please decrease K or provide more observations (in the past) for mv")}


if(last_mv != last_r_t & type=="weekly") { stop("x and mv must end in the same period")}
if(w_days_r_t != 6 & type=="weekly") { stop("the last day of x must be friday")}


mv_m<-matrix(c(rep(NA,(K+1)*N)),ncol=N)


if(type=="weekly"){
for(i in 1:N){
start_p_mv1     <- which(format(index(mv),"%Y-%m-%d")==
format(
floor_date(index(x), "week",week_start = 1)[i],"%Y-%m-%d"))
mv_m[,i]         <- mv[(start_p_mv1-K):(start_p_mv1+1)][2:(K+2)]
}

} else if (type=="monthly") {
for(i in 1:N){
start_p_mv1     <- which(format(index(mv),"%Y-%m")==format(index(x)[i],"%Y-%m"))
mv_m[,i]         <- mv[(start_p_mv1-K):(start_p_mv1)]
}
} else {
for(i in 1:N){
start_p_mv1     <- which(format(index(mv),"%Y-%m-%d")==
format(
floor_date(index(x), "quarter",week_start = 1)[i],"%Y-%m-%d"))
mv_m[,i]         <- mv[(start_p_mv1-K):(start_p_mv1+1)][2:(K+2)]
message(i)
}
}
return(mv_m)
}


#' Information Criteria
#'
#' Returns the Akaike and Bayesian information criteria.
#' @param est The output of the estimation process.
#' @return The resulting vector represents the AIC and BIC criteria.
#' @importFrom Rdpack reprompt
#' @import maxLik
#' @keywords internal

Inf_criteria<-function(est){

inf<-round(c(AIC=(2*length(stats::coef(est))-2*stats::logLik(est)),
BIC=length(stats::coef(est))*log(nrow(est$gradientObs))-2*stats::logLik(est)),6)

}


#' Loss functions
#'
#' Returns the MSE and QLIKE.
#' @param vol_est It is the estimated volatility.
#' @param vol_proxy It is the volatility proxy.
#' @return The resulting vector represents the MSE and QLIKE averages.
#' @importFrom Rdpack reprompt
#' @keywords internal

LF_f<-function(vol_est,vol_proxy){

LF_avg<-c(
"MSE(%)"=100*mean( (vol_est-vol_proxy)^2),
"QLIKE"=mean( log(vol_est) + (vol_proxy/vol_est) )
)

return(LF_avg)

}

#' Standard errors for the Quasi Maximum Likelihood estimator of the GARCH-MIDAS-based models
#'
#' Obtains the standard errors for the Quasi Maximum Likelihood (QML) estimator.
#' @param est It is the output of the maximum likelihood estimation process.
#' @return The resulting vector represents the QML standard errors.
#' @importFrom Rdpack reprompt
#' @import maxLik
#' @keywords internal

QMLE_sd<-function(est){

############### QMLE standard errors (-H)^-1 OPG (-H)^-1

H_hat_inv<- solve(-hessian(est))
OPG<-t(est$gradientObs)%*%est$gradientObs

return(diag(H_hat_inv%*%OPG%*%H_hat_inv)^0.5)
}


#' Standard errors for the Quasi Maximum Likelihood estimator of the GARCH-MIDAS-based models
#'
#' Obtains the standard errors for the Quasi Maximum Likelihood (QML) estimator.
#' @param est It is the output of the maximum likelihood estimation process.
#' @return The resulting vector represents the QML standard errors.
#' @importFrom Rdpack reprompt
#' @import maxLik
#' @keywords internal

MEM_QMLE_sd<-function(est,x,x_pred){

num_param<-ncol(est$gradientObs)

eps<-x/x_pred
var_eps<-stats::var(eps)

sd_coef	<- diag(
matrix(rep(var_eps/2,num_param^2),ncol=num_param) * -solve(hessian(est))
)^0.5 # QMLE standard errors

}

#' Print method for 'rumidas' class
#'
#' @param x An object of class 'rumidas'.
#' @param ... Further arguments passed to or from other methods.

#' @keywords internal
#' @export print.rumidas 
#' @export
print.rumidas <- function(x, ...) {

model<-x[[1]]
mat_coef<-x[[2]]

cat(
cat("\n"),
cat(paste("Model:",model),"\n"),
cat("\n"),
cat(paste("Coefficients: \n",sep="\n")),
cat(rownames(mat_coef), sep="\t", "\n"),
cat(round(mat_coef[,1],4), sep="\t", "\n"),
cat("\n"))
}

#' Summary method for 'rumidas' class
#'
#' @param object An object of class 'rumidas', that is the result of a call to \code{\link{ugmfit}} or \code{\link{umemfit}}.
#' @param ... Additional arguments affecting the summary produced.
#' @examples
#' 
#' # r_t<-sp500['2003/2010']
#' # real<-(rv5['2003/2010'])^0.5		# realized volatility
#' # fit<-umemfit(model="MEM",skew="NO",x=real)
#' # summary.rumidas(fit)
#' @importFrom utils capture.output
#' @export summary.rumidas 
#' @export
summary.rumidas <- function(object, ...) {

model<-object$model
mat_coef<-object$rob_coef_mat
Obs<-object$obs
Period<-object$period
loss<-object$loss_in_s


Period<-paste(substr(Period[1], 1, 10),"/",
substr(Period[2], 1, 10),sep="")

p_value<-mat_coef[,4]

sig<-ifelse(p_value<=0.01,"***",ifelse(p_value>0.01&p_value<=0.05,"**",
ifelse(p_value>0.05&p_value<=0.1,"*"," ")))

mat_coef<-round(mat_coef,4)

mat_coef<-cbind(mat_coef,Sig.=sig)

cat(
cat("\n"),
cat("Coefficients:\n"),
cat(utils::capture.output(mat_coef),  sep = '\n'),
cat("--- \n"),
cat("Signif. codes: 0.01 '***', 0.05 '**', 0.1 '*' \n"),
cat("\n"),
cat("Obs.:", paste(Obs, ".",sep=""), "Sample Period:", Period, "\n"),
cat("MSE(%):", paste(as.numeric(round(loss[1],6)), "; ",sep=""), 
"QLIKE:", as.numeric(round(loss[2],6)), "\n"),
cat("\n"))

}








