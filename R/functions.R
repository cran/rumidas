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
#' # K<-4 
#' # w1<-1	# by setting w1=1, only a monotonically decreasing weighting scheme is allowed
#'		#(more recent observations weigh more)
#' # w2<-5
#' beta_function(1:4,K=4,w1=1,w2=5)
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

#' Exponential Almon Lag
#'
#' Represents a tool able to accommodate various lag structures for the 
#' additional MIDAS variable observed each "low-frequency" period \eqn{t}. It can have a monotonically increasing, 
#' decreasing weighting scheme or a hump-shaped weighting scheme. As in \insertCite{ghysels_2007;textual}{rumidas},
#' here the function form uses only two parameters:
#' \deqn{\delta_k(\omega_1, \omega_2) = \frac{exp(\omega_{1}k + \omega_2 k^2)}{\sum_{k=1}^K exp(\omega_1 k + \omega_2 k^2)}.}
#' For additional details, see \insertCite{almon_1965;textual}{rumidas} and \insertCite{ghysels_2007;textual}{rumidas}.
#' @param k Lag of interest.
#' @param K Number of (lagged) realizations to consider.
#' @param w1,w2 Parameters governing the weights of each \eqn{k} lag.
#' @return The weights associated to each lag \eqn{k}, with \eqn{k=1,\cdots,K}.
#' @importFrom Rdpack reprompt
#' @references
#'  \insertAllCited{} 
#' @examples
#' # suppose to have four lags: 
#' # K<-4 # Note: the number of lags to consider
#' # w1<-1	
#' # w2<- -0.5 # by setting w2<0, the monotonically decreasing weighting scheme is used
#' exp_almon(1:4,K=4,w1=0.1,w2=-0.5)
#' @export

## Exponential Almon lag function

exp_almon<-function(k,K,w1,w2){
  j<-1:K
  num <- exp(w1*k + w2*k^2)
  den <- sum( exp(w1*j + w2*j^2) )
  exp_func <- num/den
  return(exp_func)
}

#' MIDAS variable matrix transformation
#'
#' Implements the transformation of the MIDAS variable into a matrix, whose dimension is 
#' \eqn{(K+1) \times N}, where \eqn{K} is the number of lagged realizations to consider and
#' \eqn{N} is the length of the variable \eqn{x}. 
#' @param x Variable according to which the MIDAS term has to be aligned. It must be an 'xts' object.
#' @param mv MIDAS variable, observed each period \eqn{t}. It must be an 'xts' object.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param type The frequency of the period of observations for the MIDAS variable. It can be 'weekly', 'monthly', 'quarterly' or 'yearly'.
#' @return The resulting matrix has as many rows as the number of lagged realizations (plus one) of the MIDAS variable
#' to consider, and as many columns as the length of \eqn{x}. 
#' @importFrom Rdpack reprompt
#' @importFrom xts apply.weekly
#' @importFrom xts apply.monthly
#' @importFrom xts apply.quarterly
#' @importFrom xts apply.yearly
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
#' # mv_into_mat(sp500['2002/2003-12-26'],diff(RV_weekly['/2003-12']),K=4,type="weekly")
#' # mv_into_mat(sp500['2002/2003-12-26'],diff(RV_weekly['/2005-12']),K=4,"weekly") #does not work
#'
#' # monthly frequency
#' # r_t<-sp500['2005/2010']
#' # mv_into_mat(r_t,diff(indpro),K=12,type="monthly")
#'
#' # quarterly frequency
#' # RV_quarterly_sum<-apply.quarterly(rv5,sum)
#' # RV_quarterly<-as.xts(coredata(RV_quarterly_sum),seq(as.Date("2000-04-01"), 
#' # by = "quarter", length.out = length(RV_quarterly_sum)))
#' # mv_into_mat(sp500['2004/2010'],diff(RV_quarterly),K=10,type="quarterly")
#' 
#' # yearly frequency
#' # RV_yearly_sum<-apply.yearly(rv5,sum)
#' # RV_yearly<-as.xts(coredata(RV_yearly_sum),seq(as.Date("2001-01-01"), 
#' # by = "year", length.out = length(RV_yearly_sum)))
#' # mv_into_mat(sp500['2006/2010'],diff(RV_yearly),K=2,type="yearly")
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
num_obs_k<-ifelse(type=="weekly",5,
ifelse(type=="monthly",30,
ifelse(type=="quarterly",120,365)))


w_days_r_t<- lubridate::wday(index(last(x)))


last_mv<-paste(xts::last(lubridate::year(mv)),xts::last(lubridate::month(mv)),"20",sep="-")
last_r_t<-paste(xts::last(lubridate::year(x)),xts::last(lubridate::month(x)),"20",sep="-")



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
lubridate::floor_date(index(x), "week",week_start = 1)[i],"%Y-%m-%d"))
mv_m[,i]         <- mv[(start_p_mv1-K):(start_p_mv1+1)][2:(K+2)]
}

} else if (type=="monthly") {
for(i in 1:N){
start_p_mv1     <- which(format(index(mv),"%Y-%m")==format(index(x)[i],"%Y-%m"))
mv_m[,i]         <- mv[(start_p_mv1-K):(start_p_mv1)]
}
} else if (type=="quarterly"){
for(i in 1:N){
start_p_mv1     <- which(format(index(mv),"%Y-%m-%d")==
format(
lubridate::floor_date(index(x), "quarter",week_start = 1)[i],"%Y-%m-%d"))
mv_m[,i]         <- mv[(start_p_mv1-K):(start_p_mv1+1)][2:(K+2)]
message(i)
} 

} else {
for(i in 1:N){
start_p_mv1     <- which(format(index(mv),"%Y-%m-%d")==
format(
lubridate::floor_date(index(x), "year",week_start = 1)[i],"%Y-%m-%d"))
mv_m[,i]         <- c(zoo::coredata(mv[(start_p_mv1-K):(start_p_mv1)][2:(K+1)]),0)
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


#' Standard errors for the Quasi Maximum Likelihood estimator of the MEM-based models
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

options(scipen = 999)

model<-x[[1]]
mat_coef<-x[[2]]

coef_char<-as.character(round(mat_coef[,1],4))

row_names<-gsub("\\s", " ", format(rownames(mat_coef), width=9))
coef_val<-gsub("\\s", " ", format(coef_char,width=9))

cat(
cat("\n"),
cat(paste("Model:",model),"\n"),
cat("\n"),
cat(paste("Coefficients: \n",sep="\n")),
cat(row_names, sep=" ", "\n"),
cat(coef_val, sep=" ", "\n"),
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

#' Summation function for the multi-step-ahead predictions of the GARCH--MIDAS models with the '--X' part.
#'
#' For details, see Eq. (20) of \insertCite{amendola_candila_gallo_2020;textual}{rumidas}.
#' @param COEF The sum of the parameters \eqn{alpha}, \eqn{\beta}, and \eqn{\gamma/2}, if present.
#' @param DELTA The AR coefficient.
#' @param H The length of the multi-step-ahead predictions.
#' @return The vector of the summation for each H.
#' @importFrom Rdpack reprompt
#' @references
#'  \insertAllCited{} 
#' @import maxLik
#' @keywords internal

sum_X_f<-function(COEF,DELTA,H){

res<-list()
for(j in 0:(H-2)){
res[[j+1]]<-COEF^j*DELTA^(H-j-1)
}
return(unlist(res))
}



#' Multi--step--ahead predictions of the GARCH--MIDAS--based models with and without the '--X' part.
#'
#' Calculates the multi--step--ahead predictions for the GARCH--MIDAS and DAGM models, according to the procedure
#' suggested by \insertCite{amendola_candila_gallo_2020;textual}{rumidas}.
#' @param est The estimation object as resulting by the \code{\link{ugmfit}} function
#' @param h The length of the multi-step-ahead predictions
#' @param X **optional**. The '--X' variable. NULL by default. It hat to be equal to the 'X' used in the \code{\link{ugmfit}} function
#' @return The multi-step-ahead predictions, for the following h days, starting from the last day of the 
#' chosen in-sample period adopted in the 'est' object.
#' @importFrom stats coef
#' @importFrom stats arima
#' @references
#'  \insertAllCited{} 
#' @details
#' The multi--step--ahead procedure calculates the volatility predictions keeping fixed the information set at the last
#' observation available and projecting forward the forecasts. The procedure calculates the volatility predictions conditionally
#' to the parameters estimated in the in-sample period. Therefore, the estimation object (through the \code{\link{ugmfit}} function)
#' has to be provided. For additional details, see Eq. (20) in \insertCite{amendola_candila_gallo_2020;textual}{rumidas}.
#' @examples
#' # r_t<-sp500['2008']
#' # X<-(rv5['2008'])^0.5
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly") 
#' # fit<-ugmfit(model="GMX",skew="YES",distribution="norm",r_t,mv_m,K=12,X=X)
#' ### ten days predictions
#' # multi_step_ahead_pred(fit,h=10,X)
#' @importFrom Rdpack reprompt
#' @import maxLik
#' @export

multi_step_ahead_pred<-function(est,h,X=NULL){

### check

if((est$model=="GMX"|est$model=="DAGMX")&(missing(X))) { stop(cat("#Warning:\n If the estimated model includes the 'X' term, then the 'X' variable has to be provided \n"))}

### verify if the object 'est' has the out-of-sample option

out_of_sample<-NULL
out_of_sample<-ifelse(is.null(est$est_vol_oos),0,length(est$est_vol_oos))

vol<-est$est_vol_in_s
tau_t<-est$est_lr_in_s
g_it<-vol/tau_t

g_1_t_plus_1<-xts::last(g_it)
tau_t_plus_1<-xts::last(tau_t)

coef_est<-cbind(est$rob_coef_mat[,1])
rownames(coef_est)<-rownames(est$rob_coef_mat)

alpha_hat<-coef_est[which(rownames(coef_est)=="alpha")]
beta_hat<-coef_est[which(rownames(coef_est)=="beta")]

gamma_hat<-coef_est[which(rownames(coef_est)=="gamma")]
gamma_hat<-ifelse(identical(gamma_hat, numeric(0)),0,gamma_hat)

coef_sum<-alpha_hat+gamma_hat/2+beta_hat

pred<-rep(NA,h)

pred[1]<-(g_1_t_plus_1)*tau_t_plus_1


if (est$model=="GMX"|est$model=="DAGMX"){

z_hat<-coef_est[which(rownames(coef_est)=="z")]

N<-length(vol)+out_of_sample

X_in_s<-X[1:(N-out_of_sample)]

x_last<-xts::last(X_in_s)
x_last<-zoo::coredata(x_last)
z_x_last<-z_hat*x_last
delta<-stats::coef(stats::arima(X_in_s,order=c(1,0,0),include.mean=T))[1]


suppressWarnings(
for(i in 2:h){
pred[i]<-(1+coef_sum^(i-1)*(g_1_t_plus_1-1)+
sum(sum_X_f(coef_sum,delta,i)*z_x_last))*tau_t_plus_1
}
)
} else {

for(i in 2:h){
pred[i]<-(1+coef_sum^(i-1)*(g_1_t_plus_1-1))*tau_t_plus_1
}

}

return(pred)
}



