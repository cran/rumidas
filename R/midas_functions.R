#' GARCH-MIDAS log-likelihood (with skewness)
#'
#' Obtains the log-likelihood of the GARCH-MIDAS, with an asymmetric term linked to past negative returns,
#' according to two errors' conditional distributions: Normal and Student-t. 
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param distribution The conditional density to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # conditional density of the innovations: normal
#' start_val<-c(alpha=0.01,beta=0.8,gamma=0.05,m=0,theta=0.1,w2=2)
#' r_t<-sp500['2005/2010']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' sum(GM_loglik(start_val,r_t,mv_m,K=12,distribution="norm"))
#' 
#' # conditional density of the innovations: Student-t
#' start_val<-c(alpha=0.01,beta=0.8,gamma=0.05,m=0,theta=0.1,w2=2,shape=5)
#' r_t<-sp500['2005/2010']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' sum(GM_loglik(start_val,r_t,mv_m,K=12,distribution="std"))
#' }
#' @keywords internal
#' @export

GM_loglik<-function(param,
daily_ret,
mv_m,
K,
distribution,
lag_fun="Beta"){

if(distribution=="norm"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta		<- param[5]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[6]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)


tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 


####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

ll <- as.numeric(
stats::dnorm(daily_ret, mean(daily_ret), sqrt(g_it*tau_d), log = TRUE)
)

} else {

    alpha           	<- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta		<- param[5]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[6]
	  v			<- param[7]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run


weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)



tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 


####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### variance 

h_it<-zoo::coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+zoo::coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

)

}

return(ll)
##### end
}

#' GARCH-MIDAS conditional volatility (with skewness)
#'
#' Obtains the estimated conditional volatility for the GARCH-MIDAS model, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of estimated values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#'  \insertAllCited{} 
#' @examples
#' \donttest{
#' # estimated volatility
#' est_val<-c(alpha=0.01,beta=0.8,gamma=0.05,m=2,theta=0.1,w2=2)
#' r_t<-sp500['/2010']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' head(GM_cond_vol(est_val,r_t,mv_m,K=12))
#' }
#' @keywords internal
#' @export

GM_cond_vol<-function(param,
daily_ret,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta			<- param[5]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[6]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

cond_vol<- sqrt(g_it*tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

##### end
}

#' GARCH-MIDAS (daily) long-run (with skewness)
#'
#' Obtains the estimated daily long-run volatility for the GARCH-MIDAS model, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of estimated values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#'  \insertAllCited{} 
#' @examples
#' \donttest{
#' # estimated volatility
#' est_val<-c(alpha=0.01,beta=0.8,gamma=0.05,m=2,theta=0.1,w2=2)
#' r_t<-sp500['/2010']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' head(GM_long_run_vol(est_val,r_t,mv_m,K=12))
#' }
#' @keywords internal
#' @export

GM_long_run_vol<-function(param,
daily_ret,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta			  <- param[5]
        w1			  <- ifelse(lag_fun=="Beta",1,0)
        w2			  <- param[6]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)


tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

cond_vol<- sqrt(tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

##### end
}

#' GARCH-MIDAS-X log-likelihood (with skewness)
#'
#' Obtains the log-likelihood of the GARCH-MIDAS-X, with an asymmetric term linked to past negative returns,
#' according to two errors' conditional distributions: Normal and Student-t. 
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param X Additional "X" variable, which must be an "xts" object. Morever, "X" must be observed for the same days of daily_ret.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param distribution The conditional density to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # conditional density of the innovations: normal
#' start_val<-c(alpha=0.01,beta=0.8,gamma=0.05,z=0.1,m=0,theta=0.1,w2=2)
#' r_t<-sp500['2005/2010']
#' X<-rv5['2005/2010']^0.5
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' sum(GM_X_loglik(start_val,r_t,X=X,mv_m,K=12,distribution="norm"))
#' 
#' # conditional density of the innovations: Student-t
#' start_val<-c(alpha=0.01,beta=0.8,gamma=0.05,z=0.1,m=0,theta=0.1,w2=2,shape=5)
#' r_t<-sp500['2005/2010']
#' X<-rv5['2005/2010']^0.5
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' sum(GM_X_loglik(start_val,r_t,X=X,mv_m,K=12,distribution="std"))
#' }
#' @keywords internal
#' @export

GM_X_loglik<-function(param,
daily_ret,
X,
mv_m,
K,
distribution,
lag_fun="Beta"){

if(distribution=="norm"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
	   z				  <- param[4]
        m               <- param[5]
        theta				<- param[6]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[7]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)


tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 


####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d + z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

ll <- as.numeric(
stats::dnorm(daily_ret, mean(daily_ret), sqrt(g_it*tau_d), log = TRUE)
)

} else {

    alpha           	<- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
		z				<- param[4]
        m               <- param[5]
        theta		<- param[6]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[7]
	  v			<- param[8]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run


weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)



tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d+ z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### variance 

h_it<-zoo::coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+zoo::coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

)

}

return(ll)
##### end
}

#' GARCH-MIDAS-X conditional volatility (with skewness)
#'
#' Obtains the estimated conditional volatility for the GARCH-MIDAS-X model, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of estimated values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param X Additional "X" variable, which must be an "xts" object. Morever, "X" must be observed for the same days of daily_ret.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#'  \insertAllCited{} 
#' @examples
#' \donttest{
#' # estimated volatility
#' est_val<-c(alpha=0.01,beta=0.8,gamma=0.05,z=0.1,m=2,theta=0.1,w2=2)
#' r_t<-sp500['/2010']
#' X<-rv5['/2010']^0.5
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' head(GM_X_cond_vol(est_val,r_t,X,mv_m,K=12))
#' }
#' @keywords internal
#' @export

GM_X_cond_vol<-function(param,
daily_ret,
X,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
		z				<- param[4]
        m               <- param[5]
        theta			<- param[6]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[7]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run


weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)


tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d + z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

cond_vol<- sqrt(g_it*tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

##### end
}

#' GARCH-MIDAS-X (daily) long-run (with skewness)
#'
#' Obtains the estimated daily long-run volatility for the GARCH-MIDAS-X model, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of estimated values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param X Additional "X" variable, which must be an "xts" object. Morever, "X" must be observed for the same days of daily_ret.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#'  \insertAllCited{} 
#' @examples
#' \donttest{
#' # estimated volatility
#' est_val<-c(alpha=0.01,beta=0.8,gamma=0.05,z=0.1,m=2,theta=0.1,w2=2)
#' r_t<-sp500['/2010']
#' X<-rv5['/2010']^0.5
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' head(GM_X_long_run_vol(est_val,r_t,X,mv_m,K=12))
#' }
#' @keywords internal
#' @export

GM_X_long_run_vol<-function(param,
daily_ret,
X,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        z			  <- param[4]
        m               <- param[5]
        theta			  <- param[6]
        w1			  <- ifelse(lag_fun=="Beta",1,0)
        w2			  <- param[7]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run


weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)


tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d + z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

cond_vol<- sqrt(tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

##### end
}


#' GARCH-MIDAS log-likelihood (no skewness)
#'
#' Obtains the log-likelihood of the GARCH-MIDAS, according to two errors' conditional distributions: Normal and Student-t.
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of starting values.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param distribution The conditional density to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # conditional density of the innovations: normal
#' start_val<-c(alpha=0.01,beta=0.8,m=0,theta=0.1,w2=2)
#' r_t<-sp500['2005/2010']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' sum(GM_loglik_no_skew(start_val,r_t,mv_m,K=12,distribution="norm"))
#' 
#' # conditional density of the innovations: Student-t
#' start_val<-c(alpha=0.01,beta=0.8,m=0,theta=0.1,w2=2,shape=5)
#' r_t<-sp500['2005/2010']
#' mv_m<-mv_into_mat(r_t,indpro,K=12,"monthly")
#' sum(GM_loglik_no_skew(start_val,r_t,mv_m,K=12,distribution="std"))
#' }
#' @keywords internal
#' @export

GM_loglik_no_skew<-function(param,
daily_ret,
mv_m,
K,
distribution,
lag_fun="Beta"){

if(distribution=="norm"){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta		<- param[4]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[5]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        ll                                      <- 0 

##### daily long-run


weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)


tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 


####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

ll <- as.numeric(
stats::dnorm(daily_ret, mean(daily_ret), sqrt(g_it*tau_d), log = TRUE)
)

} else {

    alpha           	<- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta		<- param[4]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[5]
	  v			<- param[6]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        ll                                      <- 0 

##### daily long-run

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 


####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### variance 

h_it<-zoo::coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+zoo::coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

)

}

return(ll)
##### end
}

#' GARCH-MIDAS conditional volatility (without skewness)
#'
#' Obtains the estimated conditional volatility for the GARCH-MIDAS model, without the skewness parameter in the short--run.
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of estimated values.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' est_val<-c(alpha=0.01,beta=0.8,m=2,theta=0.1,w2=2)
#' r_t<-sp500['/2010']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' head(GM_cond_vol_no_skew(est_val,r_t,mv_m,K=12))
#' }
#' @keywords internal
#' @export

GM_cond_vol_no_skew<-function(param,
daily_ret,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta		<- param[4]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[5]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        ll                                      <- 0 

##### daily long-run


weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)


tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 


####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

cond_vol<- sqrt(g_it*tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))
##### end
}

#' GARCH-MIDAS (daily) long-run volatility (without skewness)
#'
#' Obtains the daily long-run volatility for the GARCH-MIDAS model, without the skewness parameter in the short-run..
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of estimated values.  
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' est_val<-c(alpha=0.01,beta=0.8,m=2,theta=0.1,w2=2)
#' r_t<-sp500['/2010']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' head(GM_long_run_vol_no_skew(est_val,r_t,mv_m,K=12))
#' }
#' @keywords internal
#' @export

GM_long_run_vol_no_skew<-function(param,
daily_ret,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta		<- param[4]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[5]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        ll                                      <- 0 

##### daily long-run



weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)


tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

cond_vol<- sqrt(tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))
##### end
}

#' GARCH-MIDAS-X log-likelihood (no skewness)
#'
#' Obtains the log-likelihood of the GARCH-MIDAS-X, according to two errors' conditional distributions: Normal and Student-t.
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param X Additional "X" variable, which must be an "xts" object. Morever, "X" must be observed for the same days of daily_ret.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param distribution The conditional density to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # conditional density of the innovations: normal
#' start_val<-c(alpha=0.01,beta=0.8,z=0.1,m=0,theta=0.1,w2=2)
#' r_t<-sp500['2005/2010']
#' X<-rv5['2005/2010']^0.5
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' sum(GM_X_loglik_no_skew(start_val,r_t,X,mv_m,K=12,distribution="norm"))
#' 
#' # conditional density of the innovations: Student-t
#' start_val<-c(alpha=0.01,beta=0.8,z=0.1,m=0,theta=0.1,w2=2,shape=5)
#' r_t<-sp500['2005/2010']
#' X<-rv5['2005/2010']^0.5
#' mv_m<-mv_into_mat(r_t,indpro,K=12,"monthly")
#' sum(GM_X_loglik_no_skew(start_val,r_t,X,mv_m,K=12,distribution="std"))
#' }
#' @keywords internal
#' @export

GM_X_loglik_no_skew<-function(param,
daily_ret,
X,
mv_m,
K,
distribution,
lag_fun="Beta"){

if(distribution=="norm"){

        alpha           <- param[1]
        beta            <- param[2]
	   z				   <- param[3]
        m               <- param[4]
        theta		<- param[5]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[6]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        ll                                      <- 0 

##### daily long-run


weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)


tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 


####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d+z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

ll <- as.numeric(
stats::dnorm(daily_ret, mean(daily_ret), sqrt(g_it*tau_d), log = TRUE)
)

} else {

    alpha           	<- param[1]
        beta            <- param[2]
		z				<- param[3]
        m               <- param[4]
        theta		<- param[5]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[6]
	  v			<- param[7]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        ll                                      <- 0 

##### daily long-run

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 


####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d+z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### variance 

h_it<-zoo::coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+zoo::coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

)

}

return(ll)
##### end
}

#' GARCH-MIDAS-X conditional volatility (without skewness)
#'
#' Obtains the estimated conditional volatility for the GARCH-MIDAS-X model, without the skewness parameter in the short--run.
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of estimated values.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param X Additional "X" variable, which must be an "xts" object. Morever, "X" must be observed for the same days of daily_ret.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' est_val<-c(alpha=0.01,beta=0.8,z=0.1,m=2,theta=0.1,w2=2)
#' r_t<-sp500['/2010']
#' X<-rv5['/2010']^0.5
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' head(GM_X_cond_vol_no_skew(est_val,r_t,X,mv_m,K=12))
#' }
#' @keywords internal
#' @export

GM_X_cond_vol_no_skew<-function(param,
daily_ret,
X,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
	    z				<- param[3]
        m               <- param[4]
        theta		<- param[5]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[6]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        ll                                      <- 0 

##### daily long-run


weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)


tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 


####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d+z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

cond_vol<- sqrt(g_it*tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))
##### end
}

#' GARCH-MIDAS-X (daily) long-run volatility (without skewness)
#'
#' Obtains the daily long-run volatility for the GARCH-MIDAS-X model, without the skewness parameter in the short--run.
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of estimated values.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param X Additional "X" variable, which must be an "xts" object. Morever, "X" must be observed for the same days of daily_ret.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' est_val<-c(alpha=0.01,beta=0.8,z=0.1,m=2,theta=0.1,w2=2)
#' r_t<-sp500['/2010']
#' X<-rv5['/2010']^0.5
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' head(GM_X_long_run_vol_no_skew(est_val,r_t,X,mv_m,K=12))
#' }
#' @keywords internal
#' @export

GM_X_long_run_vol_no_skew<-function(param,
daily_ret,
X,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
		z			  <- param[3]
        m               <- param[4]
        theta		<- param[5]
        w1			<- ifelse(lag_fun=="Beta",1,0)
        w2			<- param[6]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        ll                                      <- 0 

##### daily long-run



weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)


tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d+z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

cond_vol<- sqrt(tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))
##### end
}

#' DAGM log-likelihood (with skewness)
#'
#' Obtains the log-likelihood of the DAGM, with an asymmetric term linked to past negative returns,
#' according to two errors' conditional distributions: Normal and Student-t. 
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param distribution The conditional density to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # conditional density of the innovations: normal
#' start_val<-c(0.01,0.80,0.05,0,0,1.1,0,1.1)
#' r_t<-sp500['2005/2010']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' sum(DAGM_loglik(start_val,r_t,mv_m,K=12,distribution="norm"))
#' }
#' @keywords internal
#' @export

DAGM_loglik<-function(param,
daily_ret,
mv_m,
K,
distribution,
lag_fun="Beta"){

if(distribution=="norm"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta_pos       <- param[5]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[6]
        theta_neg       <- param[7]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[8]
	   
        TT              <- length(daily_ret)

	   MV_pos		<- ifelse(mv_m>=0,mv_m,0)
         MV_neg		<- ifelse(mv_m<0,mv_m,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg   <- ifelse(zoo::coredata(daily_ret) <0,1,0)

        ll              <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### loglik

ll <- as.numeric(
stats::dnorm(daily_ret, mean(daily_ret), sqrt(g_it*tau_d), log = TRUE)
)


} else {

  	 alpha            <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta_pos       <- param[5]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[6]
        theta_neg       <- param[7]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[8]
	  v	      	<- param[9]
	   
        TT              <- length(daily_ret)

	   MV_pos		<- ifelse(mv_m>=0,mv_m,0)
         MV_neg		<- ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) <0,1,0)

        ll                                      <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}


###### variance 

h_it<-zoo::coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+zoo::coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

)

}

return(ll)
##### end
}

#' DAGM-2M log-likelihood (with skewness)
#'
#' Obtains the log-likelihood of the DAGM with two MIDAS variables, with an asymmetric term linked to past negative returns,
#' according to two errors' conditional distributions: Normal and Student-t. 
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m_1 first MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param mv_m_2 second MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K_1 Number of (lagged) realizations of the first MIDAS variable to consider.
#' @param K_2 Number of (lagged) realizations of the second MIDAS variable to consider.
#' @param distribution The conditional density to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # conditional density of the innovations: normal
#' start_val<-c(0.01,0.80,0.05,0.2,0.1,1.1,0.4,1.1,0.5,1.1,0,1.1)
#' r_t<-sp500['2005/2010']
#' mv_m_1<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' mv_m_2<-mv_into_mat(r_t,diff(indpro),K=24,"monthly")
#' sum(DAGM_2M_loglik(start_val,r_t,mv_m_1,mv_m_2,K_1=12,K_2=24,distribution="norm"))
#' }
#' @keywords internal
#' @export

DAGM_2M_loglik<-function(param,
daily_ret,
mv_m_1,
mv_m_2,
K_1,
K_2,
distribution,
lag_fun="Beta"){

if(distribution=="norm"){

        alpha           	<- param[1]
        beta            	<- param[2]
        gamma_1         	<- param[3]
        m               	<- param[4]
        theta_pos_1      	<- param[5]
        w1_pos_1          	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_1          	<- param[6]
        theta_neg_1       	<- param[7]
        w1_neg_1          	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_1          	<- param[8]
 		theta_pos_2      	<- param[9]
        w1_pos_2          	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_2          	<- param[10]
        theta_neg_2       <- param[11]
        w1_neg_2         	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_2          	<- param[12]
	   
        TT              <- length(daily_ret)

	   MV_pos_1		<- ifelse(mv_m_1>=0,mv_m_1,0)
       MV_neg_1	<- ifelse(mv_m_1<0,mv_m_1,0)

       MV_pos_2		<- ifelse(mv_m_2>=0,mv_m_2,0)
       MV_neg_2	<- ifelse(mv_m_2<0,mv_m_2,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg   <- ifelse(zoo::coredata(daily_ret) <0,1,0)

        ll              <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_pos_1,w2_pos_1))[2:(K_1+1)],0)
betas_neg_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_neg_1,w2_neg_1))[2:(K_1+1)],0)

betas_pos_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_pos_2,w2_pos_2))[2:(K_2+1)],0)
betas_neg_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_neg_2,w2_neg_2))[2:(K_2+1)],0)


tau_d                     <- exp(
m+
theta_pos_1*suppressWarnings(roll_sum(MV_pos_1, c(K_1+1),weights = betas_pos_1))[K_1+1,] + 
theta_neg_1*suppressWarnings(roll_sum(MV_neg_1, c(K_1+1),weights = betas_neg_1))[K_1+1,] + 
theta_pos_2*suppressWarnings(roll_sum(MV_pos_2, c(K_2+1),weights = betas_pos_2))[K_2+1,] + 
theta_neg_2*suppressWarnings(roll_sum(MV_neg_2, c(K_2+1),weights = betas_neg_2))[K_2+1,] 
)


####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### loglik

ll <- as.numeric(
stats::dnorm(daily_ret, mean(daily_ret), sqrt(g_it*tau_d), log = TRUE)
)


} else {

  	 alpha            	<- param[1]
        beta            	<- param[2]
        gamma_1         	<- param[3]
        m               	<- param[4]
        theta_pos_1      	<- param[5]
        w1_pos_1          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos_1          <- param[6]
        theta_neg_1       <- param[7]
        w1_neg_1          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg_1          <- param[8]
 		theta_pos_2      	<- param[9]
        w1_pos_2         	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_2          	<- param[10]
        theta_neg_2       	<- param[11]
        w1_neg_2         	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_2          	<- param[12]
	  	v	      		<- param[13]
	   
        TT              <- length(daily_ret)

	   MV_pos_1		<- ifelse(mv_m_1>=0,mv_m_1,0)
       MV_neg_1		<- ifelse(mv_m_1<0,mv_m_1,0)

       MV_pos_2		<- ifelse(mv_m_2>=0,mv_m_2,0)
       MV_neg_2		<- ifelse(mv_m_2<0,mv_m_2,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) <0,1,0)

        ll                                      <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)


betas_pos_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_pos_1,w2_pos_1))[2:(K_1+1)],0)
betas_neg_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_neg_1,w2_neg_1))[2:(K_1+1)],0)

betas_pos_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_pos_2,w2_pos_2))[2:(K_2+1)],0)
betas_neg_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_neg_2,w2_neg_2))[2:(K_2+1)],0)



tau_d                     <- exp(
m+
theta_pos_1*suppressWarnings(roll_sum(MV_pos_1, c(K_1+1),weights = betas_pos_1))[K_1+1,] + 
theta_neg_1*suppressWarnings(roll_sum(MV_neg_1, c(K_1+1),weights = betas_neg_1))[K_1+1,] + 
theta_pos_2*suppressWarnings(roll_sum(MV_pos_2, c(K_2+1),weights = betas_pos_2))[K_2+1,] + 
theta_neg_2*suppressWarnings(roll_sum(MV_neg_2, c(K_2+1),weights = betas_neg_2))[K_2+1,] 
)

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}


###### variance 

h_it<-zoo::coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+zoo::coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

)

}

return(ll)
##### end
}



#' DAGM log-likelihood (without skewness)
#'
#' Obtains the log-likelihood of the DAGM, without the asymmetric term linked to past negative returns,
#' according to two errors' conditional distributions: Normal and Student-t. 
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of starting values.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param distribution The conditional density to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # conditional density of the innovations: normal
#' start_val<-c(alpha=0.01,beta=0.80,gamma_1=0.05,m=0,theta_pos=0,w2_pos=1.1,theta_neg=0,w2_neg=1.1)
#' r_t<-sp500['2005/2010']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' sum(DAGM_loglik(start_val,r_t,mv_m,K=12,distribution="norm"))
#' 
#' # conditional density of the innovations: Student-t
#' start_val<-c(0.01,0.80,0.05,0,0,1.1,0,1.1,5)
#' r_t<-sp500['2005/2010']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' sum(DAGM_loglik(start_val,r_t,mv_m,K=12,distribution="std"))
#' }
#' @keywords internal
#' @export

DAGM_loglik_no_skew<-function(param,
daily_ret,
mv_m,
K,
distribution,
lag_fun="Beta"){

if(distribution=="norm"){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta_pos       <- param[4]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[5]
        theta_neg       <- param[6]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[7]
	   
        TT              <- length(daily_ret)

		MV_pos<-ifelse(mv_m>=0,mv_m,0)
         MV_neg<-ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        ll                                      <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### loglik

ll <- as.numeric(
stats::dnorm(daily_ret, mean(daily_ret), sqrt(g_it*tau_d), log = TRUE)
)


} else {

  	 alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta_pos       <- param[4]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[5]
        theta_neg       <- param[6]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[7]
		v			  <- param[8]
	   
        TT              <- length(daily_ret)

		MV_pos<-ifelse(mv_m>=0,mv_m,0)
         MV_neg<-ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       

        ll                                      <- 0 


###### long-run 


weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}


###### variance 

h_it<-zoo::coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+zoo::coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

)

}

return(ll)
##### end
}


#' DAGM-2M log-likelihood (without skewness)
#'
#' Obtains the log-likelihood of the DAGM with two MIDAS variables,
#' according to two errors' conditional distributions: Normal and Student-t. 
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m_1 first MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param mv_m_2 second MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K_1 Number of (lagged) realizations of the first MIDAS variable to consider.
#' @param K_2 Number of (lagged) realizations of the second MIDAS variable to consider.
#' @param distribution The conditional density to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # conditional density of the innovations: normal
#' start_val<-c(0.01,0.80,0.2,0.1,1.1,0.4,1.1,0.5,1.1,0,1.1)
#' r_t<-sp500['2005/2010']
#' mv_m_1<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' mv_m_2<-mv_into_mat(r_t,diff(indpro),K=24,"monthly")
#' sum(DAGM_2M_loglik_no_skew(start_val,r_t,mv_m_1,mv_m_2,K_1=12,K_2=24,distribution="norm"))
#' }
#' @keywords internal
#' @export

DAGM_2M_loglik_no_skew<-function(param,
daily_ret,
mv_m_1,
mv_m_2,
K_1,
K_2,
distribution,
lag_fun="Beta"){

if(distribution=="norm"){

        alpha           	<- param[1]
        beta            	<- param[2]
        m               	<- param[3]
        theta_pos_1      	<- param[4]
        w1_pos_1          	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_1          	<- param[5]
        theta_neg_1       	<- param[6]
        w1_neg_1          	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_1          	<- param[7]
 		theta_pos_2      	<- param[8]
        w1_pos_2          	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_2          	<- param[9]
        theta_neg_2       <- param[10]
        w1_neg_2         	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_2          	<- param[11]
	   
        TT              <- length(daily_ret)

	   MV_pos_1		<- ifelse(mv_m_1>=0,mv_m_1,0)
       MV_neg_1	<- ifelse(mv_m_1<0,mv_m_1,0)

       MV_pos_2		<- ifelse(mv_m_2>=0,mv_m_2,0)
       MV_neg_2	<- ifelse(mv_m_2<0,mv_m_2,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       
        ll              <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_pos_1,w2_pos_1))[2:(K_1+1)],0)
betas_neg_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_neg_1,w2_neg_1))[2:(K_1+1)],0)

betas_pos_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_pos_2,w2_pos_2))[2:(K_2+1)],0)
betas_neg_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_neg_2,w2_neg_2))[2:(K_2+1)],0)


tau_d                     <- exp(
m+
theta_pos_1*suppressWarnings(roll_sum(MV_pos_1, c(K_1+1),weights = betas_pos_1))[K_1+1,] + 
theta_neg_1*suppressWarnings(roll_sum(MV_neg_1, c(K_1+1),weights = betas_neg_1))[K_1+1,] + 
theta_pos_2*suppressWarnings(roll_sum(MV_pos_2, c(K_2+1),weights = betas_pos_2))[K_2+1,] + 
theta_neg_2*suppressWarnings(roll_sum(MV_neg_2, c(K_2+1),weights = betas_neg_2))[K_2+1,] 
)


####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### loglik

ll <- as.numeric(
stats::dnorm(daily_ret, mean(daily_ret), sqrt(g_it*tau_d), log = TRUE)
)


} else {

  	 alpha            	<- param[1]
        beta            	<- param[2]
        m               	<- param[3]
        theta_pos_1      	<- param[4]
        w1_pos_1          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos_1          <- param[5]
        theta_neg_1       <- param[6]
        w1_neg_1          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg_1          <- param[7]
 		theta_pos_2      	<- param[8]
        w1_pos_2         	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_2          	<- param[9]
        theta_neg_2       	<- param[10]
        w1_neg_2         	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_2          	<- param[11]
	  	v	      		<- param[12]
	   
        TT              <- length(daily_ret)

	   MV_pos_1		<- ifelse(mv_m_1>=0,mv_m_1,0)
       MV_neg_1		<- ifelse(mv_m_1<0,mv_m_1,0)

       MV_pos_2		<- ifelse(mv_m_2>=0,mv_m_2,0)
       MV_neg_2		<- ifelse(mv_m_2<0,mv_m_2,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        ll                                      <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)


betas_pos_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_pos_1,w2_pos_1))[2:(K_1+1)],0)
betas_neg_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_neg_1,w2_neg_1))[2:(K_1+1)],0)

betas_pos_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_pos_2,w2_pos_2))[2:(K_2+1)],0)
betas_neg_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_neg_2,w2_neg_2))[2:(K_2+1)],0)



tau_d                     <- exp(
m+
theta_pos_1*suppressWarnings(roll_sum(MV_pos_1, c(K_1+1),weights = betas_pos_1))[K_1+1,] + 
theta_neg_1*suppressWarnings(roll_sum(MV_neg_1, c(K_1+1),weights = betas_neg_1))[K_1+1,] + 
theta_pos_2*suppressWarnings(roll_sum(MV_pos_2, c(K_2+1),weights = betas_pos_2))[K_2+1,] + 
theta_neg_2*suppressWarnings(roll_sum(MV_neg_2, c(K_2+1),weights = betas_neg_2))[K_2+1,] 
)

####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}


###### variance 

h_it<-zoo::coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+zoo::coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

)

}

return(ll)
##### end
}

#' DAGM conditional volatility (with skewness)
#'
#' Obtains the conditional volatility for the DAGM, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of estimated values. It must be a eight- or nine- dimensional vector. See the examples below.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # estimated volatility 
#' est_val<-c(0.01,0.80,0.05,0,0.1,1.1,-0.3,1.1)
#' r_t<-sp500['2005/2010']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' vol<-DAGM_cond_vol(est_val,r_t,mv_m,K=12)
#' head(vol)
#' }
#' @keywords internal
#' @export

DAGM_cond_vol<-function(param,
daily_ret,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta_pos       <- param[5]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[6]
        theta_neg       <- param[7]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[8]
	   
        TT              <- length(daily_ret)

	   MV_pos		<- ifelse(mv_m>=0,mv_m,0)
         MV_neg		<- ifelse(mv_m<0,mv_m,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg   <- ifelse(zoo::coredata(daily_ret) <0,1,0)

        ll                                      <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### cond vol.

cond_vol<- sqrt(g_it*tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

return(cond_vol)

##### end
}

#' DAGM (daily) long-run volatility (with skewness)
#'
#' Obtains the daily long-run volatility for the DAGM, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of estimated values. It must be a eight- or nine- dimensional vector. See the examples below.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' est_val<-c(0.01,0.80,0.05,0,0.1,1.1,-0.3,1.1)
#' r_t<-sp500['/2010']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' head(DAGM_long_run_vol(est_val,r_t,mv_m,K=12))
#' }
#' @keywords internal
#' @export

DAGM_long_run_vol<-function(param,
daily_ret,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta_pos       <- param[5]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[6]
        theta_neg       <- param[7]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[8]
	   
        TT              <- length(daily_ret)

	   MV_pos		<- ifelse(mv_m>=0,mv_m,0)
         MV_neg		<- ifelse(mv_m<0,mv_m,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg   <- ifelse(zoo::coredata(daily_ret) <0,1,0)

        ll              <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)

tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### cond vol.

cond_vol<- sqrt(tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

return(cond_vol)

##### end
}


#' DAGM conditional volatility (no skewness)
#'
#' Obtains the conditional volatility for the DAGM.
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of estimated values. It must be a seven- or eight- dimensional vector. See the examples below.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' # est_val<-c(0.01,0.80,0,0.1,1.1,-0.3,1.1)
#' # r_t<-sp500['/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # head(DAGM_cond_vol_no_skew(est_val,r_t,mv_m,K=12))
#' @keywords internal
#' @export

DAGM_cond_vol_no_skew<-function(param,
daily_ret,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta_pos       <- param[4]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[5]
        theta_neg       <- param[6]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[7]
	   
        TT              <- length(daily_ret)

		MV_pos<-ifelse(mv_m>=0,mv_m,0)
         MV_neg<-ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        ll                                      <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### cond vol.

cond_vol<- sqrt(g_it*tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

return(cond_vol)

##### end
}

#' DAGM (daily) long-run volatility (no skewness)
#'
#' Obtains the daily long-run volatility for the DAGM.
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of estimated values. It must be a seven- or eight- dimensional vector. See the examples below.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the long-run volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' # est_val<-c(0.01,0.80,0,0.1,1.1,-0.3,1.1)
#' # r_t<-sp500['/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # head(DAGM_long_run_vol_no_skew(est_val,r_t,mv_m,K=12))
#' @keywords internal
#' @export

DAGM_long_run_vol_no_skew<-function(param,
daily_ret,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta_pos       <- param[4]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[5]
        theta_neg       <- param[6]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[7]
	   
        TT              <- length(daily_ret)

		MV_pos<-ifelse(mv_m>=0,mv_m,0)
         MV_neg<-ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        ll                                      <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### cond vol.

cond_vol<- sqrt(tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

return(cond_vol)
##### end
}


#' DAGM-2M conditional volatility (with skewness)
#'
#' Obtains the conditional volatility of the DAGM with two MIDAS variables. 
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m_1 first MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param mv_m_2 second MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K_1 Number of (lagged) realizations of the first MIDAS variable to consider.
#' @param K_2 Number of (lagged) realizations of the second MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' start_val<-c(0.01,0.80,0.05,0.2,0.1,1.1,0.4,1.1,0.5,1.1,0,1.1)
#' r_t<-sp500['2005/2010']
#' mv_m_1<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' mv_m_2<-mv_into_mat(r_t,diff(indpro),K=24,"monthly")
#' head(DAGM_2M_cond_vol(start_val,r_t,mv_m_1,mv_m_2,K_1=12,K_2=24))
#' }
#' @keywords internal
#' @export

DAGM_2M_cond_vol<-function(param,
daily_ret,
mv_m_1,
mv_m_2,
K_1,
K_2,
lag_fun="Beta"){

        alpha           	<- param[1]
        beta            	<- param[2]
        gamma_1         	<- param[3]
        m               	<- param[4]
        theta_pos_1      	<- param[5]
        w1_pos_1          	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_1          	<- param[6]
        theta_neg_1       	<- param[7]
        w1_neg_1          	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_1          	<- param[8]
 		theta_pos_2      	<- param[9]
        w1_pos_2          	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_2          	<- param[10]
        theta_neg_2       <- param[11]
        w1_neg_2         	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_2          	<- param[12]
	   
        TT              <- length(daily_ret)

	   MV_pos_1		<- ifelse(mv_m_1>=0,mv_m_1,0)
       MV_neg_1	<- ifelse(mv_m_1<0,mv_m_1,0)

       MV_pos_2		<- ifelse(mv_m_2>=0,mv_m_2,0)
       MV_neg_2	<- ifelse(mv_m_2<0,mv_m_2,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg   <- ifelse(zoo::coredata(daily_ret) <0,1,0)

       
###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_pos_1,w2_pos_1))[2:(K_1+1)],0)
betas_neg_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_neg_1,w2_neg_1))[2:(K_1+1)],0)

betas_pos_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_pos_2,w2_pos_2))[2:(K_2+1)],0)
betas_neg_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_neg_2,w2_neg_2))[2:(K_2+1)],0)


tau_d                     <- exp(
m+
theta_pos_1*suppressWarnings(roll_sum(MV_pos_1, c(K_1+1),weights = betas_pos_1))[K_1+1,] + 
theta_neg_1*suppressWarnings(roll_sum(MV_neg_1, c(K_1+1),weights = betas_neg_1))[K_1+1,] + 
theta_pos_2*suppressWarnings(roll_sum(MV_pos_2, c(K_2+1),weights = betas_pos_2))[K_2+1,] + 
theta_neg_2*suppressWarnings(roll_sum(MV_neg_2, c(K_2+1),weights = betas_neg_2))[K_2+1,] 
)


####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### cond vol.

cond_vol<- sqrt(g_it*tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

return(cond_vol)

##### end
}


#' DAGM-2M (daily) long-run volatility (with skewness)
#'
#' Obtains the long-run volatility of the DAGM with two MIDAS variables. 
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m_1 first MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param mv_m_2 second MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K_1 Number of (lagged) realizations of the first MIDAS variable to consider.
#' @param K_2 Number of (lagged) realizations of the second MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the long-run volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' start_val<-c(0.01,0.80,0.05,0.2,0.1,1.1,0.4,1.1,0.5,1.1,0,1.1)
#' r_t<-sp500['2005/2010']
#' mv_m_1<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' mv_m_2<-mv_into_mat(r_t,diff(indpro),K=24,"monthly")
#' head(DAGM_2M_long_run(start_val,r_t,mv_m_1,mv_m_2,K_1=12,K_2=24))
#' }
#' @keywords internal
#' @export

DAGM_2M_long_run<-function(param,
daily_ret,
mv_m_1,
mv_m_2,
K_1,
K_2,
lag_fun="Beta"){

        alpha           	<- param[1]
        beta            	<- param[2]
        gamma_1         	<- param[3]
        m               	<- param[4]
        theta_pos_1      	<- param[5]
        w1_pos_1          	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_1          	<- param[6]
        theta_neg_1       	<- param[7]
        w1_neg_1          	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_1          	<- param[8]
 		theta_pos_2      	<- param[9]
        w1_pos_2          	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_2          	<- param[10]
        theta_neg_2       <- param[11]
        w1_neg_2         	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_2          	<- param[12]
	   
        TT              <- length(daily_ret)

	   MV_pos_1		<- ifelse(mv_m_1>=0,mv_m_1,0)
       MV_neg_1	<- ifelse(mv_m_1<0,mv_m_1,0)

       MV_pos_2		<- ifelse(mv_m_2>=0,mv_m_2,0)
       MV_neg_2	<- ifelse(mv_m_2<0,mv_m_2,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg   <- ifelse(zoo::coredata(daily_ret) <0,1,0)

       
###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_pos_1,w2_pos_1))[2:(K_1+1)],0)
betas_neg_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_neg_1,w2_neg_1))[2:(K_1+1)],0)

betas_pos_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_pos_2,w2_pos_2))[2:(K_2+1)],0)
betas_neg_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_neg_2,w2_neg_2))[2:(K_2+1)],0)


tau_d                     <- exp(
m+
theta_pos_1*suppressWarnings(roll_sum(MV_pos_1, c(K_1+1),weights = betas_pos_1))[K_1+1,] + 
theta_neg_1*suppressWarnings(roll_sum(MV_neg_1, c(K_1+1),weights = betas_neg_1))[K_1+1,] + 
theta_pos_2*suppressWarnings(roll_sum(MV_pos_2, c(K_2+1),weights = betas_pos_2))[K_2+1,] + 
theta_neg_2*suppressWarnings(roll_sum(MV_neg_2, c(K_2+1),weights = betas_neg_2))[K_2+1,] 
)


####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### cond vol.

cond_vol<- sqrt(tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

return(cond_vol)

##### end
}

#' DAGM-2M conditional volatility (no skewness)
#'
#' Obtains the conditional volatility of the DAGM with two MIDAS variables. 
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m_1 first MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param mv_m_2 second MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K_1 Number of (lagged) realizations of the first MIDAS variable to consider.
#' @param K_2 Number of (lagged) realizations of the second MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' start_val<-c(0.01,0.80,0.2,0.1,1.1,0.4,1.1,0.5,1.1,0,1.1)
#' r_t<-sp500['2005/2010']
#' mv_m_1<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' mv_m_2<-mv_into_mat(r_t,diff(indpro),K=24,"monthly")
#' head(DAGM_2M_cond_vol_no_skew(start_val,r_t,mv_m_1,mv_m_2,K_1=12,K_2=24))
#' }
#' @keywords internal
#' @export

DAGM_2M_cond_vol_no_skew<-function(param,
daily_ret,
mv_m_1,
mv_m_2,
K_1,
K_2,
lag_fun="Beta"){

        alpha           	<- param[1]
        beta            	<- param[2]
        m               	<- param[3]
        theta_pos_1      	<- param[4]
        w1_pos_1          	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_1          	<- param[5]
        theta_neg_1       	<- param[6]
        w1_neg_1          	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_1          	<- param[7]
 		theta_pos_2      	<- param[8]
        w1_pos_2          	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_2          	<- param[9]
        theta_neg_2       <- param[10]
        w1_neg_2         	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_2          	<- param[11]
	   
        TT              <- length(daily_ret)

	   MV_pos_1		<- ifelse(mv_m_1>=0,mv_m_1,0)
       MV_neg_1	<- ifelse(mv_m_1<0,mv_m_1,0)

       MV_pos_2		<- ifelse(mv_m_2>=0,mv_m_2,0)
       MV_neg_2	<- ifelse(mv_m_2<0,mv_m_2,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       
    
###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_pos_1,w2_pos_1))[2:(K_1+1)],0)
betas_neg_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_neg_1,w2_neg_1))[2:(K_1+1)],0)

betas_pos_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_pos_2,w2_pos_2))[2:(K_2+1)],0)
betas_neg_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_neg_2,w2_neg_2))[2:(K_2+1)],0)


tau_d                     <- exp(
m+
theta_pos_1*suppressWarnings(roll_sum(MV_pos_1, c(K_1+1),weights = betas_pos_1))[K_1+1,] + 
theta_neg_1*suppressWarnings(roll_sum(MV_neg_1, c(K_1+1),weights = betas_neg_1))[K_1+1,] + 
theta_pos_2*suppressWarnings(roll_sum(MV_pos_2, c(K_2+1),weights = betas_pos_2))[K_2+1,] + 
theta_neg_2*suppressWarnings(roll_sum(MV_neg_2, c(K_2+1),weights = betas_neg_2))[K_2+1,] 
)


####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### cond vol.

cond_vol<- sqrt(g_it*tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

return(cond_vol)

##### end
}


#' DAGM-2M (daily) long-run volatility (no skewness)
#'
#' Obtains the long-run volatility of the DAGM with two MIDAS variables. 
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m_1 first MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param mv_m_2 second MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K_1 Number of (lagged) realizations of the first MIDAS variable to consider.
#' @param K_2 Number of (lagged) realizations of the second MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the long-run volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' start_val<-c(0.01,0.80,0.2,0.1,1.1,0.4,1.1,0.5,1.1,0,1.1)
#' r_t<-sp500['2005/2010']
#' mv_m_1<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' mv_m_2<-mv_into_mat(r_t,diff(indpro),K=24,"monthly")
#' head(DAGM_2M_long_run_no_skew(start_val,r_t,mv_m_1,mv_m_2,K_1=12,K_2=24))
#' }
#' @keywords internal
#' @export

DAGM_2M_long_run_no_skew<-function(param,
daily_ret,
mv_m_1,
mv_m_2,
K_1,
K_2,
lag_fun="Beta"){

        alpha           	<- param[1]
        beta            	<- param[2]
        m               	<- param[3]
        theta_pos_1      	<- param[4]
        w1_pos_1          	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_1          	<- param[5]
        theta_neg_1       	<- param[6]
        w1_neg_1          	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_1          	<- param[7]
 		theta_pos_2      	<- param[8]
        w1_pos_2          	<- ifelse(lag_fun=="Beta",1,0)
        w2_pos_2          	<- param[9]
        theta_neg_2       <- param[10]
        w1_neg_2         	<- ifelse(lag_fun=="Beta",1,0)
        w2_neg_2          	<- param[11]
	   
        TT              <- length(daily_ret)

	   MV_pos_1		<- ifelse(mv_m_1>=0,mv_m_1,0)
       MV_neg_1	<- ifelse(mv_m_1<0,mv_m_1,0)

       MV_pos_2		<- ifelse(mv_m_2>=0,mv_m_2,0)
       MV_neg_2	<- ifelse(mv_m_2<0,mv_m_2,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       
       
###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_pos_1,w2_pos_1))[2:(K_1+1)],0)
betas_neg_1<-c(rev(weight_fun(1:(K_1+1),(K_1+1),w1_neg_1,w2_neg_1))[2:(K_1+1)],0)

betas_pos_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_pos_2,w2_pos_2))[2:(K_2+1)],0)
betas_neg_2<-c(rev(weight_fun(1:(K_2+1),(K_2+1),w1_neg_2,w2_neg_2))[2:(K_2+1)],0)


tau_d                     <- exp(
m+
theta_pos_1*suppressWarnings(roll_sum(MV_pos_1, c(K_1+1),weights = betas_pos_1))[K_1+1,] + 
theta_neg_1*suppressWarnings(roll_sum(MV_neg_1, c(K_1+1),weights = betas_neg_1))[K_1+1,] + 
theta_pos_2*suppressWarnings(roll_sum(MV_pos_2, c(K_2+1),weights = betas_pos_2))[K_2+1,] + 
theta_neg_2*suppressWarnings(roll_sum(MV_neg_2, c(K_2+1),weights = betas_neg_2))[K_2+1,] 
)


####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### cond vol.

cond_vol<- sqrt(tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

return(cond_vol)

##### end
}


#' DAGM-X log-likelihood (with skewness)
#'
#' Obtains the log-likelihood of the DAGM-X, with an asymmetric term linked to past negative returns,
#' according to two errors' conditional distributions: Normal and Student-t. 
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param X Additional "X" variable, which must be an "xts" object. Morever, "X" must be observed for the same days of daily_ret.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param distribution The conditional density to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # conditional density of the innovations: normal
#' start_val<-c(0.01,0.80,0.05,0.05,0,0,1.1,0,1.1)
#' r_t<-sp500['2005/2010']
#' X<-rv5['2005/2010']^0.5
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' sum(DAGM_X_loglik(start_val,r_t,X,mv_m,K=12,distribution="norm"))
#' }
#' @keywords internal
#' @export

DAGM_X_loglik<-function(param,
daily_ret,
X,
mv_m,
K,
distribution,
lag_fun="Beta"){

if(distribution=="norm"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
	    z			  <- param[4]
        m               <- param[5]
        theta_pos       <- param[6]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[7]
        theta_neg       <- param[8]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[9]
	   
        TT              <- length(daily_ret)

	   MV_pos		<- ifelse(mv_m>=0,mv_m,0)
         MV_neg		<- ifelse(mv_m<0,mv_m,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg   <- ifelse(zoo::coredata(daily_ret) <0,1,0)

        ll              <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d + z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### loglik

ll <- as.numeric(
stats::dnorm(daily_ret, mean(daily_ret), sqrt(g_it*tau_d), log = TRUE)
)


} else {

  	 alpha            <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
		z			  <- param[4]
        m               <- param[5]
        theta_pos       <- param[6]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[7]
        theta_neg       <- param[8]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[9]
	  v	      	<- param[10]
	   
        TT              <- length(daily_ret)

	   MV_pos		<- ifelse(mv_m>=0,mv_m,0)
         MV_neg		<- ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) <0,1,0)

        ll                                      <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d + z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}


###### variance 

h_it<-zoo::coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+zoo::coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

)

}

return(ll)
##### end
}

#' DAGM-X conditional volatility (with skewness)
#'
#' Obtains the conditional volatility for the DAGM-X, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of estimated values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param X Additional "X" variable, which must be an "xts" object. Morever, "X" must be observed for the same days of daily_ret.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # estimated volatility 
#' est_val<-c(0.01,0.80,0.05,0.05,0,0.1,1.1,-0.3,1.1)
#' r_t<-sp500['2005/2010']
#' X<-rv5['2005/2010']^0.5
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' vol<-DAGM_X_cond_vol(est_val,r_t,X,mv_m,K=12)
#' head(vol)
#' }
#' @keywords internal
#' @export

DAGM_X_cond_vol<-function(param,
daily_ret,
X,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
		z			  <- param[4]
        m               <- param[5]
        theta_pos       <- param[6]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[7]
        theta_neg       <- param[8]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[9]
	   
        TT              <- length(daily_ret)

	   MV_pos		<- ifelse(mv_m>=0,mv_m,0)
         MV_neg		<- ifelse(mv_m<0,mv_m,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg   <- ifelse(zoo::coredata(daily_ret) <0,1,0)

        ll                                      <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d + z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### cond vol.

cond_vol<- sqrt(g_it*tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

return(cond_vol)
##### end
}

#' DAGM-X (daily) long-run volatility (with skewness)
#'
#' Obtains the daily long-run volatility for the DAGM-X, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of estimated values. It must be a eight- or nine- dimensional vector. See the examples below.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param X Additional "X" variable, which must be an "xts" object. Morever, "X" must be observed for the same days of daily_ret.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' est_val<-c(0.01,0.80,0.05,0.05,0,0.1,1.1,-0.3,1.1)
#' r_t<-sp500['/2010']
#' X<-rv5['/2010']^0.5
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' head(DAGM_X_long_run_vol(est_val,r_t,X,mv_m,K=12))
#' }
#' @keywords internal
#' @export

DAGM_X_long_run_vol<-function(param,
daily_ret,
X,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
		z			  <- param[4]
        m               <- param[5]
        theta_pos       <- param[6]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[7]
        theta_neg       <- param[8]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[9]
	   
        TT              <- length(daily_ret)

	   MV_pos		<- ifelse(mv_m>=0,mv_m,0)
         MV_neg		<- ifelse(mv_m<0,mv_m,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg   <- ifelse(zoo::coredata(daily_ret) <0,1,0)

        ll              <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)

tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d + z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### cond vol.

cond_vol<- sqrt(tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

return(cond_vol)
##### end
}



#' DAGM-X log-likelihood (no skewness)
#'
#' Obtains the log-likelihood of the DAGM-X, according to two errors' conditional distributions: Normal and Student-t. 
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param X Additional "X" variable, which must be an "xts" object. Morever, "X" must be observed for the same days of daily_ret.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param distribution The conditional density to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # conditional density of the innovations: normal
#' start_val<-c(0.01,0.80,0.05,0,0,1.1,0,1.1)
#' r_t<-sp500['2005/2010']
#' X<-rv5['2005/2010']^0.5
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' sum(DAGM_X_loglik_no_skew(start_val,r_t,X,mv_m,K=12,distribution="norm"))
#' }
#' @keywords internal
#' @export

DAGM_X_loglik_no_skew<-function(param,
daily_ret,
X,
mv_m,
K,
distribution,
lag_fun="Beta"){

if(distribution=="norm"){

        alpha           <- param[1]
        beta            <- param[2]
        z			  <- param[3]
        m               <- param[4]
        theta_pos       <- param[5]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[6]
        theta_neg       <- param[7]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[8]
	   
        TT              <- length(daily_ret)

	   MV_pos		<- ifelse(mv_m>=0,mv_m,0)
         MV_neg		<- ifelse(mv_m<0,mv_m,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       
        ll              <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d + z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### loglik

ll <- as.numeric(
stats::dnorm(daily_ret, mean(daily_ret), sqrt(g_it*tau_d), log = TRUE)
)


} else {

  	 alpha            <- param[1]
        beta            <- param[2]
		z			  <- param[3]
        m               <- param[4]
        theta_pos       <- param[5]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[6]
        theta_neg       <- param[7]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[8]
	  v	      	<- param[9]
	   
        TT              <- length(daily_ret)

	   MV_pos		<- ifelse(mv_m>=0,mv_m,0)
         MV_neg		<- ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        ll                                      <- 0 


###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d + z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}


###### variance 

h_it<-zoo::coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+zoo::coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

)

}

return(ll)
##### end
}

#' DAGM-X conditional volatility (no skewness)
#'
#' Obtains the conditional volatility for the DAGM-X.
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of estimated values. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param X Additional "X" variable, which must be an "xts" object. Morever, "X" must be observed for the same days of daily_ret.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # estimated volatility 
#' est_val<-c(0.01,0.80,0.05,0,0.1,1.1,-0.3,1.1)
#' r_t<-sp500['2005/2010']
#' X<-rv5['2005/2010']^0.5
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' vol<-DAGM_X_cond_vol_no_skew(est_val,r_t,X,mv_m,K=12)
#' head(vol)
#' }
#' @keywords internal
#' @export

DAGM_X_cond_vol_no_skew<-function(param,
daily_ret,
X,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
        z			  <- param[3]
        m               <- param[4]
        theta_pos       <- param[5]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[6]
        theta_neg       <- param[7]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[8]
	   
        TT              <- length(daily_ret)

	   MV_pos		<- ifelse(mv_m>=0,mv_m,0)
         MV_neg		<- ifelse(mv_m<0,mv_m,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       

###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d + z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### cond vol.

cond_vol<- sqrt(g_it*tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

return(cond_vol)
##### end
}

#' DAGM-X (daily) long-run volatility (no skewness)
#'
#' Obtains the daily long-run volatility for the DAGM-X.
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of estimated values. It must be a eight- or nine- dimensional vector. See the examples below.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param X Additional "X" variable, which must be an "xts" object. Morever, "X" must be observed for the same days of daily_ret.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' est_val<-c(0.01,0.80,0.05,0,0.1,1.1,-0.3,1.1)
#' r_t<-sp500['/2010']
#' X<-rv5['/2010']^0.5
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' head(DAGM_X_long_run_vol_no_skew(est_val,r_t,X,mv_m,K=12))
#' }
#' @keywords internal
#' @export

DAGM_X_long_run_vol_no_skew<-function(param,
daily_ret,
X,
mv_m,
K,
lag_fun="Beta"){

        alpha           <- param[1]
        beta            <- param[2]
		z			  <- param[3]
        m               <- param[4]
        theta_pos       <- param[5]
        w1_pos          <- ifelse(lag_fun=="Beta",1,0)
        w2_pos          <- param[6]
        theta_neg       <- param[7]
        w1_neg          <- ifelse(lag_fun=="Beta",1,0)
        w2_neg          <- param[8]
	   
        TT              <- length(daily_ret)

	   MV_pos		<- ifelse(mv_m>=0,mv_m,0)
         MV_neg		<- ifelse(mv_m<0,mv_m,0)

        tau_d           <- rep(NA,TT)
	   
        g_it            <- rep(1,TT)            # daily conditional variance 
		       



###### long-run 

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas_pos<-c(rev(weight_fun(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(weight_fun(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)

tau_d                     <- exp(
m+
theta_pos*suppressWarnings(roll_sum(MV_pos, c(K+1),weights = betas_pos)) + 
theta_neg*suppressWarnings(roll_sum(MV_neg, c(K+1),weights = betas_neg))
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1<-(1-alpha-beta)+
(alpha)*(daily_ret)^2/tau_d + z*X

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### cond vol.

cond_vol<- sqrt(tau_d)
cond_vol<-as.xts(cond_vol,index(daily_ret))

return(cond_vol)
##### end
}


#' MEM-MIDAS log-likelihood (with skewness parameter)
#'
#' Obtains the log-likelihood of the MEM-MIDAS, with an asymmetric term linked to past negative returns.
#' @param param Vector of starting values. 
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @param daily_ret Daily returns, which must be an "xts" object, and with the same length of x.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' start_val<-c(alpha=0.10,beta=0.8,gamma=0.1,m=0,theta=-0.16,w2=5)
#' r_t<-sp500['/2010']
#' real<-(rv5['/2010'])^0.5		# realized volatility
#' mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' sum(MEM_MIDAS_loglik(start_val,real,r_t,mv_m,K=12))
#' }
#' @keywords internal
#' @export

MEM_MIDAS_loglik<-function(param,
x,
daily_ret,
mv_m,
K){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta			  <- param[5]
        w1			  <- 1
        w2			  <- param[6]
        		
       N              <- length(daily_ret)
		x			<- zoo::coredata(x)
		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1 <- (1-alpha-beta-gamma_1/2)*sample_mean_x+
(alpha+gamma_1*daily_ret_neg)*(x)/tau_d

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

#############################

ll <- - (
log(mu_it*tau_d) + x/(mu_it*tau_d)
)


return(ll)
##### end
}

#' MEM-MIDAS log-likelihood (no skewness parameter)
#'
#' Obtains the log-likelihood of the MEM-MIDAS.
#' @param param Vector of starting values
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function
#' @param K Number of (lagged) realizations of the MIDAS variable to consider
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' start_val<-c(alpha=0.10,beta=0.8,m=0,theta=-0.16,w2=5)
#' real<-(rv5['/2010'])^0.5		# realized volatility
#' mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' sum(MEM_MIDAS_loglik_no_skew(start_val,real,mv_m,K=12))
#' }
#' @keywords internal
#' @export

MEM_MIDAS_loglik_no_skew<-function(param,
x,
mv_m,
K){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta			  <- param[4]
        w1			  <- 1
        w2			  <- param[5]
        		
       N              <- length(x)

		tau_d<-rep(NA,N)
		x			<- zoo::coredata(x)
  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1 <- (1-alpha-beta)*sample_mean_x+
(alpha)*(x)/tau_d

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

#############################

ll <- - (
log(mu_it*tau_d) + x/(mu_it*tau_d)
)


return(ll)
##### end
}

#' MEM-MIDAS-X log-likelihood (with skewness parameter)
#'
#' Obtains the log-likelihood of the MEM-MIDAS-X, with an asymmetric term linked to past negative returns and an additional
#' X part (for instance, the VIX).
#' @param param Vector of starting values
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object
#' @param daily_ret Daily returns, which must be an "xts" object, and with the same length of x
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function
#' @param K Number of (lagged) realizations of the MIDAS variable to consider
#' @param z Additional daily variable which must be an "xts" object, and with the same length of x
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' start_val<-c(alpha=0.10,beta=0.8,gamma=0.1,m=0,theta=-0.16,w2=5,delta=0.1)
#' r_t<-sp500['2010']
#' real<-(rv5['2010'])^0.5		# realized volatility
#' mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' sum(MEM_MIDAS_X_loglik(start_val,real,r_t,mv_m,K=12,z=vix['2010']))
#' }
#' @keywords internal
#' @export

MEM_MIDAS_X_loglik<-function(param,
x,
daily_ret,
mv_m,
K,
z){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta			  <- param[5]
        w1			  <- 1
        w2			  <- param[6]
	    delta			  <- param[7]	 
        		
       N              <- length(daily_ret)
		x			<- zoo::coredata(x)
		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1 <- (1-alpha-beta-gamma_1/2)*sample_mean_x+
(alpha+gamma_1*daily_ret_neg)*(x)/tau_d + delta*z

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

#############################

ll <- - (
log(mu_it*tau_d) + x/(mu_it*tau_d)
)


return(ll)
##### end
}


#' MEM-MIDAS-X one-step-ahead predictions (with skewness parameter)
#'
#' Predicts the dependent variable, usually the realized volatility, for the MEM-MIDAS-X, with an asymmetric term linked to past negative returns and
#' an additional X part (for instance, the VIX).
#' @param param Vector of starting values
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object
#' @param daily_ret Daily returns, which must be an "xts" object, and with the same length of x
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function
#' @param K Number of (lagged) realizations of the MIDAS variable to consider
#' @param z Additional daily variable which must be an "xts" object, and with the same length of x
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

MEM_MIDAS_X_pred<-function(param, 
x,
daily_ret,
mv_m,
K,
z){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta			  <- param[5]
        w1			  <- 1
        w2			  <- param[6]
	    delta			  <- param[7]	 
        		
       N              <- length(daily_ret)
		# x			<- zoo::coredata(x)
		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1 <- (1-alpha-beta-gamma_1/2)*sample_mean_x+
(alpha+gamma_1*daily_ret_neg)*(x)/tau_d + delta*z

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

####### predictions

x_hat<-mu_it*tau_d

x_hat<-as.xts(x_hat,zoo::index(x))

return(x_hat)


##### end
}

#' MEM-MIDAS-X long-run one-step-ahead predictions (with skewness parameter)
#'
#' Predicts the long-run term of the dependent variable, usually the realized volatility, for the MEM-MIDAS-X, 
#' with an asymmetric term linked to past negative returns and an additional X part (for instance, the VIX).
#' @param param Vector of starting values
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object
#' @param daily_ret Daily returns, which must be an "xts" object, and with the same length of x
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function
#' @param K Number of (lagged) realizations of the MIDAS variable to consider
#' @param z Additional daily variable which must be an "xts" object, and with the same length of x
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

MEM_MIDAS_X_lr_pred<-function(param, 
x,
daily_ret,
mv_m,
K,
z){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta			  <- param[5]
        w1			  <- 1
        w2			  <- param[6]
	    delta			  <- param[7]	 
        		
       N              <- length(daily_ret)
		# x			<- zoo::coredata(x)
		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1 <- (1-alpha-beta-gamma_1/2)*sample_mean_x+
(alpha+gamma_1*daily_ret_neg)*(x)/tau_d + delta*z

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

####### predictions

x_hat<-tau_d

x_hat<-as.xts(x_hat,zoo::index(x))

return(x_hat)

##### end
}

#' MEM-MIDAS-X log-likelihood (no skewness parameter)
#'
#' Obtains the log-likelihood of the MEM-MIDAS-X, with an additional X part (for instance, the VIX).
#' @param param Vector of starting values
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function
#' @param K Number of (lagged) realizations of the MIDAS variable to consider
#' @param z Additional daily variable which must be an "xts" object, and with the same length of x
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' start_val<-c(alpha=0.10,beta=0.8,m=0,theta=-0.16,w2=5,delta=0.1)
#' real<-(rv5['2010'])^0.5		# realized volatility
#' mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' sum(MEM_MIDAS_X_loglik_no_skew(start_val,real,mv_m,K=12,z=vix['2010']))
#' }
#' @keywords internal
#' @export

MEM_MIDAS_X_loglik_no_skew<-function(param,
x,
mv_m,
K,
z){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta			  <- param[4]
        w1			  <- 1
        w2			  <- param[5]
	    delta			  <- param[6]	 
        		
       N              <- length(x)
		x			<- zoo::coredata(x)
		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1 <- (1-alpha-beta)*sample_mean_x+
(alpha)*(x)/tau_d + delta*z

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

#############################

ll <- - (
log(mu_it*tau_d) + x/(mu_it*tau_d)
)


return(ll)
##### end
}

#' MEM-MIDAS-X one-step-ahead predictions (no skewness parameter)
#'
#' Predicts the dependent variable, usually the realized volatility, for the MEM-MIDAS-X, with an additional X part (for instance, the VIX).
#' @param param Vector of starting values
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function
#' @param K Number of (lagged) realizations of the MIDAS variable to consider
#' @param z Additional daily variable which must be an "xts" object, and with the same length of x
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

MEM_MIDAS_X_pred_no_skew<-function(param,
x,
mv_m,
K,
z){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta			  <- param[4]
        w1			  <- 1
        w2			  <- param[5]
	    delta			  <- param[6]	 
        		
       N              <- length(x)
		# x			<- zoo::coredata(x)
		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1 <- (1-alpha-beta)*sample_mean_x+
(alpha)*(x)/tau_d + delta*z

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

####### predictions

x_hat<-mu_it*tau_d

x_hat<-as.xts(x_hat,zoo::index(x))

return(x_hat)

##### end
}

#' MEM-MIDAS-X long-run one-step-ahead predictions (no skewness parameter)
#'
#' Predicts the long-run term of the dependent variable, usually the realized volatility, for the MEM-MIDAS-X, with an additional X part (for instance, the VIX).
#' @param param Vector of starting values
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function
#' @param K Number of (lagged) realizations of the MIDAS variable to consider
#' @param z Additional daily variable which must be an "xts" object, and with the same length of x
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

MEM_MIDAS_X_lr_pred_no_skew<-function(param,
x,
mv_m,
K,
z){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta			  <- param[4]
        w1			  <- 1
        w2			  <- param[5]
	    delta			  <- param[6]	 
        		
       N              <- length(x)
		# x			<- zoo::coredata(x)
		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1 <- (1-alpha-beta)*sample_mean_x+
(alpha)*(x)/tau_d + delta*z

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

####### predictions

x_hat<-tau_d

x_hat<-as.xts(x_hat,zoo::index(x))

return(x_hat)


##### end
}


#' MEM log-likelihood (with skewness parameter)
#'
#' Obtains the log-likelihood of the base MEM, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{engle_gallo_2006;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @param daily_ret Daily returns, which must be an "xts" object, and with the same length of x.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' start_val<-c(alpha=0.10,beta=0.8,gamma=0.05)
#' real<-(rv5['/2010'])^0.5		# realized volatility
#' r_t<-sp500['/2010']
#' sum(MEM_loglik(start_val,real,r_t))
#' }
#' @keywords internal
#' @export

MEM_loglik<-function(param,
x,
daily_ret){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
                		
       N              <- length(daily_ret)
			x		<- zoo::coredata(x)
		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

#######  

step_1 <- (1-alpha-beta-gamma_1/2)*sample_mean_x+
(alpha+gamma_1*daily_ret_neg)*(x)

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

#############################

ll <- - (
log(mu_it) + x/(mu_it)
)


return(ll)
##### end
}

#' MEM-X log-likelihood (with skewness parameter)
#'
#' Obtains the log-likelihood of the base MEM, with an asymmetric term linked to past negative returns and an additional
#' X part (for instance, the VIX).
#' @param param Vector of starting values. 
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @param daily_ret Daily returns, which must be an "xts" object, and with the same length of x.
#' @param z Additional daily variable which must be an "xts" object, and with the same length of x.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' start_val<-c(alpha=0.10,beta=0.8,gamma=0.05,delta=0.01)
#' real<-(rv5['/2010'])^0.5		# realized volatility
#' r_t<-sp500['/2010']
#' z<-vix['2010']
#' sum(MEM_X_loglik(start_val,real,r_t,z))
#' }
#' @keywords internal
#' @export

MEM_X_loglik<-function(param,
x,
daily_ret,
z){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        delta           <- param[4]   
     		
       N              <- length(daily_ret)
			x		<- zoo::coredata(x)
		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

#######  

step_1 <- (1-alpha-beta-gamma_1/2)*sample_mean_x+
(alpha+gamma_1*daily_ret_neg)*(x) + delta*z

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

#############################

ll <- - (
log(mu_it) + x/(mu_it)
)


return(ll)
##### end
}

#' MEM-X one-step-ahead predictions (with skewness parameter)
#'
#' Predicts the dependent variable, usually the realized volatility, for the base MEM, with an asymmetric term 
#' linked to past negative returns and an additional X part (for instance, the VIX).
#' @param param Vector of estimated values. 
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @param daily_ret Daily returns, which must be an "xts" object, and with the same length of x.
#' @param z Additional daily variable which must be an "xts" object, and with the same length of x.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

MEM_X_pred<-function(param,
x,
daily_ret,
z){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        delta           <- param[4]   
     		
       N              <- length(daily_ret)
			
		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

#######  

step_1 <- (1-alpha-beta-gamma_1/2)*sample_mean_x+
(alpha+gamma_1*daily_ret_neg)*(x) + delta*z

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

############################# predictions


x_hat<-mu_it

x_hat<-as.xts(x_hat,zoo::index(x))

return(x_hat)

##### end
}

#' MEM-X log-likelihood (no skewness parameter)
#'
#' Obtains the log-likelihood of the base MEM, with an additional X part (for instance, the VIX).
#' @param param Vector of starting values. 
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @param z Additional daily variable which must be an "xts" object, and with the same length of x.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' start_val<-c(alpha=0.10,beta=0.8,delta=0.01)
#' real<-(rv5['/2010'])^0.5		# realized volatility
#' z<-vix['2010']
#' sum(MEM_X_loglik_no_skew(start_val,real,z))
#' }
#' @keywords internal
#' @export

MEM_X_loglik_no_skew<-function(param,
x,
z){

        alpha           <- param[1]
        beta            <- param[2]
        delta           <- param[3]   
     		
       N              <- length(x)
			x		<- zoo::coredata(x)
		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        ll                                      <- 0 

#######  

step_1 <- (1-alpha-beta)*sample_mean_x+
(alpha)*(x) + delta*z

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

#############################

ll <- - (
log(mu_it) + x/(mu_it)
)


return(ll)
##### end
}

#' MEM-X one-step-ahead predictions (no skewness parameter)
#'
#' Predicts the dependent variable, usually the realized volatility, for the base MEM, with an additional X part (for instance, the VIX).
#' @param param Vector of estimated values. 
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @param z Additional daily variable which must be an "xts" object, and with the same length of x.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

MEM_X_pred_no_skew<-function(param,
x,
z){

        alpha           <- param[1]
        beta            <- param[2]
        delta           <- param[3]   
     		
       N              <- length(x)
			
		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        ll                                      <- 0 

#######  

step_1 <- (1-alpha-beta)*sample_mean_x+
(alpha)*(x) + delta*z

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

############################# predictions


x_hat<-mu_it

x_hat<-as.xts(x_hat,zoo::index(x))

return(x_hat)

##### end
}


#' MEM log-likelihood (no skewness parameter)
#'
#' Obtains the log-likelihood of the base MEM.
#' For details, see \insertCite{engle_gallo_2006;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' start_val<-c(alpha=0.10,beta=0.8)
#' real<-(rv5['/2010'])^0.5		# realized volatility
#' sum(MEM_loglik_no_skew(start_val,real))
#' }
#' @keywords internal
#' @export

MEM_loglik_no_skew<-function(param,x){

        alpha           <- param[1]
        beta            <- param[2]
              		
        N              <- length(x)

		tau_d<-rep(NA,N)
		x			<- zoo::coredata(x)
  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        ll                                      <- 0 

#######  

step_1 <- (1-alpha-beta)*sample_mean_x+
(alpha)*(x)

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

#############################

ll <- - (
log(mu_it) + x/(mu_it)
)


return(ll)
##### end
}

#' MEM one-step-ahead predictions (with skewness parameter)
#'
#' Predicts the dependent variable, usually the realized volatility, for the base MEM, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{engle_gallo_2006;textual}{rumidas}.
#' @param param Vector of estimated values. 
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @param daily_ret Daily returns, which must be an "xts" object, and with the same length of x.
#' @return The resulting vector is the one-step-ahead prediction for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' est_val<-c(alpha=0.10,beta=0.8,gamma=0.05)
#' real<-(rv5['/2010'])^0.5		# realized volatility
#' r_t<-sp500['/2010']
#' head(MEM_pred(est_val,real,r_t))
#' }
#' @keywords internal
#' @export

MEM_pred<-function(param,
x,
daily_ret){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
                		
       N              <- length(daily_ret)

		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

#######  

step_1 <- (1-alpha-beta-gamma_1/2)*sample_mean_x+
(alpha+gamma_1*daily_ret_neg)*(x)

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

####### predictions

x_hat<-mu_it

x_hat<-as.xts(x_hat,zoo::index(x))

return(x_hat)

##### end
}


#' MEM one-step-ahead predictions (no skewness parameter)
#'
#' Predicts the dependent variable, usually the realized volatility, for the base MEM.
#' For details, see \insertCite{engle_gallo_2006;textual}{rumidas}.
#' @param param Vector of starting values. 
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' est_val<-c(alpha=0.10,beta=0.8)
#' real<-(rv5['/2010'])^0.5		# realized volatility
#' head(MEM_pred_no_skew(est_val,real))
#' }
#' @keywords internal
#' @export

MEM_pred_no_skew<-function(param,x){

        alpha           <- param[1]
        beta            <- param[2]
              		
        N              <- length(x)

		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        ll                                      <- 0 

#######  

step_1 <- (1-alpha-beta)*sample_mean_x+
(alpha)*(x)

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}


####### predictions

x_hat<-mu_it

x_hat<-as.xts(x_hat,zoo::index(x))

return(x_hat)

##### end
}


#' MEM-MIDAS one-step-ahead predictions (with skewness parameter)
#'
#' Predicts the dependent variable, usually the realized volatility, for the MEM-MIDAS, with an asymmetric term linked to past negative returns.
#' @param param Vector of estimated values. It must be a six--dimensional vector. See the example below.
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @return The resulting vector is the one-step-ahead prediction for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' est_val<-c(alpha=0.10,beta=0.8,gamma=0.1,m=0,theta=-0.16,w2=5)
#' r_t<-sp500['/2010']
#' real<-rv5['/2010']^0.5 # realized volatility
#' mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' est_vol<-MEM_MIDAS_pred(est_val,real,r_t,mv_m,K=12)
#' head(est_vol)
#' }
#' @keywords internal
#' @export

MEM_MIDAS_pred<-function(param,
x,
daily_ret,
mv_m,
K){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta			  <- param[5]
        w1			  <- 1
        w2			  <- param[6]
        		
       N              <- length(daily_ret)

		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        daily_ret_neg                       <- ifelse(zoo::coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1 <- (1-alpha-beta-gamma_1/2)*sample_mean_x+
(alpha+gamma_1*daily_ret_neg)*(x)/tau_d

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

####### predictions

x_hat<-mu_it*tau_d

x_hat<-as.xts(x_hat,zoo::index(x))

return(x_hat)

##### end
}

#' MEM-MIDAS one-step-ahead predictions (no skewness parameter)
#'
#' Predicts the dependent variable, usually the realized volatility, for the MEM-MIDAS.
#' @param param Vector of starting values. 
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @return The resulting vector is the one-step-ahead prediction for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' est_val<-c(alpha=0.10,beta=0.8,m=0,theta=-0.16,w2=5)
#' real<-(rv5['/2010'])^0.5		# realized volatility
#' mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' sum(MEM_MIDAS_pred_no_skew(est_val,real,mv_m,K=12))
#' }
#' @keywords internal
#' @export

MEM_MIDAS_pred_no_skew<-function(param,
x,
mv_m,
K){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta			  <- param[4]
        w1			  <- 1
        w2			  <- param[5]
        		
       N              <- length(x)

		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1 <- (1-alpha-beta)*sample_mean_x+
(alpha)*(x)/tau_d

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

####### predictions

x_hat<-mu_it*tau_d

x_hat<-as.xts(x_hat,zoo::index(x))

return(x_hat)

##### end
}

#' MEM-MIDAS long-run one-step-ahead predictions (with skewness parameter)
#'
#' Predicts the long-run term of the dependent variable, usually the realized volatility, for the MEM-MIDAS, with an asymmetric term linked to past negative returns.
#' @param param Vector of estimated values. It must be a six--dimensional vector. See the example below.
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' est_val<-c(alpha=0.10,beta=0.8,gamma=0.1,m=0,theta=-0.16,w2=5)
#' r_t<-sp500['/2010']
#' real<-(rv5['/2010'])^0.5		# realized volatility
#' mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' est_vol<-MEM_MIDAS_pred(est_val,real,r_t,mv_m,K=12)
#' head(est_vol)
#' }
#' @keywords internal
#' @export

MEM_MIDAS_lr_pred<-function(param,
x,
daily_ret,
mv_m,
K){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta			  <- param[5]
        w1			  <- 1
        w2			  <- param[6]
        		
       N              <- length(daily_ret)

		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        daily_ret_neg                       <- ifelse(coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1 <- (1-alpha-beta-gamma_1/2)*sample_mean_x+
(alpha+gamma_1*daily_ret_neg)*(x)/tau_d

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

####### predictions

x_hat<-tau_d

x_hat<-as.xts(x_hat,zoo::index(x))

return(x_hat)

##### end
}

#' MEM-MIDAS long-run one-step-ahead predictions (no skewness parameter)
#'
#' Predicts the long-run term of the dependent variable, usually the realized volatility, for the MEM-MIDAS.
#' @param param Vector of starting values. 
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' est_val<-c(alpha=0.10,beta=0.8,m=0,theta=-0.16,w2=5)
#' real<-(rv5['/2010'])^0.5		# realized volatility
#' mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' sum(MEM_MIDAS_lr_pred_no_skew(est_val,real,mv_m,K=12))
#' }
#' @keywords internal
#' @export

MEM_MIDAS_lr_pred_no_skew<-function(param,
x,
mv_m,
K){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta			  <- param[4]
        w1			  <- 1
        w2			  <- param[5]
        		
       N              <- length(x)

		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

tau_d                      <- exp(
m+
theta*suppressWarnings(roll_sum(mv_m, c(K+1),weights = betas)) 
)

tau_d<-tau_d[(K+1),]	 

####### short-run 

step_1 <- (1-alpha-beta)*sample_mean_x+
(alpha)*(x)/tau_d

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

####### predictions

x_hat<-tau_d

x_hat<-as.xts(x_hat,zoo::index(x))

return(x_hat)

##### end
}


#' Methods for obtaining (and evaluating) a variety of GARCH-MIDAS-based models
#'
#' Estimates several GARCH-MIDAS-based models, according to two errors' conditional distributions: Normal and Student-t, and 
#' the presence of asymmetric terms in the short- and long-run components.
#' @param model Model to estimate. Valid choices are: "GM" for GARCH-MIDAS, "GMX" for GARCH-MIDAS-X, 
#' "DAGM" for Double Asymmetric GARCH-MIDAS (DAGM), "DAGM2M" for DAGM with two MIDAS variables, and "DAGMX" for DAGM-X
#' @param skew The skewness parameter to include in the short--run equation. Valid choices are: "YES" and "NO"
#' @param distribution The conditional density to use for the innovations. At the moment, valid choices are "norm" and "std", 
#' for the Normal and Student-t distribution, respectively
#' @param daily_ret Daily returns, which must be an "xts" objectparamter
#' @param X **optional** Additional "X" variable, which must be an "xts" object. Moreover, "X" must be observed for the same days of daily_ret
#' @param mv_m (first) MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function
#' @param mv_m_2 **optional** second MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function
#' @param lag_fun **optional** Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively
#' @param K Number of (lagged) realizations of the (first) MIDAS variable to consider
#' @param K_2 **optional** Number of (lagged) realizations of the (second) MIDAS variable to consider
#' @param out_of_sample **optional** A positive integer indicating the number of periods before the last to keep for out of sample forecasting
#' @param vol_proxy **optional** If present, the 'vol_proxy' is the volatility proxy used for the in-sample and out-of-sample (again, if present) evaluation. 
#' It could be the realized variance. If it left unspecified, vol_proxy is replaced by the squared daily returns
#' @param R **optional**.  A positive integer indicating the number of replications used to find the best starting values. 
#' Equal to 100 by default
#' @return \code{ugmfit} returns an object of class 'rumidas'. The function \code{\link{summary.rumidas}} 
#' can be used to print a summary of the results. Moreover, an object of class 'rumidas' is a list containing the following components:
#' \itemize{
#' 	\item model: The model used for the estimation.
#'   \item rob_coef_mat: The matrix of estimated coefficients, with the QML standard errors. 
#'    For details, see: \insertCite{Bollerslev_Wooldridge_1992;textual}{rumidas}. 
#'   \item obs: The number of daily observations used for the (in-sample) estimation.
#'     \item period: The period of the in-sample estimation.
#'   \item loglik: The value of the log-likelihood at the maximum.
#'	\item inf_criteria: The AIC and BIC information criteria.
#'	\item loss_in_s: The in-sample MSE and QLIKE averages, calculated considering the distance with respect to the volatility proxy (if provided) or the squared daily returns.
#'	\item est_in_s: The one-step-ahead volatility, for the in-sample period, that is: \eqn{\sqrt{\hat{\tau}_t \times \hat{g}_{i,t} }}. 
#'	\item est_lr_in_s: The one-step-ahead long-run volatility, for the in-sample period.
#'	\item loss_oos: The out-of-sample MSE and QLIKE averages, calculated considering the distance with respect to the volatility proxy (if provided) or the squared daily returns.
#'	\item est_oos: The one-step-ahead volatility, for the out-of-sample period, that is: \eqn{\sqrt{\hat{\tau}_t \times \hat{g}_{i,t} }}.
#'	\item est_lr_oos: The one-step-ahead long-run volatility, for the out-of-sample period.
#' }
#' @importFrom Rdpack reprompt
#' @importFrom stats time
#' @importFrom zoo coredata
#' @import maxLik
#' @seealso \code{\link{mv_into_mat}}.
#' 
#' @details
#' Function \code{ugmfit} implements the estimation and evaluation of the GARCH--MIDAS--based models, with and without the asymmetric term 
#' linked to negative lagged daily returns, according to two distributions for the error term. The general framework assumes that: 
#' \deqn{r_{i,t}= \sqrt{\tau_t \times g_{i,t} } \epsilon_{i,t},}
#' where 
#' \itemize{
#' \item \eqn{r_{i,t}} is the daily return for the \eqn{i}-th day (\eqn{i = 1, \ldots, N_t}) 
#' of the period \eqn{t} (for example, a week, a month or a quarter; \eqn{t = 1 , \ldots, T});
#' \item \eqn{\tau_{t}} is the long-run component, varying each period \eqn{t};
#' \item \eqn{g_{i,t}} is the short--run term, varying each day \eqn{i} of the period \eqn{t};
#' \item \eqn{\epsilon_{i,t}} is an \eqn{iid} error term which has a zero mean and unit variance.
#' }
#' The short--run component of the GARCH--MIDAS (parameter "model" set to "GM"), DAGM (parameter "model" set to "DAGM"), and DAGM with two
#' MIDAS variables (parameter "model" set to "DAGM2M") models is as follows.
#' When the parameter "skew" is present (set to "YES"):
#' \deqn{g_{i,t} = \left(1-\alpha-\gamma/2-\beta\right) + \left(\alpha + \gamma \cdot I_{\left(r_{i-1,t}  < 0 \right)}\right) \frac{\left(r_{i-1,t}\right)^2}{\tau_t} + \beta g_{i-1,t},}
#' where \eqn{I_{(\cdot)}} is an indicator function. 
#' The short--run component of the GARCH--MIDAS--X (parameter "model" set to "GMX") and DAGM--X (parameter "model" set to "DAGMX"),
#' when the parameter "skew" is set to "YES", is:
#' \deqn{g_{i,t} = \left(1-\alpha-\gamma/2-\beta\right) + \left(\alpha + \gamma \cdot I_{\left(r_{i-1,t}  < 0 \right)}\right) \frac{\left(r_{i-1,t}\right)^2}{\tau_t} + \beta g_{i-1,t} + z \cdot \left(X_{i-1,t}- E(X_{i-1,t})\right).}
#' When, for the models GARCH--MIDAS, GARCH--MIDAS--X, DAGM, and DAGM--X, the parameter "skew" is set to "NO", parameter \eqn{\gamma} disappears.
#' For details on the GARCH--MIDAS--X and DAGM--X models, see \insertCite{amendola_candila_gallo_2020;textual}{rumidas}.
#' The long-run component of the GARCH-MIDAS and GARCH--MIDAS--X models is:
#' \deqn{\tau_{t} = \exp \left\{ m + \theta \sum_{j=1}^K \delta_{j}(\omega) X_{t-j}\right\},}
#' where \eqn{X_{t}} is the MIDAS term and \eqn{\delta_{j}(\omega)} is the chosen weighting function, 
#' which can be the Beta (\code{\link{beta_function}}) or Exponential Almon lag (\code{\link{exp_almon}}) functions.  
#' The long-run component of the DAGM and DAGM--X models is:
#' \deqn{\tau_t   =  \exp \left(m + \theta^{+}  \sum_{k=1}^K \delta_k(\omega)^{+}  X_{t-k} I_{\left( X_{t-k} \geq 0 \right)} +  \theta^{-}  \sum_{k=1}^K \delta_k(\omega)^{-} X_{t-k} I_{\left( X_{t-k} < 0 \right)} \right).}
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # estimate a GARH-MIDAS model, without the skewness parameter
#' r_t<-sp500['2008']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly") 
#' fit<-ugmfit(model="GM",skew="NO",distribution="norm",r_t,mv_m,K=12)
#' fit
#' summary.rumidas(fit)
#' names(fit)
#' 
#' # to see the estimated coefficients with the QML standard errors:
#' fit$rob_coef_mat
#'
#' # estimate a DAGM model, with the skewness parameter, 
#' # including the volatility proxy (realized variance), and
#' # leaving the last 100 observations for the out-of-sample evaluation
#' r_t<-sp500['2002/2020']
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly") 
#' fit_2<-ugmfit(model="DAGM",skew="YES",distribution="norm",r_t,
#' mv_m,K=12,vol_proxy=rv5['2002/2020'],out_of_sample=100)
#' fit_2
#' summary.rumidas(fit_2)
#' # estimate a GM-X model, without the skewness parameter
#' r_t<-sp500['2010/2013']
#' X<-vix['2010/2013'] 
#' mv_m<-mv_into_mat(r_t,diff(indpro),K=36,"monthly") 
#' fit_3<-ugmfit(model="GMX",skew="NO",distribution="norm",r_t,mv_m,K=36,X=X)
#' summary.rumidas(fit_3)
#' }
#' @export

ugmfit<-function(
model,
skew,
distribution,
daily_ret,
mv_m,
mv_m_2=NULL,
K,
K_2=NULL,
lag_fun="Beta",
X=NULL,
out_of_sample=NULL,
vol_proxy=NULL,
R=100
){

############################# check on valid choices

if((model != "GM")&(model != "GMX")&(model != "DAGM")&(model != "DAGM2M")&(model != "DAGMX")) { stop(cat("#Warning:\n Valid choices for the parameter 'model' are 'GM','GMX','DAGM', 'DAGM2M' and 'DAGMX' \n"))}
if((model == "GMX"|model == "DAGMX")&(missing(X))) { stop(cat("#Warning:\n If the model chosen includes the 'X' term, then the 'X' variable has to be provided \n"))}
if((model == "DAGM2M")&(missing(mv_m_2)|missing(K_2))) { stop(cat("#Warning:\n If the model is the DAGM with two MIDAS variables (DAGM2M), then the 'mv_m_2' and 'K_2' parameters have to be provided \n"))}
if((skew != "YES")&(skew != "NO")) { stop(cat("#Warning:\n Valid choices for the parameter 'skew' are 'YES' and 'NO' \n"))}
if((distribution != "norm")&(distribution != "std")) { stop(cat("#Warning:\n Valid choices for the parameter 'distribution' are 'norm' and 'std' \n"))}
if((lag_fun != "Beta")&(lag_fun!= "Almon")) { stop(cat("#Warning:\n Valid choices for the parameter 'lag_fun' are 'Beta' and 'Almon' \n"))}

N<-length(daily_ret)

cond_r_t<- class(daily_ret)[1]
cond_mv_m<- class(mv_m)[1]


if(cond_r_t != "xts") { stop(
cat("#Warning:\n Parameter 'daily_ret' must be an xts object. Please provide it in the correct form \n")
)}

if(cond_mv_m != "matrix") { stop(
cat("#Warning:\n Parameter 'mv_m' must be a matrix. Please provide it in the correct form \n")
)}

if(dim(mv_m)[2] != N) { stop(
cat("#Warning:\n The columns of the matrix 'mv_m' must be equal to the length of vector 'daily_ret'. Please provide it in the correct form \n")
)}

if(!missing(mv_m_2)){
cond_mv_m_2<- class(mv_m_2)[1]

if(cond_mv_m_2 != "matrix") { stop(
cat("#Warning:\n Parameter 'mv_m_2' must be a matrix. Please provide it in the correct form \n")
)}

if(dim(mv_m_2)[2] != N) { stop(
cat("#Warning:\n The columns of the matrix 'mv_m_2' must be equal to the length of vector 'daily_ret'. Please provide it in the correct form \n")
)}

}

############## check if the vol_proxy is provided and if it has the same time span of daily_ret

if(!missing(vol_proxy)){ 
if(any(range(time(daily_ret))!=range(time(vol_proxy)))){
stop(
cat("#Warning:\n The vector 'vol_proxy' has to be observed during the same time span of 'daily_ret' \n")
)}}

############## check if the X is provided and if it has the same time span of daily_ret

if(!missing(X)){ 
if(any(range(time(daily_ret))!=range(time(X)))){
stop(
cat("#Warning:\n The vector 'X' has to be observed during the same time span of 'daily_ret' \n")
)}}

############## check if the X is provided and if it has the same length of daily_ret

if(!missing(X)){ 
if(length(X)!=length(daily_ret)){
stop(
cat("#Warning:\n The vector 'X' must have the same length of 'daily_ret' \n")
)}}

############## in sample and out of sample

if (missing(out_of_sample)){ # out_of_sample parameter is missing

r_t_in_s<-daily_ret
mv_m_in_s<-mv_m

if(!missing(mv_m_2)){
mv_m_2_in_s<-mv_m_2
}

if(missing(vol_proxy)){
vol_proxy_in_s<-r_t_in_s^2
} else {
vol_proxy_in_s<-vol_proxy
}
X_in_s<-X

} else {					# out_of_sample parameter is present

r_t_in_s<-daily_ret[1:(N-out_of_sample)]
mv_m_in_s<-mv_m[,1:(N-out_of_sample)]

if(!missing(mv_m_2)){
mv_m_2_in_s<-mv_m_2[,1:(N-out_of_sample)]
}

if(missing(vol_proxy)){
vol_proxy_in_s<-r_t_in_s^2
} else {
vol_proxy_in_s<-vol_proxy[1:(N-out_of_sample)]
}

X_in_s<-X[1:(N-out_of_sample)]
}

############################ setting of the parameters for each model

if (model=="GM"&skew=="YES"&distribution=="norm"&lag_fun=="Beta"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=6)
colnames(begin_val)<-c("alpha","beta","gamma","m","theta","w2")
begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-stats::runif(R,min=-1,max=1)
begin_val[,6]<-1.01

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(GM_loglik(begin_val[i,],
r_t_in_s,mv_m_in_s,K=K,distribution="norm",lag_fun="Beta"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0),	     		## alpha+beta+gamma/2<1
c(0,0,0,0,0,1))				## w2>1.001

ci<-c(-0.0001,-0.001,0.999,-1.001)

LOGLIK<-GM_loglik
cond_vol<-GM_cond_vol
long_run_vol<-GM_long_run_vol

} else if (model=="GM"&skew=="YES"&distribution=="std"&lag_fun=="Beta"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=7)
colnames(begin_val)<-c("alpha","beta","gamma","m","theta","w2","shape")
begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-stats::runif(R,min=-1,max=1)
begin_val[,6]<-1.01
begin_val[,7]<-stats::runif(R,min=2.01,max=10)

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(GM_loglik(begin_val[i,],
r_t_in_s,mv_m_in_s,K=K,distribution="std",lag_fun="Beta"))
}

start_val<-begin_val[which.max(which_row),]


ui<-rbind(
c(1,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0,0),	     			## alpha+beta+gamma/2<1
c(0,0,0,0,0,1,0),				 	## w2>1.001
c(0,0,0,0,0,0,1))					## shape > 2.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-2.001)

LOGLIK<-GM_loglik
cond_vol<-GM_cond_vol
long_run_vol<-GM_long_run_vol

} else if (model=="GM"&skew=="YES"&distribution=="norm"&lag_fun=="Almon"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=6)
colnames(begin_val)<-c("alpha","beta","gamma","m","theta","w2")
begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-stats::runif(R,min=-1,max=1)
begin_val[,6]<- -0.1

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(GM_loglik(begin_val[i,],
r_t_in_s,mv_m_in_s,K=K,distribution="norm",lag_fun="Almon"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0),	     		## alpha+beta+gamma/2<1
c(0,0,0,0,0,-1))				## w2<0

ci<-c(-0.0001,-0.001,0.999,0)

LOGLIK<-GM_loglik
cond_vol<-GM_cond_vol
long_run_vol<-GM_long_run_vol

} else if (model=="GM"&skew=="YES"&distribution=="std"&lag_fun=="Almon"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=7)
colnames(begin_val)<-c("alpha","beta","gamma","m","theta","w2","shape")
begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-stats::runif(R,min=-1,max=1)
begin_val[,6]<- -0.01
begin_val[,7]<-stats::runif(R,min=2.01,max=10)

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(GM_loglik(begin_val[i,],
r_t_in_s,mv_m_in_s,K=K,distribution="std",lag_fun="Almon"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0,0),	     			## alpha+beta+gamma/2<1
c(0,0,0,0,0,-1,0),				## w2<0
c(0,0,0,0,0,0,1))					## shape > 2.001

ci<-c(-0.0001,-0.001,0.999,0,-2.001)

LOGLIK<-GM_loglik
cond_vol<-GM_cond_vol
long_run_vol<-GM_long_run_vol

} else if (model=="GMX"&skew=="YES"&distribution=="norm"&lag_fun=="Beta"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.8,gamma=0.05,z=0.05,m=0,theta=0,w2=2)

ui<-rbind(
c(1,0,0,0,0,0,0),       	 	 ## alpha>0.0001
c(0,1,0,0,0,0,0),        		 ## beta>0.001
c(-1,-1,-0.5,0,0,0,0),	     ## alpha+beta+gamma/2<1
c(0,0,0,0,0,0,1))				## w2>1.001

ci<-c(-0.0001,-0.001,0.999,-1.001)

LOGLIK<-GM_X_loglik
cond_vol<-GM_X_cond_vol
long_run_vol<-GM_X_long_run_vol

} else if (model=="GMX"&skew=="YES"&distribution=="std"&lag_fun=="Beta"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.8,gamma=0.05,z=0.05,m=0,theta=0,w2=2,shape=5)


ui<-rbind(
c(1,0,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0),	     			## alpha+beta+gamma/2<1
c(0,0,0,0,0,0,1,0),				 	## w2>1.001
c(0,0,0,0,0,0,0,1))					## shape > 2.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-2.001)

LOGLIK<-GM_X_loglik
cond_vol<-GM_X_cond_vol
long_run_vol<-GM_X_long_run_vol

} else if (model=="GMX"&skew=="YES"&distribution=="norm"&lag_fun=="Almon"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.8,gamma=0.05,z=0.05,m=0,theta=0,w2=-0.1)


ui<-rbind(
c(1,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0,0),	     		## alpha+beta+gamma/2<1
c(0,0,0,0,0,0,-1))			## w2<0

ci<-c(-0.0001,-0.001,0.999,0)

LOGLIK<-GM_X_loglik
cond_vol<-GM_X_cond_vol
long_run_vol<-GM_X_long_run_vol

} else if (model=="GMX"&skew=="YES"&distribution=="std"&lag_fun=="Almon"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.8,gamma=0.05,z=0.05,m=0,theta=0,w2=-0.1,shape=5)

ui<-rbind(
c(1,0,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0),	     			## alpha+beta+gamma/2<1
c(0,0,0,0,0,0,-1,0),				## w2<0
c(0,0,0,0,0,0,0,1))					## shape > 2.001

ci<-c(-0.0001,-0.001,0.999,0,-2.001)

LOGLIK<-GM_X_loglik
cond_vol<-GM_X_cond_vol
long_run_vol<-GM_X_long_run_vol

} else if (model=="GM"&skew=="NO"&distribution=="norm"&lag_fun=="Beta"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=5)
colnames(begin_val)<-c("alpha","beta","m","theta","w2")
begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=-1,max=1)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-1.01

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(GM_loglik_no_skew(begin_val[i,],
r_t_in_s,mv_m_in_s,K=K,distribution="norm",lag_fun="Beta"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0),       	 			## alpha>0.001
c(0,1,0,0,0),        		 		## beta>0.001
c(-1,-1,0,0,0),	     				## alpha+beta<1
c(0,0,0,0,1))				 	## w2>1.001

ci<-c(-0.001,-0.001,0.999,-1.001)

LOGLIK<-GM_loglik_no_skew
cond_vol<-GM_cond_vol_no_skew
long_run_vol<-GM_long_run_vol_no_skew

} else if (model=="GM"&skew=="NO"&distribution=="std"&lag_fun=="Beta"){ 


start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=6)
colnames(begin_val)<-c("alpha","beta","m","theta","w2","shape")
begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=-1,max=1)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-1.01
begin_val[,6]<-stats::runif(R,min=2.01,max=10)

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(GM_loglik_no_skew(begin_val[i,],
r_t_in_s,mv_m_in_s,K=K,distribution="std",lag_fun="Beta"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0),       	 			## alpha>0.0001
c(0,1,0,0,0,0),        		 		## beta>0.001
c(-1,-1,0,0,0,0),	     				## alpha+beta<1
c(0,0,0,0,1,0),				 	## w2>1.001
c(0,0,0,0,0,1))					## shape>2.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-2.001)

LOGLIK<-GM_loglik_no_skew
cond_vol<-GM_cond_vol_no_skew
long_run_vol<-GM_long_run_vol_no_skew

} else if (model=="GM"&skew=="NO"&distribution=="norm"&lag_fun=="Almon"){ 


start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=5)
colnames(begin_val)<-c("alpha","beta","m","theta","w2")
begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=-1,max=1)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<- -0.1

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(GM_loglik_no_skew(begin_val[i,],
r_t_in_s,mv_m_in_s,K=K,distribution="norm",lag_fun="Almon"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0),       	 			## alpha>0.001
c(0,1,0,0,0),        		 		## beta>0.001
c(-1,-1,0,0,0),	     				## alpha+beta<1
c(0,0,0,0,-1))				 	## w2<0

ci<-c(-0.001,-0.001,0.999,0)

LOGLIK<-GM_loglik_no_skew
cond_vol<-GM_cond_vol_no_skew
long_run_vol<-GM_long_run_vol_no_skew

 
} else if (model=="GM"&skew=="NO"&distribution=="std"&lag_fun=="Almon"){ 

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=6)
colnames(begin_val)<-c("alpha","beta","m","theta","w2","shape")
begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=-1,max=1)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<- -0.01
begin_val[,6]<-stats::runif(R,min=2.01,max=10)

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(GM_loglik_no_skew(begin_val[i,],
r_t_in_s,mv_m_in_s,K=K,distribution="std",lag_fun="Almon"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0),       	 			## alpha>0.0001
c(0,1,0,0,0,0),        		 		## beta>0.001
c(-1,-1,0,0,0,0),	     				## alpha+beta<1
c(0,0,0,0,-1,0),				 	## w2<0
c(0,0,0,0,0,1))					## shape>2.001

ci<-c(-0.0001,-0.001,0.999,0,-2.001)

LOGLIK<-GM_loglik_no_skew
cond_vol<-GM_cond_vol_no_skew
long_run_vol<-GM_long_run_vol_no_skew

} else if (model=="GMX"&skew=="NO"&distribution=="norm"&lag_fun=="Beta"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.8,z=0,m=0,theta=0,w2=2)

ui<-rbind(
c(1,0,0,0,0,0),       	 	 ## alpha>0.0001
c(0,1,0,0,0,0),        		 ## beta>0.001
c(-1,-1,0,0,0,0),	             ## alpha+beta<1
#c(0,0,1,0,0,0),				## z>0
c(0,0,0,0,0,1))				## w2>1.001
#c(0,0,-1,0,0,0))				## z<z_up_bound

ci<-c(-0.0001,-0.001,0.999,-1.001)

LOGLIK<-GM_X_loglik_no_skew
cond_vol<-GM_X_cond_vol_no_skew
long_run_vol<-GM_X_long_run_vol_no_skew

} else if (model=="GMX"&skew=="NO"&distribution=="std"&lag_fun=="Beta"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.8,z=0,m=0,theta=0,w2=2,shape=5)

#z_up_bound<- -0.7/ min(zoo::coredata(X) - mean(zoo::coredata(X)))

ui<-rbind(
c(1,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0,0,0),	     			## alpha+beta<1
c(0,0,0,0,0,1,0),				 	## w2>1.001
c(0,0,0,0,0,0,1))				## shape > 2.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-2.001)

LOGLIK<-GM_X_loglik_no_skew
cond_vol<-GM_X_cond_vol_no_skew
long_run_vol<-GM_X_long_run_vol_no_skew

} else if (model=="GMX"&skew=="NO"&distribution=="norm"&lag_fun=="Almon"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.8,z=0,m=0,theta=0,w2=-0.1)

ui<-rbind(
c(1,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0,0),	     		## alpha+beta<1
c(0,0,0,0,0,-1))				## w2<0

ci<-c(-0.0001,-0.001,0.999,0)

LOGLIK<-GM_X_loglik_no_skew
cond_vol<-GM_X_cond_vol_no_skew
long_run_vol<-GM_X_long_run_vol_no_skew

} else if (model=="GMX"&skew=="NO"&distribution=="std"&lag_fun=="Almon"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.8,z=0,m=0,theta=0,w2=-0.1,shape=5)

ui<-rbind(
c(1,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0,0,0),	     			## alpha+beta<1
c(0,0,0,0,0,-1,0),				## w2<0
c(0,0,0,0,0,0,1))					## shape > 2.001

ci<-c(-0.0001,-0.001,0.999,0,-2.001)

LOGLIK<-GM_X_loglik_no_skew
cond_vol<-GM_X_cond_vol_no_skew
long_run_vol<-GM_X_long_run_vol_no_skew

} else if (model=="DAGM"&skew=="YES"&distribution=="norm"&lag_fun=="Beta"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=8)
colnames(begin_val)<-c("alpha","beta","gamma","m","theta_pos","w2_pos","theta_neg","w2_neg")
begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-stats::runif(R,min=-1,max=1)
begin_val[,6]<-1.01
begin_val[,7]<-stats::runif(R,min=-1,max=1)
begin_val[,8]<-1.01

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(DAGM_loglik(begin_val[i,],
r_t_in_s,mv_m_in_s,K=K,distribution="norm",lag_fun="Beta"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0),        		 		## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0),			 	## alpha+beta+gamma/2<1
c(0,0,0,0,0,1,0,0),				 	## w2_pos>1.001
c(0,0,0,0,0,0,0,1))				 	## w2_neg>1.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-1.001)

LOGLIK<-DAGM_loglik
cond_vol<-DAGM_cond_vol
long_run_vol<-DAGM_long_run_vol

} else if (model=="DAGM"&skew=="YES"&distribution=="std"&lag_fun=="Beta"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=9)
colnames(begin_val)<-c("alpha","beta","gamma","m","theta_pos","w2_pos",
"theta_neg","w2_neg","shape")
begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-stats::runif(R,min=-1,max=1)
begin_val[,6]<-1.01
begin_val[,7]<-stats::runif(R,min=-1,max=1)
begin_val[,8]<-1.01
begin_val[,9]<-stats::runif(R,min=2.001,max=10)

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(DAGM_loglik(begin_val[i,],
r_t_in_s,mv_m_in_s,K=K,distribution="std",lag_fun="Beta"))
}

start_val<-begin_val[which.max(which_row),]


ui<-rbind(
c(1,0,0,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0,0),			 	## alpha+beta+gamma/2<1
c(0,0,0,0,0,1,0,0,0),				 	## w2_pos>1.001
c(0,0,0,0,0,0,0,1,0),				 	## w2_neg>1.001
c(0,0,0,0,0,0,0,0,1))					## shape>2.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-1.001,-2.001)

LOGLIK<-DAGM_loglik
cond_vol<-DAGM_cond_vol
long_run_vol<-DAGM_long_run_vol

} else if (model=="DAGM"&skew=="YES"&distribution=="norm"&lag_fun=="Almon"){


start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=8)
colnames(begin_val)<-c("alpha","beta","gamma","m","theta_pos","w2_pos","theta_neg","w2_neg")
begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-stats::runif(R,min=-1,max=1)
begin_val[,6]<- -0.1
begin_val[,7]<-stats::runif(R,min=-1,max=1)
begin_val[,8]<- -0.1

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(DAGM_loglik(begin_val[i,],
r_t_in_s,mv_m_in_s,K=K,distribution="norm",lag_fun="Almon"))
}

start_val<-begin_val[which.max(which_row),]



ui<-rbind(
c(1,0,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0),        		 		## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0),			 	## alpha+beta+gamma/2<1
c(0,0,0,0,0,-1,0,0),				 	## w2_pos<0
c(0,0,0,0,0,0,0,-1))				 	## w2_neg<0

ci<-c(-0.0001,-0.001,0.999,0,0)

LOGLIK<-DAGM_loglik
cond_vol<-DAGM_cond_vol
long_run_vol<-DAGM_long_run_vol

} else if (model=="DAGM"&skew=="YES"&distribution=="std"&lag_fun=="Almon"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=9)
colnames(begin_val)<-c("alpha","beta","gamma","m","theta_pos","w2_pos",
"theta_neg","w2_neg","shape")
begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-stats::runif(R,min=-1,max=1)
begin_val[,6]<- -0.1
begin_val[,7]<-stats::runif(R,min=-1,max=1)
begin_val[,8]<- -0.1
begin_val[,9]<-stats::runif(R,min=2.001,max=10)

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(DAGM_loglik(begin_val[i,],
r_t_in_s,mv_m_in_s,K=K,distribution="std",lag_fun="Almon"))
}

start_val<-begin_val[which.max(which_row),]


ui<-rbind(
c(1,0,0,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0,0),			## alpha+beta+gamma/2<1
c(0,0,0,0,0,-1,0,0,0),				## w2_pos<0
c(0,0,0,0,0,0,0,-1,0),				## w2_neg<0
c(0,0,0,0,0,0,0,0,1))				## shape>2.001

ci<-c(-0.0001,-0.001,0.999,0,0,-2.001)

LOGLIK<-DAGM_loglik
cond_vol<-DAGM_cond_vol
long_run_vol<-DAGM_long_run_vol

} else if (model=="DAGMX"&skew=="YES"&distribution=="norm"&lag_fun=="Beta"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.90,gamma=0,z=0,m=0,theta_pos=0, 
w2_pos=2,theta_neg=0,w2_neg=2)


ui<-rbind(
c(1,0,0,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0),        		 		## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0,0),			 	## alpha+beta+gamma/2<1
c(0,0,0,0,0,0,1,0,0),				 	## w2_pos>1.001
c(0,0,0,0,0,0,0,0,1))				 	## w2_neg>1.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-1.001)

LOGLIK<-DAGM_X_loglik
cond_vol<-DAGM_X_cond_vol
long_run_vol<-DAGM_X_long_run_vol

} else if (model=="DAGMX"&skew=="YES"&distribution=="std"&lag_fun=="Beta"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.90,gamma=0,z=0,m=0,theta_pos=0, 
w2_pos=2,theta_neg=0,w2_neg=2,shape=5)


ui<-rbind(
c(1,0,0,0,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0,0,0),			 	## alpha+beta+gamma/2<1
c(0,0,0,0,0,0,1,0,0,0),				 	## w2_pos>1.001
c(0,0,0,0,0,0,0,0,1,0),				 	## w2_neg>1.001
c(0,0,0,0,0,0,0,0,0,1))					## shape>2.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-1.001,-2.001)

LOGLIK<-DAGM_X_loglik
cond_vol<-DAGM_X_cond_vol
long_run_vol<-DAGM_X_long_run_vol

} else if (model=="DAGMX"&skew=="YES"&distribution=="norm"&lag_fun=="Almon"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.90,gamma=0,z=0,m=0,theta_pos=0, 
w2_pos=-0.1,theta_neg=0,w2_neg=-0.1)

ui<-rbind(
c(1,0,0,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0),        		 		## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0,0),			 	## alpha+beta+gamma/2<1
c(0,0,0,0,0,0,-1,0,0),				 	## w2_pos<0
c(0,0,0,0,0,0,0,0,-1))				 	## w2_neg<0

ci<-c(-0.0001,-0.001,0.999,0,0)

LOGLIK<-DAGM_X_loglik
cond_vol<-DAGM_X_cond_vol
long_run_vol<-DAGM_X_long_run_vol

} else if (model=="DAGMX"&skew=="YES"&distribution=="std"&lag_fun=="Almon"){

start_val<-ui<-ci<-NULL


start_val<-c(alpha=0.01,beta=0.90,gamma=0,z=0,m=0,theta_pos=0, 
w2_pos=-0.1,theta_neg=0,w2_neg=-0.1,shape=5)

ui<-rbind(
c(1,0,0,0,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0,0,0),			## alpha+beta+gamma/2<1
c(0,0,0,0,0,0,-1,0,0,0),				## w2_pos<0
c(0,0,0,0,0,0,0,0,-1,0),				## w2_neg<0
c(0,0,0,0,0,0,0,0,0,1))				## shape>2.001

ci<-c(-0.0001,-0.001,0.999,0,0,-2.001)

LOGLIK<-DAGM_X_loglik
cond_vol<-DAGM_X_cond_vol
long_run_vol<-DAGM_X_long_run_vol

} else if (model=="DAGM"&skew=="NO"&distribution=="norm"&lag_fun=="Beta"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.90,m=0,theta_pos=0, 
w2_pos=2,theta_neg=0,w2_neg=2)

ui<-rbind(
c(1,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0),        		 		## beta>0.001
c(-1,-1,0,0,0,0,0),			 		## alpha+beta<1
c(0,0,0,0,1,0,0),				 		## w2_pos>1.001
c(0,0,0,0,0,0,1))				 		## w2_neg>1.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-1.001)

LOGLIK<-DAGM_loglik_no_skew
cond_vol<-DAGM_cond_vol_no_skew
long_run_vol<-DAGM_long_run_vol_no_skew

} else if (model=="DAGM"&skew=="NO"&distribution=="std"&lag_fun=="Beta"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.90,m=0,theta_pos=0, 
w2_pos=2,theta_neg=0,w2_neg=2,shape=5)

ui<-rbind(
c(1,0,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0,0,0,0),			 		## alpha+beta<1
c(0,0,0,0,1,0,0,0),				 	## w2_pos>1.001
c(0,0,0,0,0,0,1,0),				 	## w2_neg>1.001
c(0,0,0,0,0,0,0,1))					## shape>2.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-1.001,-2.001)

LOGLIK<-DAGM_loglik_no_skew
cond_vol<-DAGM_cond_vol_no_skew
long_run_vol<-DAGM_long_run_vol_no_skew

} else if (model=="DAGM"&skew=="NO"&distribution=="norm"&lag_fun=="Almon"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.90,m=0,theta_pos=0, 
w2_pos=-0.1,theta_neg=0,w2_neg=-0.1)

ui<-rbind(
c(1,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0),        		 		## beta>0.001
c(-1,-1,0,0,0,0,0),			 		## alpha+beta<1
c(0,0,0,0,-1,0,0),				 		## w2_pos<0
c(0,0,0,0,0,0,-1))				 		## w2_neg<0

ci<-c(-0.0001,-0.001,0.999,0,0)

LOGLIK<-DAGM_loglik_no_skew
cond_vol<-DAGM_cond_vol_no_skew
long_run_vol<-DAGM_long_run_vol_no_skew

} else if (model=="DAGM"&skew=="NO"&distribution=="std"&lag_fun=="Almon"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.90,m=0,theta_pos=0, 
w2_pos=-0.1,theta_neg=0,w2_neg=-0.1,shape=5)

ui<-rbind(
c(1,0,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0,0,0,0),			 		## alpha+beta<1
c(0,0,0,0,-1,0,0,0),				 	## w2_pos<0
c(0,0,0,0,0,0,-1,0),				 	## w2_neg<0
c(0,0,0,0,0,0,0,1))					## shape>2.001

ci<-c(-0.0001,-0.001,0.999,0,0,-2.001)

LOGLIK<-DAGM_loglik_no_skew
cond_vol<-DAGM_cond_vol_no_skew
long_run_vol<-DAGM_long_run_vol_no_skew

} else if (model=="DAGMX"&skew=="NO"&distribution=="norm"&lag_fun=="Beta"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.90,z=0,m=0,theta_pos=0, 
w2_pos=2,theta_neg=0,w2_neg=2)

ui<-rbind(
c(1,0,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0),        		 		## beta>0.001
c(-1,-1,0,0,0,0,0,0),			 	## alpha+beta<1
c(0,0,0,0,0,1,0,0),				 	## w2_pos>1.001
c(0,0,0,0,0,0,0,1))				 	## w2_neg>1.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-1.001)

LOGLIK<-DAGM_X_loglik_no_skew
cond_vol<-DAGM_X_cond_vol_no_skew
long_run_vol<-DAGM_X_long_run_vol_no_skew

} else if (model=="DAGMX"&skew=="NO"&distribution=="std"&lag_fun=="Beta"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.90,z=0,m=0,theta_pos=0, 
w2_pos=2,theta_neg=0,w2_neg=2,shape=5)

ui<-rbind(
c(1,0,0,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0,0,0,0,0),			 	## alpha+beta<1
c(0,0,0,0,0,1,0,0,0),				 	## w2_pos>1.001
c(0,0,0,0,0,0,0,1,0),				 	## w2_neg>1.001
c(0,0,0,0,0,0,0,0,1))					## shape>2.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-1.001,-2.001)

LOGLIK<-DAGM_X_loglik_no_skew
cond_vol<-DAGM_X_cond_vol_no_skew
long_run_vol<-DAGM_X_long_run_vol_no_skew

} else if (model=="DAGMX"&skew=="NO"&distribution=="norm"&lag_fun=="Almon"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.90,z=0,m=0,theta_pos=0, 
w2_pos=-0.1,theta_neg=0,w2_neg=-0.1)

ui<-rbind(
c(1,0,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0),        		 		## beta>0.001
c(-1,-1,0,0,0,0,0,0),			 	## alpha+beta<1
c(0,0,0,0,0,-1,0,0),				 	## w2_pos<0
c(0,0,0,0,0,0,0,-1))				 	## w2_neg<0

ci<-c(-0.0001,-0.001,0.999,0,0)

LOGLIK<-DAGM_X_loglik_no_skew
cond_vol<-DAGM_X_cond_vol_no_skew
long_run_vol<-DAGM_X_long_run_vol_no_skew

} else if (model=="DAGMX"&skew=="NO"&distribution=="std"&lag_fun=="Almon"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.90,z=0.1,m=0,theta_pos=0, 
w2_pos=-0.1,theta_neg=0,w2_neg=-0.1,shape=5)

ui<-rbind(
c(1,0,0,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0,0,0,0,0),			## alpha+beta<1
c(0,0,0,0,0,-1,0,0,0),				## w2_pos<0
c(0,0,0,0,0,0,0,-1,0),				## w2_neg<0
c(0,0,0,0,0,0,0,0,1))				## shape>2.001


ci<-c(-0.0001,-0.001,0.999,0,0,-2.001)

LOGLIK<-DAGM_X_loglik_no_skew
cond_vol<-DAGM_X_cond_vol_no_skew
long_run_vol<-DAGM_X_long_run_vol_no_skew

} else if (model=="DAGM2M"&skew=="YES"&distribution=="norm"&lag_fun=="Beta"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=12)
colnames(begin_val)<-c("alpha","beta","gamma","m","theta_pos_1","w2_pos_1","theta_neg_1","w2_neg_1",
"theta_pos_2","w2_pos_2","theta_neg_2","w2_neg_2"
)
begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-stats::runif(R,min=-1,max=1)
begin_val[,6]<-1.01
begin_val[,7]<-stats::runif(R,min=-1,max=1)
begin_val[,8]<-1.01
begin_val[,9]<-stats::runif(R,min=-1,max=1)
begin_val[,10]<-1.01
begin_val[,11]<-stats::runif(R,min=-1,max=1)
begin_val[,12]<-1.01

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(DAGM_2M_loglik(begin_val[i,],
r_t_in_s,mv_m_1=mv_m_in_s,mv_m_2=mv_m_2_in_s,K_1=K,K_2=K_2,distribution="norm",lag_fun="Beta"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0,0,0,0),        		 		## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0,0,0,0,0),			 	## alpha+beta+gamma/2<1
c(0,0,0,0,0,1,0,0,0,0,0,0),				 	## w2_pos_1>1.001
c(0,0,0,0,0,0,0,1,0,0,0,0),				 	## w2_neg_1>1.001
c(0,0,0,0,0,0,0,0,0,1,0,0),				 	## w2_pos_2>1.001
c(0,0,0,0,0,0,0,0,0,0,0,1))				 	## w2_neg_2>1.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-1.001,-1.001,-1.001)

LOGLIK<-DAGM_2M_loglik
cond_vol<-DAGM_2M_cond_vol
long_run_vol<-DAGM_2M_long_run

} else if (model=="DAGM2M"&skew=="YES"&distribution=="std"&lag_fun=="Beta"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=13)
colnames(begin_val)<-c("alpha","beta","gamma","m","theta_pos_1","w2_pos_1","theta_neg_1","w2_neg_1",
"theta_pos_2","w2_pos_2","theta_neg_2","w2_neg_2","shape")

begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-stats::runif(R,min=-1,max=1)
begin_val[,6]<-1.01
begin_val[,7]<-stats::runif(R,min=-1,max=1)
begin_val[,8]<-1.01
begin_val[,9]<-stats::runif(R,min=-1,max=1)
begin_val[,10]<-1.01
begin_val[,11]<-stats::runif(R,min=-1,max=1)
begin_val[,12]<-1.01
begin_val[,13]<-stats::runif(R,min=2.001,max=10)

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(DAGM_2M_loglik(begin_val[i,],
r_t_in_s,mv_m_1=mv_m_in_s,mv_m_2=mv_m_2_in_s,K_1=K,K_2=K_2,distribution="std",lag_fun="Beta"))
}

start_val<-begin_val[which.max(which_row),]


ui<-rbind(
c(1,0,0,0,0,0,0,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0,0,0,0,0,0),			 	## alpha+beta+gamma/2<1
c(0,0,0,0,0,1,0,0,0,0,0,0,0),				 	## w2_pos_1>1.001
c(0,0,0,0,0,0,0,1,0,0,0,0,0),				 	## w2_neg_1>1.001
c(0,0,0,0,0,1,0,0,0,1,0,0,0),				 	## w2_pos_2>1.001
c(0,0,0,0,0,0,0,1,0,0,0,1,0),				 	## w2_neg_2>1.001
c(0,0,0,0,0,0,0,0,0,0,0,0,1))					## shape>2.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-1.001,-1.001,-1.001,-2.001)

LOGLIK<-DAGM_2M_loglik
cond_vol<-DAGM_2M_cond_vol
long_run_vol<-DAGM_2M_long_run


} else if (model=="DAGM2M"&skew=="YES"&distribution=="norm"&lag_fun=="Almon"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=12)
colnames(begin_val)<-c("alpha","beta","gamma","m","theta_pos_1","w2_pos_1","theta_neg_1","w2_neg_1",
"theta_pos_2","w2_pos_2","theta_neg_2","w2_neg_2")

begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-stats::runif(R,min=-1,max=1)
begin_val[,6]<- -0.1
begin_val[,7]<-stats::runif(R,min=-1,max=1)
begin_val[,8]<- -0.1
begin_val[,9]<-stats::runif(R,min=-1,max=1)
begin_val[,10]<- -0.1
begin_val[,11]<-stats::runif(R,min=-1,max=1)
begin_val[,12]<- -0.1
which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(DAGM_2M_loglik(begin_val[i,],
r_t_in_s,mv_m_1=mv_m_in_s,mv_m_2=mv_m_2_in_s,K_1=K,K_2=K_2,distribution="norm",lag_fun="Almon"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0,0,0,0),        		 		## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0,0,0,0,0),			 	## alpha+beta+gamma/2<1
c(0,0,0,0,0,-1,0,0,0,0,0,0),				 	## w2_pos_1<0
c(0,0,0,0,0,0,0,-1,0,0,0,0),				 	## w2_neg_1<0
c(0,0,0,0,0,0,0,0,0,-1,0,0),				 	## w2_pos_2<0
c(0,0,0,0,0,0,0,0,0,0,0,-1))				 	## w2_neg_2<0

ci<-c(-0.0001,-0.001,0.999,0,0,0,0)

LOGLIK<-DAGM_2M_loglik
cond_vol<-DAGM_2M_cond_vol
long_run_vol<-DAGM_2M_long_run

} else if (model=="DAGM2M"&skew=="YES"&distribution=="std"&lag_fun=="Almon"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=13)
colnames(begin_val)<-c("alpha","beta","gamma","m","theta_pos_1","w2_pos_1","theta_neg_1","w2_neg_1",
                       "theta_pos_2","w2_pos_2","theta_neg_2","w2_neg_2","shape")


begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-stats::runif(R,min=-1,max=1)
begin_val[,6]<- -0.1
begin_val[,7]<-stats::runif(R,min=-1,max=1)
begin_val[,8]<- -0.1
begin_val[,9]<-stats::runif(R,min=-1,max=1)
begin_val[,10]<- -0.1
begin_val[,11]<-stats::runif(R,min=-1,max=1)
begin_val[,12]<- -0.1
begin_val[,13]<-stats::runif(R,min=2.001,max=10)

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(DAGM_2M_loglik(begin_val[i,],
r_t_in_s,mv_m_1=mv_m_in_s,mv_m_2=mv_m_2_in_s,K_1=K,K_2=K_2,distribution="std",lag_fun="Almon"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0,0,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0,0,0,0,0),        		 		## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0,0,0,0,0,0),			 	## alpha+beta+gamma/2<1
c(0,0,0,0,0,-1,0,0,0,0,0,0,0),				 	## w2_pos_1<0
c(0,0,0,0,0,0,0,-1,0,0,0,0,0),				 	## w2_neg_1<0
c(0,0,0,0,0,0,0,0,0,-1,0,0,0),				 	## w2_pos_2<0
c(0,0,0,0,0,0,0,0,0,0,0,-1,0),				 	## w2_neg_2<0
c(0,0,0,0,0,0,0,0,0,0,0,0,1))				## shape>2.001

ci<-c(-0.0001,-0.001,0.999,0,0,0,0,-2.001)

LOGLIK<-DAGM_2M_loglik
cond_vol<-DAGM_2M_cond_vol
long_run_vol<-DAGM_2M_long_run

} else if (model=="DAGM2M"&skew=="NO"&distribution=="norm"&lag_fun=="Beta"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=11)
colnames(begin_val)<-c("alpha","beta","m","theta_pos_1","w2_pos_1","theta_neg_1","w2_neg_1",
"theta_pos_2","w2_pos_2","theta_neg_2","w2_neg_2")

begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=-1,max=1)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-1.01
begin_val[,6]<-stats::runif(R,min=-1,max=1)
begin_val[,7]<-1.01
begin_val[,8]<-stats::runif(R,min=-1,max=1)
begin_val[,9]<-1.01
begin_val[,10]<-stats::runif(R,min=-1,max=1)
begin_val[,11]<-1.01

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(DAGM_2M_loglik_no_skew(begin_val[i,],
r_t_in_s,mv_m_1=mv_m_in_s,mv_m_2=mv_m_2_in_s,K_1=K,K_2=K_2,distribution="norm",lag_fun="Beta"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0,0,0,0,0,0,0),			 	## alpha+beta<1
c(0,0,0,0,1,0,0,0,0,0,0),				 	## w2_pos_1>1.001
c(0,0,0,0,0,0,1,0,0,0,0),				 	## w2_neg_1>1.001
c(0,0,0,0,0,0,0,0,1,0,0),				 	## w2_pos_2>1.001
c(0,0,0,0,0,0,0,0,0,0,1))				 	## w2_neg_2>1.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-1.001,-1.001,-1.001)

LOGLIK<-DAGM_2M_loglik_no_skew
cond_vol<-DAGM_2M_cond_vol_no_skew
long_run_vol<-DAGM_2M_long_run_no_skew


} else if (model=="DAGM2M"&skew=="NO"&distribution=="std"&lag_fun=="Beta"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=12)
colnames(begin_val)<-c("alpha","beta","m","theta_pos_1","w2_pos_1","theta_neg_1","w2_neg_1",
"theta_pos_2","w2_pos_2","theta_neg_2","w2_neg_2","shape")

begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=-1,max=1)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-1.01
begin_val[,6]<-stats::runif(R,min=-1,max=1)
begin_val[,7]<-1.01
begin_val[,8]<-stats::runif(R,min=-1,max=1)
begin_val[,9]<-1.01
begin_val[,10]<-stats::runif(R,min=-1,max=1)
begin_val[,11]<-1.01
begin_val[,12]<-stats::runif(R,min=2.001,max=10)

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(DAGM_2M_loglik_no_skew(begin_val[i,],
r_t_in_s,mv_m_1=mv_m_in_s,mv_m_2=mv_m_2_in_s,K_1=K,K_2=K_2,distribution="std",lag_fun="Beta"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0,0,0,0,0,0,0),       	 	 		## alpha>0.0001
c(0,1,0,0,0,0,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0,0,0,0,0,0,0,0),			 		## alpha+beta<1
c(0,0,0,0,1,0,0,0,0,0,0,0),				 	## w2_pos_1>1.001
c(0,0,0,0,0,0,1,0,0,0,0,0),				 	## w2_neg_1>1.001
c(0,0,0,0,0,0,0,0,1,0,0,0),				 	## w2_pos_2>1.001
c(0,0,0,0,0,0,0,0,0,0,1,0),				 	## w2_neg_2>1.001
c(0,0,0,0,0,0,0,0,0,0,0,1))				 	## shape>2

ci<-c(-0.0001,-0.001,0.999,-1.001,-1.001,-1.001,-1.001,-2.001)

LOGLIK<-DAGM_2M_loglik_no_skew
cond_vol<-DAGM_2M_cond_vol_no_skew
long_run_vol<-DAGM_2M_long_run_no_skew


} else if (model=="DAGM2M"&skew=="NO"&distribution=="norm"&lag_fun=="Almon"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=11)
colnames(begin_val)<-c("alpha","beta","m","theta_pos_1","w2_pos_1","theta_neg_1","w2_neg_1",
"theta_pos_2","w2_pos_2","theta_neg_2","w2_neg_2")

begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=-1,max=1)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<- -0.1
begin_val[,6]<-stats::runif(R,min=-1,max=1)
begin_val[,7]<- -0.1
begin_val[,8]<-stats::runif(R,min=-1,max=1)
begin_val[,9]<- -0.1
begin_val[,10]<-stats::runif(R,min=-1,max=1)
begin_val[,11]<- -0.1

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(DAGM_2M_loglik_no_skew(begin_val[i,],
r_t_in_s,mv_m_1=mv_m_in_s,mv_m_2=mv_m_2_in_s,K_1=K,K_2=K_2,distribution="norm",lag_fun="Almon"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0,0,0,0,0,0,0),			 	## alpha+beta<1
c(0,0,0,0,-1,0,0,0,0,0,0),				 	## w2_pos_1<0
c(0,0,0,0,0,0,-1,0,0,0,0),				 	## w2_neg_1<0
c(0,0,0,0,0,0,0,0,-1,0,0),				 	## w2_pos_2<0
c(0,0,0,0,0,0,0,0,0,0,-1))				 	## w2_neg_2<0

ci<-c(-0.0001,-0.001,0.999,0,0,0,0)

LOGLIK<-DAGM_2M_loglik_no_skew
cond_vol<-DAGM_2M_cond_vol_no_skew
long_run_vol<-DAGM_2M_long_run_no_skew

} else if (model=="DAGM2M"&skew=="NO"&distribution=="std"&lag_fun=="Almon"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=12)
colnames(begin_val)<-c("alpha","beta","m","theta_pos_1","w2_pos_1","theta_neg_1","w2_neg_1",
"theta_pos_2","w2_pos_2","theta_neg_2","w2_neg_2","shape")

begin_val[,1]<-stats::runif(R,min=0.001,max=0.095)
begin_val[,2]<-stats::runif(R,min=0.6,max=0.8)
begin_val[,3]<-stats::runif(R,min=-1,max=1)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<- -1.01
begin_val[,6]<-stats::runif(R,min=-1,max=1)
begin_val[,7]<- -1.01
begin_val[,8]<-stats::runif(R,min=-1,max=1)
begin_val[,9]<- -1.01
begin_val[,10]<-stats::runif(R,min=-1,max=1)
begin_val[,11]<- -1.01
begin_val[,12]<-stats::runif(R,min=2.001,max=10)

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(DAGM_2M_loglik_no_skew(begin_val[i,],
r_t_in_s,mv_m_1=mv_m_in_s,mv_m_2=mv_m_2_in_s,K_1=K,K_2=K_2,distribution="std",lag_fun="Almon"))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0,0,0,0,0,0,0),       	 	 		## alpha>0.0001
c(0,1,0,0,0,0,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0,0,0,0,0,0,0,0),			 		## alpha+beta<1
c(0,0,0,0,-1,0,0,0,0,0,0,0),				 	## w2_pos_1<0
c(0,0,0,0,0,0,-1,0,0,0,0,0),				 	## w2_neg_1<0
c(0,0,0,0,0,0,0,0,-1,0,0,0),				 	## w2_pos_2<0
c(0,0,0,0,0,0,0,0,0,0,-1,0),				 	## w2_neg_2<0
c(0,0,0,0,0,0,0,0,0,0,0,1))				 	## shape>2

ci<-c(-0.0001,-0.001,0.999,0,0,0,0,-2.001)

LOGLIK<-DAGM_2M_loglik_no_skew
cond_vol<-DAGM_2M_cond_vol_no_skew
long_run_vol<-DAGM_2M_long_run_no_skew


}

####################################### begin estimation

r_t_in_s_est<-zoo::coredata(r_t_in_s)
if(!missing(X)){ 
X_in_s_est<-zoo::coredata(X_in_s) - mean(zoo::coredata(X_in_s))
} else {
X_in_s_est<-X_in_s
}

if(!missing(X)){ 
est<-suppressWarnings(maxLik(
logLik=LOGLIK,
start=start_val,
daily_ret=r_t_in_s_est,
X=X_in_s_est,
mv_m=mv_m_in_s,
K=K,
lag_fun=lag_fun,
distribution=distribution,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))
} else if (model=="DAGM2M") {
est<-suppressWarnings(maxLik(
logLik=LOGLIK,
start=start_val,
daily_ret=r_t_in_s_est,
mv_m_1=mv_m_in_s,
mv_m_2=mv_m_2_in_s,
K_1=K,
K_2=K_2,
lag_fun=lag_fun,
distribution=distribution,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))
} else {
est<-suppressWarnings(maxLik(
logLik=LOGLIK,
start=start_val,
daily_ret=r_t_in_s_est,
mv_m=mv_m_in_s,
K=K,
lag_fun=lag_fun,
distribution=distribution,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))
}

N_coef<-length(stats::coef(est))

mat_coef<-data.frame(rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef))
colnames(mat_coef)<-c("Estimate","Std. Error","t value","Pr(>|t|)")

rownames(mat_coef)<-names(stats::coef(est))

mat_coef[,1]<-round(stats::coef(est),6)
mat_coef[,2]<-round(QMLE_sd(est),6)
mat_coef[,3]<-round(stats::coef(est)/QMLE_sd(est),6)
mat_coef[,4]<-round(apply(rbind(stats::coef(est)/QMLE_sd(est)),1,function(x) 2*(1-stats::pnorm(abs(x)))),6)

## change the order of variables
if (skew=="YES") {

mat_coef2<-mat_coef

mat_coef2[2,]<-mat_coef[3,] #gamma
rownames(mat_coef2)[3]<-"temp"
rownames(mat_coef2)[2]<-"gamma"
mat_coef2[3,]<-mat_coef[2,] #beta
rownames(mat_coef2)[3]<-"beta"
mat_coef<-mat_coef2
}

######## in-sample and out-of-sample estimation and evaluation

est_coef<-stats::coef(est)

if (distribution=="std"){

est_coef<-est_coef[-N_coef]

}

if(!missing(X)){ 
vol_est<-(cond_vol(est_coef,r_t_in_s,X_in_s,mv_m_in_s,K=K,lag_fun=lag_fun))^2
lr_vol<-long_run_vol(est_coef,r_t_in_s,X_in_s,mv_m_in_s,K=K,lag_fun=lag_fun)
} else if (model=="DAGM2M") {
vol_est<-(cond_vol(est_coef,r_t_in_s,mv_m_1=mv_m_in_s,mv_m_2=mv_m_2_in_s,
K_1=K,K_2=K_2,lag_fun=lag_fun))^2
lr_vol<-long_run_vol(est_coef,r_t_in_s,mv_m_1=mv_m_in_s,mv_m_2=mv_m_2_in_s,
K_1=K,K_2=K_2,lag_fun=lag_fun)
} else {
vol_est<-(cond_vol(est_coef,r_t_in_s,mv_m_in_s,K=K,lag_fun=lag_fun))^2
lr_vol<-long_run_vol(est_coef,r_t_in_s,mv_m_in_s,K=K,lag_fun=lag_fun)
}


if (missing(out_of_sample)){

N_est<-length(vol_est)

vol_est_lf<-zoo::coredata(vol_est)
vol_proxy_lf<-zoo::coredata(vol_proxy_in_s)

res<-list(
model=model,
rob_coef_mat=mat_coef,
obs=N,
period=range(stats::time(r_t_in_s)),
loglik=as.numeric(stats::logLik(est)),
inf_criteria=Inf_criteria(est),
loss_in_s=LF_f(vol_est_lf,vol_proxy_lf),
est_vol_in_s=vol_est^0.5,
est_lr_in_s=lr_vol)

} else {

r_t_oos<-daily_ret[(N-out_of_sample+1):N]
mv_m_oos<-mv_m[,(N-out_of_sample+1):(N)]

if (model=="DAGM2M"){
mv_m_oos_2<-mv_m_2[,(N-out_of_sample+1):(N)]
}

if(missing(vol_proxy)){
vol_proxy_oos<-r_t_oos^2
} else {
vol_proxy_oos<-vol_proxy[(N-out_of_sample+1):N]
}

N_est<-length(vol_est)

vol_est_lf<-zoo::coredata(vol_est)
vol_proxy_lf<-zoo::coredata(vol_proxy_in_s)
vol_proxy_oos_lf<-zoo::coredata(vol_proxy_oos)

if(!missing(X)){ 
X_oos<-X-mean(X)
X_oos_f<-X_oos[(N-out_of_sample+1):N]
vol_est_oos<-(cond_vol(est_coef,r_t_oos,X_oos_f,mv_m_oos,K=K,lag_fun=lag_fun))^2
LR_oos<-(long_run_vol(est_coef,r_t_oos,X_oos_f,mv_m_oos,K=K,lag_fun=lag_fun))
} else if (model=="DAGM2M") {
vol_est_oos<-(cond_vol(est_coef,r_t_oos,mv_m_1=mv_m_oos,mv_m_2=mv_m_oos_2,
K_1=K,K_2=K_2,lag_fun=lag_fun))^2
LR_oos<-(long_run_vol(est_coef,r_t_oos,mv_m_1=mv_m_oos,mv_m_2=mv_m_oos_2,
K_1=K,K_2=K_2,lag_fun=lag_fun))
} else {
vol_est_oos<-(cond_vol(est_coef,r_t_oos,mv_m_oos,K=K,lag_fun=lag_fun))^2
LR_oos<-(long_run_vol(est_coef,r_t_oos,mv_m_oos,K=K,lag_fun=lag_fun))
}

vol_est_oos_lf<-zoo::coredata(vol_est_oos)

res<-list(
model=model,
rob_coef_mat=mat_coef,
obs=length(r_t_in_s),
period=range(time(r_t_in_s)),
loglik=as.numeric(stats::logLik(est)),
inf_criteria=Inf_criteria(est),
loss_in_s=LF_f(vol_est_lf,vol_proxy_lf),
est_vol_in_s=vol_est^0.5,
est_lr_in_s=lr_vol,
loss_oos=LF_f(vol_est_oos_lf,vol_proxy_oos_lf),
est_vol_oos=vol_est_oos^0.5,
est_lr_oos=LR_oos)

}

class(res)<-c("rumidas")
return(res)
print.rumidas(res)


}


#' Methods for obtaining (and evaluating) a variety of MEM(-MIDAS)-based models
#'
#' Estimates several MEM and MEM-MIDAS-based models. For details, see \insertCite{amendola2024doubly;textual}{rumidas}
#' @param model Model to estimate. Valid choices are: "MEMMIDAS" for MEM-MIDAS, "MEM" for base MEM, "MEMX" for the base MEM with an X term, 
#' "MEMMIDASX" for the MEM-MIDAS-X
#' @param skew The skewness parameter linked to lagged daily returns. Valid choices are: "YES" and "NO"
#' @param x Dependent variable to predict. Usually the realized volatility. It must be positive and "xts" object
#' @param daily_ret **optional**. Daily returns, which must be an "xts" object. NULL by default
#' @param mv_m **optional**. MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function. NULL by default
#' @param K **optional**. Number of (lagged) realizations of the MIDAS variable to consider. NULL by default
#' for the Beta and Exponential Almon lag functions, respectively
#' @param z  **optional**. Additional daily variable, which must be an "xts" object, and with the same length of x. NULL by default
#' @param out_of_sample **optional**. A positive integer indicating the number of periods before the last to keep for out of sample forecasting
#' @param R **optional**.  A positive integer indicating the number of replications used to find the best starting values. 
#' Equal to 100 by default
#' @return \code{umemfit} returns an object of class 'rumidas'. The function \code{\link{summary.rumidas}} 
#' can be used to print a summary of the results. Moreover, an object of class 'rumidas' is a list containing the following components:
#' \itemize{
#' 	\item model: The model used for the estimation.
#'   \item rob_coef_mat: The matrix of estimated coefficients, with the QML standard errors. 
#'   \item obs: The number of daily observations used for the (in-sample) estimation.
#'   \item period: The period of the in-sample estimation.
#'   \item loglik: The value of the log-likelihood at the maximum.
#'	\item inf_criteria: The AIC and BIC information criteria.
#'	\item loss_in_s: The in-sample MSE and QLIKE averages, calculated considering the distance with respect to the dependent variable.
#'	\item est_in_s: The in-sample predicted dependent variable.
#'	\item est_lr_in_s: The in-sample predicted long-run component (if present) of the dependent variable.
#'	\item loss_oos: The out-of-sample MSE and QLIKE averages, calculated considering the distance with respect to the dependent variable.
#'	\item est_oos: The out-of-sample predicted dependent variable.
#'	\item est_lr_oos: The out-of-sample predicted long-run component (if present) of the dependent variable.
#' }
#' @importFrom Rdpack reprompt
#' @importFrom stats time
#' @importFrom zoo coredata
#' @import maxLik
#' @seealso \code{\link{mv_into_mat}}.
#'
#' @details
#' Function \code{umemfit} implements the estimation and evaluation of the MEM, MEM-MIDAS MEM-X and MEM-MIDAS-X models, 
#' with and without the asymmetric term linked to negative lagged daily returns. The general framework assumes that: 
#' \deqn{x_{i,t}= \mu_{i,t}\epsilon_{i,t} = \tau_{t} \xi_{i,t} \epsilon_{i,t},}
#' where 
#' \itemize{
#' \item \eqn{x_{i,t}} is a time series coming from a non-negative discrete time process for the \eqn{i}-th day (\eqn{i = 1, \ldots, N_t}) 
#' of the period \eqn{t} (for example, a week, a month or a quarter; \eqn{t = 1 , \ldots, T});
#' \item \eqn{\tau_{t}} is the long-run component, determining the average level of the conditional mean, varying each period \eqn{t};
#' \item \eqn{\xi_{i,t}} is a factor centered around one, labelled as the short--run term, which plays the role of dumping or amplifying \eqn{\tau_{i,t}};
#' \item \eqn{\epsilon_{i,t}} is an \eqn{iid} error term which, conditionally on the information set, has a unit mean, an unknown variance, and a probability density function 
#' defined over a non-negative support.
#' }
#' The short--run component of the MEM-MIDAS-X is:
#' \deqn{\xi_{i,t}=(1-\alpha-\gamma/2-\beta) + \left(\alpha +  \gamma \cdot {I}_{\left(r_{i-1,t}  < 0 \right)}\right) \frac{x_{i-1,t}}{\tau_t} + \beta \xi_{i-1,t} + \delta \left(Z_{i-1,t}-E(Z)\right),}
#' where \eqn{I_{(\cdot)}} is an indicator function, \eqn{r_{i,t}} is the daily return of the day \eqn{i} of the period \eqn{t} and \eqn{Z} is
#' an additional X term (for instance, the VIX). When the X part is absent, then the parameter \eqn{\delta} cancels.
#' The long-run component of the MEM-MIDAS and MEM-MIDAS-X is:
#' \deqn{\tau_{t} = \exp \left\{ m + \theta \sum_{k=1}^K \delta_{k}(\omega) X_{t-k}\right\},}
#' where \eqn{X_{t}} is the MIDAS term. When the "skew" parameter is set to "NO", \eqn{\gamma} disappears.  
#' The MEM and MEM-X models do not have the long- and short-run components. Therefore, they directly evolve according to \eqn{\mu_{i,t}}.
#' When the "skew" and X parameters are present, the MEM-X is:
#' \deqn{\mu_{i,t}= \left(1-\alpha - \gamma / 2 - \beta   \right)\mu + (\alpha + \gamma I_{\left(r_{i-1,t}  < 0 \right)}) x_{i-1,t} +   \beta \mu_{i-1,t}+\delta \left(Z_{i-1,t}-E(Z)\right),}
#' where \eqn{\mu=E(x_{i,t})}. When the "skew" parameter is set to "NO", in the previous equation \eqn{\gamma} cancels.
#' Finally, when the additional X part is not present, then we have the MEM model, where \eqn{\delta} disappears.
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' # estimate the base MEM, without the asymmetric term linked to negative lagged returns
#' real<-(rv5['2003/2010'])^0.5		# realized volatility
#' fit<-umemfit(model="MEM",skew="NO",x=real)
#' fit
#' summary.rumidas(fit)
#' # to see the estimated coefficients with the QML standard errors:
#' fit$rob_coef_mat
#'
#' # All the other elements of fit are:
#' names(fit)
#'
#' # estimate the MEM-MIDAS, with the asymmetric term linked to negative lagged returns,
#' # leaving the last 200 observations for the out-of-sample analysis
#' r_t<-sp500['2003/2010']
#' real<-(rv5['2003/2010'])^0.5		# realized volatility
#' mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' fit_2<-umemfit(model="MEMMIDAS",skew="YES",x=real,daily_ret=r_t,mv_m=mv_m,K=12,out_of_sample=200)
#' fit_2
#' summary.rumidas(fit_2)
#' # to see the estimated coefficients with the QML standard errors:
#' fit_2$rob_coef_mat
#' }
#' @export

umemfit<-function(
model,
skew,
x,
daily_ret=NULL,
mv_m=NULL,
K=NULL,
z=NULL,
out_of_sample=NULL,
R=100
){

############################# check on valid choices

if((model != "MEM")&(model != "MEMMIDAS")&(model != "MEMX")&(model != "MEMMIDASX")) { stop(cat("#Warning:\n Valid choices for the parameter 'model' are currently 'MEM', 'MEMX', 'MEMMIDAS' and 'MEMMIDASX' \n"))}
if((skew != "YES")&(skew != "NO")) { stop(cat("#Warning:\n Valid choices for the parameter 'skew' are 'YES' and 'NO' \n"))}

################################### checks on x variable

N<-length(x)
cond_x<- class(x)[1]
if(all(x<0)) { stop(
cat("#Warning:\n Parameter 'x' must have all entries positive \n")
)}
if(cond_x != "xts") { stop(
cat("#Warning:\n Parameter 'x' must be an 'xts' object \n")
)}

################################### checks on z variable

if(!is.null(z)&length(z)!=N){ stop(
cat("#Warning:\n Parameter 'z' must have the same length of 'x' \n")
)}

################################### checks on the skew parameter

if(skew=="YES"&missing(daily_ret)) { stop(
cat("#Warning:\n Parameter 'daily_ret' is missing. Please provide it \n")
)}

################################### checks on the MEM-MIDAS setting

#### if the model to estimate is the MEM-MIDAS or MEM-MIDAS-X with skew parameter
if ((model=="MEMMIDAS"|model=="MEMMIDASX")&skew=="YES") { 

N_r_t<-length(daily_ret)

cond_r_t<- class(daily_ret)[1]

cond_mv_m<- class(mv_m)[1]

if(cond_r_t != "xts") { stop(
cat("#Warning:\n Vector 'daily_ret' must be an xts object. Please provide it in the correct form \n")
)}
if(N != N_r_t) { stop(
cat("#Warning:\n Vectors 'x' and 'daily_ret' must have the same length \n")
)}

if(cond_mv_m != "matrix") { stop(
cat("#Warning:\n Parameter 'mv_m' must be a matrix. Please provide it in the correct form \n")
)}
if(dim(mv_m)[2] != N) { stop(
cat("#Warning:\n The columns of the matrix 'mv_m' must be equal to the length of vector 'x'. Please provide it in the correct form \n")
)}

}					# end

#### if the model to estimate is the MEM-MIDAS without the skew parameter

if ((model=="MEMMIDAS"|model=="MEMMIDASX")&skew=="NO") { 

cond_mv_m<- class(mv_m)[1]

if(cond_mv_m != "matrix") { stop(
cat("#Warning:\n Parameter 'mv_m' must be a matrix. Please provide it in the correct form \n")
)}
if(dim(mv_m)[2] != N) { stop(
cat("#Warning:\n The columns of the parameter 'mv_m' must be equal to the length of parameter 'x'. Please provide it in the correct form \n")
)}

}

####################################### configuration of the in-sample and out-of-sample variables

if (missing(out_of_sample)){ # out_of_sample parameter is missing

x_in_s<-x
r_t_in_s<-daily_ret
mv_m_in_s<-mv_m
z_in_s<-z


} else {					# out_of_sample parameter is present

x_in_s<-x[1:(N-out_of_sample)]
r_t_in_s<-daily_ret[1:(N-out_of_sample)]
mv_m_in_s<-mv_m[,1:(N-out_of_sample)]
z_in_s<-z[1:(N-out_of_sample)]

}



############################ setting of the parameters for each model

if (model=="MEMMIDAS"&skew=="YES") {

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=6)
colnames(begin_val)<-c("alpha","beta","gamma","m","theta","w2")
begin_val[,1:3]<-stats::runif(R*3)
begin_val[,1:3]<-begin_val[,1:3]/apply(begin_val[,1:3],1,sum)
begin_val[,3]<-begin_val[,3]/2
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-stats::runif(R,min=-1,max=1)
begin_val[,6]<-1.01

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(MEM_MIDAS_loglik(begin_val[i,],
x_in_s,r_t_in_s,mv_m_in_s,K=K))
}

start_val<-begin_val[which.max(which_row),]

ui<-rbind(
c(1,0,0,0,0,0),       	 	 	## alpha>0.001
c(0,1,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0),	     		## alpha+beta+gamma/2<1
c(0,0,0,0,0,1))					## w2>1.001

ci<-c(-0.001,-0.001,0.999,-1.001)

LOGLIK<-MEM_MIDAS_loglik
cond_vol<-MEM_MIDAS_pred
long_run_vol<-MEM_MIDAS_lr_pred

} else if (model=="MEMMIDAS"&skew=="NO"){

start_val<-begin_val<-ui<-ci<-NULL

begin_val<-matrix(NA,nrow=R,ncol=5)
colnames(begin_val)<-c("alpha","beta","m","theta","w2")
begin_val[,1:2]<-stats::runif(R*2)
begin_val[,1:2]<-begin_val[,1:2]/apply(begin_val[,1:2],1,sum)
begin_val[,1]<-begin_val[,1]/2
begin_val[,3]<-stats::runif(R,min=-1,max=1)
begin_val[,4]<-stats::runif(R,min=-1,max=1)
begin_val[,5]<-1.01

which_row<-rep(NA,R)

for(i in 1:R){
which_row[i]<-sum(MEM_MIDAS_loglik_no_skew(begin_val[i,],
x_in_s,mv_m_in_s,K=K))
}

start_val<-begin_val[which.max(which_row),]


ui<-rbind(
c(1,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0),	     		## alpha+beta<1
c(0,0,0,0,1))					## w2>1.001

ci<-c(-0.0001,-0.001,0.999,-1.001)

LOGLIK<-MEM_MIDAS_loglik_no_skew
cond_vol<-MEM_MIDAS_pred_no_skew
long_run_vol<-MEM_MIDAS_lr_pred_no_skew

} else if (model=="MEM"&skew=="YES"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.28,beta=0.63,gamma=0.05)

ui<-rbind(
c(1,0,0),       	 	 	 ## alpha>0.0001
c(0,1,0),        		 	 ## beta>0.001
c(-1,-1,-0.5))	     	 ## alpha+beta+gamma/2<1

ci<-c(-0.0001,-0.001,0.999)


LOGLIK<-MEM_loglik
cond_vol<-MEM_pred

} else if (model=="MEM"&skew=="NO"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.28,beta=0.63)

ui<-rbind(
c(1,0),       	 	 	## alpha>0.0001
c(0,1),        		 	## beta>0.001
c(-1,-1))	     			## alpha+beta<1					

ci<-c(-0.0001,-0.001,0.999)

LOGLIK<-MEM_loglik_no_skew
cond_vol<-MEM_pred_no_skew

} else if (model=="MEMX"&skew=="YES"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.28,beta=0.63,gamma=0.05,delta=0.01)

ui<-rbind(
c(1,0,0,0),       	 	 	 ## alpha>0.0001
c(0,1,0,0),        		 	 ## beta>0.001
c(-1,-1,-0.5,0))	     	 ## alpha+beta+gamma/2<1

ci<-c(-0.0001,-0.001,0.999)


LOGLIK<-MEM_X_loglik
cond_vol<-MEM_X_pred

} else if (model=="MEMX"&skew=="NO"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.28,beta=0.63,delta=0.01)

ui<-rbind(
c(1,0,0),       	 	 	 ## alpha>0.0001
c(0,1,0),        		 	 ## beta>0.001
c(-1,-1,0))	     	 ## alpha+beta<1

ci<-c(-0.0001,-0.001,0.999)


LOGLIK<-MEM_X_loglik_no_skew
cond_vol<-MEM_X_pred_no_skew

} else if (model=="MEMMIDASX"&skew=="YES") {

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.28,beta=0.63,gamma=0.11,m=0,theta=0,w2=5,delta=0.01)

ui<-rbind(
c(1,0,0,0,0,0,0),       	 	 	## alpha>0.001
c(0,1,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0,0),	     		## alpha+beta+gamma/2<1
c(0,0,0,0,0,1,0))					## w2>1.001

ci<-c(-0.001,-0.001,0.999,-1.001)

LOGLIK<-MEM_MIDAS_X_loglik
cond_vol<-MEM_MIDAS_X_pred
long_run_vol<-MEM_MIDAS_X_lr_pred

} else if (model=="MEMMIDASX"&skew=="NO"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.28,beta=0.63,m=0,theta=0,w2=5,delta=0.01)

ui<-rbind(
c(1,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0,0),	     		## alpha+beta<1
c(0,0,0,0,1,0))					## w2_pos>1.001

ci<-c(-0.0001,-0.001,0.999,-1.001)

LOGLIK<-MEM_MIDAS_X_loglik_no_skew
cond_vol<-MEM_MIDAS_X_pred_no_skew
long_run_vol<-MEM_MIDAS_X_lr_pred_no_skew

}


########################### begin estimate

r_t_in_s_est<-zoo::coredata(r_t_in_s)
x_in_s_est<-zoo::coredata(x_in_s) 

if (model=="MEMMIDAS"&skew=="YES"){

est<-maxLik(
logLik=LOGLIK,
start=start_val,
x=x_in_s_est,
daily_ret=r_t_in_s_est,
mv_m=mv_m_in_s,
K=K,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS")

} else if (model=="MEMMIDAS"&skew=="NO"){

est<-maxLik(
logLik=LOGLIK,
start=start_val,
x=x_in_s_est,
mv_m=mv_m_in_s,
K=K,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS")

} else if (model=="MEM"&skew=="NO"){

est<-maxLik(
logLik=LOGLIK,
start=start_val,
x=x_in_s_est,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS")

} else if (model=="MEM"&skew=="YES"){

est<-maxLik(
logLik=LOGLIK,
start=start_val,
x=x_in_s_est,
daily_ret=r_t_in_s_est,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS")

} else if (model=="MEMX"&skew=="YES"){

z_in_s_est<-zoo::coredata(z_in_s) - mean(zoo::coredata(z_in_s))

est<-suppressWarnings(maxLik(
logLik=LOGLIK,
start=start_val,
x=x_in_s_est,
z=z_in_s_est,
daily_ret=r_t_in_s_est,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))

} else if (model=="MEMX"&skew=="NO"){

z_in_s_est<-zoo::coredata(z_in_s) - mean(zoo::coredata(z_in_s))

est<-suppressWarnings(maxLik(
logLik=LOGLIK,
start=start_val,
x=x_in_s_est,
z=z_in_s_est,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))

} else if (model=="MEMMIDASX"&skew=="YES"){

z_in_s_est<-zoo::coredata(z_in_s) - mean(zoo::coredata(z_in_s))

est<-suppressWarnings(maxLik(
logLik=LOGLIK,
start=start_val,
x=x_in_s_est,
daily_ret=r_t_in_s_est,
mv_m=mv_m_in_s,
K=K,
z=z_in_s_est,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))

} else if (model=="MEMMIDASX"&skew=="NO"){

z_in_s_est<-zoo::coredata(z_in_s) - mean(zoo::coredata(z_in_s))

est<-suppressWarnings(maxLik(
logLik=LOGLIK,
start=start_val,
x=x_in_s_est,
mv_m=mv_m_in_s,
K=K,
z=z_in_s_est,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))

}

##############################################


N_coef<-length(stats::coef(est))

mat_coef<-data.frame(rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef))
colnames(mat_coef)<-c("Estimate","Std. Error","t value","Pr(>|t|)")

rownames(mat_coef)<-names(stats::coef(est))

est_coef<-stats::coef(est)


if (model=="MEMMIDAS"&skew=="YES"){

lr_vol<-long_run_vol(est_coef,x_in_s,r_t_in_s,mv_m_in_s,K=K)
vol_est<-cond_vol(est_coef,x_in_s,r_t_in_s,mv_m_in_s,K=K)

} else if (model=="MEMMIDAS"&skew=="NO"){

lr_vol<-long_run_vol(est_coef,x_in_s,mv_m_in_s,K=K)
vol_est<-cond_vol(est_coef,x_in_s,mv_m_in_s,K=K)

} else if (model=="MEMMIDASX"&skew=="YES"){

lr_vol<-long_run_vol(est_coef,x_in_s,r_t_in_s,mv_m_in_s,K=K,z_in_s)
vol_est<-cond_vol(est_coef,x_in_s,r_t_in_s,mv_m_in_s,K=K,z_in_s)

} else if (model=="MEMMIDASX"&skew=="NO"){

lr_vol<-long_run_vol(est_coef,x_in_s,mv_m_in_s,K=K,z_in_s)
vol_est<-cond_vol(est_coef,x_in_s,mv_m_in_s,K=K,z_in_s)

} else if (model=="MEM"&skew=="NO") {

lr_vol<-c("no long-run component in the MEM model")
vol_est<-cond_vol(est_coef,x_in_s)

} else if (model=="MEM"&skew=="YES") {

lr_vol<-c("no long-run component in the MEM model")
vol_est<-cond_vol(est_coef,x_in_s,r_t_in_s)

} else if (model=="MEMX"&skew=="YES") {

lr_vol<-c("no long-run component in the MEM-X model")
vol_est<-cond_vol(est_coef,x_in_s,r_t_in_s,z_in_s)

} else if (model=="MEMX"&skew=="NO") {

lr_vol<-c("no long-run component in the MEM-X model")
vol_est<-cond_vol(est_coef,x_in_s,z_in_s)

}

#std_est<-MEM_QMLE_sd(est,x_in_s,vol_est)
std_est<-QMLE_sd(est)

mat_coef[,1]<-round(stats::coef(est),6)
mat_coef[,2]<-round(std_est,6)
mat_coef[,3]<-round(stats::coef(est)/std_est,6)
mat_coef[,4]<-round(apply(rbind(stats::coef(est)/std_est),1,function(xx) 2*(1-stats::pnorm(abs(xx)))),6)

## change the order of variables
if (skew=="YES") {

mat_coef2<-mat_coef

mat_coef2[2,]<-mat_coef[3,] #gamma
rownames(mat_coef2)[3]<-"temp"
rownames(mat_coef2)[2]<-"gamma"
mat_coef2[3,]<-mat_coef[2,] #beta
rownames(mat_coef2)[3]<-"beta"
mat_coef<-mat_coef2
}


######## in-sample and out-of-sample estimation and evaluation


if (missing(out_of_sample)){

res<-list(
model=model,
rob_coef_mat=mat_coef,
obs=N,
period=range(stats::time(x)),
loglik=as.numeric(stats::logLik(est)),
inf_criteria=Inf_criteria(est),
loss_in_s=LF_f(vol_est,x_in_s),
est_in_s=vol_est,
est_lr_in_s=lr_vol)

} else {

x_oos<-x[(N-out_of_sample+1):N]
r_t_oos<-daily_ret[(N-out_of_sample+1):N]
mv_m_oos<-mv_m[,(N-out_of_sample+1):(N)]

if (!is.null(z)){
z_oos<-z - mean(z)
z_oos<-z_oos[(N-out_of_sample+1):N]
}

if (model=="MEMMIDAS"&skew=="YES"){

LR_oos<-long_run_vol(est_coef,x_oos,r_t_oos,mv_m_oos,K=K)
vol_est_oos<-cond_vol(est_coef,x_oos,r_t_oos,mv_m_oos,K=K)

} else if (model=="MEMMIDAS"&skew=="NO"){

LR_oos<-long_run_vol(est_coef,x_oos,mv_m_oos,K=K)
vol_est_oos<-cond_vol(est_coef,x_oos,mv_m_oos,K=K)

} else if (model=="MEMMIDASX"&skew=="YES"){

LR_oos<-long_run_vol(est_coef,x_oos,r_t_oos,mv_m_oos,K=K,z_oos)
vol_est_oos<-cond_vol(est_coef,x_oos,r_t_oos,mv_m_oos,K=K,z_oos)

} else if (model=="MEMMIDASX"&skew=="NO"){

LR_oos<-long_run_vol(est_coef,x_oos,mv_m_oos,K=K,z_oos)
vol_est_oos<-cond_vol(est_coef,x_oos,mv_m_oos,K=K,z_oos)

} else if (model=="MEM"&skew=="NO") {

LR_oos<-c("no long-run component in the MEM model")
vol_est_oos<-cond_vol(est_coef,x_oos)

} else if (model=="MEM"&skew=="YES") {

LR_oos<-c("no long-run component in the MEM model")
vol_est_oos<-cond_vol(est_coef,x_oos,r_t_oos)

} else if (model=="MEMX"&skew=="YES") {

LR_oos<-c("no long-run component in the MEM-X model")
vol_est_oos<-cond_vol(est_coef,x_oos,r_t_oos,z_oos)

} else if (model=="MEMX"&skew=="NO") {

LR_oos<-c("no long-run component in the MEM-X model")
vol_est_oos<-cond_vol(est_coef,x_oos,z_oos)

}

res<-list(
model=model,
rob_coef_mat=mat_coef,
obs=N,
period=range(stats::time(x)),
loglik=as.numeric(stats::logLik(est)),
inf_criteria=Inf_criteria(est),
loss_in_s=LF_f(vol_est,x_in_s),
est_in_s=vol_est,
est_lr_in_s=lr_vol,
loss_oos=LF_f(vol_est_oos,x_oos),
est_oos=vol_est_oos,
est_lr_oos=LR_oos
)
}

class(res)<-c("rumidas")
return(res)
print.rumidas(res)

}





