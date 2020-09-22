#' GARCH-MIDAS log-likelihood (with skewness)
#'
#' Obtains the log-likelihood of the GARCH-MIDAS, with an asymmetric term linked to past negative returns,
#' according to two errors' conditional distributions: Normal and Student-t. 
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of starting values. It must be a six- or seven- dimensional vector. See the examples below.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param distribution The conditional density to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' # conditional density of the innovations: normal
#' # start_val<-c(alpha=0.01,beta=0.8,gamma=0.05,m=0,theta=0.1,w2=2)
#' # r_t<-sp500['2005/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # sum(GM_loglik(start_val,r_t,mv_m,K=12,distribution="norm"))
#' 
#' # conditional density of the innovations: Student-t
#' # start_val<-c(alpha=0.01,beta=0.8,gamma=0.05,m=0,theta=0.1,w2=2,shape=5)
#' # r_t<-sp500['2005/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # sum(GM_loglik(start_val,r_t,mv_m,K=12,distribution="std"))
#' @keywords internal
#' @export

GM_loglik<-function(param,
daily_ret,
mv_m,
K,
distribution){

if(distribution=="norm"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta			<- param[5]
        w1			<- 1
        w2			<- param[6]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
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
        theta			<- param[5]
        w1			<- 1
        w2			<- param[6]
		v			<- param[7]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
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

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

 for(i in 2:TT){
  g_it[i]      <- sum(step_1[i-1],beta*g_it[i-1],na.rm=T)
			}

###### variance 

h_it<-coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

)

}

return(ll)
##### end
}

#' GARCH-MIDAS conditional volatility (with skewness)
#'
#' Obtains the estimated conditional volatility for the GARCH-MIDAS model, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of estimated values. It must be a six- or seven- dimensional vector. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#'  \insertAllCited{} 
#' @examples
#' # estimated volatility
#' # est_val<-c(alpha=0.01,beta=0.8,gamma=0.05,m=2,theta=0.1,w2=2)
#' # r_t<-sp500['/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # head(GM_cond_vol(est_val,r_t,mv_m,K=12))
#' @keywords internal
#' @export

GM_cond_vol<-function(param,
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
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
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
#' @param param Vector of estimated values. It must be a six- or seven- dimensional vector. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#'  \insertAllCited{} 
#' @examples
#' # estimated volatility
#' # est_val<-c(alpha=0.01,beta=0.8,gamma=0.05,m=2,theta=0.1,w2=2)
#' # r_t<-sp500['/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # head(GM_long_run_vol(est_val,r_t,mv_m,K=12))
#' @keywords internal
#' @export

GM_long_run_vol<-function(param,
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
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
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

step_1<-(1-alpha-beta-gamma_1/2)+
(alpha+gamma_1*daily_ret_neg)*(daily_ret)^2/tau_d

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
#' @param param Vector of starting values. It must be a six- or seven- dimensional vector. See the examples below.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param distribution The conditional dentisity to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' # conditional density of the innovations: normal
#' # start_val<-c(alpha=0.01,beta=0.8,m=0,theta=0.1,w2=2)
#' # r_t<-sp500['2005/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # sum(GM_loglik_no_skew(start_val,r_t,mv_m,K=12,distribution="norm"))
#' 
#' # conditional density of the innovations: Student-t
#' # start_val<-c(alpha=0.01,beta=0.8,m=0,theta=0.1,w2=2,shape=5)
#' # r_t<-sp500['2005/2010']
#' # mv_m<-mv_into_mat(r_t,indpro,K=12,"monthly")
#' # sum(GM_loglik_no_skew(start_val,r_t,mv_m,K=12,distribution="std"))
#' @keywords internal
#' @export

GM_loglik_no_skew<-function(param,
daily_ret,
mv_m,
K,
distribution){

if(distribution=="norm"){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta			<- param[4]
        w1			<- 1
        w2			<- param[5]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

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
        theta			<- param[4]
        w1			<- 1
        w2			<- param[5]
		v			<- param[6]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

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

h_it<-coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

)

}

return(ll)
##### end
}

#' GARCH-MIDAS conditional volatility (without skewness)
#'
#' Obtains the estimated conditional volatility for the GARCH-MIDAS model, without the skewness parameter in the short-run..
#' For details, see \insertCite{engle_ghysels_sohn_2013;textual}{rumidas} and \insertCite{conrad_lock_2015;textual}{rumidas}.
#' @param param Vector of estimated values. It must be a six- or seven- dimensional vector. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' # est_val<-c(alpha=0.01,beta=0.8,m=2,theta=0.1,w2=2)
#' # r_t<-sp500['/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # head(GM_cond_vol_no_skew(est_val,r_t,mv_m,K=12))
#' @keywords internal
#' @export

GM_cond_vol_no_skew<-function(param,
daily_ret,
mv_m,
K){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta			  <- param[4]
        w1			  <- 1
        w2			  <- param[5]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

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
#' @param param Vector of estimated values. It must be a six- or seven- dimensional vector. 
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' # est_val<-c(alpha=0.01,beta=0.8,m=2,theta=0.1,w2=2)
#' # r_t<-sp500['/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # head(GM_long_run_vol_no_skew(est_val,r_t,mv_m,K=12))
#' @keywords internal
#' @export

GM_long_run_vol_no_skew<-function(param,
daily_ret,
mv_m,
K){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta			  <- param[4]
        w1			  <- 1
        w2			  <- param[5]
         
		
       TT              <- length(daily_ret)
	

		tau_d<-rep(NA,TT)

        g_it                    <- rep(1,TT)            # daily conditional variance 
        #g_it[[1]]        		<- 0
		       
        ll                                      <- 0 

##### daily long-run

betas<-c(rev(beta_function(1:(K+1),(K+1),w1,w2))[2:(K+1)],0)

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

#' DAGM log-likelihood (with skewness)
#'
#' Obtains the log-likelihood of the DAGM, with an asymmetric term linked to past negative returns,
#' according to two errors' conditional distributions: Normal and Student-t. 
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of starting values. It must be a eight- or nine- dimensional vector. See the examples below.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param distribution The conditional dentisity to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' # conditional density of the innovations: normal
#' # start_val<-c(0.01,0.80,0.05,0,0,1.1,0,1.1)
#' # r_t<-sp500['2005/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # sum(DAGM_loglik(start_val,r_t,mv_m,K=12,distribution="norm"))
#' 
#' @keywords internal
#' @export

DAGM_loglik<-function(param,
daily_ret,
mv_m,
K,
distribution){

if(distribution=="norm"){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta_pos       <- param[5]
        w1_pos          <- 1
        w2_pos          <- param[6]
        theta_neg       <- param[7]
        w1_neg          <- 1
        w2_neg          <- param[8]
	   
        TT              <- length(daily_ret)

		MV_pos<-ifelse(mv_m>=0,mv_m,0)
         MV_neg<-ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg                       <- ifelse(coredata(daily_ret) <0,1,0)

        ll                                      <- 0 


###### long-run 

betas_pos<-c(rev(beta_function(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(beta_function(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


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

  	 alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta_pos       <- param[5]
        w1_pos          <- 1
        w2_pos          <- param[6]
        theta_neg       <- param[7]
        w1_neg          <- 1
        w2_neg          <- param[8]
		v			  <- param[9]
	   
        TT              <- length(daily_ret)

		MV_pos<-ifelse(mv_m>=0,mv_m,0)
         MV_neg<-ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg                       <- ifelse(coredata(daily_ret) <0,1,0)

        ll                                      <- 0 


###### long-run 

betas_pos<-c(rev(beta_function(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(beta_function(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


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

h_it<-coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

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
#' @param param Vector of starting values. It must be a seven- or eight- dimensional vector. See the examples below.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param distribution The conditional dentisity to use for the innovations. At the moment, valid choices are "norm" and "std", for the Normal 
#' and Student-t distributions.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' # conditional density of the innovations: normal
#' # start_val<-c(alpha=0.01,beta=0.80,gamma_1=0.05,m=0,theta_pos=0,w2_pos=1.1,theta_neg=0,w2_neg=1.1)
#' # r_t<-sp500['2005/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # sum(DAGM_loglik(start_val,r_t,mv_m,K=12,distribution="norm"))
#' 
#' # conditional density of the innovations: Student-t
#' # start_val<-c(0.01,0.80,0.05,0,0,1.1,0,1.1,5)
#' # r_t<-sp500['2005/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # sum(DAGM_loglik(start_val,r_t,mv_m,K=12,distribution="std"))
#' @keywords internal
#' @export

DAGM_loglik_no_skew<-function(param,
daily_ret,
mv_m,
K,
distribution){

if(distribution=="norm"){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta_pos       <- param[4]
        w1_pos          <- 1
        w2_pos          <- param[5]
        theta_neg       <- param[6]
        w1_neg          <- 1
        w2_neg          <- param[7]
	   
        TT              <- length(daily_ret)

		MV_pos<-ifelse(mv_m>=0,mv_m,0)
         MV_neg<-ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        ll                                      <- 0 


###### long-run 

betas_pos<-c(rev(beta_function(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(beta_function(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


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
        w1_pos          <- 1
        w2_pos          <- param[5]
        theta_neg       <- param[6]
        w1_neg          <- 1
        w2_neg          <- param[7]
		v			  <- param[8]
	   
        TT              <- length(daily_ret)

		MV_pos<-ifelse(mv_m>=0,mv_m,0)
         MV_neg<-ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       

        ll                                      <- 0 


###### long-run 

betas_pos<-c(rev(beta_function(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(beta_function(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


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

h_it<-coredata(g_it*tau_d)

###### loglik

ll <- log(

gamma((v+1)/2)*gamma(v/2)^-1*((v-2)*h_it)^-0.5*
(1+coredata(daily_ret)^2*h_it^-1*(v-2)^-1)^(-(v+1)/2)

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
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' # estimated volatility 
#' # est_val<-c(0.01,0.80,0.05,0,0.1,1.1,-0.3,1.1)
#' # r_t<-sp500['2005/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # vol<-DAGM_cond_vol(est_val,r_t,mv_m,K=12)
#' # head(vol)
#' @keywords internal
#' @export

DAGM_cond_vol<-function(param,
daily_ret,
mv_m,
K){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta_pos       <- param[5]
        w1_pos          <- 1
        w2_pos          <- param[6]
        theta_neg       <- param[7]
        w1_neg          <- 1
        w2_neg          <- param[8]
	   
        TT              <- length(daily_ret)

		MV_pos<-ifelse(mv_m>=0,mv_m,0)
         MV_neg<-ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg                       <- ifelse(coredata(daily_ret) <0,1,0)

        ll                                      <- 0 


###### long-run 

betas_pos<-c(rev(beta_function(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(beta_function(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


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
#' @return The resulting vector is an "xts" object representing the conditional volatility.
#' @importFrom Rdpack reprompt
#' @import roll
#' @seealso \code{\link{mv_into_mat}}.
#' @references
#' \insertAllCited{} 
#' @examples
#' # est_val<-c(0.01,0.80,0.05,0,0.1,1.1,-0.3,1.1)
#' # r_t<-sp500['/2010']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly")
#' # head(DAGM_long_run_vol(est_val,r_t,mv_m,K=12))
#' @keywords internal
#' @export

DAGM_long_run_vol<-function(param,
daily_ret,
mv_m,
K){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
        m               <- param[4]
        theta_pos       <- param[5]
        w1_pos          <- 1
        w2_pos          <- param[6]
        theta_neg       <- param[7]
        w1_neg          <- 1
        w2_neg          <- param[8]
	   
        TT              <- length(daily_ret)

		MV_pos<-ifelse(mv_m>=0,mv_m,0)
         MV_neg<-ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        daily_ret_neg                       <- ifelse(coredata(daily_ret) <0,1,0)

        ll                                      <- 0 


###### long-run 

betas_pos<-c(rev(beta_function(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(beta_function(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


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
#' Obtains the conditional volatility for the DAGM, without the asymmetric term linked to past negative returns.
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of estimated values. It must be a seven- or eight- dimensional vector. See the examples below.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
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
K){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta_pos       <- param[4]
        w1_pos          <- 1
        w2_pos          <- param[5]
        theta_neg       <- param[6]
        w1_neg          <- 1
        w2_neg          <- param[7]
	   
        TT              <- length(daily_ret)

		MV_pos<-ifelse(mv_m>=0,mv_m,0)
         MV_neg<-ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        ll                                      <- 0 


###### long-run 

betas_pos<-c(rev(beta_function(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(beta_function(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


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
#' Obtains the daily long-run volatility for the DAGM, without the asymmetric term linked to past negative returns.
#' For details, see \insertCite{amendola_candila_gallo:2019;textual}{rumidas}.
#' @param param Vector of estimated values. It must be a seven- or eight- dimensional vector. See the examples below.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
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
#' # head(DAGM_long_run_vol_no_skew(est_val,r_t,mv_m,K=12))
#' @keywords internal
#' @export

DAGM_long_run_vol_no_skew<-function(param,
daily_ret,
mv_m,
K){

        alpha           <- param[1]
        beta            <- param[2]
        m               <- param[3]
        theta_pos       <- param[4]
        w1_pos          <- 1
        w2_pos          <- param[5]
        theta_neg       <- param[6]
        w1_neg          <- 1
        w2_neg          <- param[7]
	   
        TT              <- length(daily_ret)

		MV_pos<-ifelse(mv_m>=0,mv_m,0)
         MV_neg<-ifelse(mv_m<0,mv_m,0)

        tau_d                   <- rep(NA,TT)
	   
        g_it                    <- rep(1,TT)            # daily conditional variance 
		       
        ll                                      <- 0 


###### long-run 

betas_pos<-c(rev(beta_function(1:(K+1),(K+1),w1_pos,w2_pos))[2:(K+1)],0)
betas_neg<-c(rev(beta_function(1:(K+1),(K+1),w1_neg,w2_neg))[2:(K+1)],0)


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

#' MEM-MIDAS log-likelihood (with skewness parameter)
#'
#' Obtains the log-likelihood of the MEM-MIDAS, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{amendola2020doubly;textual}{rumidas} and \insertCite{engle_gallo_2006;textual}{rumidas}.
#' @param param Vector of starting values. It must be a six--dimensional vector. See the example below.
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
#' # start_val<-c(alpha=0.10,beta=0.8,gamma=0.1,m=0,theta=-0.16,w2=5)
#' # r_t<-sp500['/2010']
#' # real<-(rv5['/2010'])^0.5		# realized volatility
#' # mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' # sum(MEM_MIDAS_loglik(start_val,real,r_t,mv_m,K=12))
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

#############################

ll <- - (
log(mu_it*tau_d) + x/(mu_it*tau_d)
)


return(coredata(ll))
##### fine
}

#' MEM-MIDAS log-likelihood (no skewness parameter)
#'
#' Obtains the log-likelihood of the MEM-MIDAS.
#' For details, see \insertCite{amendola2020doubly;textual}{rumidas} and \insertCite{engle_gallo_2006;textual}{rumidas}.
#' @param param Vector of starting values. It must be a five--dimensional vector. See the example below.
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
#' # start_val<-c(alpha=0.10,beta=0.8,m=0,theta=-0.16,w2=5)
#' # real<-(rv5['/2010'])^0.5		# realized volatility
#' # mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' # sum(MEM_MIDAS_loglik_no_skew(start_val,real,mv_m,K=12))
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


return(coredata(ll))
##### fine
}

#' MEM log-likelihood (with skewness parameter)
#'
#' Obtains the log-likelihood of the base MEM, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{engle_gallo_2006;textual}{rumidas}.
#' @param param Vector of starting values. It must be a two--dimensional vector. See the example below.
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @param daily_ret Daily returns, which must be an "xts" object, and with the same length of x.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @examples
#' # start_val<-c(alpha=0.10,beta=0.8,gamma=0.05)
#' # real<-(rv5['/2010'])^0.5		# realized volatility
#' # r_t<-sp500['/2010']
#' # sum(MEM_loglik(start_val,real,r_t))
#' @keywords internal
#' @export

MEM_loglik<-function(param,
x,
daily_ret){

        alpha           <- param[1]
        beta            <- param[2]
        gamma_1         <- param[3]
                		
       N              <- length(daily_ret)

		tau_d<-rep(NA,N)

  sample_mean_x		<-		mean(x)
  mu_it               <- 		rep(sample_mean_x,N)  
       
        daily_ret_neg                       <- ifelse(coredata(daily_ret) < 0,1,0)

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


return(coredata(ll))
##### fine
}

#' MEM log-likelihood (no skewness parameter)
#'
#' Obtains the log-likelihood of the base MEM.
#' For details, see \insertCite{engle_gallo_2006;textual}{rumidas}.
#' @param param Vector of starting values. It must be a two--dimensional vector. See the example below.
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @examples
#' # start_val<-c(alpha=0.10,beta=0.8)
#' # real<-(rv5['/2010'])^0.5		# realized volatility
#' # sum(MEM_loglik_no_skew(start_val,real))
#' @keywords internal
#' @export

MEM_loglik_no_skew<-function(param,x){

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

#############################

ll <- - (
log(mu_it) + x/(mu_it)
)


return(coredata(ll))
##### fine
}

#' MEM one-step-ahead predictions (with skewness parameter)
#'
#' Predicts the dependent variable, usually the realized volatility, for the base MEM, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{engle_gallo_2006;textual}{rumidas}.
#' @param param Vector of starting values. It must be a three--dimensional vector. See the example below.
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @param daily_ret Daily returns, which must be an "xts" object, and with the same length of x.
#' @return The resulting vector is the one-step-ahead prediction for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @examples
#' # est_val<-c(alpha=0.10,beta=0.8,gamma=0.05)
#' # real<-(rv5['/2010'])^0.5		# realized volatility
#' # r_t<-sp500['/2010']
#' # head(MEM_pred(est_val,real,r_t))
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
       
        daily_ret_neg                       <- ifelse(coredata(daily_ret) < 0,1,0)

        ll                                      <- 0 

#######  

step_1 <- (1-alpha-beta-gamma_1/2)*sample_mean_x+
(alpha+gamma_1*daily_ret_neg)*(x)

 for(i in 2:N){
  mu_it[i]      <- sum(step_1[i-1],beta*mu_it[i-1],na.rm=T)
			}

####### predictions

x_hat<-mu_it

x_hat<-as.xts(x_hat,index(x))

return(x_hat)

##### fine
}


#' MEM one-step-ahead predictions (no skewness parameter)
#'
#' Predicts the dependent variable, usually the realized volatility, for the base MEM.
#' For details, see \insertCite{engle_gallo_2006;textual}{rumidas}.
#' @param param Vector of starting values. It must be a two--dimensional vector. See the example below.
#' @param x Dependent variable, usually the realized volatility. It must be positive and "xts" object.
#' @return The resulting vector is the log-likelihood value for each \eqn{i,t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @examples
#' # est_val<-c(alpha=0.10,beta=0.8)
#' # real<-(rv5['/2010'])^0.5		# realized volatility
#' # head(MEM_pred_no_skew(est_val,real))
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

x_hat<-as.xts(x_hat,index(x))

return(x_hat)

##### fine
}


#' MEM-MIDAS one-step-ahead predictions (with skewness parameter)
#'
#' Predicts the dependent variable, usually the realized volatility, for the MEM-MIDAS, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{amendola2020doubly;textual}{rumidas} and \insertCite{engle_gallo_2006;textual}{rumidas}.
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
#' # est_val<-c(alpha=0.10,beta=0.8,gamma=0.1,m=0,theta=-0.16,w2=5)
#' # r_t<-sp500['/2010']
#' # real<-rv5['/2010']^0.5 # realized volatility
#' # mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' # est_vol<-MEM_MIDAS_pred(est_val,real,r_t,mv_m,K=12)
#' # head(est_vol)
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

x_hat<-mu_it*tau_d

x_hat<-as.xts(x_hat,index(x))

return(x_hat)

##### fine
}

#' MEM-MIDAS one-step-ahead predictions (no skewness parameter)
#'
#' Predicts the dependent variable, usually the realized volatility, for the MEM-MIDAS.
#' For details, see \insertCite{amendola2020doubly;textual}{rumidas} and \insertCite{engle_gallo_2006;textual}{rumidas}.
#' @param param Vector of starting values. It must be a five--dimensional vector. See the example below.
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
#' # est_val<-c(alpha=0.10,beta=0.8,m=0,theta=-0.16,w2=5)
#' # real<-(rv5['/2010'])^0.5		# realized volatility
#' # mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' # sum(MEM_MIDAS_pred_no_skew(est_val,real,mv_m,K=12))
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

x_hat<-as.xts(x_hat,index(x))

return(x_hat)

##### fine
}

#' MEM-MIDAS long-run one-step-ahead predictions (with skewness parameter)
#'
#' Predicts the long-run term of the dependent variable, usually the realized volatility, for the MEM-MIDAS, with an asymmetric term linked to past negative returns.
#' For details, see \insertCite{amendola2020doubly;textual}{rumidas} and \insertCite{engle_gallo_2006;textual}{rumidas}.
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
#' # est_val<-c(alpha=0.10,beta=0.8,gamma=0.1,m=0,theta=-0.16,w2=5)
#' # r_t<-sp500['/2010']
#' # real<-(rv5['/2010'])^0.5		# realized volatility
#' # mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' # est_vol<-MEM_MIDAS_pred(est_val,real,r_t,mv_m,K=12)
#' # head(est_vol)
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

x_hat<-as.xts(x_hat,index(x))

return(x_hat)

##### fine
}

#' MEM-MIDAS long-run one-step-ahead predictions (no skewness parameter)
#'
#' Predicts the long-run term of the dependent variable, usually the realized volatility, for the MEM-MIDAS.
#' For details, see \insertCite{amendola2020doubly;textual}{rumidas} and \insertCite{engle_gallo_2006;textual}{rumidas}.
#' @param param Vector of starting values. It must be a five--dimensional vector. See the example below.
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
#' # est_val<-c(alpha=0.10,beta=0.8,m=0,theta=-0.16,w2=5)
#' # real<-(rv5['/2010'])^0.5		# realized volatility
#' # mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' # sum(MEM_MIDAS_lr_pred_no_skew(est_val,real,mv_m,K=12))
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

x_hat<-as.xts(x_hat,index(x))

return(x_hat)

##### fine
}


#' Methods for obtaining (and evaluating) a variety of GARCH-MIDAS-based models
#'
#' Estimates several GARCH-MIDAS-based models, according to two errors' conditional distributions: Normal and Student-t, and 
#' the presence of asymmetric terms in the short- and long-run components.
#' @param model Model to estimate. Valid choices are: "GM" for GARCH-MIDAS, "DAGM" for Double Asymmetric GARCH-MIDAS.
#' @param skew The skewness parameter to include in the short-run equation. Valid choices are: "YES" and "NO".
#' @param distribution The conditional density to use for the innovations. At the moment, valid choices are "norm" and "std", 
#' for the Normal and Student-t distribution, respectively.
#' @param daily_ret Daily returns, which must be an "xts" object.
#' @param mv_m MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function.
#' @param K Number of (lagged) realizations of the MIDAS variable to consider.
#' @param out_of_sample **optional**. A positive integer indicating the number of periods before the last to keep for out of sample forecasting.
#' @param vol_proxy **optional**. If present, the vol_proxy is the volatilty proxy used for the in-sample and out-of-sample (again, if present) evaluation. 
#' It could be the realized variance. If it left unspecified, vol_proxy is replaced by the squared daily returns. 
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
#'	\item est_in_s: The in-sample predicted volatility, that is: \eqn{\sqrt{\hat{\tau}_t \times \hat{g}_{i,t} }}.
#'	\item est_lr_in_s: The in-sample predicted long-run component of the dependent variable.
#'	\item loss_oos: The out-of-sample MSE and QLIKE averages, calculated considering the distance with respect to the volatility proxy (if provided) or the squared daily returns.
#'	\item est_oos: The out-of-sample predicted dependent variable, that is: \eqn{\sqrt{\hat{\tau}_t \times \hat{g}_{i,t} }}.
#'	\item est_lr_oos: The out-of-sample predicted long-run component of the dependent variable.
#' }
#' @importFrom Rdpack reprompt
#' @importFrom stats time
#' @importFrom zoo coredata
#' @import maxLik
#' @seealso \code{\link{mv_into_mat}}.
#' 
#' @details
#' Function \code{ugmfit} implements the estimation and evaluation of the GARCH-MIDAS-based models, with and without the asymmetric term 
#' linked to negative lagged daily returns, according to two distributions for the error term. The general framework assumes that: 
#' \deqn{r_{i,t}= \sqrt{\tau_t \times g_{i,t} } \epsilon_{i,t},}
#' where 
#' \itemize{
#' \item \eqn{r_{i,t}} is the daily return for the \eqn{i}-th day (\eqn{i = 1, \ldots, N_t}) 
#' of the period \eqn{t} (for example, a week, a month or a quarter; \eqn{t = 1 , \ldots, T});
#' \item \eqn{\tau_{t}} is the long-run component, varying each period \eqn{t};
#' \item \eqn{g_{i,t}} is the short-run term, varying each day \eqn{i} of the period \eqn{t};
#' \item \eqn{\epsilon_{i,t}} is an \eqn{iid} error term which has a zero mean and unit variance.
#' }
#' The short-run component of the GARCH-MIDAS (parameter "model" set to "GM") and DAGM, when the parameter "skew" is "YES", is:
#' \deqn{g_{i,t} = \left(1-\alpha-\gamma/2-\beta\right) + \left(\alpha + \gamma \cdot I_{\left(r_{i-1,t}  < 0 \right)}\right) \frac{\left(r_{i-1,t}\right)^2}{\tau_t} + \beta g_{i-1,t},}
#' where \eqn{I_{(\cdot)}} is an indicator function. When, for both the models, the parameter "skew" is set to "NO", \eqn{\gamma} disappears.
#' The long-run component of the GARCH-MIDAS is:
#' \deqn{\tau_{t} = \exp \left\{ m + \theta \sum_{j=1}^K \delta_{j}(\omega) X_{t-j}\right\},}
#' where \eqn{X_{t}} is the MIDAS term.  
#' The long-run component of the DAGM model is:
#' \deqn{\tau_t   =  \exp \left( m + \theta^{+}  \sum_{k=1}^K \delta_k(\omega)^{+}  X_{t-k} I_{\left( X_{t-k} \geq 0 \right)} +  \theta^{-}  \sum_{k=1}^K \delta_k(\omega)^{-} X_{t-k} I_{\left( X_{t-k} < 0 \right)} \right)}
#'
#' @examples
#'
#' # estimate a GARH-MIDAS model, without the skewness parameter
#' # r_t<-sp500['2008']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly") 
#' # fit<-ugmfit(model="GM",skew="NO",distribution="norm",r_t,mv_m,K=12)
#' # fit
#' # summary.rumidas(fit)
#' # names(fit)
#' 
#' # to see the estimated coefficients with the QML standard errors:
#' # fit$rob_coef_mat
#'
#' # estimate a DAGM model, with the skewness parameter, 
#' # including the volatility proxy (realized variance), and
#' # leaving the last 100 observations for the out-of-sample evaluation
#' # r_t<-sp500['2002/2020']
#' # mv_m<-mv_into_mat(r_t,diff(indpro),K=12,"monthly") 
#' # fit_2<-ugmfit(model="DAGM",skew="YES",distribution="norm",r_t,
#' # mv_m,K=12,vol_proxy=rv5['2002/2020'],out_of_sample=100)
#' # fit_2
#' # summary.rumidas(fit_2)
#' @export

ugmfit<-function(
model,
skew,
distribution,
daily_ret,
mv_m,
K,
out_of_sample=NULL,
vol_proxy=NULL
){

############################# check on valid choices

if((model != "GM")&(model != "DAGM")) { stop(cat("#Warning:\n Valid choices for the parameter 'model' are currently 'GM' and 'DAGM' \n"))}
if((skew != "YES")&(skew != "NO")) { stop(cat("#Warning:\n Valid choices for the parameter 'skew' are 'YES' and 'NO' \n"))}

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


############## check if the vol_proxy is provided and if it has the same time span of daily_ret

if(!missing(vol_proxy)){ 
if(any(range(time(daily_ret))!=range(time(vol_proxy)))){
stop(
cat("#Warning:\n The vector 'vol_proxy' has to be observed during the same time span of 'daily_ret' \n")
)}}

if (missing(out_of_sample)){ # out_of_sample parameter is missing

r_t_in_s<-daily_ret
mv_m_in_s<-mv_m
vol_proxy_in_s<-ifelse(missing(vol_proxy),r_t_in_s^2,vol_proxy)


} else {					# out_of_sample parameter is present

r_t_in_s<-daily_ret[1:(N-out_of_sample)]
mv_m_in_s<-mv_m[,1:(N-out_of_sample)]
vol_proxy_in_s<-ifelse(missing(vol_proxy),r_t_in_s^2,vol_proxy[1:(N-out_of_sample)])

}

############################ setting of the parameters for each model

if (model=="GM"&skew=="YES"&distribution=="norm"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.8,gamma=0.05,m=0,theta=0,w2=2)

ui<-rbind(
c(1,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0),	     		## alpha+beta+gamma/2<1
c(0,0,0,0,0,1))				 	## w2_pos>1.001

ci<-c(-0.0001,-0.001,0.999,-1.001)

LOGLIK<-GM_loglik
cond_vol<-GM_cond_vol
long_run_vol<-GM_long_run_vol



} else if (model=="GM"&skew=="YES"&distribution=="std"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.8,gamma=0.05,m=0,theta=0,w2=2,shape=5)

ui<-rbind(
c(1,0,0,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0,0),	     	## alpha+beta+gamma/2<1
c(0,0,0,0,0,1,0),				 	## w2_pos>1.001
c(0,0,0,0,0,0,1))					## shape > 2.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-2.001)

LOGLIK<-GM_loglik
cond_vol<-GM_cond_vol
long_run_vol<-GM_long_run_vol

} else if (model=="GM"&skew=="NO"&distribution=="norm"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.8,m=0,theta=0,w2=2)

ui<-rbind(
c(1,0,0,0,0),       	 			## alpha>0.001
c(0,1,0,0,0),        		 		## beta>0.001
c(-1,-1,0,0,0),	     			## alpha+beta<1
c(0,0,0,0,1))				 		## w2_pos>1.001

ci<-c(-0.001,-0.001,0.999,-1.001)

LOGLIK<-GM_loglik_no_skew
cond_vol<-GM_cond_vol_no_skew
long_run_vol<-GM_long_run_vol_no_skew


} else if (model=="GM"&skew=="NO"&distribution=="std"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.8,m=0,theta=0,w2=2,shape=5)

ui<-rbind(
c(1,0,0,0,0,0),       	 			## alpha>0.0001
c(0,1,0,0,0,0),        		 		## beta>0.001
c(-1,-1,0,0,0,0),	     				## alpha+beta<1
c(0,0,0,0,1,0),				 		## w2_pos>1.001
c(0,0,0,0,0,1))						## shape>2.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-2.001)

LOGLIK<-GM_loglik_no_skew
cond_vol<-GM_cond_vol_no_skew
long_run_vol<-GM_long_run_vol_no_skew

}

else if (model=="DAGM"&skew=="YES"&distribution=="norm"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.90,gamma=0,m=0,theta_pos=0, 
w2_pos=2,theta_neg=0,w2_neg=2)

ui<-rbind(
c(1,0,0,0,0,0,0,0),       	 	 	 	## alpha>0.0001
c(0,1,0,0,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0,0,0),			 	## alpha+beta+gamma/2<1
c(0,0,0,0,0,1,0,0),				 	## w2_pos>1.001
c(0,0,0,0,0,0,0,1))				 	## w2_neg>1.001

ci<-c(-0.0001,-0.001,0.999,-1.001,-1.001)

LOGLIK<-DAGM_loglik
cond_vol<-DAGM_cond_vol
long_run_vol<-DAGM_long_run_vol

}

else if (model=="DAGM"&skew=="YES"&distribution=="std"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.01,beta=0.90,gamma=0,m=0,theta_pos=0, 
w2_pos=2,theta_neg=0,w2_neg=2,shape=5)

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

} else if (model=="DAGM"&skew=="NO"&distribution=="norm"){

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

} else if (model=="DAGM"&skew=="NO"&distribution=="std"){

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

}

r_t_in_s_est<-zoo::coredata(r_t_in_s)

est<-suppressWarnings(maxLik(
logLik=LOGLIK,
start=start_val,
daily_ret=r_t_in_s_est,
mv_m=mv_m_in_s,
K=K,
distribution=distribution,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS"))


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

vol_est<-(cond_vol(est_coef,r_t_in_s,mv_m_in_s,K=K))^2
lr_vol<-long_run_vol(est_coef,r_t_in_s,mv_m_in_s,K=K)


if (missing(out_of_sample)){

res<-list(
model=model,
rob_coef_mat=mat_coef,
obs=N,
period=range(time(r_t_in_s)),
loglik=as.numeric(stats::logLik(est)),
inf_criteria=Inf_criteria(est),
loss_in_s=LF_f(vol_est,vol_proxy_in_s),
est_vol_in_s=vol_est^0.5,
est_lr_in_s=lr_vol)

} else {

r_t_oos<-daily_ret[(N-out_of_sample+1):N]
mv_m_oos<-mv_m[,(N-out_of_sample+1):(N)]
vol_proxy_oos<-ifelse(missing(vol_proxy),r_t_oos^2,vol_proxy[(N-out_of_sample+1):N])

vol_est_oos<-(cond_vol(est_coef,r_t_oos,mv_m_oos,K=K))^2
LR_oos<-(long_run_vol(est_coef,r_t_oos,mv_m_oos,K=K))


res<-list(
model=model,
rob_coef_mat=mat_coef,
obs=N,
period=range(time(r_t_in_s)),
loglik=as.numeric(stats::logLik(est)),
inf_criteria=Inf_criteria(est),
loss_in_s=LF_f(vol_est,vol_proxy_in_s),
est_vol_in_s=vol_est^0.5,
est_lr_in_s=lr_vol,
loss_oos=LF_f(vol_est_oos,vol_proxy_oos),
est_vol_oos=vol_est_oos^0.5,
est_lr_oos=LR_oos)

}


class(res)<-c("rumidas")

return(res)
print.rumidas(res)


}


#' Methods for obtaining (and evaluating) a variety of MEM(-MIDAS)-based models
#'
#' Estimates several MEM-MIDAS-based models.
#' @param model Model to estimate. Valid choices are: "MEMMIDAS" for MEM-MIDAS, "MEM" for base MEM.
#' @param skew The skewness parameter linked to lagged daily returns. Valid choices are: "YES" and "NO".
#' @param x Dependent variable to predict. Usually the realized volatility. It must be positive and "xts" object.
#' @param daily_ret **optional**. Daily returns, which must be an "xts" object. NULL by default.
#' @param mv_m **optional**. MIDAS variable already transformed into a matrix, through \code{\link{mv_into_mat}} function. NULL by default.
#' @param K **optional**. Number of (lagged) realizations of the MIDAS variable to consider. NULL by default.
#' @param out_of_sample **optional**. A positive integer indicating the number of periods before the last to keep for out of sample forecasting.
#' @return \code{umemfit} returns an object of class 'rumidas'. The function \code{\link{summary.rumidas}} 
#' can be used to print a summary of the results. Moreover, an object of class 'rumidas' is a list containing the following components:
#' \itemize{
#' 	\item model: The model used for the estimation.
#'   \item rob_coef_mat: The matrix of estimated coefficients, with the QML standard errors. 
#'    For details, see: \insertCite{Bollerslev_Wooldridge_1992;textual}{rumidas}, \insertCite{engle_gallo_2006;textual}{rumidas}, and  \insertCite{amendola2020doubly;textual}{rumidas}. 
#'   \item obs: The number of daily observations used for the (in-sample) estimation.
#'     \item period: The period of the in-sample estimation.
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
#' Function \code{umemfit} implements the estimation and evaluation of the MEM and MEM-MIDAS models, with and without the asymmetric term 
#' linked to negative lagged daily returns. The general framework assumes that: 
#' \deqn{x_{i,t}= \mu_{i,t}\epsilon_{i,t} = \tau_{t} \xi_{i,t} \epsilon_{i,t},}
#' where 
#' \itemize{
#' \item \eqn{x_{i,t}} is a time series coming from a non-negative discrete time process for the \eqn{i}-th day (\eqn{i = 1, \ldots, N_t}) 
#' of the period \eqn{t} (for example, a week, a month or a quarter; \eqn{t = 1 , \ldots, T});
#' \item \eqn{\tau_{t}} is the long-run component, determining the average level of the conditional mean, varying each period \eqn{t};
#' \item \eqn{\xi_{i,t}} is a factor centered around one, labelled as the short-run term, which plays the role of dumping or amplifying \eqn{\tau_{i,t}};
#' \item \eqn{\epsilon_{i,t}} is an \eqn{iid} error term which, conditionally on the information set, has a unit mean, an unknown variance, and a probability density function 
#' defined over a non-negative support.
#' }
#' The short-run component of the MEM-MIDAS is:
#' \deqn{\xi_{i,t}=(1-\alpha-\gamma/2-\beta) + \left(\alpha +  \gamma \cdot {I}_{\left(r_{i-1,t}  < 0 \right)}\right) \frac{x_{i-1,t}}{\tau_t} + \beta \xi_{i-1,t},}
#' where \eqn{I_{(\cdot)}} is an indicator function and \eqn{r_{i,t}} is the daily return of the day \eqn{i} of the period \eqn{t}.
#' The long-run component of the MEM-MIDAS is:
#' \deqn{\tau_{t} = \exp \left\{ m + \theta \sum_{k=1}^K \delta_{k}(\omega) X_{t-k}\right\},}
#' where \eqn{X_{t}} is the MIDAS term. When the "skew" parameter is set to "NO", \eqn{\gamma} disappears.  
#' The MEM model does not have the difference between the long- and short-run components. Therefore, it directly evolves according to \eqn{\mu_{i,t}}.
#' When the "skew" parameter is present:
#' \deqn{\mu_{i,t}= \left(1-\alpha - \gamma_{1} / 2 - \beta   \right)\mu + (\alpha + \gamma I_{\left(r_{i-1,t}  < 0 \right)}) x_{i-1,t} +   \beta \mu_{i-1,t},}
#' where \eqn{\mu=E(x_{i,t})}. When the "skew" parameter is set to "NO", in the previous equation \eqn{\gamma} cancels.
#' @examples
#' 
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
#' # r_t<-sp500['2003/2010']
#' # real<-(rv5['2003/2010'])^0.5		# realized volatility
#' # mv_m<-mv_into_mat(real,diff(indpro),K=12,"monthly")
#' # fit_2<-umemfit(model="MEMMIDAS",skew="YES",x=real,daily_ret=r_t,mv_m=mv_m,K=12,out_of_sample=200)
#' # fit_2
#' # summary.rumidas(fit_2)
#' # to see the estimated coefficients with the QML standard errors:
#' # fit_2$rob_coef_mat
#' 
#' @export

umemfit<-function(
model,
skew,
x,
daily_ret=NULL,
mv_m=NULL,
K=NULL,
out_of_sample=NULL
){

############################# check on valid choices

if((model != "MEM")&(model != "MEMMIDAS")) { stop(cat("#Warning:\n Valid choices for the parameter 'model' are currently 'MEM' and 'MEMMIDAS' \n"))}
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

################################### checks on the skew parameter

if(skew=="YES"&missing(daily_ret)) { stop(
cat("#Warning:\n Parameter 'daily_ret' is missing. Please provide it \n")
)}

################################### checks on the MEM-MIDAS setting

#### if the model to estimate is the MEM-MIDAS with skew parameter
if (model=="MEMMIDAS"&skew=="YES") { 

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

if (model=="MEMMIDAS"&skew=="NO") { 

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


} else {					# out_of_sample parameter is present

x_in_s<-x[1:(N-out_of_sample)]
r_t_in_s<-daily_ret[1:(N-out_of_sample)]
mv_m_in_s<-mv_m[,1:(N-out_of_sample)]

}



############################ setting of the parameters for each model

if (model=="MEMMIDAS"&skew=="YES") {

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.28,beta=0.63,gamma=0.11,m=0,theta=0,w2=5)

ui<-rbind(
c(1,0,0,0,0,0),       	 	 	## alpha>0.001
c(0,1,0,0,0,0),        		 	## beta>0.001
c(-1,-1,-0.5,0,0,0),	     		## alpha+beta+gamma/2<1
c(0,0,0,0,0,1))					## w2_pos>1.001

ci<-c(-0.001,-0.001,0.999,-1.001)

LOGLIK<-MEM_MIDAS_loglik
cond_vol<-MEM_MIDAS_pred
long_run_vol<-MEM_MIDAS_lr_pred

} else if (model=="MEMMIDAS"&skew=="NO"){

start_val<-ui<-ci<-NULL

start_val<-c(alpha=0.28,beta=0.63,m=0,theta=0,w2=5)

ui<-rbind(
c(1,0,0,0,0),       	 	 	## alpha>0.0001
c(0,1,0,0,0),        		 	## beta>0.001
c(-1,-1,0,0,0),	     		## alpha+beta<1
c(0,0,0,0,1))					## w2_pos>1.001

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

} else if (model=="MEM"){

est<-maxLik(
logLik=LOGLIK,
start=start_val,
x=x_in_s_est,
daily_ret=r_t_in_s_est,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS")

}



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

} else if (model=="MEM"&skew=="NO") {

lr_vol<-c("no long-run component in the MEM model")
vol_est<-cond_vol(est_coef,x_in_s)

} else if (model=="MEM"&skew=="YES") {

lr_vol<-c("no long-run component in the MEM model")
vol_est<-cond_vol(est_coef,x_in_s,r_t_in_s)

}

std_est<-MEM_QMLE_sd(est,x_in_s,vol_est)

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


if (model=="MEMMIDAS"&skew=="YES"){

LR_oos<-long_run_vol(est_coef,x_oos,r_t_oos,mv_m_oos,K=K)
vol_est_oos<-cond_vol(est_coef,x_oos,r_t_oos,mv_m_oos,K=K)

} else if (model=="MEMMIDAS"&skew=="NO"){

LR_oos<-long_run_vol(est_coef,x_oos,mv_m_oos,K=K)
vol_est_oos<-cond_vol(est_coef,x_oos,mv_m_oos,K=K)

} else if (model=="MEM"&skew=="NO") {

LR_oos<-c("no long-run component in the MEM model")
vol_est_oos<-cond_vol(est_coef,x_oos)

} else if (model=="MEM"&skew=="YES") {

LR_oos<-c("no long-run component in the MEM model")
vol_est_oos<-cond_vol(est_coef,x_oos,r_t_oos)

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





