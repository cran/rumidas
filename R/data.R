#' S&P 500 open-to-close daily log-returns
#'
#' Daily data on S&P 500 collected from the realized library of
#' the Oxford-Man Institute \insertCite{heber_2009}{rumidas}. 
#'
#' @docType data
#'
#' @usage data(sp500)
#'
#' @format An object of class \code{"xts"}.
#'
#' @keywords datasets
#'
#' @importFrom Rdpack reprompt
#' @import xts
#' @import highfrequency
#' @import tseries
#' @references
#' \insertAllCited{} 
#'
#' @source Realized library of the \href{https://realized.oxford-man.ox.ac.uk/data/download}{Oxford-Man Institute}
#'
#' @examples
#' head(sp500)
#' summary(sp500)
#' plot(sp500)
"sp500"

#' S&P 500 realized variance at 5-minutes
#'
#' Daily data on the realized variance of the S&P 500 collected from the realized library of
#' the Oxford-Man Institute \insertCite{heber_2009}{rumidas}. The realized variance has been calculated
#' using intradaily intervals of five minutes \insertCite{andersen_boll_1998}{rumidas}.
#'
#' @docType data
#'
#' @usage data(rv5)
#'
#' @format An object of class \code{"xts"}.
#'
#' @keywords datasets
#'
#' @importFrom Rdpack reprompt
#' @import xts
#' @import highfrequency
#' @import tseries
#' @references
#' \insertAllCited{} 
#'
#' @source Realized library of the \href{https://realized.oxford-man.ox.ac.uk/data/download}{Oxford-Man Institute}
#'
#' @examples
#' head(rv5)
#' summary(rv5)
#' plot(rv5)
"rv5"

#' Monthly U.S. Industrial Production
#'
#' Monthly data on the U.S. Industrial Production index (IP, index 2012=100, seasonally adjusted) collected from the 
#' Federal Reserve Economic Data (FRED) archive. The IP has been used as MIDAS term in different contributions  
#' (see, for instance, \insertCite{engle_ghysels_sohn_2013;textual}{rumidas}, \insertCite{conrad_lock_2015;textual}{rumidas}, and
#' \insertCite{amendola_candila_scognamillo_2017;textual}{rumidas}).
#'
#' @docType data
#'
#' @usage data(indpro)
#'
#' @format An object of class \code{"xts"}.
#'
#' @keywords datasets
#'
#' @importFrom Rdpack reprompt
#' @import xts
#' @import highfrequency
#' @import tseries
#' @references
#' \insertAllCited{} 
#'
#' @source Archive of the Federal Reserve Economic Data \href{https://fred.stlouisfed.org/series/INDPRO}{(FRED)}
#'
#' @examples
#' head(indpro)
#' summary(indpro)
#' plot(indpro)
"indpro"

#' VIX daily data 
#'
#' Daily data on VIX collected from Yahoo Finance site. The VIX data 
#' have been de-annualized and multiplied by \eqn{100^{-1}}.
#'
#' @docType data
#'
#' @usage data(vix)
#'
#' @format An object of class \code{"xts"}.
#'
#' @keywords datasets
#'
#' @importFrom Rdpack reprompt
#' @import xts
#' @import highfrequency
#' @import tseries
#'
#' @source \href{https://finance.yahoo.com/}{Yahoo Finance}
#'
#' @examples
#' head(vix)
#' summary(vix)
#' plot(vix)
"vix"


