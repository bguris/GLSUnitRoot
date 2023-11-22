#' Su and Nguyen(2013) GLS detrending in Sollis nonlinear unit root test function
#'
#' This function allows you to make Su and Nguyen(2013) GLS detrending method for the unit root test procedure developed by Kruse 2009.
#' @param data_name series name,
#' @param case if demeaned data 1 if detrended data 2,
#' @param lags maximum lag
#' @param lsm lag selection methods if 1 AIC, if 2 BIC, if 3 t-stat significance
#' @return "Model" Estimated model
#' @return "Selected lag" the lag order
#' @return "Test Statistic" the value of the test statistic
#' @references
#' Jen-Je Su & Jeremy K. Nguyen (2013) GLS detrending in Sollis nonlinear unit root tests, Applied Economics Letters, 20:13, 1259-1262
#'
#'
#' Burak Guris, R Uygulamalı Dogrusal Olmayan Zaman Serileri Analizi, DER Yayinevi, 2020.
#' @keywords GLS nonlinear unit root test
#' @export
#' @importFrom stats lm coef
#' @importFrom NonlinearTSA Sollis2009_Unit_Root
#' @examples
#'
#' x <- rnorm(1000)
#' Su_Nguyen_2013(x, case = 2, lags = 3, lsm = 2)
#'
#'
#' y <- cumsum(rnorm(1000))
#' Su_Nguyen_2013(y, 1, 12, 1)
#'
#'
#'
Su_Nguyen_2013 <- function(data_name,case,lags,lsm){

seri_adi = data_name 
if (case == 1){
  sayi = 10
  a = 1-(sayi/length(seri_adi))
  yf = seri_adi-a*(c(NA,seri_adi[1:(length(seri_adi)-1)]))
  yf[1] = seri_adi[1]
  xf = rep(0,length(seri_adi))
  xf[1:length(seri_adi)] = 1 - a
  xf[1] = 1
  model1 = lm(yf~xf-1)
  hata = seri_adi-coef(model1)
  Sollis2009_Unit_Root(hata,1,lags,lsm)
} else {
  sayi = 18
  a = 1-(sayi/length(seri_adi))
  yf = seri_adi-a*(c(NA,seri_adi[1:(length(seri_adi)-1)]))
  yf[1] = seri_adi[1]
  xf = rep(0,length(seri_adi))
  xf[1:length(seri_adi)] = 1 - a
  xf[1] = 1
  as = 0:(length(seri_adi)-1)
  zf = rep(0,length(seri_adi))
  zf[1:length(seri_adi)] = as-a*(c(NA,as[1:(length(as)-1)]))
  zf[1] = 0
  model1 = lm(yf~xf+zf-1)
  hata = seri_adi-coef(model1)[1]-(coef(model1)[2]*as)
  Sollis2009_Unit_Root(hata,1,lags,lsm)
  
}

}