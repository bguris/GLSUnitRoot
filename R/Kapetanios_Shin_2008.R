#' Kapetanios and Shin(2008) GLS detrending in KSS(2003) nonlinear unit root test function
#'
#' This function allows you to make Kapetanios and Shin(2008) GLS detrending method for the unit root test procedure developed by Kapetanios, Shin and Snell(2003).
#' @param data_name series name,
#' @param case if demeaned data 1 if detrended data 2,
#' @param lags maximum lag
#' @param lsm lag selection methods if 1 AIC, if 2 BIC, if 3 t-stat significance
#' @return "Model" Estimated model
#' @return "Selected lag" the lag order
#' @return "Test Statistic" the value of the test statistic
#' @references
#' Kapetanios, G., Shin, Y., & Snell, A. (2003). Testing for a unit root in the nonlinear STAR framework. Journal of econometrics, 112(2), 359-379.
#'
#'
#' Burak Guris, R Uygulamalı Dogrusal Olmayan Zaman Serileri Analizi, DER Yayinevi, 2020.
#' @keywords GLS nonlinear unit root test
#' @export
#' @importFrom stats lm coef
#' @importFrom NonlinearTSA KSS_Unit_Root
#' @examples
#'
#' x <- rnorm(1000)
#' Kapetanios_Shin_2008(x, case = 1, lags = 6, lsm =1)
#'
#'
#' y <- cumsum(rnorm(1000))
#' Kapetanios_Shin_2008(y, 2, 12, 2)
#'
#'
#'
Kapetanios_Shin_2008 <- function(data_name,case,lags,lsm){ 
seri_adi = data_name
  
if (case == 1){
  sayi = 0
  a = 1-(sayi/length(seri_adi))
  yf = seri_adi-a*(c(NA,seri_adi[1:(length(seri_adi)-1)]))
  yf[1] = seri_adi[1]
  xf = rep(0,length(seri_adi))
  xf[1:length(seri_adi)] = 1 - a
  xf[1] = 1
  model1 = lm(yf~xf-1)
  hata = seri_adi-coef(model1)
  KSS_Unit_Root(hata,1,lags,lsm)
} else {
  sayi = -17.5
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
  KSS_Unit_Root(hata,1,lags,lsm)
  
}

}