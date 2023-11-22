#' Elliott, Rothenberg ve Stock(1996) GLS detrending in ADF unit root test function
#'
#' This function allows you to make Elliott, Rothenberg ve Stock(1996) GLS detrending method for the unit root test procedure developed by Dickey and Fuller(1981). 
#' @param data_name series name,
#' @param case if demeaned data 1 if detrended data 2,
#' @param max_lag maximum lag
#' @param lsm lag selection methods if 1 AIC, if 2 BIC
#' @return "Model" Estimated model
#' @return "Selected lag" the lag order
#' @return "Test Statistic" the value of the test statistic
#' @references
#' Elliott, G., T. J. Rothenberg, and J. H. Stock. 1996. Efficient tests for an autoregressive unit root. Econometrica 64 (4):813–36.
#'
#'
#' Burak Guris, R Uygulamalı Dogrusal Olmayan Zaman Serileri Analizi, DER Yayinevi, 2020.
#' @keywords GLS nonlinear unit root test
#' @export
#' @importFrom stats lm coef embed coefficients AIC BIC
#' @importFrom NonlinearTSA Kruse_Unit_Root
#' @examples
#'
#' x <- rnorm(1000)
#' GLS_ADF(x, case = 1, lags = 6, lsm =1)
#'
#'
#' y <- cumsum(rnorm(1000))
#' GLS_ADF(y, 1, 3, 3)
#'
#'
#'

GLS_ADF <- function(data_name,case,max_lag,lsm){
seri_adi = data_name

AICs = NULL
BICs = NULL
if (case == 1){
  for(gecikme in 1:max_lag){
    sayi = 7
    a = 1-(sayi/length(seri_adi))
    yf = seri_adi-a*(c(NA,seri_adi[1:(length(seri_adi)-1)]))
    yf[1] = seri_adi[1]
    xf = rep(0,length(seri_adi))
    xf[1:length(seri_adi)] = 1 - a
    xf[1] = 1
    model1 = lm(yf~xf-1)
    hata = seri_adi-coef(model1)
    gfark = embed(diff(hata),(gecikme+1))[,2:(gecikme+1)]
    modelson = lm(diff(hata)[(gecikme+1):length(diff(hata))]~hata[(gecikme+1):(length(diff(hata)))]+gfark-1)
    modelson0 = lm(diff(hata)~hata[1:(length(diff(hata)))]-1)
    AICs[gecikme + 1] = AIC(modelson)
    BICs[gecikme + 1] = BIC(modelson)
    AICs[1] = AIC(modelson0)
    BICs[1] = BIC(modelson0)
  }
  if (lsm == 1) {
    uygun_lag = which.min(AICs) - 1
  } else if (lsm == 2) {
    uygun_lag = which.min(BICs) - 1
  }
  if (uygun_lag == 0) {
    modeluygunson = lm(diff(hata)~hata[1:(length(diff(hata)))]-1)
    son <- summary(lm(diff(hata)~hata[1:(length(diff(hata)))]-1))$coefficients[1,3]
  } else {
    guygunfark = embed(diff(hata),(uygun_lag+1))[,2:(uygun_lag+1)]
    modeluygunson = lm(diff(hata)[(uygun_lag+1):length(diff(hata))]~hata[(uygun_lag+1):(length(diff(hata)))]+guygunfark-1)
    son <- summary(lm(diff(hata)[(uygun_lag+1):length(diff(hata))]~hata[(uygun_lag+1):(length(diff(hata)))]+guygunfark-1))$coefficients[1,3]
  }
} else {
  for(gecikme in 1:max_lag){
    sayi = 13.5
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
    gfark = embed(diff(hata),(gecikme+1))[,2:(gecikme+1)]
    modelson = lm(diff(hata)[(gecikme+1):length(diff(hata))]~hata[(gecikme+1):(length(diff(hata)))]+gfark-1)
    modelson0 = lm(diff(hata)[(gecikme+1):length(diff(hata))]~hata[(gecikme+1):(length(diff(hata)))]-1)
    AICs[gecikme + 1] = AIC(modelson)
    BICs[gecikme + 1] = BIC(modelson)
    AICs[1] = AIC(modelson0)
    BICs[1] = BIC(modelson0)
  }
  if (lsm == 1) {
    uygun_lag = which.min(AICs) - 1
  } else if (lsm == 2) {
    uygun_lag = which.min(BICs) - 1
  }
  if (uygun_lag == 0) {
    modeluygunson = lm(diff(hata)~hata[1:(length(diff(hata)))]-1)
    son <- summary(lm(diff(hata)~hata[1:(length(diff(hata)))]-1))$coefficients[1,3]
  } else {
    guygunfark = embed(diff(hata),(uygun_lag+1))[,2:(uygun_lag+1)]
    modeluygunson = lm(diff(hata)[(uygun_lag+1):length(diff(hata))]~hata[(uygun_lag+1):(length(diff(hata)))]+guygunfark-1)
    son <- summary(lm(diff(hata)[(uygun_lag+1):length(diff(hata))]~hata[(uygun_lag+1):(length(diff(hata)))]+guygunfark-1))$coefficients[1,3]
    }
}
my_list <- list(Model = summary(modeluygunson), `Selected lag` = uygun_lag, `Test Statistic` = son)
return(my_list)

}
