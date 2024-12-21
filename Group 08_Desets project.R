
###############################################################################
# Project:Decomposing and forecasting seasonal time series by DeSeaTS with a  #
# detailed description of the methods, algorithms and forecasting approaches  #
#                                                                             #
###############################################################################


library(fracdiff)
library(ggplot2)
library(zoo)
library(deseats)
library(TSA)

library(tidyr)



trend_ols <- function(yt, order_poly = 1) {
  n <- length(yt)
  t <- 1:n
  df <- data.frame(
    Obs = yt,
    Time = t
  )
  form <- as.formula(
    paste0("Obs ~ ", paste0("I(t^", 1:order_poly, ")", collapse = " + "))
  )
  est <- lm(form, data = df)
  est
}

S_trigo <- function(yt) {
  P <- frequency(yt)
  n <- length(yt)
  t <- 1:n
  trigo_vars <- vector(mode = "list", length = P - 1)
  lambda1 <- 2 * pi / P
  for (i in 1:(trunc((P - 1) / 2))) {
    trigo_vars[[2 * i - 1]] <- cos(i * lambda1 * t)
    trigo_vars[[2 * i]] <- sin(i * lambda1 * t)
  }
  if (P %% 2 == 0) {
    trigo_vars[[P - 1]] <- cos(pi * t)
  }
  trigo_vars[["yt"]] <- c(yt)
  df <- as.data.frame(trigo_vars)
  est <- lm(yt ~ . - 1, data = df)
  est
}
  


data <- read.csv("RSNSRN.csv") 
dim(data)

sales <- ts(data$RSNSRN, start = c(1992, 1), freq = 12)

plot_fig_1 <- autoplot.zoo(sales) + 
  xlab("Year")+
  ylab("Millions of Dollars Sales  ")+
  ggtitle(" Advance Retail Sales: Nonstore Retailers in USA from 1992 to 2024")+
  xlim(1992, 2024+10/12)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

plot_fig_1




sales_log <- ts(log(data$RSNSRN), start =  c(1992, 1), freq = 12)

plot_fig_2 <- autoplot.zoo(sales_log) + 
  xlab("Year")+
  ylab("Log-transformed Millions of dollar")+
  ggtitle("Log-transformed Nonstore Retailers Sales in USA from 1992 to 2024")+
  xlim(1992, 2024+10/12)+theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

plot_fig_2

#############################################################################

periodogram(sales, main = "Periodogram of the US Nonstore Retail Sales ")

#############################################################################


est_OLS <- trend_ols(c(sales_log), order_poly = 3)
trend_OLS <- est_OLS$fitted.values
detrended_OLS <- sales_log - trend_OLS
season_trigo <- S_trigo(detrended_OLS)$fitted.values
TS <- cbind(
  "Observations" = sales_log, 
  "Trend" = trend_OLS, 
  "Seasonality" = season_trigo + mean(sales_log)
)
autoplot.zoo(TS, facets = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  xlab("Year") +
  ylab("Log Transformed Millions of dollar") +
  ggtitle("Log-sales with cubic trend and (shifted) trigonometric seasonality") +
  scale_color_manual(name = "Series", values = c("grey60", "red", "blue"))



##################################################################################



smoothing_options <- set_options(order_poly = 3)
est_deseats <- deseats(sales_log, smoothing_options)

bwidth(est_deseats)

deseats::autoplot(est_deseats, which = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  xlab("Year") +
  ylab("Log of Millions of dollar") +
  ggtitle("Log-transformed Nonstore Retail Sales with trend and seasonality according to DeSeaTS")

###############################################################################



# Deseasonalize using DeSeaTS
deseason_deseats <- exp(deseasonalize(est_deseats))

# Prepare data for plotting
deseason_deseats_data <- data.frame(
  "Year" = c(time(sales)),
  "Observations" = c(sales),
  "DeSeaTS" = c(deseason_deseats)
)

deseason_deseats_long <- pivot_longer(deseason_deseats_data, 
                                      cols = "DeSeaTS", 
                                      names_to = "Method", 
                                      values_to = "Value")

# Plotting for DeSeaTS
plot_deseats <- ggplot(deseason_deseats_long, aes(x = Year, y = Observations)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  geom_line(aes(color = "1")) +
  geom_line(aes(x = Year, y = Value, color = "2"), linewidth = 0.7) +
  scale_color_manual(name = "Series", values = c("1" = "grey50", "2" = "red"),
                     labels = c("1" = "NSA", "2" = "SA")) +
  ylab("Millions of USD") +
  ggtitle("Observed US-Retail sales with DeSeaTS Seasonally Adjusted Version")

plot_deseats




###############################################################################


library(tidyr)

deseason_deseats <- exp(deseasonalize(est_deseats))
deseason_deseats_1 <- data.frame(
  "Year" = c(time(sales)),
  "Observations" = c(sales),
  "DeSeaTS" = c(deseason_deseats)
)

parOld <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), cex = 1, mar = c(4, 4, 3, 2) + 0.1)
periodogram(deseason_deseats_2$DeSeaTS, main = "(A) DeSeaTS")
par(parOld)


res_deseats <- residuals(est_deseats)
deseats <- cbind(
  "DeSeaTS" = res_deseats
)
plot_statio <- autoplot.zoo(res_deseats) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ylab("Log of millions of USD") +
  ggtitle("Stationarized Log-transformed Advance Retail Sales: Nonstore Retailers according to deseats methods") +
  facet_free()
plot_statio

parOld <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), cex = 1, mar = c(4, 4, 3, 2) + 0.1)

periodogram(res_all[, "DeSeaTS"], main = "(B) DeSeaTS")
par(parOld)





###############################################################################


library(forecast)
library(tseries)

acf_res_deseats <- ggAcf(c(residuals(est_deseats))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("ACF of the DeSeaTS residuals (additive model)")
acf_res_deseats

model <- s_semiarma(sales_log, set_options(order_poly = 3))
model@par_model

jarque.bera.test(model@par_model$residuals)

set.seed(123)
fc <- predict(model, n.ahead = 24, method = "boot")
fc@pred
fc@interv

deseats::autoplot(fc) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  xlab("Year") +
  ylab("Log of millions of dollar") +
  ggtitle("Log-Nonstore Retail Sales with point and interval forecasts") +
  coord_cartesian(xlim = c(2015, 2027), ylim =c(10.5, 12.5))

fc_retransf <- expo(fc)

fc_retransf@pred
fc_retransf@interv

deseats::autoplot(fc_retransf) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  xlab("Year") +
  ylab("Millions of dollar ") +
  ggtitle("Nonstore Retail Sales with point and interval forecasts") +
  coord_cartesian(xlim = c(2015, 2027), ylim =c(8701, 211148))



################################ END #########################################

