bohman <- function(t, tau = 0.25){
  boh <- (1-t/tau)*cos(pi*t/tau)+sin(pi*t/tau)/pi
  boh <- ifelse(t>=tau, 0, boh)
}