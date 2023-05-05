
#' Title
#'
#' @param x 
#' @param model 
#' @param fmax 
#' @param fmin 
#'
#' @return
#' @export
#'
#' @examples

EI <- function(x,
               model,
               fmax,
               fmin
){
  if(is.null(nrow(x))) x <- matrix(x, nrow = 1)
  
  pred <- laGP::predGPsep(model,  x, lite=TRUE, nonug = TRUE)
  predmean <- pred$mean
  predsd <- sqrt(pred$s2)
  
  
  
  meas1 <- (predmean-fmax)/predsd
  res1 <- (predmean - fmax)*pnorm(meas1) + predsd *dnorm(meas1)
  
  meas2 <- (fmin-predmean)/predsd
  res2 <- (fmin-predmean)*pnorm(meas2) + predsd *dnorm(meas2)
  
  
  return(res1+res2)
}

#' Title
#'
#' @param X1 
#' @param X2 
#' @param theta 
#' @param correlation 
#'
#' @return
#' @export
#'
#' @examples
covMatrix <- function(X1,
                      X2 = NULL,
                      theta,
                      correlation = "gauss"){
  
  if(is.null(ncol(X1))) X1 <- matrix(X1, ncol=1)
  if(is.null(X2)) X2 <- X1
  d <- ncol(X1)
  
  
  
  if(correlation == "gauss"){
    D <- 0
    
    for(i in 1:d){
      D <- D + laGP::distance(X1[,i], X2[,i])/theta[i]
    }
    out <- exp(-D)
    
  }else{
    stop("correlation has to be either 'gauss', 'matern3_2'/'matern_3_2' or 'matern5_2'/'matern_5_2'")
  }
  
  
  return(out)
}

#' Title
#'
#' @param X 
#' @param y 
#' @param GP 
#' @param f 
#'
#' @return
#' @export
#'
#' @examples
plotGP <- function(X, y, GP, f){
  
  x <- sort(c(seq(0,1, length = 100), X))
  pred <- laGP::predGPsep(GP, matrix(x), lite = TRUE, nonug = TRUE)
  
  q1 <- pred$mean + qnorm(0.05, 0, sqrt(pred$s2))
  q2 <- pred$mean + qnorm(0.95, 0, sqrt(pred$s2))
  fx <- f(x, 0)
  
  plot(x, fx, type = "l", col = "#1d3583",lwd = 2, xlab = "x", ylab = "f(x)", ylim = c(min(min(q1), min(fx))-0.1, max(max(q2), max(fx)) + 0.5))
  
  polygon(c(rev(sort(x)), (sort(x))),
          c(rev(q1), q2),
          col = adjustcolor("#099d02", alpha.f = 0.15), border = NA)
  lines(x, q1)
  lines(x, q2)
  points(X, y, pch = 20, cex = 2)
  points(x, f(x, 0), type = "l", col = "#1d3583",lwd = 2, xlab = "x", ylab = "f(x)")
  lines(x, pred$mean, col = "#099d02", lty = 2, lwd = 2)
  legend("top", legend = c("Actual f(x)", "Fitted values"), lty = 1:2, col = c("#1d3583", "#099d02"), lwd = 2)
}

#' Title
#'
#' @param gpi 
#' @param new.p 
#'
#' @return
#' @export
#'
#' @examples
create_hist <- function(gpi, new.p){
  K<- 10000
  x.test <- matrix(seq(0,1, length.out = 200), ncol = 1)
  pred <- predGPsep(gpi, x.test, nonug = TRUE)
  counts <- xdist(x.test, pred$mean, pred$Sigma,addmean = rep(0, 200), sdobs = 0, maxmin = "both", K = 10000)
  test.x <- as.matrix(expand.grid(x.test, x.test))
  counts.x <- rbind(counts, test.x)
  
  x <- counts.x[,2]
  y <- counts.x[,1]
  
  legendtitle <- list(yref='paper',xref="paper",y=c(1.3, -0.1, 0.5),x=c(1.4, 0.5, -0.1), text=c(TeX("p_*(\\chi_{min}, \\chi_{max})"),TeX("\\chi_{min}"), TeX("\\chi_{max}")),showarrow=F, font = list( size = 17), textangle = c(0,0,270))
  
  s <- subplot(
    plot_ly(x = counts[,2], type = "histogram",  marker = list(color = "#23c686",
                                                               line = list(color = "#099d02",
                                                                           width = 2)), nbinsx = 100) ,
    plotly_empty(type = "scatter", mode = "markers"),
    
    plot_ly(x = x, y = y, type = "histogram2dcontour", colors = hcl.colors(99, "Greens", rev = T), histnorm = "probability") %>%
      layout(legend=list(title=list(text='<b> Trend </b>')),
             xaxis = list(title = TeX("X_{min}"),
                          zerolinecolor = '#ffff',
                          zerolinewidth = 2,
                          gridcolor = 'ffff'),
             yaxis = list(title = TeX("X_{max}"),
                          zerolinecolor = '#ffff',
                          zerolinewidth = 2,
                          gridcolor = 'ffff'),
             annotations=legendtitle
             ,shapes = list( list(
               type = "line",
               y0 = 0,
               y1 = 1,
               yref = "paper",
               x0 = new.p,
               x1 = new.p,
               line = list(color = "red", dash="dot")
             ),list(
               type = "line",
               x0 = 0,
               x1 = 1,
               yref = "paper",
               y0 = new.p,
               y1 = new.p,
               line = list(color = "red", dash="dot")
             ))
      ) %>%
      config(mathjax = "cdn") ,
    
    plot_ly(y = counts[,1], type = "histogram", marker = list(color = "#23c686",
                                                              line = list(color = "#099d02",
                                                                          width = 1))) ,
    
    nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2), margin = 0,
    shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
  )
  fig <- layout(s, showlegend = FALSE) %>%
    config(mathjax = "cdn")
  
  fig
}


#' Title
#'
#' @param gpi 
#' @param xdiscrete 
#' @param theta 
#' @param X 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
create_entropies <- function(gpi, xdiscrete, theta, X, y){
  
  pred <- laGP::predGPsep(gpi, xdiscrete, nonug = TRUE)
  Sigma <- pred$Sigma
  sigma <- NULL
  K <- 2000
  M <- nrow(xdiscrete)
  entropies <- numeric(M)
  direction <- "both"
  for(j in 1:M){
    
    prednew <- GPupdate(X = xdiscrete, X.cur = X, y.cur = y, X.add = xdiscrete[j,], theta = theta, mean.cur = pred$mean, Sigma = Sigma, sigma = sigma)
    
    simres <- as.data.frame(xdist(xdiscrete, pred$mean, prednew$Sigma, prednew$addmean, Sigma[j,j], "both", K))
    colnames(simres) <- paste0("x", 1:2)
    simres <- simres %>%
      dplyr::group_by(.dots = colnames(simres)) %>%
      dplyr::summarise(counts = dplyr::n()/K, .groups = 'drop')
    
    entropies[j] <- -sum(simres$counts*log(simres$counts))
  }
  entropies
}

#' Title
#'
#' @param X 
#' @param X.cur 
#' @param y.cur 
#' @param X.add 
#' @param y.add 
#' @param theta 
#' @param mean.cur 
#' @param Sigma 
#' @param sigma 
#' @param mean.add 
#' @param correlation 
#'
#' @return
#' @export
#'
#' @examples
GPupdate <- function(X,
                     X.cur,
                     y.cur = NULL,
                     X.add,
                     y.add = NULL,
                     theta,
                     mean.cur,
                     Sigma,
                     sigma = NULL,
                     mean.add = NULL,
                     correlation = "gauss"){
  
  if(is.null(nrow(X.add))) X.add <- matrix(X.add, nrow = 1)
  g <- theta[length(theta)]
  
  if(is.null(sigma)){
    if(is.null(y.cur)){
      stop("y.cur must be given if sigma is NULL")
    }else{
      sigma <- drop(t(y.cur) %*%  solve(covMatrix(X.cur, X.cur, theta, correlation) + diag(g, nrow = nrow(X.cur))) %*% y.cur / (length(y.cur)))
    }
    
  }
  
  SxX <-  (covMatrix(X, X.add, theta, correlation) -
             covMatrix(X, X.cur, theta, correlation)%*%
             solve( covMatrix(X.cur, X.cur, theta, correlation)+ diag(g, nrow = nrow(X.cur)) ) %*%
             covMatrix(X.cur, X.add, theta, correlation))
  
  SXX <-  solve( covMatrix(X.add, X.add, theta, correlation)  -
                   covMatrix(X.add, X.cur, theta, correlation)%*%
                   solve(covMatrix(X.cur, X.cur, theta, correlation) + diag(g, nrow = nrow(X.cur)) ) %*%
                   covMatrix(X.cur, X.add, theta, correlation)  + g)
  
  addmean <-  SxX %*% SXX
  newSigma <- Sigma - sigma*(SxX %*%SXX%*% t(SxX))
  res <- list(addmean = addmean, Sigma = newSigma)
  
  if(!is.null(y.add) & !is.null(mean.add)){
    newmean <- mean.cur + addmean *(y.add - mean.add)
    list$mean <- newmean
  }
  return(res)
}
