#' @include create_senna.R
#'
NULL

#' Binary regionation along curve axis
#'
#' Performs binary regionation by assigning each spatial spot to one of two regions based on its relative position to the curve axis (CA).  
#' Also computes the C–S distance (distance from curve axis) for each spot.  
#' Only supports full curves (not trimmed, except for islet curves).
#'
#' @importFrom parallel mclapply
#' @importFrom BiocGenerics rbind
#' @importFrom methods is
#' @param senna A `SENNA` object with a generated curve axis
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A `SENNA` object with updated `Spatial` coordinates including region labels and C–S distance
#' @export

TissueRegionation <- function(senna, cores = 1L){

  if(!is(senna, "SENNA")) stop("INVALID SENNA OBJECT")
  if(any(is.na(senna@CurveAxis[["type"]]))) stop("Curve does not generated yet.")
  if("trimmed" %in% senna@CurveAxis[["type"]]) {
    if(!"islet" %in% senna@CurveAxis[["type"]]){
      stop("To run regionation, use full curve instead of trimmed one.")
    }
  }

  fun <- senna@CurveAxis[["fun"]]
  query <- senna@Coord[["Spatial"]]
  cell_ind <- 1:nrow(query)


  if("spline" %in% senna@CurveAxis[["type"]]) {
    if("islet" %in% senna@CurveAxis[["type"]]) {
      dfun <- trderpolcoef(fun)
      brk <- length(dfun$t)
      exd <- senna@CurveAxis[["exd"]]

      coef1 <- trpol_expan(fun[["x.coef"]])
      coef2 <- trpol_expan(fun[["y.coef"]])
      der1 <- trpol_expan(dfun[["x.deriv"]])
      der2 <- trpol_expan(dfun[["y.deriv"]])

      per_cell <- parallel::mclapply(
        cell_ind, function(i){
          spot <- query[i, ]
          cell <- sign_a_celli(spot, brk, coef1, der1, coef2, der2, exd)
          return(cell)
        }, mc.cores = cores)

      res <- do.call(BiocGenerics::rbind, per_cell)
    }

    else{
      dfun <- derpolcoef(fun)
      brk <- length(dfun$t)


      coef1 <- pol_expan(fun[["x.coef"]])
      coef2 <- pol_expan(fun[["y.coef"]])
      der1 <- pol_expan(dfun[["x.deriv"]])
      der2 <- pol_expan(dfun[["y.deriv"]])


      per_cell <- parallel::mclapply(
        cell_ind, function(i){
          spot <- query[i, ]
          cell <- sign_a_cell(spot, brk, coef1, der1, coef2, der2)
          return(cell)
        }, mc.cores = cores)

      res <- do.call(BiocGenerics::rbind, per_cell)
    }
  }

  else if(senna@CurveAxis[["type"]] == "straight" & nrow(senna@CurveAxis[["knots"]]) != 2) {
    knots <- senna@CurveAxis[["knots"]]
    coef1 <- pol_expan(fun[["x.coef"]])
    coef2 <- pol_expan(fun[["y.coef"]])

    brk <- fun[["t"]]
    brk <- brk[3:(length(brk) - 2)]

    if_flip <- BiocGenerics::cbind(fun[["x.coef"]][2:(nrow(fun[["x.coef"]]) - 1), 1],
                                   fun[["y.coef"]][2:(nrow(fun[["y.coef"]]) - 1), 1])

    fliped <- sapply(X = 1:(nrow(if_flip) - 1),
                     FUN = function(i){
                       val <- ifelse(all(if_flip[i,] %*% if_flip[i+1,] < 0), TRUE, FALSE)
                     })

    ind <- brk[fliped]

    per_cell <- parallel::mclapply(
      cell_ind, function(i){
        spot <- query[i, ]
        cell <- sign_a_cell_lin(spot, brk, coef1, coef2, ind, if_flip, knots)}, mc.cores = cores)

    res <- do.call(BiocGenerics::rbind, per_cell)
  }

  else{
    dfun <- derpolcoef(fun)
    brk <- length(dfun$t)
    
    coef1 <- pol_expan(fun[["x.coef"]])
    coef2 <- pol_expan(fun[["y.coef"]])
    der1 <- pol_expan(dfun[["x.deriv"]])
    der2 <- pol_expan(dfun[["y.deriv"]])


    per_cell <- parallel::mclapply(
      cell_ind, function(i){
        spot <- query[i, ]
        cell <- sign_a_cell(spot, brk, coef1, der1, coef2, der2)
        return(cell)
      }, mc.cores = cores)

    res <- do.call(BiocGenerics::rbind, per_cell)
  }

  senna@Coord[["Spatial"]] <- res
  return(senna)
}



#' Classify spatial spot relative to spline curve
#'
#' Classifies a spatial spot into a binary region (1 or -1) based on its position relative to a spline curve axis (CA).  
#' Uses the sign of the cross product between the curve tangent and the vector from curve to spot.  
#' Also negates the C–S distance if the spot lies in region -1.
#'
#' @importFrom pracma polyval
#' @importFrom pracma ceil
#' @importFrom dplyr mutate
#' @param df A data frame row representing a single spatial spot with curve parameter `t`, `X1`, `X2`
#' @param brk Number of curve segments (internal spline breaks)
#' @param xcoef Polynomial coefficients for x(t)
#' @param xder Derivatives of x(t)
#' @param ycoef Polynomial coefficients for y(t)
#' @param yder Derivatives of y(t)
#' @return The input row with additional columns: `region` (1 or -1) and updated `distance` (signed)
#' @export

sign_a_cell <- function(df, brk, xcoef, xder, ycoef, yder){
  t <- df[["t"]]
  x0 <- df[["X1"]]
  y0 <- df[["X2"]]

  if(t < 1 + 1e-5){
    xt <- pracma::polyval(xcoef[[1]], t)
    yt <- pracma::polyval(ycoef[[1]], t)
    xdt <- pracma::polyval(xder[[1]], t)
    ydt <- pracma::polyval(yder[[1]], t)
  }else if(t >= brk - 1e-5){
    xt <- pracma::polyval(xcoef[[brk + 1]], t)
    yt <- pracma::polyval(ycoef[[brk + 1]], t)
    xdt <- pracma::polyval(xder[[brk + 1]], t)
    ydt <- pracma::polyval(yder[[brk + 1]], t)
  }else {
    if(floor(t) == floor(t + 1e-5)){
      xt <- pracma::polyval(xcoef[[pracma::ceil(t + 1e-5)]], t)
      yt <- pracma::polyval(ycoef[[pracma::ceil(t + 1e-5)]], t)
      xdt <- pracma::polyval(xder[[pracma::ceil(t + 1e-5)]], t)
      ydt <- pracma::polyval(yder[[pracma::ceil(t + 1e-5)]], t)
    } else{
      xt <- pracma::polyval(xcoef[[pracma::ceil(t - 1e-5)]], t)
      yt <- pracma::polyval(ycoef[[pracma::ceil(t - 1e-5)]], t)
      xdt <- pracma::polyval(xder[[pracma::ceil(t - 1e-5)]], t)
      ydt <- pracma::polyval(yder[[pracma::ceil(t - 1e-5)]], t)
    }

  }


  test <- xdt * (yt - y0) - ydt * (xt - x0)

  if(test >= 0){
    df <- dplyr::mutate(df,
                        region = 1L)
  } else {
    df <- dplyr::mutate(df,
                        distance = - distance,
                        region = -1L)
  }

  return(df)
}



#' Classify spot relative to straight curve axis
#'
#' Assigns region label (1 or -1) based on the spot’s position relative to a piecewise linear curve axis.  
#' Uses tangent direction and, if needed, cross product to handle flipped segments.
#'
#' @importFrom pracma polyval
#' @importFrom pracma ceil
#' @importFrom pracma cross
#' @importFrom dplyr mutate
#' @param df A row with spot coordinates and curve parameter
#' @param brk Break points of the curve
#' @param xcoef Coefficients for x(t)
#' @param ycoef Coefficients for y(t)
#' @param ind Indices of segments with tangent direction change
#' @param if_flip Tangent vectors per segment
#' @param knot.df Knot coordinates used for cross-product decision
#' @export

sign_a_cell_lin <- function(df, brk, xcoef, ycoef, ind, if_flip, knot.df){
  t <- df[["t"]]
  x0 <- df[["X1"]]
  y0 <- df[["X2"]]

  if(! t %in% ind){
    if(t < 1) {
      xt <- pracma::polyval(xcoef[[1]], t)
      yt <- pracma::polyval(ycoef[[1]], t)
      xdt <- if_flip[1, 1]
      ydt <- if_flip[1, 2]
    }

    else if(t >= max(brk)) {
      xt <- pracma::polyval(xcoef[[max(brk) + 1]], t)
      yt <- pracma::polyval(ycoef[[max(brk) + 1]], t)
      xdt <- if_flip[nrow(if_flip), 1]
      ydt <- if_flip[nrow(if_flip), 2]
    }

    else {
      xt <- pracma::polyval(xcoef[[pracma::ceil(t)]], t)
      yt <- pracma::polyval(ycoef[[pracma::ceil(t)]], t)
      xdt <- if_flip[floor(t), 1]
      ydt <- if_flip[floor(t), 2]
    }

    test <- xdt * (yt - y0) - ydt * (xt - x0)

    if(test >= 0){
      df <- dplyr::mutate(df,
                          region = 1L)
    } else {
      df <- dplyr::mutate(df,
                          distance = - distance,
                          region = -1L)
    }

    return(df)


  } else {
    v1 <- c(if_flip[(t - 1), ], 0)
    v2 <- c(if_flip[t, ], 0)

    if(sum(pracma::cross(v1, v2)) >= 0){
      df <- dplyr::mutate(df,
                          region = 1L)
    } else {
      df <- dplyr::mutate(df,
                          distance = - distance,
                          region = -1L)
    }

    return(df)

  }
}



#' Classify spot relative to islet curve axis
#'
#' Assigns region label (1 or -1) based on the spot’s position relative to a closed spline curve axis (islet).  
#' For edge segments, uses extrapolated tangent vectors and cross-product logic to determine region.
#'
#' @importFrom pracma polyval
#' @importFrom dplyr mutate
#' @param df A row with spot coordinates and curve parameter
#' @param brk Number of internal curve segments
#' @param xcoef Coefficients for x(t)
#' @param xder Derivatives of x(t)
#' @param ycoef Coefficients for y(t)
#' @param yder Derivatives of y(t)
#' @param exd List containing tangent vectors at start and end of closed curve (`vec1`, `vec2`)
#' @return The input row with region label and possibly signed distance
#' @export

sign_a_celli <- function(df, brk, xcoef, xder, ycoef, yder, exd){
  t <- df[["t"]]
  x0 <- df[["X1"]]
  y0 <- df[["X2"]]

  if(t > 1 & t < brk){
    xt <- pracma::polyval(xcoef[[floor(t)]], t)
    yt <- pracma::polyval(ycoef[[floor(t)]], t)
    xdt <- pracma::polyval(xder[[floor(t)]], t)
    ydt <- pracma::polyval(yder[[floor(t)]], t)

    test <- xdt * (yt - y0) - ydt * (xt - x0)
  } else{
    if(exd[["vec1"]] %*% exd[["vec2"]] > 0){
      xt <- pracma::polyval(xcoef[[1]], 1)
      yt <- pracma::polyval(ycoef[[1]], 1)
      xdt <- pracma::polyval(xder[[1]], 1)
      ydt <- pracma::polyval(yder[[1]], 1)

      test <- xdt * (yt - y0) - ydt * (xt - x0)
    } else{
      xt <- pracma::polyval(xcoef[[1]], 1)
      yt <- pracma::polyval(ycoef[[1]], 1)
      xdt1 <- exd[["vec1"]][[1]]
      ydt1 <- exd[["vec1"]][[2]]
      xdt2 <- exd[["vec2"]][[1]]
      ydt2 <- exd[["vec2"]][[2]]
      
      test1 <- xdt1 * (yt - y0) - ydt1 * (xt - x0)
      test2 <- xdt2 * (yt - y0) - ydt2 * (xt - x0)
      
      v1 <- c(exd[["vec1"]], 0)
      v2 <- c(exd[["vec2"]], 0)
      
      if(sum(pracma::cross(v1, v2)) <= 0){
        if(test1 >= 0 | test2 >= 0) test <- 1L
        else test <- -1L
      } else {
        if(test1 >= 0 & test2 >= 0) test <- 1L
        else test <- -1L
      }
    }}

  if(test >= 0){
    df <- dplyr::mutate(df,
                        region = 1L)

  } else {
    df <- dplyr::mutate(df,
                        distance = - distance,
                        region = -1L)
  }
  return(df)
}


