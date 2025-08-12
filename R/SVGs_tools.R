#' @include create_senna.R
#'
NULL

#' Generic for regionation-based SVG analysis
#'
#' Defines a generic for identifying or visualizing region-based spatially variable genes (SVGs).
#'
#' @seealso \code{\link[=RegionSVGs,SENNA-method]{RegionSVGs for SENNA}}, 
#'          \code{\link[=RegionSVGs,mSENNA-method]{RegionSVGs for mSENNA}}
#'
#' @importFrom methods setGeneric
#' @param senna A SENNA or mSENNA object
#' @param ... Additional arguments
#' @return A plot showing region-based spatially variable genes (SVGs).
#' @export

setGeneric("RegionSVGs", function(senna, ...) {
  standardGeneric("RegionSVGs")
})



#' SVGs for Regionation Analysis
#'
#' Identify region-based spatially variable genes (SVGs) using linear regression on the C–S distance.  
#' Optionally includes spots filtered by `ActiveIdent()` and can regress out the curve parameter.
#'
#' @importFrom methods setMethod
#' @importFrom stats complete.cases
#' @importFrom parallel mclapply
#' @importFrom stats lm
#' @importFrom tibble tibble
#' @importFrom BiocGenerics rbind
#' @importFrom stats p.adjust
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @param senna A SENNA object
#' @param FDR_level FDR threshold (default = 0.01)
#' @param grad_cutoff Threshold for the absolute value of regression coefficients
#' @param active If `TRUE`, only use spots selected by `ActiveIdent()`
#' @param direction If `1` or `-1`, only include spots in that region. Default is `NULL`.
#' @param regress_out_cc If `TRUE`, regress out the curve parameter `t`
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A SENNA object with regionation SVGs stored in `Gene$R.SVGs`
#' @exportMethod RegionSVGs
#' @seealso \code{\link[=RegionSVGs]{RegionSVGs generic}}, 
#'          \code{\link[=RegionSVGs,mSENNA-method]{RegionSVGs for mSENNA}}


setMethod("RegionSVGs",
          signature(senna = "SENNA"),
          function(senna,
                   FDR_level = 0.01,
                   grad_cutoff = 0.1,
                   active = FALSE,
                   direction = NULL,
                   regress_out_cc = FALSE,
                   cores = 1L){

            coord_mat <- senna@Coord[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]), ]
            gene_mat <- senna@Gene[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]), ]
            if(active){
              coord_mat <- coord_mat[rownames(coord_mat) %in% senna@Coord[["Activated"]][["cells"]], ]
              gene_mat <- gene_mat[rownames(gene_mat) %in% senna@Coord[["Activated"]][["cells"]], ]
            }
            if(!is.null(direction)){
              coord_mat <- coord_mat[coord_mat[["region"]] == direction, ]
              gene_mat <- gene_mat[rownames(gene_mat) %in% rownames(coord_mat), ]
            }

            t_scale = max(coord_mat$t) - min(coord_mat$t)
            t_center = min(coord_mat$t)
            d_scale = max(coord_mat$distance) - min(coord_mat$distance)
            d_center = min(coord_mat$distance)
            coord_mat <- dplyr::mutate(coord_mat,
                                       t = (t - t_center) / t_scale,
                                       distance = (distance - d_center) / d_scale)

            if(regress_out_cc) {

              singular_feat <- names(which(colSums(gene_mat) == 0))
              gene_mat <- gene_mat[, which(colSums(gene_mat) != 0)]
              whole_features <- colnames(gene_mat)

              raw_pval <- parallel::mclapply(whole_features, function(feat){
                by_feature <- summary(stats::lm(gene_mat[ , feat] ~ coord_mat$t + coord_mat$distance))$coefficient
                p_feature <- by_feature[3, 4]
                b_feature <- by_feature[3, 1]
                se_feature <- by_feature[3, 2]
                return(tibble::tibble(Gene = feat,
                                      p = p_feature,
                                      Coefficients = b_feature,
                                      SE_coef = se_feature))
              }, mc.cores = cores)

              raw_pval <- do.call(BiocGenerics::rbind, raw_pval)
            }
            else if(!regress_out_cc) {

              singular_feat <- names(which(colSums(gene_mat) == 0))
              gene_mat <- gene_mat[, which(colSums(gene_mat) != 0)]
              whole_features <- colnames(gene_mat)

              raw_pval <- parallel::mclapply(whole_features, function(feat){
                by_feature <- summary(stats::lm(gene_mat[ , feat] ~ coord_mat$distance))$coefficient
                p_feature <- by_feature[2, 4]
                b_feature <- by_feature[2, 1]
                se_feature <- by_feature[2, 2]
                return(tibble::tibble(Gene = feat,
                                      p = p_feature,
                                      Coefficients = b_feature))
              }, mc.cores = cores)

              raw_pval <- do.call(BiocGenerics::rbind, raw_pval)
            }


            test_feat <- raw_pval %>%
              dplyr::mutate(adj.p = stats::p.adjust(p, method = "BH"))

            if(!is.null(direction)){
              test_feat <- dplyr::mutate(test_feat,
                                         Coefficients = Coefficients * direction)
            }

            test_pos <- subset(test_feat, Coefficients >= 0) %>%
              dplyr::filter(adj.p < FDR_level) %>%
              dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
              dplyr::arrange(log(adj.p)) %>%
              dplyr::pull(Gene)

            test_neg <- subset(test_feat, Coefficients < 0) %>%
              dplyr::filter(adj.p < FDR_level) %>%
              dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
              dplyr::arrange(log(adj.p)) %>%
              dplyr::pull(Gene)

            senna@Gene[["R.SVGs"]] <- list(Report = test_feat,
                                           Variable_gene = list(positive = test_pos,
                                                                negative = test_neg),
                                           Filtered_gene = singular_feat)
            return(senna)
          })



#' SVGs for Regionation Analysis
#'
#' Identify region-based spatially variable genes (SVGs) from multiple SENNA objects.  
#' Performs regression on C–S distance (optionally adjusting for curve parameter and batch) and aggregates results into an `msR` object.
#'
#' @importFrom methods setMethod
#' @importFrom stats complete.cases
#' @importFrom parallel mclapply
#' @importFrom stats lm
#' @importFrom tibble tibble
#' @importFrom BiocGenerics rbind
#' @importFrom stats p.adjust
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @param senna A `mSENNA` object
#' @param FDR_level FDR threshold (default = 0.01)
#' @param grad_cutoff Threshold for the absolute value of regression coefficients
#' @param active If `TRUE`, only use spots selected by `ActiveIdent()`
#' @param direction A numeric vector specifying which region to include per SENNA object (e.g., c(1, -1, NULL))
#' @param regress_out_cc If `TRUE`, regress out the curve parameter `t`
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A `mSENNA` object with updated `msR` slot containing regionation SVG results
#' @exportMethod RegionSVGs
#' @seealso \code{\link[=RegionSVGs]{RegionSVGs generic}}, 
#'          \code{\link[=RegionSVGs,SENNA-method]{RegionSVGs for SENNA}}


setMethod("RegionSVGs",
          signature(senna = "mSENNA"),
          function(senna,
                   FDR_level = 0.01,
                   grad_cutoff = 0.1,
                   active = FALSE,
                   direction = NULL,
                   regress_out_cc = FALSE,
                   cores = 1L){
            slen <- length(senna@SENNA)
            vlen <- 1:slen

            cmat <-
              parallel::mclapply(
                vlen,
                function(v){
                  msen <- senna@SENNA[[v]]
                  coord_mat <- msen@Coord[["Spatial"]][stats::complete.cases(msen@Coord[["Spatial"]]), ]

                  if(active){
                    idx <- rownames(coord_mat) %in% msen@Coord[["Activated"]][["cells"]]
                    coord_mat <- coord_mat[idx, ]
                  }
                  if(!is.null(direction)){
                    coord_mat <- coord_mat[coord_mat[["region"]] == direction[[v]], ]
                    coord_mat <- dplyr::mutate(coord_mat,
                                               distance = distance * region)
                  }

                  t_scale = max(coord_mat$t) - min(coord_mat$t)
                  t_center = min(coord_mat$t)
                  d_scale = max(coord_mat$distance) - min(coord_mat$distance)
                  d_center = min(coord_mat$distance)
                  coord_mat <- dplyr::mutate(coord_mat,
                                             t = (t - t_center) / t_scale,
                                             distance = (distance - d_center) / d_scale,
                                             bat = v)

                  return(coord_mat)
                }, mc.cores = cores)
            cmat <- do.call(rbind, cmat)

            gmat <-
              parallel::mclapply(
                vlen,
                function(v){
                  msen <- senna@SENNA[[v]]
                  gene_mat <- msen@Gene[["Spatial"]][stats::complete.cases(msen@Coord[["Spatial"]]), ]
                  idx <- rownames(dplyr::filter(cmat, bat == v))
                  gene_mat <- gene_mat[idx, ]

                  return(gene_mat)
                }, mc.cores = cores)

            gmat <- do.call(rbind, gmat)

            if(regress_out_cc) {

              singular_feat <- names(which(colSums(gmat) == 0))
              gmat <- gmat[, which(colSums(gmat) != 0)]
              whole_features <- colnames(gmat)

              raw_pval <- parallel::mclapply(
                whole_features,
                function(feat){
                  regdat <- tibble::as_tibble(
                    merge(dplyr::select(gmat, tidyselect::all_of(feat)),
                          dplyr::select(cmat, tidyselect::all_of(c("t", "distance", "bat"))),
                          by = 0))
                  regdat <- dplyr::mutate(
                    regdat,
                    gcount = regdat[[feat]])

                  svgmod <- stats::lm(gcount ~ t + distance + factor(bat),
                                      data = regdat)
                  by_feature <- summary(svgmod)$coefficients
                  p_feature <- by_feature[3, 4]
                  b_feature <- by_feature[3, 1]
                  se_feature <- by_feature[3, 2]
                  return(tibble::tibble(Gene = feat,
                                        p = p_feature,
                                        Coefficients = b_feature,
                                        SE_coef = se_feature))
                }, mc.cores = cores)
              raw_pval <- do.call(BiocGenerics::rbind, raw_pval)
            }
            else if(!regress_out_cc) {

              singular_feat <- names(which(colSums(gmat) == 0))
              gmat <- gmat[, which(colSums(gmat) != 0)]
              whole_features <- colnames(gmat)

              raw_pval <- parallel::mclapply(
                whole_features,
                function(feat){
                  regdat <- tibble::as_tibble(
                    merge(dplyr::select(gmat, tidyselect::all_of(feat)),
                          dplyr::select(cmat, tidyselect::all_of(c("t", "distance", "bat"))),
                          by = 0))
                  regdat <- dplyr::mutate(
                    regdat,
                    gcount = regdat[[feat]])

                  svgmod <- stats::lm(gcount ~ distance + factor(bat),
                                      data = regdat)
                  by_feature <- summary(svgmod)$coefficients
                  p_feature <- by_feature[2, 4]
                  b_feature <- by_feature[2, 1]
                  se_feature <- by_feature[2, 2]
                  return(tibble::tibble(Gene = feat,
                                        p = p_feature,
                                        Coefficients = b_feature,
                                        SE_coef = se_feature))
                }, mc.cores = cores)
              raw_pval <- do.call(BiocGenerics::rbind, raw_pval)
            }

            test_feat <- raw_pval %>%
              dplyr::mutate(adj.p = stats::p.adjust(p, method = "BH"))

            test_pos <- subset(test_feat, Coefficients >= 0) %>%
              dplyr::filter(adj.p < FDR_level) %>%
              dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
              dplyr::arrange(adj.p) %>%
              dplyr::pull(Gene)

            test_neg <- subset(test_feat, Coefficients < 0) %>%
              dplyr::filter(adj.p < FDR_level) %>%
              dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
              dplyr::arrange(adj.p) %>%
              dplyr::pull(Gene)

            res <- methods::new(
              "msR",
              Report = test_feat,
              Variable_gene = list(
                positive = test_pos,
                negative = test_neg
              ),
              Filtered_gene = singular_feat
            )
            
            if(is.null(senna@msR)) senna@msR <- res
            else senna@msR <- c(list(senna@msR), res)
              
            return(senna)
          })



#' Generics for Progression Analysis
#'
#' Define generic function `ProgSVGs()`.
#'
#' @seealso \code{\link[=ProgSVGs,SENNA-method]{ProgSVGs for SENNA}}, 
#'          \code{\link[=ProgSVGs,mSENNA-method]{ProgSVGs for mSENNA}}
#'
#' @importFrom methods setGeneric
#' @param senna A SENNA or mSENNA object
#' @param ... Additional arguments
#' @return A plot with progression SVGs.
#' @export


setGeneric("ProgSVGs", function(senna, ...) {
  standardGeneric("ProgSVGs")
})



#' Internal generic for weighted progression analysis (Gaussian)
#'
#' Define generic function `psvgs_gwei()`. Used internally by `ProgSVGs()`.
#'
#' @importFrom methods setGeneric
#' @param senna A SENNA object
#' @param ... Additional arguments
#' @export

setGeneric("psvgs_gwei", function(senna, ...) {
  standardGeneric("psvgs_gwei")
})



#' Internal generic for weighted progression analysis (Box)
#'
#' Define generic function `psvgs_box()`. Used internally by `ProgSVGs()`.
#'
#' @importFrom methods setGeneric
#' @param senna A SENNA object
#' @param ... Additional arguments
#' @export

setGeneric("psvgs_box", function(senna, ...) {
  standardGeneric("psvgs_box")
})



#' SVGs for Progression Analysis
#'
#' Identify progression-based spatially variable genes (SVGs) using weighted regression on the C–S distance.
#' Supports Gaussian or box-shaped weighting and optional spot selection via `ActiveIdent()`.
#'
#' @importFrom methods setMethod
#' @param senna A SENNA object
#' @param weight `"gaussian"` or `"box"` weighting scheme
#' @param FDR_level FDR threshold
#' @param grad_cutoff Threshold for the absolute value of regression coefficients
#' @param active If `TRUE`, only include cells selected by `ActiveIdent()`
#' @param interval Size of interval for sliding window (used in box or Gaussian weighting)
#' @param intensity Intensity of Gaussian weights
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A SENNA object with progression SVGs stored in `Gene$P.SVGs`
#' @exportMethod ProgSVGs
#' @seealso \code{\link[=ProgSVGs]{ProgSVGs generic}}, 
#'          \code{\link[=ProgSVGs,mSENNA-method]{ProgSVGs for mSENNA}}

setMethod("ProgSVGs",
          signature(senna = "SENNA"),
          function(senna,
                   weight = "box",
                   FDR_level = 0.01,
                   grad_cutoff = 0.1,
                   active = FALSE,
                   interval = sqrt(2),
                   intensity = 1,
                   cores = 1L){

            if(!weight %in% c("gaussian", "box")) {
              stop("weight should be one of the 'gaussian' or 'box'.")
            } else if(weight == "gaussian"){
              res <- psvgs_gwei(senna,
                                FDR_level,
                                grad_cutoff,
                                active,
                                interval,
                                intensity,
                                cores)
            }

            else if(weight == "box"){
              res <- psvgs_box(senna,
                               FDR_level,
                               grad_cutoff,
                               active,
                               interval,
                               cores)
            }

            return(res)})



#' SVGs for Progression Analysis
#'
#' Identify progression-based spatially variable genes (SVGs) across multiple SENNA objects,
#' using weighted regression on the C–S distance. Supports Gaussian or box-shaped weighting schemes
#' and batch correction via factor modeling.
#'
#' @importFrom methods setMethod
#' @param senna A mSENNA object.
#' @param weight `"gaussian"` or `"box"` weighting scheme
#' @param intervals A numeric vector of intervals for progression analysis
#' @param FDR_level FDR threshold (default = 0.01)
#' @param grad_cutoff Threshold for the absolute value of regression coefficients
#' @param intensity Intensity of Gaussian weights
#' @param active If `TRUE`, only include cells selected by `ActiveIdent()`
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A `msR` (multi-SENNA report) object
#' @exportMethod ProgSVGs
#' @seealso \code{\link[=ProgSVGs]{ProgSVGs generic}}, 
#'          \code{\link[=ProgSVGs,SENNA-method]{ProgSVGs for SENNA}}

setMethod("ProgSVGs",
          signature(senna = "mSENNA"),
          function(senna,
                   intervals,
                   weight = "box",
                   FDR_level = 0.01,
                   grad_cutoff = 0.1,
                   intensity = 1,
                   active = FALSE,
                   cores = 1L){

            if(!weight %in% c("gaussian", "box")) {
              stop("weight should be one of the 'gaussian' or 'box'.")
            }

            else if(weight == "gaussian"){
              res <- psvgs_gwei(senna,
                                FDR_level,
                                grad_cutoff,
                                active,
                                intervals,
                                intensity,
                                cores)
            }

            else if(weight == "box"){
              res <- psvgs_box(senna,
                               FDR_level,
                               grad_cutoff,
                               active,
                               intervals,
                               cores)
            }

            if(is.null(senna@msR)) senna@msR <- res
            else senna@msR <- c(list(senna@msR), res)
              
            return(senna)})



#' Progression SVGs with Gaussian Weights
#'
#' Identify progression-related spatially variable genes (SVGs) using weighted linear regression.  
#' Weights are assigned to each spot based on a Gaussian kernel centered at the curve axis.
#'
#' @importFrom methods setMethod
#' @importFrom stats complete.cases
#' @importFrom parallel mclapply
#' @importFrom stats lm
#' @importFrom stats dnorm
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom BiocGenerics rbind
#' @importFrom stats p.adjust
#' @importFrom stats integrate
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @importFrom tidyselect all_of
#' @param senna A SENNA object
#' @param FDR_level False discovery rate threshold for detecting SVGs
#' @param grad_cutoff Threshold for the absolute value of regression coefficients
#' @param active If `TRUE`, only include spots selected by `ActiveIdent()`
#' @param interval Spots with distance greater than this value are excluded from regression
#' @param intensity Bandwidth factor for the Gaussian kernel (higher = flatter weights)
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A SENNA object with progression SVGs stored in `Gene$P.SVGs`
#' @exportMethod psvgs_gwei
#' @seealso \code{\link[=ProgSVGs]{ProgSVGs generic}},
#'          \code{\link[=ProgSVGs,SENNA-method]{ProgSVGs for SENNA}},
#'          \code{\link[=ProgSVGs,mSENNA-method]{ProgSVGs for mSENNA}}

setMethod("psvgs_gwei",
          signature = list(senna = "SENNA"),
          function(senna,
                   FDR_level,
                   grad_cutoff,
                   active,
                   interval,
                   intensity,
                   cores){
            coord_mat <- senna@Coord[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]), ]
            gene_mat <- senna@Gene[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]), ]
            if(active){
              idx <- rownames(coord_mat) %in% senna@Coord[["Activated"]][["cells"]]
              coord_mat <- coord_mat[idx, ]
              gene_mat <- gene_mat[idx, ]
            }
            
            fx <- senna@CurveAxis[["fun"]][["x.coef"]]
            fy <- senna@CurveAxis[["fun"]][["y.coef"]]
            knots <- senna@CurveAxis[["knots"]]
            kn <- senna@CurveAxis[["fun"]][["t"]]

            if("spline" %in% senna@CurveAxis[["type"]]){
              if("trimmed" %in% senna@CurveAxis[["type"]]) {
                ll <- lapply(kn[-length(kn)],
                             function(k){
                               qf <- function(x){
                                 sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                        (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                               }
                               return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                             })
                ll <- c(0, unlist(ll))
                cll <- cumsum(ll)
                tv <- coord_mat[["t"]]
                tprime <- 
                  ifelse(is.nan(tv),
                         tv,
                         {
                           fidx <- floor(tv)
                           fidx[fidx == 0] <- 1L
                           fidx[fidx > length(cll)] <- length(cll)
                           
                           (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                             cll[length(cll)]
                         })
                coord_mat[["t"]] <- tprime
              } else {
                kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                tv <- coord_mat[["t"]]
                tmin <- min(tv); tmax <- max(tv)
                llw <- c(sqrt(((fx[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fx[1,] %*% c(1, 1, 1, 1)))^2 +
                                ((fy[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fy[1,] %*% c(1, 1, 1, 1)))^2))
                lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fx[mk,] %*% c(0, 0, 0, 1)))^2 +
                                ((fy[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fy[mk,] %*% c(0, 0, 0, 1)))^2))
                ll <- lapply(2:(mk-1),
                             function(k){
                               qf <- function(x){
                                 sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                        (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                               }
                               return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                             })
                ll <- c(0, llw, unlist(ll), lhi)
                cll <- cumsum(ll)
                tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                tv <- ifelse(tv > (max(kn) - 1L), 
                             (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                             tv)
                tv <- tv + 1L
                tprime <- {
                  fidx <- floor(tv)
                  fidx[fidx == 0] <- 1L
                  fidx[fidx >= length(cll)] <- length(cll) - 1L
                  (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                    cll[length(cll)]
                }
                coord_mat[["t"]] <- tprime
              }
            } else if("straight" %in% senna@CurveAxis[["type"]]){
              if("trimmed" %in% senna@CurveAxis[["type"]]) {
                ll <- lapply(kn[-length(kn)],
                             function(k){
                               dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                               dist <- sqrt(sum(dist))
                               return(dist)
                             })
                ll <- c(0, unlist(ll))
                cll <- cumsum(ll)
                tv <- coord_mat[["t"]]
                tprime <- 
                  ifelse(is.nan(tv),
                         tv,
                         {
                           fidx <- floor(tv)
                           fidx[fidx == 0] <- 1L
                           fidx[fidx > length(cll)] <- length(cll)
                           
                           (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                             cll[length(cll)]
                         })
                coord_mat[["t"]] <- tprime
              } else {
                kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                tv <- coord_mat[["t"]]
                tmin <- min(tv); tmax <- max(tv)
                llw <- c(sqrt(((fx[1,] %*% c(tmin, 1)) - (fx[1,] %*% c(1, 1)))^2 +
                                ((fy[1,] %*% c(tmin, 1)) - (fy[1,] %*% c(1, 1)))^2))
                lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1), 1)) - (fx[mk,] %*% c(0, 1)))^2 +
                                ((fy[mk,] %*% c((tmax%%1), 1)) - (fy[mk,] %*% c(0, 1)))^2))
                ll <- lapply(2:(mk-2),
                             function(k){
                               dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                               dist <- sqrt(sum(dist))
                               return(dist)
                             })
                ll <- c(0, llw, unlist(ll), lhi)
                cll <- cumsum(ll)
                tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                tv <- ifelse(tv > (max(kn) - 1L), 
                             (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                             tv)
                tv <- tv + 1L
                tprime <- {
                  fidx <- floor(tv)
                  fidx[fidx == 0] <- 1L
                  fidx[fidx >= length(cll)] <- length(cll) - 1L
                  
                  (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                    cll[length(cll)]
                }
                coord_mat[["t"]] <- tprime
              }
            }

            singular_feat <- names(which(colSums(gene_mat) == 0))
            gene_mat <- gene_mat[, which(colSums(gene_mat) != 0)]
            whole_features <- colnames(gene_mat)

            raw_pval <- parallel::mclapply(
              whole_features, function(feat){
                regdat <- tibble::as_tibble(merge(dplyr::select(gene_mat, tidyselect::all_of(feat)), coord_mat, by = 0))
                regdat <- regdat  %>%
                  dplyr::mutate(gcount = regdat[[feat]],
                                gwei = stats::dnorm(distance) / stats::dnorm(0) / intensity)
                svgmod <- stats::lm(gcount ~ t,
                                    weights = gwei,
                                    data = dplyr::filter(regdat, distance <= interval))

                by_feature <- summary(svgmod)$coefficients
                p_feature <- by_feature[2, 4]
                b_feature <- by_feature[2, 1]
                se_feature <- by_feature[2, 2]
                return(tibble::tibble(Gene = feat,
                                      p = p_feature,
                                      Coefficients = b_feature,
                                      SE_coef = se_feature))
              }, mc.cores = cores)

            raw_pval <- do.call(rbind, raw_pval)
            test_feat <- raw_pval %>%
              dplyr::mutate(adj.p = stats::p.adjust(p, method = "BH"))

            test_pos <- subset(test_feat, Coefficients >= 0) %>%
              dplyr::filter(adj.p < FDR_level) %>%
              dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
              dplyr::arrange(adj.p) %>%
              dplyr::pull(Gene)

            test_neg <- subset(test_feat, Coefficients < 0) %>%
              dplyr::filter(adj.p < FDR_level) %>%
              dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
              dplyr::arrange(adj.p) %>%
              dplyr::pull(Gene)

            senna@Gene[["P.SVGs"]] <- list(Report = test_feat,
                                           Variable_gene = list(positive = test_pos,
                                                                negative = test_neg),
                                           Filtered_gene = singular_feat)

            return(senna)
          })




#' Progression SVGs with zero-one weights
#'
#' Identify spatially variable genes (SVGs) along the curve axis using a box-shaped (zero-one) weighting scheme.  
#' Only spots within a specified distance from the curve are included in the regression.
#'
#' @importFrom methods setMethod
#' @importFrom stats complete.cases
#' @importFrom parallel mclapply
#' @importFrom stats lm
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom BiocGenerics rbind
#' @importFrom stats p.adjust
#' @importFrom stats integrate
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @importFrom tidyselect all_of
#' @param senna A SENNA object
#' @param FDR_level False discovery rate threshold
#' @param grad_cutoff Threshold for the absolute value of regression coefficients
#' @param active If `TRUE`, only include spots selected by `ActiveIdent()`
#' @param interval Distance threshold for inclusion (box width)
#' @param cores Number of cores for parallel computation (not supported on Windows)
#' @return A SENNA object with progression SVGs stored in `Gene$P.SVGs`
#' @exportMethod psvgs_box

setMethod("psvgs_box",
          signature = list(senna = "SENNA"),
          function(senna,
                   FDR_level,
                   grad_cutoff,
                   active,
                   interval,
                   cores){
            coord_mat <- senna@Coord[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]), ]
            gene_mat <- senna@Gene[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]), ]
            if(active){
              idx <- rownames(coord_mat) %in% senna@Coord[["Activated"]][["cells"]]
              coord_mat <- coord_mat[idx, ]
              gene_mat <- gene_mat[idx, ]
            }

            fx <- senna@CurveAxis[["fun"]][["x.coef"]]
            fy <- senna@CurveAxis[["fun"]][["y.coef"]]
            knots <- senna@CurveAxis[["knots"]]
            kn <- senna@CurveAxis[["fun"]][["t"]]
            
            if("spline" %in% senna@CurveAxis[["type"]]){
              if("trimmed" %in% senna@CurveAxis[["type"]]) {
                ll <- lapply(kn[-length(kn)],
                             function(k){
                               qf <- function(x){
                                 sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                        (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                               }
                               return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                             })
                ll <- c(0, unlist(ll))
                cll <- cumsum(ll)
                tv <- coord_mat[["t"]]
                tprime <- 
                  ifelse(is.nan(tv),
                         tv,
                         {
                           fidx <- floor(tv)
                           fidx[fidx == 0] <- 1L
                           fidx[fidx > length(cll)] <- length(cll)
                           
                           (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                             cll[length(cll)]
                         })
              } else {
                kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                tv <- coord_mat[["t"]]
                tmin <- min(tv); tmax <- max(tv)
                llw <- c(sqrt(((fx[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fx[1,] %*% c(1, 1, 1, 1)))^2 +
                                ((fy[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fy[1,] %*% c(1, 1, 1, 1)))^2))
                lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fx[mk,] %*% c(0, 0, 0, 1)))^2 +
                                ((fy[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fy[mk,] %*% c(0, 0, 0, 1)))^2))
                ll <- lapply(2:(mk-1),
                             function(k){
                               qf <- function(x){
                                 sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                        (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                               }
                               return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                             })
                ll <- c(0, llw, unlist(ll), lhi)
                cll <- cumsum(ll)
                tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                tv <- ifelse(tv > (max(kn) - 1L), 
                             (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                             tv)
                tv <- tv + 1L
                tprime <- {
                  fidx <- floor(tv)
                  fidx[fidx == 0] <- 1L
                  fidx[fidx >= length(cll)] <- length(cll) - 1L
                  (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                    cll[length(cll)]
                }
              }
            } else if( "straight"%in% senna@CurveAxis[["type"]]){
              if("trimmed" %in% senna@CurveAxis[["type"]]) {
                ll <- lapply(kn[-length(kn)],
                             function(k){
                               dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                               dist <- sqrt(sum(dist))
                               return(dist)
                             })
                ll <- c(0, unlist(ll))
                cll <- cumsum(ll)
                tv <- coord_mat[["t"]]
                tprime <- 
                  ifelse(is.nan(tv),
                         tv,
                         {
                           fidx <- floor(tv)
                           fidx[fidx == 0] <- 1L
                           fidx[fidx > length(cll)] <- length(cll)
                           
                           (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                             cll[length(cll)]
                         })
              } else {
                kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                tv <- coord_mat[["t"]]
                tmin <- min(tv); tmax <- max(tv)
                llw <- c(sqrt(((fx[1,] %*% c(tmin, 1)) - (fx[1,] %*% c(1, 1)))^2 +
                                ((fy[1,] %*% c(tmin, 1)) - (fy[1,] %*% c(1, 1)))^2))
                lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1), 1)) - (fx[mk,] %*% c(0, 1)))^2 +
                                ((fy[mk,] %*% c((tmax%%1), 1)) - (fy[mk,] %*% c(0, 1)))^2))
                ll <- lapply(2:(mk-2),
                             function(k){
                               dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                               dist <- sqrt(sum(dist))
                               return(dist)
                             })
                ll <- c(0, llw, unlist(ll), lhi)
                cll <- cumsum(ll)
                tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                tv <- ifelse(tv > (max(kn) - 1L), 
                             (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                             tv)
                tv <- tv + 1L
                tprime <- {
                  fidx <- floor(tv)
                  fidx[fidx == 0] <- 1L
                  fidx[fidx >= length(cll)] <- length(cll) - 1L
                  
                  (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                    cll[length(cll)]
                }
              }
            }
            
            coord_mat <- dplyr::mutate(coord_mat,
                                       t = tprime)

            singular_feat <- names(which(colSums(gene_mat) == 0))
            gene_mat <- gene_mat[, which(colSums(gene_mat) != 0)]
            whole_features <- colnames(gene_mat)

            raw_pval <- parallel::mclapply(
              whole_features, function(feat){
                regdat <- tibble::as_tibble(merge(dplyr::select(gene_mat, tidyselect::all_of(feat)), coord_mat, by = 0))
                regdat <- regdat  %>%
                  dplyr::mutate(gcount = regdat[[feat]],
                                box = ifelse(distance <= interval, 1, 0))

                svgmod <- stats::lm(gcount ~ t,
                                    data = dplyr::filter(regdat, box == 1))

                by_feature <- summary(svgmod)$coefficients
                p_feature <- by_feature[2, 4]
                b_feature <- by_feature[2, 1]
                se_feature <- by_feature[2, 2]
                return(tibble::tibble(Gene = feat,
                                      p = p_feature,
                                      Coefficients = b_feature,
                                      SE_coef = se_feature))
              }, mc.cores = 1L)

            raw_pval <- do.call(rbind, raw_pval)
            test_feat <- raw_pval %>%
              dplyr::mutate(adj.p = stats::p.adjust(p, method = "BH"))

            test_pos <- subset(test_feat, Coefficients >= 0) %>%
              dplyr::filter(adj.p < FDR_level) %>%
              dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
              dplyr::arrange(adj.p) %>%
              dplyr::pull(Gene)

            test_neg <- subset(test_feat, Coefficients < 0) %>%
              dplyr::filter(adj.p < FDR_level) %>%
              dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
              dplyr::arrange(adj.p) %>%
              dplyr::pull(Gene)

            senna@Gene[["P.SVGs"]] <- list(Report = test_feat,
                                           Variable_gene = list(positive = test_pos,
                                                                negative = test_neg),
                                           Filtered_gene = singular_feat)
            return(senna)
          })



#' Multi-progression SVGs with Gaussian weights
#'
#' Identify spatially variable genes (SVGs) from multiple SENNA objects using Gaussian-weighted regression on the curve axis.  
#' The method assigns a separate regression model for each gene, weighted by a Gaussian kernel centered at each spot.  
#' Distance cutoffs and kernel intensities are defined per SENNA object.
#'
#' @importFrom methods setMethod
#' @importFrom stats complete.cases
#' @importFrom parallel mclapply
#' @importFrom stats lm
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom BiocGenerics rbind
#' @importFrom stats p.adjust
#' @importFrom stats integrate
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @importFrom tidyselect all_of
#'
#' @param senna A \code{mSENNA} object containing multiple SENNA instances.
#' @param FDR_level False discovery rate threshold.
#' @param grad_cutoff Threshold for absolute value of regression coefficients (gradients).
#' @param active If \code{TRUE}, only use activated spots defined by \code{ActiveIdent()}.
#' @param intervals A numeric vector specifying the maximum distance for each SENNA object to include in the regression.
#' @param intensity A numeric value scaling the Gaussian kernel.
#' @param cores Number of CPU cores to use (not supported on Windows).
#'
#' @return An \code{msR} object containing:  
#' \itemize{
#'   \item \code{Report}: Full regression statistics for each gene
#'   \item \code{Variable_gene}: A list of positively and negatively correlated genes
#'   \item \code{Filtered_gene}: Genes excluded from testing due to low expression
#' }
#'
#' @exportMethod psvgs_gwei

setMethod("psvgs_gwei",
          signature = list(senna = "mSENNA"),
          function(senna,
                   FDR_level,
                   grad_cutoff,
                   active,
                   intervals,
                   intensity,
                   cores){

            slen <- length(senna@SENNA)
            vlen <- 1:slen

            cmat <-
              parallel::mclapply(
                vlen,
                function(v){
                  msen <- senna@SENNA[[v]]
                  coord_mat <- msen@Coord[["Spatial"]][stats::complete.cases(msen@Coord[["Spatial"]]), ]

                  if(active){
                    idx <- rownames(coord_mat) %in% msen@Coord[["Activated"]][["cells"]]
                    coord_mat <- coord_mat[idx, ]
                  }
                  fx <- msen@CurveAxis[["fun"]][["x.coef"]]
                  fy <- msen@CurveAxis[["fun"]][["y.coef"]]
                  knots <- msen@CurveAxis[["knots"]]
                  kn <- msen@CurveAxis[["fun"]][["t"]]
                  
                  if("spline" %in% msen@CurveAxis[["type"]]){
                    if("trimmed" %in% msen@CurveAxis[["type"]]) {
                      ll <- lapply(kn[-length(kn)],
                                   function(k){
                                     qf <- function(x){
                                       sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                              (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                                     }
                                     return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                                   })
                      ll <- c(0, unlist(ll))
                      cll <- cumsum(ll)
                      tv <- coord_mat[["t"]]
                      tprime <- 
                        ifelse(is.nan(tv),
                               tv,
                               {
                                 fidx <- floor(tv)
                                 fidx[fidx == 0] <- 1L
                                 fidx[fidx > length(cll)] <- length(cll)
                                 
                                 (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                                   cll[length(cll)]
                               })
                    } else {
                      kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                      tv <- coord_mat[["t"]]
                      tmin <- min(tv); tmax <- max(tv)
                      llw <- c(sqrt(((fx[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fx[1,] %*% c(1, 1, 1, 1)))^2 +
                                      ((fy[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fy[1,] %*% c(1, 1, 1, 1)))^2))
                      lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fx[mk,] %*% c(0, 0, 0, 1)))^2 +
                                      ((fy[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fy[mk,] %*% c(0, 0, 0, 1)))^2))
                      ll <- lapply(2:(mk-1),
                                   function(k){
                                     qf <- function(x){
                                       sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                              (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                                     }
                                     return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                                   })
                      ll <- c(0, llw, unlist(ll), lhi)
                      cll <- cumsum(ll)
                      tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                      tv <- ifelse(tv > (max(kn) - 1L), 
                                   (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                                   tv)
                      tv <- tv + 1L
                      tprime <- {
                        fidx <- floor(tv)
                        fidx[fidx == 0] <- 1L
                        fidx[fidx >= length(cll)] <- length(cll) - 1L
                        (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                          cll[length(cll)]
                      }
                    }
                  } else if( "straight"%in% msen@CurveAxis[["type"]]){
                    if("trimmed" %in% msen@CurveAxis[["type"]]) {
                      ll <- lapply(kn[-length(kn)],
                                   function(k){
                                     dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                                     dist <- sqrt(sum(dist))
                                     return(dist)
                                   })
                      ll <- c(0, unlist(ll))
                      cll <- cumsum(ll)
                      tv <- coord_mat[["t"]]
                      tprime <- 
                        ifelse(is.nan(tv),
                               tv,
                               {
                                 fidx <- floor(tv)
                                 fidx[fidx == 0] <- 1L
                                 fidx[fidx > length(cll)] <- length(cll)
                                 
                                 (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                                   cll[length(cll)]
                               })
                    } else {
                      kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                      tv <- coord_mat[["t"]]
                      tmin <- min(tv); tmax <- max(tv)
                      llw <- c(sqrt(((fx[1,] %*% c(tmin, 1)) - (fx[1,] %*% c(1, 1)))^2 +
                                      ((fy[1,] %*% c(tmin, 1)) - (fy[1,] %*% c(1, 1)))^2))
                      lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1), 1)) - (fx[mk,] %*% c(0, 1)))^2 +
                                      ((fy[mk,] %*% c((tmax%%1), 1)) - (fy[mk,] %*% c(0, 1)))^2))
                      ll <- lapply(2:(mk-2),
                                   function(k){
                                     dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                                     dist <- sqrt(sum(dist))
                                     return(dist)
                                   })
                      ll <- c(0, llw, unlist(ll), lhi)
                      cll <- cumsum(ll)
                      tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                      tv <- ifelse(tv > (max(kn) - 1L), 
                                   (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                                   tv)
                      tv <- tv + 1L
                      tprime <- {
                        fidx <- floor(tv)
                        fidx[fidx == 0] <- 1L
                        fidx[fidx >= length(cll)] <- length(cll) - 1L
                        
                        (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                          cll[length(cll)]
                      }
                    }
                  }
                  
                  coord_mat <- dplyr::mutate(coord_mat,
                                             t = tprime,
                                             bat = v)

                  return(coord_mat)
                }, mc.cores = cores)
            cmat <- do.call(rbind, cmat)

            gmat <-
              parallel::mclapply(
                vlen,
                function(v){
                  msen <- senna@SENNA[[v]]
                  gene_mat <- msen@Gene[["Spatial"]][stats::complete.cases(msen@Coord[["Spatial"]]), ]
                  if(active){
                    idx <- rownames(gene_mat) %in% msen@Coord[["Activated"]][["cells"]]
                    gene_mat <- gene_mat[idx, ]
                  }

                  return(gene_mat)
                }, mc.cores = cores)

            gmat <- do.call(rbind, gmat)


            singular_feat <- names(which(colSums(gmat) == 0))
            gmat <- gmat[, which(colSums(gmat) != 0)]
            whole_features <- colnames(gmat)

            raw_pval <- parallel::mclapply(
              whole_features,
              function(feat){
                regdat <- tibble::as_tibble(
                  merge(dplyr::select(gmat, tidyselect::all_of(feat)),
                        dplyr::select(cmat, tidyselect::all_of(c("t", "distance", "bat"))),
                        by = 0))
                regdat <- dplyr::mutate(
                  regdat,
                  gcount = regdat[[feat]],
                  gwei = stats::dnorm(distance) / stats::dnorm(0) / intensity)
                regdat <- dplyr::filter(regdat, distance <= intervals[bat])

                svgmod <- stats::lm(gcount ~ t + factor(bat),
                                    weights = gwei,
                                    data = regdat)

                by_feature <- summary(svgmod)$coefficients
                p_feature <- by_feature[2, 4]
                b_feature <- by_feature[2, 1]
                se_feature <- by_feature[2, 2]
                return(tibble::tibble(Gene = feat,
                                      p = p_feature,
                                      Coefficients = b_feature,
                                      SE_coef = se_feature))
              },
              mc.cores = cores)


            raw_pval <- do.call(rbind, raw_pval)

            test_feat <- raw_pval %>%
              dplyr::mutate(adj.p = stats::p.adjust(p, method = "BH"))

            test_pos <- subset(test_feat, Coefficients >= 0) %>%
              dplyr::filter(adj.p < FDR_level) %>%
              dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
              dplyr::arrange(adj.p) %>%
              dplyr::pull(Gene)

            test_neg <- subset(test_feat, Coefficients < 0) %>%
              dplyr::filter(adj.p < FDR_level) %>%
              dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
              dplyr::arrange(adj.p) %>%
              dplyr::pull(Gene)

            res <- methods::new(
              "msR",
              Report = test_feat,
              Variable_gene = list(
                positive = test_pos,
                negative = test_neg
              ),
              Filtered_gene = singular_feat
            )

            return(res)
          })



#' Multi-progression SVGs with zero-one weights
#'
#' Identify spatially variable genes (SVGs) along multiple curve axes using linear regression.  
#' Binary (box) weights are applied based on C–S distance thresholds for each axis.  
#' Used internally by `ProgSVGs()`.
#'
#' @importFrom methods setMethod
#' @importFrom stats complete.cases
#' @importFrom parallel mclapply
#' @importFrom stats lm
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom BiocGenerics rbind
#' @importFrom stats p.adjust
#' @importFrom stats integrate
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @importFrom tidyselect all_of
#'
#' @param senna A `mSENNA` object
#' @param FDR_level FDR threshold for significance
#' @param grad_cutoff Threshold for absolute value of regression coefficients
#' @param active If `TRUE`, only include cells selected by `ActiveIdent()`
#' @param intervals A numeric vector of C–S distance thresholds for each axis
#' @param cores Number of cores to use for parallel computation (not supported on Windows)
#'
#' @return An `msR` object containing:
#' \itemize{
#'   \item \code{Report}: A tibble with regression results and adjusted p-values
#'   \item \code{Variable_gene}: A list of positive and negative SVGs
#'   \item \code{Filtered_gene}: Genes excluded due to zero expression
#' }
#'
#' @exportMethod psvgs_box

setMethod("psvgs_box",
          signature = list(senna = "mSENNA"),
          function(senna,
                   FDR_level,
                   grad_cutoff,
                   active,
                   intervals,
                   cores){

            slen <- length(senna@SENNA)
            vlen <- 1:slen

            cmat <-
              parallel::mclapply(
                vlen,
                function(v){
                  msen <- senna@SENNA[[v]]
                  coord_mat <- msen@Coord[["Spatial"]][stats::complete.cases(msen@Coord[["Spatial"]]), ]

                  if(active){
                    idx <- rownames(coord_mat) %in% msen@Coord[["Activated"]][["cells"]]
                    coord_mat <- coord_mat[idx, ]
                  }
                  fx <- msen@CurveAxis[["fun"]][["x.coef"]]
                  fy <- msen@CurveAxis[["fun"]][["y.coef"]]
                  knots <- msen@CurveAxis[["knots"]]
                  kn <- msen@CurveAxis[["fun"]][["t"]]
                  
                  if("spline" %in% msen@CurveAxis[["type"]]){
                    if("trimmed" %in% msen@CurveAxis[["type"]]) {
                      ll <- lapply(kn[-length(kn)],
                                   function(k){
                                     qf <- function(x){
                                       sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                              (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                                     }
                                     return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                                   })
                      ll <- c(0, unlist(ll))
                      cll <- cumsum(ll)
                      tv <- coord_mat[["t"]]
                      tprime <- 
                        ifelse(is.nan(tv),
                               tv,
                               {
                                 fidx <- floor(tv)
                                 fidx[fidx == 0] <- 1L
                                 fidx[fidx > length(cll)] <- length(cll)
                                 
                                 (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                                   cll[length(cll)]
                               })
                    } else {
                      kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                      tv <- coord_mat[["t"]]
                      tmin <- min(tv); tmax <- max(tv)
                      llw <- c(sqrt(((fx[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fx[1,] %*% c(1, 1, 1, 1)))^2 +
                                      ((fy[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fy[1,] %*% c(1, 1, 1, 1)))^2))
                      lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fx[mk,] %*% c(0, 0, 0, 1)))^2 +
                                      ((fy[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fy[mk,] %*% c(0, 0, 0, 1)))^2))
                      ll <- lapply(2:(mk-1),
                                   function(k){
                                     qf <- function(x){
                                       sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                              (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                                     }
                                     return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                                   })
                      ll <- c(0, llw, unlist(ll), lhi)
                      cll <- cumsum(ll)
                      tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                      tv <- ifelse(tv > (max(kn) - 1L), 
                                   (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                                   tv)
                      tv <- tv + 1L
                      tprime <- {
                        fidx <- floor(tv)
                        fidx[fidx == 0] <- 1L
                        fidx[fidx >= length(cll)] <- length(cll) - 1L
                        (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                          cll[length(cll)]
                      }
                    }
                  } else if( "straight"%in% msen@CurveAxis[["type"]]){
                    if("trimmed" %in% msen@CurveAxis[["type"]]) {
                      ll <- lapply(kn[-length(kn)],
                                   function(k){
                                     dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                                     dist <- sqrt(sum(dist))
                                     return(dist)
                                   })
                      ll <- c(0, unlist(ll))
                      cll <- cumsum(ll)
                      tv <- coord_mat[["t"]]
                      tprime <- 
                        ifelse(is.nan(tv),
                               tv,
                               {
                                 fidx <- floor(tv)
                                 fidx[fidx == 0] <- 1L
                                 fidx[fidx > length(cll)] <- length(cll)
                                 
                                 (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                                   cll[length(cll)]
                               })
                    } else {
                      kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                      tv <- coord_mat[["t"]]
                      tmin <- min(tv); tmax <- max(tv)
                      llw <- c(sqrt(((fx[1,] %*% c(tmin, 1)) - (fx[1,] %*% c(1, 1)))^2 +
                                      ((fy[1,] %*% c(tmin, 1)) - (fy[1,] %*% c(1, 1)))^2))
                      lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1), 1)) - (fx[mk,] %*% c(0, 1)))^2 +
                                      ((fy[mk,] %*% c((tmax%%1), 1)) - (fy[mk,] %*% c(0, 1)))^2))
                      ll <- lapply(2:(mk-2),
                                   function(k){
                                     dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                                     dist <- sqrt(sum(dist))
                                     return(dist)
                                   })
                      ll <- c(0, llw, unlist(ll), lhi)
                      cll <- cumsum(ll)
                      tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                      tv <- ifelse(tv > (max(kn) - 1L), 
                                   (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                                   tv)
                      tv <- tv + 1L
                      tprime <- {
                        fidx <- floor(tv)
                        fidx[fidx == 0] <- 1L
                        fidx[fidx >= length(cll)] <- length(cll) - 1L
                        
                        (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                          cll[length(cll)]
                      }
                    }
                  }
                  
                  coord_mat <- dplyr::mutate(coord_mat,
                                             t = tprime,
                                             bat = v)
                  return(coord_mat)
                }, mc.cores = cores)
            cmat <- do.call(rbind, cmat)

            gmat <-
              parallel::mclapply(
                vlen,
                function(v){
                  msen <- senna@SENNA[[v]]
                  gene_mat <- msen@Gene[["Spatial"]][stats::complete.cases(msen@Coord[["Spatial"]]), ]
                  if(active){
                    idx <- rownames(gene_mat) %in% msen@Coord[["Activated"]][["cells"]]
                    gene_mat <- gene_mat[idx, ]
                  }

                  return(gene_mat)
                }, mc.cores = cores)

            gmat <- do.call(rbind, gmat)


            singular_feat <- names(which(colSums(gmat) == 0))
            gmat <- gmat[, which(colSums(gmat) != 0)]
            whole_features <- colnames(gmat)

            raw_pval <- parallel::mclapply(
              whole_features,
              function(feat){
                regdat <- tibble::as_tibble(
                  merge(dplyr::select(gmat, tidyselect::all_of(feat)),
                        dplyr::select(cmat, tidyselect::all_of(c("t", "distance", "bat"))),
                        by = 0))
                regdat <- dplyr::mutate(
                  regdat,
                  gcount = regdat[[feat]])
                regdat <- dplyr::filter(regdat, distance <= intervals[bat])

                svgmod <- stats::lm(gcount ~ t + factor(bat),
                                    data = regdat)

                by_feature <- summary(svgmod)$coefficients
                p_feature <- by_feature[2, 4]
                b_feature <- by_feature[2, 1]
                se_feature <- by_feature[2, 2]
                return(tibble::tibble(Gene = feat,
                                      p = p_feature,
                                      Coefficients = b_feature,
                                      SE_coef = se_feature))
              },
              mc.cores = cores)


            raw_pval <- do.call(rbind, raw_pval)

            test_feat <- raw_pval %>%
              dplyr::mutate(adj.p = stats::p.adjust(p, method = "BH"))

            test_pos <- subset(test_feat, Coefficients >= 0) %>%
              dplyr::filter(adj.p < FDR_level) %>%
              dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
              dplyr::arrange(adj.p) %>%
              
              dplyr::pull(Gene)

            test_neg <- subset(test_feat, Coefficients < 0) %>%
              dplyr::filter(adj.p < FDR_level) %>%
              dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
              dplyr::arrange(adj.p) %>%
              dplyr::pull(Gene)

            res <- methods::new(
              "msR",
              Report = test_feat,
              Variable_gene = list(
                positive = test_pos,
                negative = test_neg
              ),
              Filtered_gene = singular_feat
            )

            return(res)
          })







#' Identifying SVGs with Non-linear Patterns
#'
#' Identify spatially variable genes (SVGs) that exhibit peak-like (non-monotonic) expression patterns  
#' along either progression or regionation axes. Delegates to `midpprog()` or `midpregi()` internally.
#'
#' @importFrom methods setGeneric
#'
#' @param senna A `SENNA` or `mSENNA` object
#' @param model Either `"progression"` (or `"prog"`) or `"regionation"` (or `"regio"`)
#' @param peak_position A numeric value (or vector) indicating the peak location(s)
#' @param gene_list A character vector of gene names to search for peak patterns
#' @param ... Additional arguments passed to `ProgSVGs()` or `RegionSVGs()`
#' @return A `SENNA` or `msR` object with the peak SVGs report
#' @seealso \code{\link[=midpprog]{midpprog}}, \code{\link[=midpregi]{midpregi}}
#' @export


PeakSVGs <- function(senna,
                     model,
                     peak_position = 0.5,
                     gene_list = NULL,
                     ...) {
  if(model %in% c("progression", "prog"))
    senna <- midpprog(senna, pk = peak_position, subj = gene_list, ...)
  else if(model %in% c("regionation", "regio"))
    senna <- midpregi(senna, pk = peak_position, subj = gene_list, ...)
  else stop('`model` should be either `progression` (or `prog`) or `regionation` (or `regio`).')

  return(senna)
}



#' Generics for Progression Peak Detection Model
#'
#' Define the generic function `midpprog()` for identifying genes with peak-like expression  
#' patterns along the progression axis.
#'
#' @importFrom methods setGeneric
#'
#' @param senna A `SENNA` or `mSENNA` object
#' @param ... Additional arguments passed to method implementations
#'
#' @export

setGeneric("midpprog", function(senna, ...) {
  standardGeneric("midpprog")
})



#' Detect SVGs with peak expression at specific positions along the curve
#'
#' Identify spatially variable genes (SVGs) whose expression levels peak at a specified location (`pk`) along the curve axis. 
#' For each peak position, a custom transformation of the curve parameter `t` is applied, and gene-wise regression is performed using either Gaussian or box weights.
#'
#' @importFrom methods setMethod
#' @importFrom stats complete.cases
#' @importFrom parallel mclapply
#' @importFrom stats lm
#' @importFrom stats dnorm
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom BiocGenerics rbind
#' @importFrom stats p.adjust
#' @importFrom stats integrate
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @importFrom tidyselect all_of
#'
#' @param senna A SENNA object
#' @param pk A numeric value (or vector) specifying the peak location along the normalized curve parameter
#' @param subj A character vector of gene names to be tested; if NULL, all genes are used
#' @param weight Type of weighting kernel to apply: `"gaussian"` (default) or `"box"`
#' @param FDR_level False discovery rate threshold used to filter significant genes
#' @param grad_cutoff Threshold on the absolute value of regression coefficients
#' @param active If `TRUE`, restricts analysis to spots selected by `ActiveIdent()`
#' @param interval Radius of the neighborhood to consider in box or Gaussian weighting (used only when `weight = "box"`)
#' @param intensity Scaling factor for Gaussian weights (used only when `weight = "gaussian"`)
#' @param cores Number of CPU cores to use for parallel processing (not supported on Windows)
#'
#' @return A SENNA object with results stored in `Gene$PeakSVGs`
#' @exportMethod midpprog

setMethod("midpprog",
          signature(senna = "SENNA"),
          function(senna,
                   pk,
                   subj,
                   weight = "box",
                   FDR_level = 0.01,
                   grad_cutoff = 0.1,
                   active = FALSE,
                   interval = sqrt(2),
                   intensity = 1,
                   cores = 1L){

            coord_mat <- senna@Coord[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]), ]
            gene_mat <- senna@Gene[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]), ]

            if(active){
              idx <- rownames(coord_mat) %in% senna@Coord[["Activated"]][["cells"]]
              coord_mat <- coord_mat[idx, ]
              gene_mat <- gene_mat[idx, ]
            }

            if(!is.null(subj)) gene_mat <- gene_mat[, subj]

            singular_feat <- names(which(colSums(gene_mat) == 0))
            gene_mat <- gene_mat[, which(colSums(gene_mat) != 0)]
            whole_features <- colnames(gene_mat)
            
            fx <- senna@CurveAxis[["fun"]][["x.coef"]]
            fy <- senna@CurveAxis[["fun"]][["y.coef"]]
            knots <- senna@CurveAxis[["knots"]]
            kn <- senna@CurveAxis[["fun"]][["t"]]

            peaklist <- lapply(pk, function(p){
              if("spline" %in% senna@CurveAxis[["type"]]){
                if("trimmed" %in% senna@CurveAxis[["type"]]) {
                  ll <- lapply(kn[-length(kn)],
                               function(k){
                                 qf <- function(x){
                                   sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                          (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                                 }
                                 return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                               })
                  ll <- c(0, unlist(ll))
                  cll <- cumsum(ll)
                  tv <- coord_mat[["t"]]
                  tprime <- 
                    ifelse(is.nan(tv),
                           tv,
                           {
                             fidx <- floor(tv)
                             fidx[fidx == 0] <- 1L
                             fidx[fidx > length(cll)] <- length(cll)
                             
                             (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                               cll[length(cll)]
                           })
                } else {
                  kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                  tv <- coord_mat[["t"]]
                  tmin <- min(tv); tmax <- max(tv)
                  llw <- c(sqrt(((fx[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fx[1,] %*% c(1, 1, 1, 1)))^2 +
                                  ((fy[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fy[1,] %*% c(1, 1, 1, 1)))^2))
                  lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fx[mk,] %*% c(0, 0, 0, 1)))^2 +
                                  ((fy[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fy[mk,] %*% c(0, 0, 0, 1)))^2))
                  ll <- lapply(2:(mk-1),
                               function(k){
                                 qf <- function(x){
                                   sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                          (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                                 }
                                 return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                               })
                  ll <- c(0, llw, unlist(ll), lhi)
                  cll <- cumsum(ll)
                  tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                  tv <- ifelse(tv > (max(kn) - 1L), 
                               (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                               tv)
                  tv <- tv + 1L
                  tprime <- {
                    fidx <- floor(tv)
                    fidx[fidx == 0] <- 1L
                    fidx[fidx >= length(cll)] <- length(cll) - 1L
                    (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                      cll[length(cll)]
                  }
                }
              } else if( "straight"%in% senna@CurveAxis[["type"]]){
                if("trimmed" %in% senna@CurveAxis[["type"]]) {
                  ll <- lapply(kn[-length(kn)],
                               function(k){
                                 dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                                 dist <- sqrt(sum(dist))
                                 return(dist)
                               })
                  ll <- c(0, unlist(ll))
                  cll <- cumsum(ll)
                  tv <- coord_mat[["t"]]
                  tprime <- 
                    ifelse(is.nan(tv),
                           tv,
                           {
                             fidx <- floor(tv)
                             fidx[fidx == 0] <- 1L
                             fidx[fidx > length(cll)] <- length(cll)
                             
                             (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                               cll[length(cll)]
                           })
                } else {
                  kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                  tv <- coord_mat[["t"]]
                  tmin <- min(tv); tmax <- max(tv)
                  llw <- c(sqrt(((fx[1,] %*% c(tmin, 1)) - (fx[1,] %*% c(1, 1)))^2 +
                                  ((fy[1,] %*% c(tmin, 1)) - (fy[1,] %*% c(1, 1)))^2))
                  lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1), 1)) - (fx[mk,] %*% c(0, 1)))^2 +
                                  ((fy[mk,] %*% c((tmax%%1), 1)) - (fy[mk,] %*% c(0, 1)))^2))
                  ll <- lapply(2:(mk-2),
                               function(k){
                                 dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                                 dist <- sqrt(sum(dist))
                                 return(dist)
                               })
                  ll <- c(0, llw, unlist(ll), lhi)
                  cll <- cumsum(ll)
                  tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                  tv <- ifelse(tv > (max(kn) - 1L), 
                               (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                               tv)
                  tv <- tv + 1L
                  tprime <- {
                    fidx <- floor(tv)
                    fidx[fidx == 0] <- 1L
                    fidx[fidx >= length(cll)] <- length(cll) - 1L
                    
                    (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                      cll[length(cll)]
                  }
                }
              }
              tpk <- ifelse(tprime <= p,
                            tprime,
                            (1 - tprime) / (1 - p))

              coord_mat <- dplyr::mutate(coord_mat,
                                         t = tpk)

              if(weight == "gaussian") wei <- stats::dnorm(coord_mat$distance) / stats::dnorm(0) / intensity
              else if (weight == "box") wei <- ifelse(coord_mat$distance <= interval, 1, 0)

              raw_pval <- parallel::mclapply(
                whole_features, function(feat){
                  regdat <- tibble::as_tibble(merge(dplyr::select(gene_mat, tidyselect::all_of(feat)), coord_mat, by = 0))
                  regdat <- regdat  %>%
                    dplyr::mutate(gcount = regdat[[feat]],
                                  wei = wei)
                  svgmod <- stats::lm(gcount ~ t,
                                      weights = wei,
                                      data = dplyr::filter(regdat, distance <= interval))

                  by_feature <- summary(svgmod)$coefficients
                  p_feature <- by_feature[2, 4]
                  b_feature <- by_feature[2, 1]
                  se_feature <- by_feature[2, 2]
                  return(tibble::tibble(Gene = feat,
                                        p = p_feature,
                                        Coefficients = b_feature,
                                        SE_coef = se_feature))
                }, mc.cores = cores)

              raw_pval <- do.call(rbind, raw_pval)
              test_feat <- raw_pval %>%
                dplyr::mutate(adj.p = stats::p.adjust(p, method = "BH"))

              test_pos <- subset(test_feat, Coefficients >= 0) %>%
                dplyr::filter(adj.p < FDR_level) %>%
                dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
                dplyr::arrange(adj.p) %>%
                dplyr::pull(Gene)

              test_neg <- subset(test_feat, Coefficients < 0) %>%
                dplyr::filter(adj.p < FDR_level) %>%
                dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
                dplyr::arrange(adj.p) %>%
                dplyr::pull(Gene)

              oup <- list(Peak = p,
                          Report = test_feat,
                          Variable_gene = list(positive = test_pos,
                                               negative = test_neg),
                          Filtered_gene = singular_feat)

              return(oup)
            })

            senna@Gene[["PeakSVGs"]] <- list(Reports = peaklist,
                                             type = "Progression")

            return(senna)
            })



#' Detect peak-like SVGs across multiple SENNA curves
#'
#' Identify spatially variable genes (SVGs) that peak at a specified location (`pk`) along multiple curve axes in a `mSENNA` object.
#' The method performs piecewise transformation of the curve parameter `t` and runs weighted regression across batches.
#'
#' @importFrom methods setMethod
#' @importFrom stats complete.cases
#' @importFrom parallel mclapply
#' @importFrom stats lm
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom BiocGenerics rbind
#' @importFrom stats p.adjust
#' @importFrom stats integrate
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @importFrom tidyselect all_of
#'
#' @param senna A `mSENNA` object
#' @param pk A numeric value or vector specifying peak locations along the normalized curve axis
#' @param subj A character vector of gene names to test; if `NULL`, all genes are tested
#' @param intervals A numeric vector of distance thresholds (one per curve) used to filter nearby spots
#' @param weight Type of weight to apply: `"gaussian"` or `"box"`
#' @param FDR_level False discovery rate threshold for selecting SVGs
#' @param grad_cutoff Minimum absolute regression coefficient to retain SVGs
#' @param intensity Scaling factor for Gaussian kernel (used only when `weight = "gaussian"`)
#' @param active If `TRUE`, use only cells activated by `ActiveIdent()`
#' @param cores Number of cores to use for parallel computation (not supported on Windows)
#'
#' @return A `mSENNA` object with peak SVG reports stored in `@msR`
#' @exportMethod midpprog

#'
setMethod("midpprog",
          signature(senna = "mSENNA"),
          function(senna,
                   pk,
                   subj,
                   intervals,
                   weight = "box",
                   FDR_level = 0.01,
                   grad_cutoff = 0.1,
                   intensity = 1,
                   active = FALSE,
                   cores = 1L){

            slen <- length(senna@SENNA)
            vlen <- 1:slen

            cmat <-
              parallel::mclapply(
                vlen,
                function(v){
                  msen <- senna@SENNA[[v]]
                  coord_mat <- msen@Coord[["Spatial"]][stats::complete.cases(msen@Coord[["Spatial"]]), ]

                  if(active){
                    idx <- rownames(coord_mat) %in% msen@Coord[["Activated"]][["cells"]]
                    coord_mat <- coord_mat[idx, ]
                  }
                  fx <- msen@CurveAxis[["fun"]][["x.coef"]]
                  fy <- msen@CurveAxis[["fun"]][["y.coef"]]
                  knots <- msen@CurveAxis[["knots"]]
                  kn <- msen@CurveAxis[["fun"]][["t"]]
                  
                  if("spline" %in% msen@CurveAxis[["type"]]){
                    if("trimmed" %in% msen@CurveAxis[["type"]]) {
                      ll <- lapply(kn[-length(kn)],
                                   function(k){
                                     qf <- function(x){
                                       sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                              (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                                     }
                                     return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                                   })
                      ll <- c(0, unlist(ll))
                      cll <- cumsum(ll)
                      tv <- coord_mat[["t"]]
                      tprime <- 
                        ifelse(is.nan(tv),
                               tv,
                               {
                                 fidx <- floor(tv)
                                 fidx[fidx == 0] <- 1L
                                 fidx[fidx > length(cll)] <- length(cll)
                                 
                                 (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                                   cll[length(cll)]
                               })
                    } else {
                      kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                      tv <- coord_mat[["t"]]
                      tmin <- min(tv); tmax <- max(tv)
                      llw <- c(sqrt(((fx[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fx[1,] %*% c(1, 1, 1, 1)))^2 +
                                      ((fy[1,] %*% c(tmin^3, tmin^2, tmin, 1)) - (fy[1,] %*% c(1, 1, 1, 1)))^2))
                      lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fx[mk,] %*% c(0, 0, 0, 1)))^2 +
                                      ((fy[mk,] %*% c((tmax%%1)^3, (tmax%%1)^2, (tmax%%1), 1)) - (fy[mk,] %*% c(0, 0, 0, 1)))^2))
                      ll <- lapply(2:(mk-1),
                                   function(k){
                                     qf <- function(x){
                                       sqrt((3*fx[k,1]*x^2 + 2*fx[k,2]*x+fx[k,3])^2 +
                                              (3*fy[k,1]*x^2 + 2*fy[k,2]*x+fy[k,3])^2)
                                     }
                                     return(stats::integrate(qf, lower = 0, upper = 1)[["value"]])
                                   })
                      ll <- c(0, llw, unlist(ll), lhi)
                      cll <- cumsum(ll)
                      tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                      tv <- ifelse(tv > (max(kn) - 1L), 
                                   (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                                   tv)
                      tv <- tv + 1L
                      tprime <- {
                        fidx <- floor(tv)
                        fidx[fidx == 0] <- 1L
                        fidx[fidx >= length(cll)] <- length(cll) - 1L
                        (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                          cll[length(cll)]
                      }
                    }
                  } else if( "straight"%in% msen@CurveAxis[["type"]]){
                    if("trimmed" %in% msen@CurveAxis[["type"]]) {
                      ll <- lapply(kn[-length(kn)],
                                   function(k){
                                     dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                                     dist <- sqrt(sum(dist))
                                     return(dist)
                                   })
                      ll <- c(0, unlist(ll))
                      cll <- cumsum(ll)
                      tv <- coord_mat[["t"]]
                      tprime <- 
                        ifelse(is.nan(tv),
                               tv,
                               {
                                 fidx <- floor(tv)
                                 fidx[fidx == 0] <- 1L
                                 fidx[fidx > length(cll)] <- length(cll)
                                 
                                 (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                                   cll[length(cll)]
                               })
                    } else {
                      kn <- seq_along(kn); kn <- kn[-length(kn)]; mk <- max(kn)
                      tv <- coord_mat[["t"]]
                      tmin <- min(tv); tmax <- max(tv)
                      llw <- c(sqrt(((fx[1,] %*% c(tmin, 1)) - (fx[1,] %*% c(1, 1)))^2 +
                                      ((fy[1,] %*% c(tmin, 1)) - (fy[1,] %*% c(1, 1)))^2))
                      lhi <- c(sqrt(((fx[mk,] %*% c((tmax%%1), 1)) - (fx[mk,] %*% c(0, 1)))^2 +
                                      ((fy[mk,] %*% c((tmax%%1), 1)) - (fy[mk,] %*% c(0, 1)))^2))
                      ll <- lapply(2:(mk-2),
                                   function(k){
                                     dist <- (knots[k, 1:2] - knots[k+1, 1:2])^2
                                     dist <- sqrt(sum(dist))
                                     return(dist)
                                   })
                      ll <- c(0, llw, unlist(ll), lhi)
                      cll <- cumsum(ll)
                      tv <- ifelse(tv <= 1L, (tv - tmin) / (1 - tmin), tv)
                      tv <- ifelse(tv > (max(kn) - 1L), 
                                   (tv - (max(kn) - 1L)) / (tmax - (max(kn) - 1L)) + (max(kn) - 1L), 
                                   tv)
                      tv <- tv + 1L
                      tprime <- {
                        fidx <- floor(tv)
                        fidx[fidx == 0] <- 1L
                        fidx[fidx >= length(cll)] <- length(cll) - 1L
                        
                        (cll[fidx] + (tv - fidx) * ll[fidx + 1]) / 
                          cll[length(cll)]
                      }
                    }
                  }
                  coord_mat <- dplyr::mutate(coord_mat,
                                             t = tprime,
                                             bat = v)

                  return(coord_mat)
                }, mc.cores = cores)
            cmat <- do.call(rbind, cmat)

            if(weight == "gaussian") wei <- stats::dnorm(cmat$distance) / stats::dnorm(0) / intensity
            else if (weight == "box") wei <- ifelse(cmat$distance <= intervals[cmat$bat], 1, 0)

            gmat <-
              parallel::mclapply(
                vlen,
                function(v){
                  msen <- senna@SENNA[[v]]
                  gene_mat <- msen@Gene[["Spatial"]][stats::complete.cases(msen@Coord[["Spatial"]]), ]
                  if(active){
                    idx <- rownames(gene_mat) %in% msen@Coord[["Activated"]][["cells"]]
                    gene_mat <- gene_mat[idx, ]
                  }

                  return(gene_mat)
                }, mc.cores = cores)

            gmat <- do.call(rbind, gmat)

            if(!is.null(subj)) gmat <- gmat[, subj]

            singular_feat <- names(which(colSums(gmat) == 0))
            gmat <- gmat[, which(colSums(gmat) != 0)]
            whole_features <- colnames(gmat)

            peaklist <- lapply(pk, function(k){
              tpk <- cmat$t
              tpk <- ifelse(tpk <= k,
                            tpk,
                            (1 - tpk) / (1 - k))

              cmat <- dplyr::mutate(cmat, t = tpk)

              raw_pval <- parallel::mclapply(
                whole_features,
                function(feat){
                  regdat <- tibble::as_tibble(
                    merge(dplyr::select(gmat, tidyselect::all_of(feat)),
                          dplyr::select(cmat, tidyselect::all_of(c("t", "distance", "bat"))),
                          by = 0))
                  regdat <- dplyr::mutate(
                    regdat,
                    gcount = regdat[[feat]],
                    gwei = stats::dnorm(distance) / stats::dnorm(0))
                  regdat <- dplyr::filter(regdat, distance <= intervals[bat])

                  svgmod <- stats::lm(gcount ~ t + factor(bat),
                                      weights = gwei,
                                      data = regdat)

                  by_feature <- summary(svgmod)$coefficients
                  p_feature <- by_feature[2, 4]
                  b_feature <- by_feature[2, 1]
                  se_feature <- by_feature[2, 2]
                  return(tibble::tibble(Gene = feat,
                                        p = p_feature,
                                        Coefficients = b_feature,
                                        SE_coef = se_feature))
                },
                mc.cores = cores)


              raw_pval <- do.call(rbind, raw_pval)

              test_feat <- raw_pval %>%
                dplyr::mutate(adj.p = stats::p.adjust(p, method = "BH"))

              test_pos <- subset(test_feat, Coefficients >= 0) %>%
                dplyr::filter(adj.p < FDR_level) %>%
                dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
                dplyr::arrange(adj.p) %>%
                dplyr::pull(Gene)

              test_neg <- subset(test_feat, Coefficients < 0) %>%
                dplyr::filter(adj.p < FDR_level) %>%
                dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
                dplyr::arrange(adj.p) %>%
                dplyr::pull(Gene)

              oup <- methods::new(
                "msR",
                Report = test_feat,
                Variable_gene = list(
                  positive = test_pos,
                  negative = test_neg,
                  peak = k,
                  type = "progression"),
                Filtered_gene = singular_feat)
              
              return(oup)
            })
            
            if(is.null(senna@msR)) senna@msR <- peaklist
            else senna@msR <- c(list(senna@msR), peaklist)
            
            return(senna)
          })



#' Generics for Regionation Peak Detection Model
#'
#' Define generic function `midpregi()`, which identifies genes with peak-like  
#' expression patterns along the regionation axis.
#'
#' @importFrom methods setGeneric
#'
#' @param senna A `SENNA` or `mSENNA` object
#' @param ... Additional arguments
#'
#' @export

setGeneric("midpregi", function(senna, ...) {
  standardGeneric("midpregi")
})



#' Regionation Peak Detection, SENNA-method
#'
#' Identify genes with peak-like spatial patterns along the C–S distance  
#' in regionation analysis, optionally controlling for curve parameter.
#'
#' @importFrom methods setMethod
#' @importFrom stats complete.cases
#' @importFrom parallel mclapply
#' @importFrom stats lm
#' @importFrom tibble tibble
#' @importFrom BiocGenerics rbind
#' @importFrom stats p.adjust
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#'
#' @param senna A `SENNA` object
#' @param pk A numeric value or vector indicating the peak position(s)
#' @param subj A character vector of genes to test
#' @param FDR_level FDR threshold (default = 0.01)
#' @param grad_cutoff Threshold for the absolute value of regression coefficients
#' @param active If `TRUE`, use only spots selected by `ActiveIdent()`
#' @param direction If `1` or `-1`, subset region 1 or 0 respectively. Default is `NULL`.
#' @param regress_out_cc If `TRUE`, regress out the curve parameter `t`
#' @param cores Number of cores for parallel computation (not supported on Windows)
#'
#' @return A `SENNA` object with peak SVGs stored in `Gene$PeakSVGs`
#' @exportMethod midpregi

setMethod("midpregi",
          signature(senna = "SENNA"),
          function(senna,
                   pk,
                   subj,
                   FDR_level = 0.01,
                   grad_cutoff = 0.1,
                   active = FALSE,
                   direction = NULL,
                   regress_out_cc = FALSE,
                   cores = 1L){

            coord_mat <- senna@Coord[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]), ]
            gene_mat <- senna@Gene[["Spatial"]][stats::complete.cases(senna@Coord[["Spatial"]]), ]

            if(active){
              coord_mat <- coord_mat[rownames(coord_mat) %in% senna@Coord[["Activated"]][["cells"]], ]
              gene_mat <- gene_mat[rownames(gene_mat) %in% senna@Coord[["Activated"]][["cells"]], ]
            }

            if(!is.null(direction)){
              coord_mat <- coord_mat[coord_mat[["region"]] == direction, ]
              gene_mat <- gene_mat[rownames(gene_mat) %in% rownames(coord_mat), ]
            }

            singular_feat <- names(which(colSums(gene_mat) == 0))
            gene_mat <- gene_mat[, which(colSums(gene_mat) != 0)]
            whole_features <- colnames(gene_mat)

            t_scale = max(coord_mat$t) - min(coord_mat$t)
            t_center = min(coord_mat$t)
            d_scale = max(coord_mat$distance) - min(coord_mat$distance)
            d_center = min(coord_mat$distance)

            peaklist <- lapply(pk, function(k){
              dpk <- (distance - d_center) / d_scale
              dpk <- ifelse(dpk <= k,
                            dpk / k,
                            (k - dpk) / (1 - k) + 1)

              coord_mat <- dplyr::mutate(coord_mat,
                                         t = (t - t_center) / t_scale,
                                         distance = dpk)

              if(regress_out_cc) {
                singular_feat <- names(which(colSums(gene_mat) == 0))
                gene_mat <- gene_mat[, which(colSums(gene_mat) != 0)]
                whole_features <- colnames(gene_mat)

                raw_pval <- parallel::mclapply(whole_features, function(feat){
                  by_feature <- summary(stats::lm(gene_mat[ , feat] ~ coord_mat$t + coord_mat$distance))$coefficient
                  p_feature <- by_feature[3, 4]
                  b_feature <- by_feature[3, 1]
                  se_feature <- by_feature[3, 2]
                  return(tibble::tibble(Gene = feat,
                                        p = p_feature,
                                        Coefficients = b_feature,
                                        SE_coef = se_feature))
                }, mc.cores = cores)

                raw_pval <- do.call(BiocGenerics::rbind, raw_pval)
              }

              else if(!regress_out_cc) {
                singular_feat <- names(which(colSums(gene_mat) == 0))
                gene_mat <- gene_mat[, which(colSums(gene_mat) != 0)]
                whole_features <- colnames(gene_mat)

                raw_pval <- parallel::mclapply(whole_features, function(feat){
                  by_feature <- summary(stats::lm(gene_mat[ , feat] ~ coord_mat$distance))$coefficient
                  p_feature <- by_feature[2, 4]
                  b_feature <- by_feature[2, 1]
                  se_feature <- by_feature[2, 2]
                  return(tibble::tibble(Gene = feat,
                                        p = p_feature,
                                        Coefficients = b_feature,
                                        SE_coef = se_feature))
                }, mc.cores = cores)

                raw_pval <- do.call(BiocGenerics::rbind, raw_pval)
              }

              test_feat <- raw_pval %>%
                dplyr::mutate(adj.p = stats::p.adjust(p, method = "BH"))
              oup <- list(Peak = k,
                          report = test_feat)
              if(!is.null(direction)){
                test_feat <- dplyr::mutate(test_feat,
                                           Coefficients = Coefficients * direction)}

              test_pos <- subset(test_feat, Coefficients >= 0) %>%
                dplyr::filter(adj.p < FDR_level) %>%
                dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
                dplyr::arrange(log(adj.p)) %>%
                dplyr::pull(Gene)

              test_neg <- subset(test_feat, Coefficients < 0) %>%
                dplyr::filter(adj.p < FDR_level) %>%
                dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
                dplyr::arrange(log(adj.p)) %>%
                dplyr::pull(Gene)

              oup <- list(Peak = k,
                          Report = test_feat,
                          Variable_gene = list(positive = test_pos,
                                               negative = test_neg),
                          Filtered_gene = singular_feat)

              return(oup)
            })

            senna@Gene[["PeakSVGs"]] <- list(Reports = peaklist,
                                             type = "Regionation")

            return(senna)
          })



#' Regionation Peak Detection, mSENNA-method
#'
#' Identify peak-like spatially variable genes (SVGs) in regionation analysis  
#' across multiple curve axes. Supports curve parameter regression and per-axis direction control.
#'
#' @importFrom methods setMethod
#' @importFrom methods new
#' @importFrom stats complete.cases
#' @importFrom parallel mclapply
#' @importFrom stats lm
#' @importFrom tibble tibble
#' @importFrom BiocGenerics rbind
#' @importFrom stats p.adjust
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#'
#' @param senna A `mSENNA` object
#' @param pk A numeric value or vector indicating the peak position(s)
#' @param subj A character vector of genes to test
#' @param FDR_level FDR threshold (default = 0.01)
#' @param grad_cutoff Threshold for the absolute value of regression coefficients
#' @param active If `TRUE`, use only spots selected by `ActiveIdent()`
#' @param direction A numeric vector indicating region direction (1 or -1) per axis. Default is `NULL`.
#' @param regress_out_cc If `TRUE`, regress out the curve parameter `t`
#' @param cores Number of cores for parallel computation (not supported on Windows)
#'
#' @return A `mSENNA` object with peak SVGs stored in `@msR` slot
#' @exportMethod midpregi

setMethod("midpregi",
          signature(senna = "mSENNA"),
          function(senna,
                   pk,
                   subj,
                   FDR_level = 0.01,
                   grad_cutoff = 0.1,
                   active = FALSE,
                   direction = NULL,
                   regress_out_cc = FALSE,
                   cores = 1L){

            slen <- length(senna@SENNA)
            vlen <- 1:slen

            cmat <-
              parallel::mclapply(
                vlen,
                function(v){
                  msen <- senna@SENNA[[v]]
                  coord_mat <- msen@Coord[["Spatial"]][stats::complete.cases(msen@Coord[["Spatial"]]), ]

                  if(active){
                    idx <- rownames(coord_mat) %in% msen@Coord[["Activated"]][["cells"]]
                    coord_mat <- coord_mat[idx, ]
                  }
                  if(!is.null(direction)){
                    coord_mat <- coord_mat[coord_mat[["region"]] == direction[[v]], ]
                    coord_mat <- dplyr::mutate(coord_mat,
                                               distance = distance * region)
                  }

                  t_scale = max(coord_mat$t) - min(coord_mat$t)
                  t_center = min(coord_mat$t)
                  d_scale = max(coord_mat$distance) - min(coord_mat$distance)
                  d_center = min(coord_mat$distance)
                  coord_mat <- dplyr::mutate(coord_mat,
                                             t = (t - t_center) / t_scale,
                                             distance = (distance - d_center) / d_scale,
                                             bat = v)

                  return(coord_mat)
                }, mc.cores = cores)
            cmat <- do.call(rbind, cmat)

            gmat <-
              parallel::mclapply(
                vlen,
                function(v){
                  msen <- senna@SENNA[[v]]
                  gene_mat <- msen@Gene[["Spatial"]][stats::complete.cases(msen@Coord[["Spatial"]]), ]
                  idx <- rownames(dplyr::filter(cmat, bat == v))
                  gene_mat <- gene_mat[idx, ]

                  return(gene_mat)
                }, mc.cores = cores)

            gmat <- do.call(rbind, gmat)

            if(!is.null(subj)) gmat <- gmat[, subj]

            singular_feat <- names(which(colSums(gmat) == 0))
            gene_mat <- gene_mat[, which(colSums(gmat) != 0)]
            whole_features <- colnames(gmat)

            peaklist <- lapply(pk, function(k){
              dpk <- cmat$distance
              dpk <- ifelse(dpk <= k,
                            dpk / k,
                            (k - dpk) / (1 - k) + 1)

              cmat <- dplyr::mutate(cmat, distance = dpk)

              if(regress_out_cc) {
                singular_feat <- names(which(colSums(gmat) == 0))
                gmat <- gmat[, which(colSums(gmat) != 0)]
                whole_features <- colnames(gmat)

                raw_pval <- parallel::mclapply(
                  whole_features,
                  function(feat){
                    regdat <- tibble::as_tibble(
                      merge(dplyr::select(gmat, tidyselect::all_of(feat)),
                            dplyr::select(cmat, tidyselect::all_of(c("t", "distance", "bat"))),
                            by = 0))
                    regdat <- dplyr::mutate(
                      regdat,
                      gcount = regdat[[feat]])

                    svgmod <- stats::lm(gcount ~ t + distance + factor(bat),
                                        data = regdat)
                    by_feature <- summary(svgmod)$coefficients
                    p_feature <- by_feature[3, 4]
                    b_feature <- by_feature[3, 1]
                    se_feature <- by_feature[3, 2]
                    return(tibble::tibble(Gene = feat,
                                          p = p_feature,
                                          Coefficients = b_feature,
                                          SE_coef = se_feature))
                  }, mc.cores = cores)
                raw_pval <- do.call(BiocGenerics::rbind, raw_pval)
              }

              else if(!regress_out_cc) {
                singular_feat <- names(which(colSums(gmat) == 0))
                gmat <- gmat[, which(colSums(gmat) != 0)]
                whole_features <- colnames(gmat)

                raw_pval <- parallel::mclapply(
                  whole_features,
                  function(feat){
                    regdat <- tibble::as_tibble(
                      merge(dplyr::select(gmat, tidyselect::all_of(feat)),
                            dplyr::select(cmat, tidyselect::all_of(c("t", "distance", "bat"))),
                            by = 0))
                    regdat <- dplyr::mutate(
                      regdat,
                      gcount = regdat[[feat]])

                    svgmod <- stats::lm(gcount ~ distance + factor(bat),
                                        data = regdat)
                    by_feature <- summary(svgmod)$coefficients
                    p_feature <- by_feature[2, 4]
                    b_feature <- by_feature[2, 1]
                    se_feature <- by_feature[2, 2]
                    return(tibble::tibble(Gene = feat,
                                          p = p_feature,
                                          Coefficients = b_feature,
                                          SE_coef = se_feature))
                  }, mc.cores = cores)
                raw_pval <- do.call(BiocGenerics::rbind, raw_pval)
              }

              test_feat <- raw_pval %>%
                dplyr::mutate(adj.p = stats::p.adjust(p, method = "BH"))

              test_pos <- subset(test_feat, Coefficients >= 0) %>%
                dplyr::filter(adj.p < FDR_level) %>%
                dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
                dplyr::arrange(adj.p) %>%
                dplyr::pull(Gene)

              test_neg <- subset(test_feat, Coefficients < 0) %>%
                dplyr::filter(adj.p < FDR_level) %>%
                dplyr::filter(abs(Coefficients) >= grad_cutoff) %>%
                dplyr::arrange(adj.p) %>%
                dplyr::pull(Gene)


              oup <- methods::new(
                "msR",
                Report = test_feat,
                Variable_gene = list(
                  positive = test_pos,
                  negative = test_neg,
                  peak = k,
                  type = "regionation"),
                Filtered_gene = singular_feat)
              
            return(oup)})
            
            
            if(is.null(senna@msR)) senna@msR <- peaklist
            else senna@msR <- c(list(senna@msR), peaklist)
            
            return(senna)
          })
