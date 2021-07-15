


cv.loop_lasso <- function(blocks, blocks_cov, points, points_cov, spp, det, nsurv, k, ni, nt, nb, nc, lambda){
  nobs.block <- nrow(block.test)
  nk <- 1/k
  leavout <- sample(rep(1:(nobs.block/(nobs.block*nk)), ceiling(nobs.block*nk))[1:nobs.block])
  
  score <- sapply(1:k, function(j){
    test <- which(leavout == j)
    b.t <- blocks[-test, ]
    b.cv <- rbind(b.t, blocks[test, ])
    
    bcov.cv <- rbind(block_cov[-test, ], block_cov[test, ])
    
    pts.t <- points[which(points$block %in% b.t$blockID), ]
    pts.test <- points[which(!(points$block %in% b.t$blockID)), ]
    pts.cv <- rbind(pts.t, pts.test)
    
    spp_pt.cv <- spp[match(pts.cv$ss, spp$StationKey), ]
    det.cv <- det[match(pts.cv$ss, rownames(det)), , 1:3]
    points_cov.cv <- points_cov[match(pts.cv$ss, rownames(points_cov)), , ]
    nsurv.cv <- nsurv[match(pts.cv$ss, nsurv$StationKey), ]
    
    all.equal(spp_pt.cv$StationKey, rownames(det.cv), rownames(points_cov.cv), pts.cv$ss, nsurv.cv$StationKey)
    
    params <- c("beta.block", "beta.point", "beta.p", "l.score") 
    zst <- apply(spp_pt.cv[, -1], 1, max, na.rm = T)  
    inits <- function() {list(z = zst[1:nrow(pts.t)], beta.block = rnorm(ncol(block_cov.cv)), beta.point = rnorm(ncol(point_cov.cv)), beta.p = rnorm(dim(det.cv)[3]))}    
    
    data <- list(ntrain.block = nrow(b.t), ntrain.site = nrow(pts.t), ntest.block = length(test), 
                 ntest.site = nrow(pts.cv) - nrow(pts.t), y = as.matrix(spp_pt.cv[, -1]), cov.block = block_cov.cv, 
                 cov.point = points_cov.cv, cov.p = det.cv, n.beta.block = ncol(bcov.cv), 
                 n.beta.point = ncol(points_cov.cv), n.beta.p = dim(det.cv)[3], n.surv = nsurv.cv$nsurv, 
                 n.block = nrow(bcov.cv), block = as.numeric(as.factor(rownames(bcov.cv))), 
                 point.block.train = as.numeric(as.factor(pts.t$block)), 
                 point.block.test = as.numeric(as.factor(pts.test$block)), indicator = zst, lambda = lambda)
    
    system.time({
      out <- jags(data = data, inits = inits, parameters.to.save =  params_1, 
                  model.file =  "1_scripts/model_scripts/block_occupancy_1_cv_lasso.txt", 
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
    })
    return(out$mean$l.score)
  })
  return(score)
}

blocks = block.test
blocks_cov = block_cov.test
points = pt.test
points_cov = point_cov
spp = spp_pt
det = det
nsurv = bg_nsurv_veg
k = 5
ni <- 1000; nt <- 1; nb <- 500; nc <- 3 

log.lambda=seq(0.1,5,length=100)
lambda=exp(log.lambda)
lambda = lambda[1]


cv.loop_lasso(blocks = block.test, blocks_cov = block_cov.test, points = pt.test, points_cov = point_cov, spp = spp_pt, det = det, 
        nsurv = bg_nsurv_veg, k = 5, ni = 1000, nt = 1, nb = 500, nc = 3, lambda = lambda)



