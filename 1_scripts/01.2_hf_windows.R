
rm(list=ls(all=TRUE)) 
gc()

library(tidyr)
library(stringr)
library(data.table)
library(lubridate)
library(maptools)
library(dplyr)
library(sqldf)
library(reshape2)


load("0_data/processed/big-grid-veg-point.Rdata")
load("0_data/processed/big-grid-point-data.Rdata")
load("0_data/processed/big-grid-veghf-grids.Rdata")


veghf_150_2016 <- read.csv(file = "0_data/raw/10_vegHF150mscale_2016_stationgridassignment.csv", header = TRUE)
veghf_150_2017 <- read.csv(file = "0_data/raw/10_vegHF150mscale_2017_stationgridassignment.csv", header = TRUE)


veghf_300_2016 <- read.csv(file = "0_data/raw/10_vegHF36000m2scale_2016_stationgridassignment.csv", header = TRUE)
veghf_300_2017 <- read.csv(file = "0_data/raw/10_vegHF36000m2scale_2017_stationgridassignment.csv", header = TRUE)


gridnum <- sqldf("SELECT b.StationKey, b.Gridnum, a.g2x2, a.g3x3, a.g4x4, a.g5x5 FROM veghf_300_2016 a JOIN bg_ss_veg b ON a.SS = b.SS")


statnum <- data.frame(ss = gridnum$StationKey, grid = gridnum$Gridnum, point = sapply(1:nrow(gridnum), function(x) as.numeric(substr(gridnum$StationKey[x], 9, 11)))) %>% arrange(grid, point)


sum.point <- function(vars, veg_ss, rad, name){
  p <- list()
  for(i in 1:length(vars)){
    r <- list()
    t <- vars[i] 
    s <- 1
    for(j in 1:ncol(veg_ss)){
      c <- colnames(veg_ss)[j] 
      if(substr(c, 1, nchar(t)) == vars[i]){
        r[[(s)]] <- j
      }else{
        next
      }
      s <- s + 1
    }
    p[[i]] <- r
  }
  
  mat <- matrix(0, nrow(veg_ss), length(vars))
  for(i in 1:ncol(mat)){
    d <- as.matrix(veg_ss[, unlist(p[[i]])])
    if(ncol(d) == 1){
      mat[, i] <- d/rad
    }else{
      mat[, i] <- rowSums(d)/rad
    }
  }
  colnames(mat) <- vars
  rownames(mat) <- veg_ss$StationKey
  
  prop <- as.matrix(rowSums(mat))
  rownames(prop) <- veg_ss$StationKey
  colnames(prop) <- paste0("prop.", name)
  
  return(prop)
  
}


sum.point.scale <- function(varlist, veg_hf, ss, rad, namelist){
  
  veg_ss <- sqldf("SELECT b.StationKey, a.* FROM veg_hf a JOIN ss b ON a.SS = b.SS")
  
  d.var <- list()
  for(i in 1:length(varlist)){
    d.var[[i]] <- sum.point(varlist[[i]], veg_ss, rad, namelist[[i]])
    
  }
  t <- as.data.frame(do.call(cbind, d.var)) 
  t$total.hf <- rowSums(t)
  t.1 <- t
  for(i in 1:nrow(t)){
    for(j in 1:ncol(t)){
      t.1[i, j] <- ifelse(t[i, j] > 1, 1, t[i, j])
    }
  }  
  return(t.1)
}


seismic <- c("SeismicLineWide", "SeismicLineNarrow")
widelin <- c("Pipeline", "TransmissionLine", "RoadTrailVegetated", "RoadVegetatedVerge", "RailVegetatedVerge", "RoadHardSurface", "RailHardSurface")
wells <- c("WellSite")
industry <- c("IndustrialSiteRural", "MineSite", "BorrowpitsDugoutsSumps", "OtherDisturbedVegetation")


hf150_2016 <- sum.point.scale(varlist = list(seismic, widelin, wells, industry), veg_hf = veghf_150_2016, ss = bg_ss_veg, rad = pi*150^2, namelist = c("seismic", "widelin", "wells", "industry"))
hf150_2017 <- sum.point.scale(varlist = list(seismic, widelin, wells, industry), veg_hf = veghf_150_2017, ss = bg_ss_veg, rad = pi*150^2, namelist = c("seismic", "widelin", "wells", "industry"))

hf300_2016 <- sum.point.scale(varlist = list(seismic, widelin, wells, industry), veg_hf = veghf_300_2016, ss = bg_ss_veg, rad = 600^2, namelist = c("seismic", "widelin", "wells", "industry"))
head(hf300_2016)

hf300_2017 <- sum.point.scale(varlist = list(seismic, widelin, wells, industry), veg_hf = veghf_300_2017, ss = bg_ss_veg, rad = 600^2, namelist = c("seismic", "widelin", "wells", "industry"))
head(hf300_2016)

save(hf150_2016, hf150_2017, hf300_2016, hf300_2017, file = "0_data/processed/big-grid-hf-point.Rdata")


# PUt it together into an analysis package
# Get the ovenbird data to use in matching vegetation years to survey years
oven_point <- bg_pa[, c(1:2, which(colnames(bg_pa) =="OVEN"))]

oven_pt <- oven_point %>% spread(survey, OVEN)
oven_pt <- oven_pt[match(bg_nsurv_veg$StationKey, oven_pt$StationKey), 1:5]
head(oven_pt)
all.equal(oven_pt$StationKey, bg_nsurv_veg$StationKey)


# Get the year in which each survey was done
data.year <- sqldf::sqldf("SELECT DISTINCT StationKey, Year FROM bg_maxyear_veg;")
data.year <- data.year[match(oven_pt$StationKey, data.year$StationKey), ]
all.equal(data.year$StationKey, oven_pt$StationKey)


# Write the function to match the veg  and hf data to the survey years
veg_year <- function(i, v2016, v2017, hf2016, hf2017){
  if(data.year$Year[i] <= 2016){
    v <- data.frame(c(v2016[which(rownames(v2016) == data.year$StationKey[i]), ], hf2016[which(rownames(hf2016) == data.year$StationKey[i]), ]))
  }else{
    v <- data.frame(c(v2017[which(rownames(v2017) == data.year$StationKey[i]), ], hf2017[which(rownames(hf2017) == data.year$StationKey[i]), ]))
  }
  return(unlist(v))
}

veg_year(1, veg_300_2016, veg_300_2017, hf300_2016, hf300_2017)

# Make the 300 meter dataset
veghf_300 <- data.frame(t(sapply(1:nrow(data.year), veg_year, v2016 = veg_300_2016, v2017 = veg_300_2017, hf2016 = hf300_2016, hf2017 = hf300_2017)))
rownames(veghf_300) <- data.year$StationKey
all.equal(rownames(veghf_300), oven_pt$StationKey)


veghf_150 <- data.frame(t(sapply(1:nrow(data.year), veg_year, v2016 = veg_150_2016, v2017 = veg_150_2017, hf2016 = hf150_2016, hf2017 = hf150_2017)))
rownames(veghf_150) <- data.year$StationKey
all.equal(rownames(veghf_150), oven_pt$StationKey)



# Get the survey variables and put into a site x survey x variable array
det.var.raw <- data.frame(bg_maxyear_veg[, c("StationKey", "survey")], int = rep(1, nrow(bg_maxyear_veg)), jday = scale(bg_maxyear_veg$julian), tssr = scale(bg_maxyear_veg$tssr), jday2 = scale(bg_maxyear_veg$julian^2))

det.var.scale <- det.var.raw[-which(det.var.raw$survey > 4), ]

det.var.m <- reshape2::melt(det.var.scale, id = c("StationKey", "survey"))

det.var_veg<- acast(det.var.m, StationKey ~ survey ~ variable)
det.var_veg <- det.var_veg[match(bg_nsurv_veg$StationKey, rownames(det.var_veg)), , ]
all.equal(rownames(det.var_veg), bg_nsurv_veg$StationKey, oven_pt$StationKey, rownames(veghf_300))

all.equal(data.year$StationKey, rownames(det.var_veg))



save(bg_count_veg, bg_maxyear_veg, bg_nsurv_veg, bg_pa, bg_ss_veg, data.year, veghf_300, veghf_150, det.var_veg, file = "0_data/processed/big-grid-point-data.Rdata")


# Check to see that models come out as expected
oven_count <- bg_count_veg[, c(1,2, which(colnames(bg_count_veg) == "OVEN"))]

oven_ct <- oven_count %>% spread(survey, OVEN)
oven_ct <- oven_ct[match(bg_nsurv_veg$StationKey, oven_ct$StationKey), 1:5]
head(oven_ct)
all.equal(oven_ct$StationKey, bg_nsurv_veg$StationKey)


oven.try <- apply(oven_ct[, -1], 1, max, na.rm = TRUE)

summary(glm(oven.try ~ 1, family = "poisson"))
summary(glm(oven.try ~ veghf_300$prop.decid, family = "poisson"))
summary(glm(oven.try ~ veghf_300$wt_mean_age + veghf_300$prop.decid, family = "poisson"))



oven.try2 <- apply(oven_pt[, -1], 1, max, na.rm = TRUE)
summary(glm(oven.try2 ~ 1, family = "binomial"))
summary(glm(oven.try2 ~ veghf_300$prop.decid, family = "binomial"))
summary(glm(oven.try2 ~ veghf_300$wt_mean_age + veghf_300$prop.decid, family = "binomial"))
summary(glm(oven.try2 ~ veghf_300$wt_mean_age + veghf_300$prop.decid + veghf_300$total.hf, family = "binomial"))

cawa_point <- bg_pa[, c(1:2, which(colnames(bg_pa) =="CAWA"))]

cawa_pt <- cawa_point %>% spread(survey, CAWA)
cawa_pt <- cawa_pt[match(bg_nsurv_veg$StationKey, cawa_pt$StationKey), 1:5]
cawa.try <- apply(cawa_pt[, -1], 1, max, na.rm = TRUE)
summary(glm(cawa.try ~ 1, family = "binomial"))
summary(glm(cawa.try ~ veghf_300$prop.decid, family = "binomial"))
summary(glm(cawa.try ~ veghf_300$wt_mean_age + veghf_300$prop.decid, family = "binomial"))
summary(glm(cawa.try ~ veghf_300$wt_mean_age + veghf_300$prop.decid + veghf_300$total.hf, family = "binomial"))

# Function to create a moving window over a sequentially numbered grid of points, where each row is a point and the row number reflects it's position in a matrix numbered vertically from top to bottom and then left to right. 

# Returns a list, where each list item represents a window, and consists of the row numbers for points within that window. 

# x.dim is the x-dimension of the matrix (number of columns)

# y.dim is the y-dimension of the matrix (number of rows)

# win.size is the dimension of the square window (number of rows x number of columns): ex. c(3, 3) for a square window or c(3, 4) for a regtangular window
win.list <- function(x.dim, y.dim, win.size){
  t.y <- y.dim - win.size[2] + 1 # The number of windows fit going down the y dimension
  t.x <- x.dim - win.size[1] + 1 # The number of windows going across the x dimension
  ts <- 1:t.y # Numbered sequence of windows in the y dimension
  tx <- 1:t.x
  
  p <- do.call(c, lapply(1:t.x, function(i){
    s.y <- ((1:win.size[1]) - 1)*y.dim + (tx[i] - 1)*y.dim
    d <- lapply(1:t.y, function(j){
      v <- ts[j]
      u <- unlist(lapply(1:win.size[1], function(x) (v:(v + win.size[2] - 1)) + s.y[x]))
      return(u)
    })
    return(d)
  }))
  return(p)
}


block.size <- list(c(2, 2), c(3, 3), c(4, 4), c(5, 5))

window.blocks <- lapply(1:length(block.size), function(x) win.list(x.dim = 10, y.dim = 10, win.size = block.size[[x]]))

point.block <- lapply(1:length(window.blocks), function(j){
  e <- do.call(rbind, lapply(1:nlevels(as.factor(statnum$grid)), function(i){
    a <- statnum[which(statnum$grid == grids[i]), ]
    
    c <- do.call(rbind, lapply(1:length(window.blocks[[j]]), function(x){
      b <- a[match(window.blocks[[j]][[x]], a$point), ]
      b$block <- paste0(b$grid, "-", x)
      return(na.omit(b))
    }))
    return(c)
  }))
})



veg_hf_sum <- function(veg_origin, veg, hf, ss, window.blocks, rad, statnum, gridsize, point.block){
  veg_ss<- sqldf("SELECT b.StationKey, a.* FROM veg_origin a JOIN ss b ON a.SS = b.SS")
  
  veg_ss <- veg_ss[, 1:which(colnames(veg_ss) == "CCSpruce4")]
  
  
  age <- c(4.5, 10, 20, 40, 60, 80, 100, 120, 140, 160) 
  age.des <- c("R", 1:9) 
  
  p <- list() 
  for(i in 1:length(age.des)){ 
    t <- list() 
    s <- 1 
    for(j in 1:ncol(veg_ss)){
      
      c <- colnames(veg_ss)[j] 
      if(substr(c, nchar(c), nchar(c)) == age.des[i]){
        t[[(s)]] <- j
      }else{
        next
      } 
      s <- s + 1 
    } 
    p[[i]] <- t 
  }
  
  age.mat <- data.frame(matrix(0, nrow(veg_ss), length(age)))
  for(i in 1:ncol(age.mat)){
    age.mat[, i] <- rowSums(veg_ss[, unlist(p[[i]])])
  }
  colnames(age.mat) <- age.des
  rownames(age.mat) <- veg_ss$StationKey
  
  blocklist <- list()
  grids <- as.numeric(levels(as.factor(statnum$grid)))
  for(j in 1:length(window.blocks)){
    
    e <- point.block[[j]]
    
    vblock <- data.frame(ss = e$ss, block = e$block, veg[match(e$ss, rownames(veg)), ])
    
    # Sum land cover types in grid 
    d <- ((vblock[-c(1:3)]*rad) %>% group_by(e$block) %>% summarize_at(colnames(veg[-c(1:3)]), list(sum = sum)))
    
    
    
    colnames(d)[-1] <- substr(colnames(d)[-1], 1, nchar(colnames(d)[-1]) - 4) # Remove "_sum" from column names
    
    
    age.bl <- data.frame(ss = e$ss, age.mat[match(e$ss, rownames(age.mat)), ])
    age.block <- ((age.bl[, -1]) %>% group_by(e$block) %>% summarize_at(colnames(age.bl[, -1]), list(sum = sum)))[, -1]/d$total.for
    wt_mean_age <- as.matrix(age.block) %*% age
    rownames(wt_mean_age) <- d$`e$block`
    colnames(wt_mean_age) <- "wt_mean_age"
    
    d.3 <- data.frame(wt_mean_age, d[, -1]/(rad*gridsize[j]))
    
    hf.bl <- data.frame(ss = e$ss, hf[match(e$ss, rownames(hf)), ])
    d.4 <- ((hf.bl[, -1]) %>% group_by(e$block) %>% summarize_at(colnames(hf.bl[, -1]), list(sum = sum)))
    colnames(d.4)[-1] <- substr(colnames(d.4)[-1], 1, nchar(colnames(d.4)[-1]) - 4) # Remove "_sum" from column names
    d.4.1 <- d.4[, -1]/(rad*gridsize[j])
    rownames(d.4.1) <- d.4$`e$block`
    
    t <- data.frame(blockID = d$`e$block`) 
    t1 <- data.frame(bID = e$block, gnum = e$grid)
    gnum <-  sqldf("SELECT DISTINCT a.blockID, b.gnum FROM t a JOIN t1 b ON a.blockID = b.bID")
    
    blocklist[[j]] <- cbind(gnum, d.3, d.4.1)
    ifelse(all.equal(rownames(d.3), gnum$blockID, rownames(d.4.1), blocklist[[j]]$blockID), next, break)
    
  } 
  return(blocklist)
}



veg_hf_sum <- function(veg_origin, veg, hf, ss, window.blocks, rad, statnum, gridsize)
  

blocks2016 <- veg_hf_sum(veg_origin = veghf_300_2016, veg = veg_300_2016, hf = hf300_2016, window.blocks = window.blocks, statnum = statnum, rad = 600^2, ss = bg_ss_veg, gridsize = c(4, 9, 16, 25), point.block = point.block)
blocks2017 <- veg_hf_sum(veg_origin = veghf_300_2017, veg = veg_300_2017, hf = hf300_2017, window.blocks = window.blocks, statnum = statnum, rad = 600^2, ss = bg_ss_veg, gridsize = c(4, 9, 16, 25), point.block = point.block)



grid_year <- sqldf("SELECT DISTINCT b.Gridnum, MAX(a.Year) grid_year FROM bg_maxyear_veg a JOIN gridnum b ON a.StationKey = b.StationKey GROUP BY Gridnum")
head(grid_year)
length(which(grid_year$`COUNT(Year)` > 1))

a <- blocks2016
b <- blocks2017


grids_veghf <- list()
for(i in 1:length(gridsize)){
  d <- a[[i]]$gnum
  e <- do.call(rbind, lapply(1:length(d), function(x) if(grid_year[which(grid_year$Gridnum == d[x]), "grid_year"] <= 2016){a[[i]][x, ]}else{b[[i]][x, ]}))
  grids_veghf[[i]] <- e
}


save(grids2016, grids2017, grids_veghf, gridID, point.block, file = "0_data/processed/big-grid-veghf-blocks.Rdata")

