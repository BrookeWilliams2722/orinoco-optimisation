
library(Matrix)
options(stringsAsFactors = FALSE)

load(file="Data.RData")

# load the coordinates data
load(file="pu.df.RData")


# View(Data)

# Data$Current.land.use.num codes:
# Livestock = 1
# Oil palm = 2
# Forestry = 3
# Rice = 4
# Soy = 5
# Natural = 6
# (Urban = 8)

####set up vectors

# FULL DATA SET
#create sparse matrix from species data 


# remove urban areas:
Data <- Data[-which(Data$Current.land.use.num == 8),]
dim(Data)

# remove planning units with zero for water or carbon
Data <- Data[-which(Data$Carbon == 0 | Data$H2Oyr.1 == 0),]
dim(Data)


ncol(Data)

# extract species dataset and remove species with no cells
startcol <- which(names(Data) == "Bird_Agamia_agami")
endcol <- which(names(Data) == "rep_Micrurus_isozonus")
sppmat <- Data[,c(startcol:endcol)]
# remove species with no cells
cs <- colSums(sppmat)
recs <- which(cs == 0)
if (length(recs > 0)) sppmat <- sppmat[,-recs]


# total current habitat area for each species (each cell is 1 km2, so units km2):
# includes all PUs
sppartot <- colSums(sppmat)

# total area of study area - with no PUs removed
artot <- dim(Data)[1]

# PROTECTED AREAS
currpa.recs <- which(Data$Current.land.use == 'Protected area')
# total current habitat area within protected areas for each species:
spppatot <- colSums(sppmat[currpa.recs,])
# calculate baseline extinction risk - no need to do this, the baseline local ER is 0
# i.e. baseline.er <- 1 - (sppartot / sppartot)^erisk.power

check.nspp <- length(spppatot)

# get the coordinates of the protected areas
if (make.maps){
	if (!exists("pacoords")){
		message("Acquiring protected area coordinates...")
		recs <- rep(0, length(currpa.recs))
		for (i in 1:length(currpa.recs)) recs[i] <- which(pu.df$id == Data$Id[currpa.recs[i]])
		pacoords <- pu.df[recs,c(2, 3)]
	}
}


# remove protected areas from the planning units as they cannot become unprotected
# so there are no decisions to make:
Data <- Data[-currpa.recs,]

# ALLOCATE EXISTING AG CELLS AND REMOVE THEM
prof.existing <- rep(0, 5)
cb.loss.existing <- 0
h2o.loss.existing <- 0

prod.area.targets <- prod.area.targets.t0
keep.existing.n <- 0

if (keep.existing.ag){
	message("Excluding existing agriculture areas...")
	message("Initial production area targets:")
	print(prod.area.targets)
	message(paste0("Total: ", sum(prod.area.targets)))


	recs <- which(Data$Current.land.use.num < 6)

	keep.existing <- list()

	if (keep.target.only){
		for (i in 1:5){
			ids <- Data$Id[which(Data$Current.land.use.num == i)]
			if (length(ids) <= prod.area.targets[i]){
				keep.existing[[i]] <- ids
			} else {
				recs2 <- which(Data$Current.land.use.num == i)
				if (i == 1)	ord <- order(Data$Livestock_profit[recs2], decreasing=TRUE)
				if (i == 2)	ord <- order(Data$Palm_profit[recs2], decreasing=TRUE)
				if (i == 3)	ord <- order(Data$Forestry_profit[recs2], decreasing=TRUE)
				if (i == 4)	ord <- order(Data$Rice_profit[recs2], decreasing=TRUE)
				if (i == 5)	ord <- order(Data$Soy_profit[recs2], decreasing=TRUE)
				keep.existing[[i]] <- ids[ord[1:prod.area.targets[i]]]
			}
		}
		recs <- which(Data$Id %in% keep.existing[[1]] | Data$Id %in% keep.existing[[2]] | 
			Data$Id %in% keep.existing[[3]] | Data$Id %in% keep.existing[[4]] | Data$Id %in% keep.existing[[5]])
	} else {
		keep.existing[[1]] <- Data$Id[which(Data$Current.land.use.num == 1)]
		keep.existing[[2]] <- Data$Id[which(Data$Current.land.use.num == 2)]
		keep.existing[[3]] <- Data$Id[which(Data$Current.land.use.num == 3)]
		keep.existing[[4]] <- Data$Id[which(Data$Current.land.use.num == 4)]
		keep.existing[[5]] <- Data$Id[which(Data$Current.land.use.num == 5)]
	}


	prof.existing <- c(sum(Data$Livestock_profit[which(Data$Id %in% keep.existing[[1]])]), 
		sum(Data$Palm_profit[which(Data$Id %in% keep.existing[[2]])]), 
		sum(Data$Forestry_profit[which(Data$Id %in% keep.existing[[3]])]), 
		sum(Data$Rice_profit[which(Data$Id %in% keep.existing[[4]])]), 
		sum(Data$Soy_profit[which(Data$Id %in% keep.existing[[5]])]))
	
	prof.existing.CI.u <- c(sum(Data$Livestock_profitU[which(Data$Id %in% keep.existing[[1]])]), 
	                   sum(Data$Palm_profitU[which(Data$Id %in% keep.existing[[2]])]), 
	                   sum(Data$Forestry_profitU[which(Data$Id %in% keep.existing[[3]])]), 
	                   sum(Data$Rice_profitU[which(Data$Id %in% keep.existing[[4]])]), 
	                   sum(Data$Soy_profitU[which(Data$Id %in% keep.existing[[5]])]))
	
	prof.existing.CI.l <- c(sum(Data$Livestock_profitL[which(Data$Id %in% keep.existing[[1]])]), 
	                        sum(Data$Palm_profitL[which(Data$Id %in% keep.existing[[2]])]), 
	                        sum(Data$Forestry_profitL[which(Data$Id %in% keep.existing[[3]])]), 
	                        sum(Data$Rice_profitL[which(Data$Id %in% keep.existing[[4]])]), 
	                        sum(Data$Soy_profitL[which(Data$Id %in% keep.existing[[5]])]))

	

	cb.loss.existing <- sum(Data$Carbon[recs])
	cb.loss.existing.er.u <- sum(Data$CarbonerU[recs])
	cb.loss.existing.er.l <- sum(Data$CarbonerL[recs])
	h2o.loss.existing <- sum(Data$H2Oyr.1[recs])
	h2o.loss.existing.er.u <- sum(Data$H2OerU[recs])
	h2o.loss.existing.er.l <- sum(Data$H2OerL[recs])

	
	
	Data <- Data[-recs,]

	# this works because the area of each PU is exactly 1km:
	prod.area.targets[1] <- max(0, prod.area.targets[1] - length(keep.existing[[1]]))
	prod.area.targets[2] <- max(0, prod.area.targets[2] - length(keep.existing[[2]]))
	prod.area.targets[3] <- max(0, prod.area.targets[3] - length(keep.existing[[3]]))
	prod.area.targets[4] <- max(0, prod.area.targets[4] - length(keep.existing[[4]]))
	prod.area.targets[5] <- max(0, prod.area.targets[5] - length(keep.existing[[5]]))

	keep.existing.n <- length(recs)

	message("Revised production area targets:")
	print(prod.area.targets)
	message(paste0("Total area pre-zoned: ", length(recs)))
	message(paste0("Total target area (original): ", sum(prod.area.targets.t0)))
	message(paste0("Total target area (revised): ", sum(prod.area.targets)))
	message(paste0("Total target area (revised + pre-zoned): ", sum(prod.area.targets) + length(recs)))
} else {
  
  cb.loss.existing <- 0
  cb.loss.existing.er.u <- 0
  cb.loss.existing.er.l <- 0
  h2o.loss.existing <- 0
  h2o.loss.existing.er.u <- 0
  h2o.loss.existing.er.l <- 0
  prof.existing.CI.u <- c(0,0,0,0,0)
  prof.existing.CI.l <- c(0,0,0,0,0)
  
}

dim(Data)


# REDUCED DATASET
# need to reset variables so that all data objects match the correct order of PUs
# the sppmat now only refers to the subset of PUs that are available for zoning decisions
sppmat <- as.matrix(Data[c(startcol:endcol)])
# remove species with no cells
cs <- colSums(sppmat)
recs <- which(cs == 0)
if (length(recs > 0)) sppmat <- sppmat[,-recs]

if (!dim(sppmat)[2] == check.nspp){
	message("ERROR: Species counts have changed")
	stop()
}

# remove species data from Data object
Data <- Data[, c(1:(startcol - 1))]
head(Data)



# make sure the ID numbers are recorded so that the results can be tied back to 
# the spatial dataset
pid <- Data$Id

# variable for number of active planning units
np <- dim(Data)[1]

# production value vector
pv <- c(Data$Livestock_profit, Data$Palm_profit, Data$Forestry_profit, Data$Rice_profit, 
             Data$Soy_profit)
summary(pv)
pv.mean <- mean(pv[1:(np*5)])
pv.sd <- sd(pv[1:(np*5)])

# water value vector (0 for all non-native veg types)
h2o <- Data$H2Oyr.1
summary(h2o)
h2o.mean <- mean(Data$H2Oyr.1)
h2o.sd <- sd(Data$H2Oyr.1)

# carbon value vector (0 for all non-native veg types)
cb <- Data$Carbon
summary(cb)
cb.mean <- mean(Data$Carbon)
cb.sd <- sd(Data$Carbon)


# species value vector

# local extinction risk is estimated as a function of the original area of habitat for each species
# in the study area, and the existing area that falls within either a protected area or a PU
# designated as native vegetation

# r = 1 − (x/A0)^z

# determine which cells are currently native veg
# must assume that areas that have not been pre-zoned as agriculture still count as native veg (or 
# would be converted back to native veg)
nveg.t0 <- rep(1, np)
# nveg.t0[which(Data$Current.land.use.num == 6)] <- 1

sp <- species.calc(nveg.t0, sppmat, spppatot, sppartot, a=erisk.power)
sp.mean <- mean(sp)
sp.sd <- sd(sp)

# standardise all objective function variables
pv.std <- (pv - pv.mean) / pv.sd
h2o.std <- (h2o - h2o.mean) / h2o.sd
cb.std <- (cb - cb.mean) / cb.sd
sp.std <- (sp - sp.mean) / sp.sd

# construct sparse constraints matrix:
# empty vectors to populatate:
row <- rep(0, 1E7)
col <- rep(0, 1E7)
z <- rep(0, 1E7)
idx <- 1

# Constraint 1: must meet the area target for each iteration of optimisation
# row 1, z value is area (1 km2 for each PU)
row[1:(np*5)] <- rep(1, np*5)
col[1:(np*5)] <- c(1:(np*5)) 
z[1:(np*5)] <- rep(1, np*5)
idx <- np*5 + 1


# Constraint set 2: must not exceed area target for each production type

# livestock:
row[idx:(idx + np - 1)] <- rep(2, np)
col[idx:(idx + np - 1)] <- c(1:np) 
z[idx:(idx + np - 1)] <- rep(1, np)
idx <- idx + np

# palm:
row[idx:(idx + np - 1)] <- rep(3, np)
col[idx:(idx + np - 1)] <- c((np + 1):(np*2)) 
z[idx:(idx + np - 1)] <- rep(1, np)
idx <- idx + np

# forestry:
row[idx:(idx + np - 1)] <- rep(4, np)
col[idx:(idx + np - 1)] <- c((np*2 + 1):(np*3)) 
z[idx:(idx + np - 1)] <- rep(1, np)
idx <- idx + np

# rice:
row[idx:(idx + np - 1)] <- rep(5, np)
col[idx:(idx + np - 1)] <- c((np*3 + 1):(np*4)) 
z[idx:(idx + np - 1)] <- rep(1, np)
idx <- idx + np

# soy:
row[idx:(idx + np - 1)] <- rep(6, np)
col[idx:(idx + np - 1)] <- c((np*4 + 1):(np*5)) 
z[idx:(idx + np - 1)] <- rep(1, np)
idx <- idx + np

# Constraint set 3: ecoregion representation


if (retention){
  et.df <- read.csv('Ecosystem_targets_new.csv', header=F)
  
} else { et.df <- read.csv('Ecosystem_targets_new_no.csv', header=F)
  }

#et.df <- read.csv('Ecosystem_targets_new.csv', header=F)

# remove any ecosystems with 0 area target (not relevant)
recs <- which(et.df[,2] == 0)
if (length(recs) > 0) et.df <- et.df[-recs,]


#
ecosystem.targets <- rep(0, dim(et.df)[1])
row.idx <- 7


for (i in 1:dim(et.df)[1]){
	recs <- which(Data$Ecosystem == et.df[i, 1])
	if (length(recs) > 0){
		for (j in 1:5){
			row[idx:(idx + length(recs) - 1)] <- rep(row.idx, length(recs))
			col[idx:(idx + length(recs) - 1)] <- recs + ((j - 1) * np)
			z[idx:(idx + length(recs) - 1)] <- rep(1, length(recs))
			idx <- idx + length(recs)
		}
		row.idx <- row.idx + 1

		# check that it is possible for the target to be met:
		if (length(recs) < et.df[i,2]){
			message(paste0("Warning: target is infeasible for ecosystem ", et.df[i,1]), ", reducing target...")
		}

		# the target needs to be set not as the area that can be allocated, but as the maximum
		# area of each ecosystem that can be converted. The following code works only because the 
		# area of each planning unit is 1 km2:
		ecosystem.targets[i] <- max(0, length(recs) - et.df[i,2]) 

	}
}


# Constraint set 5: ensure that each PU is assigned to only one land cover type
for (i in 1:np){
	row[idx:(idx + 4)] <- rep(row.idx, 5)
	col[idx:(idx + 4)] <- (c(0:4) * np) + i
	z[idx:(idx + 4)] <- rep(1, 5)
	idx <- idx + 5
	row.idx <- row.idx + 1
}


# set up rhs and sense vectors

# set the target vector for the areas to be zones in each iteration of the algorithm.
# note that the first one *must be large enough to accommodate all of the cells that are 
# currently agriculture* - this can lead to some tricky bookkeeping
# M_q <- rep(0, niters)
# M_q[1] <- length(recs.current.ag)
# M_q[2:niters] <- trunc(seq(length(recs.current.ag), sum(prod.area.targets), length=niters))[2:niters]

M_q <- trunc(seq(0, sum(prod.area.targets), length=(niters+1)))[-1]

# create vectors
# sense <- c("<=", rep("<=", 5), rep("<=", length(ecosystem.targets)), ">=", rep("<=", np))
# rhs <- c(M_q[1], prod.area.targets, ecosystem.targets, length(recs.current.ag), rep(1, np))

sense <- c("<=", rep("<=", 5), rep("<=", length(ecosystem.targets)), rep("<=", np))
rhs <- c(M_q[1], prod.area.targets, ecosystem.targets, rep(1, np))


# need to initialise this vector with zeros for the optimisation code:
prev.x <- rep(0, np * 5)


# create a coordinates object that matches the order of the PUs to help with visualisation
# this is slow, but it is only run once
if (make.maps){
	if (!exists("coords")){
		message("Preparing mapping data...")
		recs <- rep(0, np)
		for (i in 1:np)	recs[i] <- which(pu.df$id == Data$Id[i])
		coords <- pu.df[recs,c(2, 3)]

		if (keep.existing.ag){
			npk <- length(keep.existing[[1]]) + length(keep.existing[[2]]) + length(keep.existing[[3]]) + 
			length(keep.existing[[4]]) + length(keep.existing[[5]])
			recs <- rep(0, npk)
			idx2 <- 1
			for (j in 1:5){
				if (length(keep.existing[[j]]) > 0){
					for (i in 1:length(keep.existing[[j]])){
						recs[idx2] <- which(pu.df$id == keep.existing[[j]][i])
						idx2 <- idx2 + 1
					}	
				}
			}
			coords.existing <- cbind(pu.df[recs,c(2, 3)], c(rep(1, length(keep.existing[[1]])), rep(2, length(keep.existing[[2]])), rep(3, length(keep.existing[[3]])), 
				rep(4, length(keep.existing[[4]])), rep(5, length(keep.existing[[5]]))))
			names(coords.existing) <- c("x", "y", "class")
		}

		rm(pu.df)
	}
}


message("Finished initialisation")




