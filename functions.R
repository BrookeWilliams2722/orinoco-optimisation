


# this function usese the decision vector, x, to calculate the current value of each PU to 
# each species
species.calc <- function(nveg, sppmat, spppatot, sppartot, a=0.25){

	# for each spexies, how much habitat area is provided by nveg
	hab.area <- as.vector(matrix(nveg, nrow=1) %*% sppmat)

	# extrinction risk:
	er <- as.vector(1 - ((spppatot + hab.area) / sppartot)^a)

	# delta extinction risk if 1 ha area of habitat is lost
	der <- as.vector(1 - ((spppatot + hab.area - 0.01) / sppartot)^a)

	# calculate slope of tangent to er curve
	delta.er <- (der - er) / 0.01
	# a precautionary measure:
	if (min(delta.er) < 0) delta.er[which(delta.er < 0)] <- 0

	# for each PU, calculate the sum of delta extinction risks
	spp <- apply(sppmat, MAR=1, FUN=function(x){return(delta.er %*% x)})
	if (min(spp) < 0) spp[which(spp < 0)] <- 0

	return(spp)
}




objective.function <- function(v, h, c, s, w){
  y <- (w[1] * v) - (w[2] * ((w[3] * rep(h, 5)) + (w[4] * rep(c, 5)) + (w[5] * rep(s, 5))))
  return(y + abs(min(y)) + 1)
}


# function for calculating benefits from result vector
calculate_benefits <- function(x, m, w, step){
# x is a binary result vector x

  liv.prof <- ((x[1:np] %*% Data$Livestock_profit) + prof.existing[1]) / 10^6
  palm.prof <- ((x[(np+1):(2*np)] %*% Data$Palm_profit) + prof.existing[2]) / 10^6
  for.prof <- ((x[((2*np)+1):(3*np)] %*% Data$Forestry_profit) + prof.existing[3]) / 10^6
  rice.prof <- ((x[((3*np)+1):(4*np)] %*% Data$Rice_profit) + prof.existing[4]) / 10^6
  soy.prof <- ((x[((4*np)+1):(5*np)] %*% Data$Soy_profit) + prof.existing[5]) / 10^6
  
  ag <- x[1:np] + x[(np+1):(2*np)] + x[((2*np)+1):(3*np)] + 
    x[((3*np)+1):(4*np)] + x[((4*np)+1):(5*np)]

  # also divide carbon and water by 10^6 to get plot-friendly values
  carbon.loss <- ((ag %*% Data$Carbon) + cb.loss.existing) / 10^6
  water.loss <- ((ag %*% Data$H2Oyr.1) + h2o.loss.existing) / 10^6
  
  prop.ag <- (sum(x) + keep.existing.n) / artot # (np + keep.existing.n)

  # extinction risk gain is calcaulated as the change in extinction risk
  nveg <- rep(0, np)
  nveg[which(colSums(matrix(x, ncol=np, nrow=5, byrow=TRUE)) == 0)] <- 1
  hab.area <- matrix(nveg, nrow=1) %*% sppmat
  er <- as.vector(1 - ((spppatot + hab.area) / sppartot)^erisk.power)
  er.gain <- sum(er) / length(er)
  ex.avoid <- sum(er)
  ex.avoidsd <- sqrt(sum(er * (1 - er)))

  total.prof <- liv.prof + palm.prof + for.prof + rice.prof + soy.prof
  results.df <- data.frame(m, w, step, liv.prof, palm.prof, for.prof, rice.prof, soy.prof, 
     total.prof, carbon.loss, water.loss, er.gain, ex.avoid, ex.avoidsd, prop.ag)
  return(results.df)

}


plot.benefits <- function(){

	pnt.cex <- 0.25

	# biodiversity
	cols=c("yellow", "green4")
	colpal1 <- colorRampPalette(cols)
	cr <- colpal1(100)
	cols <- rep("", length(sp))

	brks <- quantile(sp, prob=seq(0, 1, len=101))
	for (j in 1:(length(brks)-1)){
		cols[which(sp >= brks[j] & sp <= brks[j+1])] <- cr[j]
	}

	png(file=paste0(prefix, "/", "map_bd_value_tm_", m, "_time0.png"), height=1000, width=1000)
	par(mar=c(0,0,0,0))
	plot(coords, type="n", xlab="", ylab="", axes=F)
	points(coords, pch=15, col=cols, cex=pnt.cex)
	dev.off()

	# carbon
	cols=c("yellow", "red")
	colpal1 <- colorRampPalette(cols)
	cr <- colpal1(100)
	cols <- rep("", length(cb))

	brks <- quantile(cb, prob=seq(0, 1, len=101))
	for (j in 1:(length(brks)-1)){
		cols[which(cb >= brks[j] & cb <= brks[j+1])] <- cr[j]
	}

	png(file=paste0(prefix, "/", "map_cb_value_tm_", m, "_time0.png"), height=1000, width=1000)
	par(mar=c(0,0,0,0))
	plot(coords, type="n", xlab="", ylab="", axes=F)
	points(coords, pch=15, col=cols, cex=pnt.cex)
	dev.off()

	# water
	cols=c("yellow", "blue")
	colpal1 <- colorRampPalette(cols)
	cr <- colpal1(100)
	cols <- rep("", length(h2o))

	brks <- quantile(h2o, prob=seq(0, 1, len=101))
	for (j in 1:(length(brks)-1)){
		cols[which(h2o >= brks[j] & h2o <= brks[j+1])] <- cr[j]
	}

	png(file=paste0(prefix, "/", "map_h20_value_tm_", m, "_time0.png"), height=1000, width=1000)
	par(mar=c(0,0,0,0))
	plot(coords, type="n", xlab="", ylab="", axes=F)
	points(coords, pch=15, col=cols, cex=pnt.cex)
	dev.off()
}


plot.objectivefunction <- function(){

	pnt.cex <- 0.25

	# carbon
	cols=c("red", "yellow", "cyan", "blue")
	colpal1 <- colorRampPalette(cols)
	cr <- colpal1(100)
	cols <- rep("", length(of))
	labels <- c("livestock", "palm", "forestry", "rice", "soy")

	brks <- quantile(of, prob=seq(0, 1, len=101))
	for (j in 1:(length(brks)-1)){
		cols[which(of >= brks[j] & of <= brks[j+1])] <- cr[j]
	}

	png(file=paste0(prefix, "/", "map_of_value_tm_", m, "_w_", w, "_i_", i, ".png"), height=2000, width=3000)
	par(mar=c(0,0,0,0), mfrow=c(2, 3))
	for (i in 1:5){
		plot(coords, type="n", xlab="", ylab="", axes=F)
		points(coords, pch=15, col=cols[c(1:np) + ((i-1)*np)], cex=pnt.cex)
		text(min(coords[,1]), max(coords[,2]), labels[i], pos=4, cex=5)
    }

    # plot the allocated zones
    plot(coords, type="n", xlab="", ylab="", axes=F)
    points(coords, pch=15, col="grey80", cex=pnt.cex)
    # plot zones
    cols <- c("chartreuse", "red", "green4", "darkorange", "purple", "steelblue")
    for (i in 1:5){
       recs <- which(result$x[c(1:np) + ((i-1)*np)] == 1)
       if (length(recs) > 0) points(coords[recs,], pch=15, col=cols[i], cex=pnt.cex)
    }
    legend("topleft", pch=15, cex=5, col=cols, legend=labels, bty="n")
	dev.off()

}





