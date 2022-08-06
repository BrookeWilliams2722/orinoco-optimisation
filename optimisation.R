#This code was developed by Dr Hawthorne Beyer (h.beyer@uq.edu.au) and Brooke Williams (brooke.williams@uq.edu.au) 
#from the University of Queensland

#By using this code you agree to inform the authors of any publications, applications for funding, funding acquired and all other applications associated with its use
#By using this code you agree to acknowledge the intellectual property of the authors in all published work, applications and dealings with its use

#The correct citation is:
#Williams B A, Grantham H S, Watson J E M, Alvarez S J, Simmonds J S, Rogéliz C A, Da Silva M, Forero Medina G, Etter A, Nogales J, Walschburger T, Hyman G, Beyer H L. (2020) Minimising the loss of biodiversity and ecosystem services in an intact landscape under risk of rapid agricultural development. Environmental Research Letters, 15, 014001. https://doi.org/10.1088/1748-9326/ab5ff7  

#For enquiries please email brooke.williams@uq.edu.au 

#The optimisation is solved using gurobi which must be licenced and installed https://www.gurobi.com/ 
#To run this code you need the following files: pu.df.r, initialise_data.r, functions.r, Final_sheet.csv, Ecosystem_targets_new or Ecosystem_targets_new_no

#setwd("C:/Users/Set_working_directory_to_current_folder")
Data <- read.csv("Final_sheet.csv")
save(Data, file="Data.RData")


options(stringsAsFactors = FALSE)
library(gurobi)

# load the functions used by this code
source("functions.R")

### set global variables
erisk.power <- 0.25
niters <- 20
pnt.cex <- 0.4
#This will create a folder to save your outputs
prefix <- "output"
if (!dir.exists(prefix)) dir.create(prefix)


#These are your expansion targets in km2, the order is Livestock, Oil palm, Forestry, Rice and Soy. 
#The example provided here covers the entire study region
prod.area.targets.max <- c(51025, 51025, 51025, 51025, 51025)

#You can use the multipliers if you want to look at proportions of the study area to be 
#converted - it will reduce your areal expansion targets by each value. 
prod.area.multipliers <- c(0.15,0.2,0.3,0.4)

# Use this if you want to create images of your solutions
# if you set make.maps = TRUE then when you first run the initialise_data.r code, it will
# load up the data needed to make the maps of the zones - this take about 5 minutes, but
# it is only run once
#make.maps <- FALSE
make.maps <- TRUE
make.diagnostic.maps <- FALSE


# this variable controls whether all existing ag cells are kept in their existing zone
# or are allocated to other zones
# it also reduces the production area targets to a minimum of 0 after accounting for
# these existing cells
# also checks whether ecosystem targes are feasible after removing these cells and adjusts
# them if not
keep.existing.ag <- TRUE
keep.target.only <- TRUE


# this variable controls whether ecosystem retention targets are used
retention <- TRUE



### RUN OPTIMISATION ###############################################################

# max. w1 v x + w2 [ w3 h x + w4 c x - w5 s x ]
# s.t. b_j x >= T_j
# 	e_j x >= E_j
# v = prod value
# h = water value
# c = carbon
# s = species metric
# b = commodity type
# e = ecosystem area



####Weight Matrix###
####This is the weight matrix where you define how much you care about each objective#######
# the first two values in each row should always sum to 1
# the last three values in each row should always sum to 1
#The order is 1. Overall agricultural production value 2. Overall biodiversity (species) and ecosystem services
#3. Water retention 4. Carbon sequestration and 5. Species persistance

#You can use it to generate a single optimal solution
weightm <- matrix(0, nrow=1, ncol=5)
weightm[1,] <- c(0.5,0.5,1/3,1/3,1/3)

#Or to generate trade-off curves and evaluate a variety of scenarios against each other
#This example evaluates the trade-off between Agricultural production value and environmental objectives, 
#when you weight all three objectives equally
# weightm <- matrix(0, nrow=11, ncol=5)
# 
# weightm[1,] <- c(1,0,1/3,1/3,1/3)
# weightm[2,] <- c(0.9,0.1,1/3,1/3,1/3)
# weightm[3,] <- c(0.8,0.2,1/3,1/3,1/3)
# weightm[4,] <- c(0.7,0.3,1/3,1/3,1/3)
# weightm[5,] <- c(0.6,0.4,1/3,1/3,1/3)
# weightm[6,] <- c(0.5,0.5,1/3,1/3,1/3)
# weightm[7,] <- c(0.4,0.6,1/3,1/3,1/3)
# weightm[8,] <- c(0.3,0.7,1/3,1/3,1/3)
# weightm[9,] <- c(0.2,0.8,1/3,1/3,1/3)
# weightm[10,] <- c(0.1,0.9,1/3,1/3,1/3)
# weightm[11,] <- c(0,1,1/3,1/3,1/3)


if (exists("benefits.df")) rm(benefits.df)

for (m in 1:length(prod.area.multipliers)){
  prod.area.targets.t0 <- trunc(prod.area.targets.max * prod.area.multipliers[m])
  # must recreate the coords object every time that the production target areas change
  if (exists("coords")) rm(coords)

  for (w in 1:dim(weightm)[1]){
    
    # ALWAYS reset everything at the beginning of a new optimisation - this is an important
    # safety precaution
    source('initialise_data.r')

    if (make.maps & make.diagnostic.maps) plot.benefits()

    # initialise the beenfits dataframe for the first time with pre-solution values:
    if (!exists("benefits.df")) benefits.df <- calculate_benefits(rep(0, 5*np), 0, 0, 0)
           
    for (i in 1:niters){
     
        # calculate / recalculate objective function
        of <- objective.function(pv.std, h2o.std, cb.std, sp.std, weightm[w,])


        # make the sparse constraint matrix
        #constr <- sparseMatrix(i=row[1:(idx-1)], j=col[1:(idx-1)], x=z[1:(idx-1)])
        constr <- sparseMatrix(i=row[1:(idx-1)], j=col[1:(idx-1)], x=z[1:(idx-1)])
        
        # update the first objective
        rhs[1] <- M_q[i]
        
        #Bcreate sense vector - this will change if PU's change
        #sense <- as.character(rep(">=", length=165603))
        #rhs <- rep(1, length=165603)
        
        # set up Gurobi model
        model <- list()
        model$obj <- of
        model$modelsense <- "max"
        model$vtype <- "B"
        model$A <- constr
        model$rhs <- rhs
        model$sense <- sense
        params <- list(Threads=4) #Presolve=2, MIPGap=0.005,
        result <- gurobi(model,params)
        
        if (make.maps & make.diagnostic.maps) plot.objectivefunction()

        message(paste0("Iteration ", i, ": ", sum(result$x), " of ", rhs[1], " allocated"))
        
        benefits.df <- rbind(benefits.df, calculate_benefits(result$x, m, w, i))

        # must add constraints to ensure that all areas that were zoned in previous iterations
        # are zoned in exactly the same way in the next iteration
        if (i == 1){
           recs <- which(result$x == 1)
           row[idx:(idx + length(recs) - 1)] <- rep(row.idx, length(recs))
           col[idx:(idx + length(recs) - 1)] <- recs
           z[idx:(idx + length(recs) - 1)] <- rep(1, length(recs))
           idx <- idx + length(recs)
           rhs <- c(rhs, sum(result$x))
           sense <- c(sense, ">=")
        
        } else {
           recs <- which(result$x == 1 & prev.x == 0)
           row[idx:(idx + length(recs) - 1)] <- rep(row.idx, length(recs))
           col[idx:(idx + length(recs) - 1)] <- recs
           z[idx:(idx + length(recs) - 1)] <- rep(1, length(recs))
           idx <- idx + length(recs)
           rhs[length(rhs)] <- sum(result$x)
        }
        
        
        prev.x <- result$x
        
        if (i < niters){
           # recalculate species vector
           # which cells are native veg?
           nveg <- rep(0, np)
           nveg[which(colSums(matrix(prev.x, ncol=np, nrow=5, byrow=TRUE)) == 0)] <- 1
           sp <- species.calc(nveg, sppmat, spppatot, sppartot, a=erisk.power)
           sp.mean <- mean(sp)
           sp.sd <- sd(sp)
           sp.std <- (sp - sp.mean) / sp.sd

           # save
           save(result, file=paste0(prefix, "/", prefix, "_tm_", m, "_w_", w, "_result.RData"))
           
        }
     
     }
     
     # save
     save(result, file=paste0(prefix, "/", prefix, "_tm_", m, "_w_", w, "_result.RData"))
     
     if (make.maps){
        # map
        png(file=paste0(prefix, "/", "map_zones_", prefix, "_tm_", m, "_w_", w, ".png"), height=800, width=800)
        par(mar=c(0,0,0,0))
        plot(coords, pch=15, col="grey80", cex=pnt.cex, xlab="", ylab="", axes=F)
        points(pacoords, pch=15, col="darkgreen", cex=pnt.cex)

        # plot zones
        cols <- c("chartreuse", "red", "goldenrod1", "dodgerblue", "maroon2", "green4")
        for (i in 1:5){
           recs <- which(result$x[c(1:np) + ((i-1)*np)] == 1)
           if (length(recs) > 0) points(coords[recs,], pch=15, col=cols[i], cex=pnt.cex)
        }

        # plot pre-existing zones
        if (keep.existing.ag){
           for (i in 1:5){
              recs <- which(coords.existing$class == i)
              if (length(recs) > 0) points(coords.existing[recs, c(1,2)], pch=15, col=cols[i], cex=pnt.cex)
           }
        }

        # natural areas
        nveg <- rep(0, np)
        nveg[which(colSums(matrix(result$x, ncol=np, nrow=5, byrow=TRUE)) == 0)] <- 1
        recs <- which(nveg == 1)
        if (length(recs) > 0) points(coords[recs,], pch=15, col=cols[6], cex=pnt.cex)

        legend("topleft", pch=15, col=c(cols, "darkgreen"), legend=c("livestock", "palm", "forestry", "rice", 
          "soy", "natural", "protected areas"), bty="n", cex=2)

        dev.off()
        
        ####create data frame with the selected ag land uses for if you want to plot a solution in arcmap###
        Livestockselection <- result$x[1:(length(result$x)/5)]
        Palmselection <- result$x[((length(result$x)/5)+1):((length(result$x)/5)*2)]
        Forestryselection <- result$x[(((length(result$x)/5)*2)+1):((length(result$x)/5)*3)]
        Riceselection <- result$x[(((length(result$x)/5)*3)+1):((length(result$x)/5)*4)]
        Soyselection <- result$x[(((length(result$x)/5)*4)+1):(length(result$x))]
        
        Solution <- data.frame(pid, Livestockselection, Palmselection, Forestryselection, Riceselection, Soyselection)
        write.csv(Solution, file =paste0(prefix, "/", prefix, "_tm_", m, "_w_", w, "_result.csv") )
        write.csv(ecosystem.targets, file =paste0(prefix, "/", prefix, "_tm_", m, "_w_", w, "_ET.csv") )
        
     }
           
  }

}

save(benefits.df, file=paste0(prefix, "/", prefix, "_benefits.df.RData"))
write.csv(benefits.df, file=paste0(prefix, "/", prefix, "_benefits.df.csv"))


### POST PROCESSING ##################################################################
##Calculating CI

liv.prof <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
palm.prof <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
for.prof <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
rice.prof <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
soy.prof <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))

liv.prof.u <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
liv.prof.l <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
palm.prof.u <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
palm.prof.l <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
for.prof.u <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
for.prof.l <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
rice.prof.u <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
rice.prof.l <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
soy.prof.u <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
soy.prof.l <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))

carbon.loss <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
carbon.loss.er.u <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
carbon.loss.er.l <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
water.loss <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
water.loss.er.u <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
water.loss.er.l <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))

er.gain <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
ex.avoid <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))
ex.avoidsd <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))

prop.ag <- matrix(0, nrow=dim(weightm)[1], ncol=length(prod.area.multipliers))

# safest to rerun the initialisation script before postprocessing

make.maps <- FALSE

for (m in 1:length(prod.area.multipliers)){
  prod.area.targets.t0 <- trunc(prod.area.targets.max * prod.area.multipliers[m])
  if (exists("coords")) rm(coords)
  source('initialise_data.r')
  
  for (w in 1:dim(weightm)[1]){
    
    load(paste0(prefix, "/", prefix, "_tm_", m, "_w_", w, "_result.RData"))
    
    # profitability divided by 10^6 to convert to millions
    liv.prof[w, m] <- ((result$x[1:np] %*% Data$Livestock_profit) + prof.existing[1]) / 10^6
    palm.prof[w, m] <- ((result$x[(np+1):(2*np)] %*% Data$Palm_profit) + prof.existing[2]) / 10^6
    for.prof[w, m] <- ((result$x[((2*np)+1):(3*np)] %*% Data$Forestry_profit) + prof.existing[3]) / 10^6
    rice.prof[w, m] <- ((result$x[((3*np)+1):(4*np)] %*% Data$Rice_profit) + prof.existing[4]) / 10^6
    soy.prof[w, m] <- ((result$x[((4*np)+1):(5*np)] %*% Data$Soy_profit) + prof.existing[5]) / 10^6

    ##Create CI for profit##
    liv.prof.u[w, m] <- ((result$x[1:np] %*% Data$Livestock_profitU) + prof.existing.CI.u[1]) / 10^6
    liv.prof.l[w, m] <- ((result$x[1:np] %*% Data$Livestock_profitL) + prof.existing.CI.l[1]) / 10^6
    palm.prof.u[w, m] <- ((result$x[(np+1):(2*np)] %*% Data$Palm_profitU) + prof.existing.CI.u[2]) / 10^6
    palm.prof.l[w, m] <- ((result$x[(np+1):(2*np)] %*% Data$Palm_profitL) + prof.existing.CI.l[2]) / 10^6
    for.prof.u[w, m] <- ((result$x[((2*np)+1):(3*np)] %*% Data$Forestry_profitU) + prof.existing.CI.u[3]) / 10^6
    for.prof.l[w, m] <- ((result$x[((2*np)+1):(3*np)] %*% Data$Forestry_profitL) + prof.existing.CI.l[3]) / 10^6
    rice.prof.u[w, m] <- ((result$x[((3*np)+1):(4*np)] %*% Data$Rice_profitU) + prof.existing.CI.u[4]) / 10^6
    rice.prof.l[w, m] <- ((result$x[((3*np)+1):(4*np)] %*% Data$Rice_profitL) + prof.existing.CI.l[4]) / 10^6
    soy.prof.u[w, m] <- ((result$x[((4*np)+1):(5*np)] %*% Data$Soy_profitU) + prof.existing.CI.u[5]) / 10^6
    soy.prof.l[w, m] <- ((result$x[((4*np)+1):(5*np)] %*% Data$Soy_profitL) + prof.existing.CI.l[5]) / 10^6
      
      
    ag <- result$x[1:np] + result$x[(np+1):(2*np)] + result$x[((2*np)+1):(3*np)] + 
      result$x[((3*np)+1):(4*np)] + result$x[((4*np)+1):(5*np)]

    # also divide carbon and water by 10^6 to get plot-friendly values
    carbon.loss[w, m] <- ((ag %*% Data$Carbon) + cb.loss.existing) / 10^6
    carbon.loss.er.u [w, m] <- ((ag %*% Data$CarbonerU) + cb.loss.existing.er.u) / 10^6
    carbon.loss.er.l [w, m] <- ((ag %*% Data$CarbonerL) + cb.loss.existing.er.l) / 10^6
    water.loss[w, m] <- ((ag %*% Data$H2Oyr.1) + h2o.loss.existing) / 10^6
    water.loss.er.u [w, m] <- ((ag %*% Data$H2OerU) + h2o.loss.existing.er.u) / 10^6
    water.loss.er.l [w, m] <- ((ag %*% Data$H2OerL) + h2o.loss.existing.er.l) / 10^6

    prop.ag[w, m] <- (sum(result$x) + keep.existing.n) / artot 

    # extinction risk gain is calcaulated as the change in extinction risk
    # baseline.er is already calculated in the initialisation code
    nveg <- rep(0, np)
    nveg[which(colSums(matrix(result$x, ncol=np, nrow=5, byrow=TRUE)) == 0)] <- 1
    hab.area <- matrix(nveg, nrow=1) %*% sppmat
    er <- as.vector(1 - ((spppatot + hab.area) / sppartot)^erisk.power)
    # er.gain[w, m] <- sum(baseline.er - er) / length(baseline.er)
    # ex.avoid[w, m] <- sum(baseline.er - er)
    # ex.avoidsd[w, m] <- sqrt(sum(baseline.er * (1 - baseline.er)) + sum(er * (1 - er)))
    er.gain[w, m] <- sum(er) / length(er)
    ex.avoid[w, m] <- sum(er)
    ex.avoidsd[w, m] <- sqrt(sum(er * (1 - er)))
  }
}

####For species, take the SD and change it into 90% CI#########
ex.avoidCI.u <- ex.avoid + (1.645*(ex.avoidsd/sqrt(145)))
ex.avoidCI.l <- ex.avoid - (1.645*(ex.avoidsd/sqrt(145)))

###Add profit vectors to create one profit value#########
total.prof <- as.vector(liv.prof) + as.vector(palm.prof) + as.vector(for.prof) + as.vector(rice.prof) + 
  as.vector(soy.prof)

total.profCI.u <- as.vector(liv.prof.u) + as.vector(palm.prof.u) + as.vector(for.prof.u) + as.vector(rice.prof.u) +
  as.vector(soy.prof.u)

total.profCI.l <- as.vector(liv.prof.l) + as.vector(palm.prof.l) + as.vector(for.prof.l) + as.vector(rice.prof.l) +
  as.vector(soy.prof.l)

results.df <- data.frame(rep(1:length(prod.area.multipliers), each=dim(weightm)[1]), rep(c(1:dim(weightm)[1]), length(prod.area.multipliers)), as.vector(liv.prof), 
  as.vector(palm.prof), as.vector(for.prof), as.vector(rice.prof), as.vector(soy.prof), as.vector(total.prof), as.vector(total.profCI.u), as.vector(total.profCI.l),
  as.vector(carbon.loss), as.vector(carbon.loss.er.u), as.vector(carbon.loss.er.l), as.vector(water.loss), as.vector(water.loss.er.u), as.vector(water.loss.er.l), as.vector(er.gain), as.vector(ex.avoid), as.vector(ex.avoidsd), as.vector(ex.avoidCI.u), as.vector(ex.avoidCI.l), 
  as.vector(prop.ag))



names(results.df) <- c("m", "weightm", "liv.prof", "palm.prof", "for.prof", "rice.prof", "soy.prof", "total.prof", "total.profCI.u", "total.profCI.l",
  "carbon.loss", "carbon.loss.er.u", "carbon.loss.er.l", "water.loss", "water.loss.er.u", "water.loss.er.l", "er.gain", "ex.avoid", "ex.avoidsd","ex.avoidCI.u", "ex.avoidCI.l", "prop.ag")

save(results.df, file=paste0(prefix, "/results.df_", prefix, ".RData"))
#######This will export the results#########
write.csv(results.df, file =paste0(prefix, "/", prefix, "results_.df.csv"))
          

### If you used a variety of weights to evaluate the trade-off curves
### PLOT THE TRADEOFF CURVES
pdf(file=paste0(prefix, "/plot_tradeoffs_", prefix, ".pdf"), width=4, height=12)
par(mfrow=c(3,1), mar=c(4.2, 5, 1, 7), mgp=c(3, 1, 0), bty = 'l')
#cols=c("red", "green4", "blue", "black", "orange", "purple")
cols=colorRampPalette(c("darkmagenta", "cyan3", "chartreuse1"))(n = m)

plot(results.df$total.prof, results.df$ex.avoid, type="n", col= "green", xlab=expression(Production~value~(million~USD~yr^{-1})), 
ylab=expression(Number~of~extinctions))
#ylab="Expected local extinctions avoided (no.species)"

for (m in 1:length(prod.area.multipliers)){
  recs <- which(results.df$m == m)
  lines(results.df$total.prof[recs], results.df$ex.avoid[recs], col=cols[m], cex=2, lwd=1)
  #If you want to plot points rather than lines
  #points(results.df$total.prof[recs], results.df$ex.avoid[recs], pch=15, col=cols[m], cex=1)
}

plot(results.df$total.prof, results.df$carbon.loss, type="n", col= "green", xlab=expression(Production~value~(million~USD~yr^{-1})), 
ylab=expression(Carbon~loss~(million~t)))

for (m in 1:length(prod.area.multipliers)){
  recs <- which(results.df$m == m)
  lines(results.df$total.prof[recs], results.df$carbon.loss[recs], col=cols[m], lwd=1)
  #points(results.df$total.prof[recs], results.df$carbon.loss[recs], pch=15, col=cols[m], cex=1)
}

#This plots the gradient legend, locations of elements here tend to change so may need adjusting
x<- 1:10; y<- 1:15; z<- outer( x,y,"+") 
zr = c(min(round(results.df$prop.ag, digits = 2)*100), max(round(results.df$prop.ag, digits = 2)*100))
library(fields)
image.plot(x,y,z, zlim = zr, legend.shrink=0.9, legend.width = 0.4, col=colorRampPalette(c("darkmagenta", "cyan3", "chartreuse1"))(n = 100), add = TRUE, legend.only=TRUE) 
text(22000, 0.85, pos=2, labels = 'Percentage of landscape converted', xpd = NA, srt = -90)
#Alternative legend type
#legend("topleft", legend = prod.area.multipliers, bty = "n", col =  cols, lty = c(1, 1, 1), lwd = 3)

plot(results.df$total.prof, results.df$water.loss, type="n", col= "green", xlab=expression(Production~value~(million~USD~yr^{-1})), 
ylab=expression(Water~loss~(trillion~l~yr^{-1})))

for (m in 1:length(prod.area.multipliers)){
  recs <- which(results.df$m == m)
  lines(results.df$total.prof[recs], results.df$water.loss[recs], col=cols[m], lwd=1)
  #points(results.df$total.prof[recs], results.df$water.loss[recs], pch=15, col=cols[m], cex=1)
}


dev.off()
