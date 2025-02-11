

rm(list =ls())
gc()

# load libraries that we need
library(ggplot2)
library(gridExtra)
library(grid)

# load simulation functions
source ("~/PathFinder/SimulationFunctions_2024dec.R")

# functions list to be used in simulation
Functions <- list(GrowthModels = GrowthModels, RegFunction= RegFunction, PostRegFunction = PostRegFunction, Deep_learning_model = Deep_learning_model)

# growth models parameters
MP <- ModelsAndParameters[[1]]

# Species codes
SP <- c(1, 2, 3)

# the number of 5-year growth projection intervals
NP = 20

# Start the PDF device
pdf("SR16 samples growth projection_newbiomass_biocheck1.pdf", width = 6, height = 8)
#png("SR16 samples growth projection_nppClass.png", width = 10*72, height = 8*72, res = 72)

# load initial data for simulations (you have 20 parts)
for(j in seq(1, 20, by=2)) { # Increment by 2 to process two datasets per iteration
  
  # Initialize an empty list to store the plots
  plots <- list()
  
  # Process two datasets per iteration
  for(k in 0:1) {
    # load initial data for simulations
    simDAT <- readRDS(paste0("~/PathFinder/simDAT", j+k, ".rds"))
    simDAT <- simDAT[!is.na(npp) & npp>0 & Species %in% c(1,2,3), ]
    simDAT <- simDAT[, Species:= as.numeric(paste(Species))]
    simDAT <- simDAT[, hl1:= hl1/10]
    simDAT[npp_class %in% c(2, 4, 7), ]
    
    
    # Group by Species and npp_class, then sample
    simDAT <- simDAT[, .SD[sample(.N, min(.N, round(.N*0.001)), replace = FALSE)], by = .(Species, npp_class)]
    
    # simulation run
    PixSim(Data =simDAT,
           Np = NP, 
           nSpecies = SP,
           functions = Functions,
           WriteOut = TRUE,
           LocalFldr = Fold,
           ModelsAndParameters = MP,
           Reg_Data=Reg_Data)
    
    # Check the results
    Results <- list.files(Fold, full.names = TRUE)
    
    figDAT <- data.table(plotID = character(), npp_class = numeric(), Species = character(),
      biomass1 = numeric(), vol1 = numeric(), ba1 = numeric(), qmd1 = numeric(),
      hl1 = numeric(),Interval = integer())
    
    # Loop over the list and bind the data frames
    for (i in seq_along(Results)) {
      temp <- fst::read_fst(Results[i], as.data.table = TRUE)  # Read the data file
      
      
      # For the first iteration, initialize values
      if (i == 1) {
        
        temp$Interval <- 0
        
        Species <- temp$Species
        plotID <- temp$plotID
        npp_class <- temp$npp_class
        # Bind only the essential columns to figDAT
        figDAT <- rbind(figDAT, temp[, c("plotID", "npp_class", "Species", "biomass1", "vol1", "ba1", "qmd1", "hl1", "Interval")])
        
        # Save these columns for later use in other iterations
        
        plotID <- temp$plotID
        npp_class <- temp$npp_class
        sp <- temp$Species
      } else {
        temp$Interval <- i - 1  # Set the correct interval for subsequent iterations
        temp$plotID <- plotID   # Retain the plotID from the first iteration
        temp$npp_class <- npp_class  # Retain npp_class from the first iteration
        temp$Species <- Species           # Retain species classification
        
        # Append data from the subsequent results
        figDAT <- rbind(figDAT, temp)
      }
      
      rm(temp)
      gc() # force garbage collection to free up memory
    }
    
    #### nnew pp clesses 
    # Create new NPP_class and set order
    figDAT[, NPP_class := factor(ifelse(npp_class %in% c(1,2), "low",
                                        ifelse(npp_class %in% c(3,4), "medium", "high")),
                                 levels = c("low", "medium", "high"), ordered = TRUE)]
    figDAT[Interval %in% c(0, 5, 10, 15, 20), ]
    # create plot
    # Define color palette
    my_colors <- c("low" = "blue", "medium" = "green", "high" = "red") # change the colors as needed
    setDT(figDAT)
    
    # Reshape to long format using melt
    figDAT_long <- melt(figDAT, 
                        id.vars = c("Interval", "plotID", "NPP_class", "Species"), 
                        measure.vars = c("vol1", "biomass1", "hl1", "qmd1", "ba1"),
                        variable.name = "variable", 
                        value.name = "value")
    
    # Create the plot with free y-axis scales for each variable
    plot <- ggplot(figDAT_long, aes(x = Interval, y = value)) +
      geom_line(aes(group = plotID), colour = "gray", alpha = 0.1) +  # Lines for each plotID
      geom_smooth(aes(colour = NPP_class), method = "lm", formula = y ~ poly(x, 3), se = F, linetype = 1) +
      facet_grid(variable ~ Species, scales = "free_y") +  # Free y-axis scales for each facet
      xlab("Simulation Intervals") +
      ylab("Value") +  # Generic y-axis label that will be overridden by facet labels
      theme_bw(base_size = 10) +
      #coord_cartesian(ylim = c(0, 2000)) +  # You may adjust ylim for each facet individually
      theme(
        strip.text = element_text(size = 11, colour = "black", angle = 0, face = "bold"),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "top"
      ) +
      scale_colour_manual(values = my_colors) +  # Adjust the color palette as needed
      ggtitle(paste0("Data: simDAT", j + k))  # Title that includes data identifier
    
    
    # Add the plot to the list
    plots[[k+1]] <- plot
    
    # remove large objects to free up memory
    rm(simDAT, Results, figDAT, plot)
    gc() # force garbage collection to free up memory
  }
  
  # After collecting the plots for the iteration, print each plot on a new page
  for (p in plots) {
    print(p)  # Directly print the plot to the PDF
  }
}

# Close the PDF device
dev.off()
#dev.off()
