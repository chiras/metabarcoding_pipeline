library(circlize)
library(viridis)

create_circos_plot <- function(physeq) {
  # Transpose the OTU matrix for plotting
  otu_matrix = otu_table(physeq)

  # Replace any NA values in the OTU matrix with 0
  otu_matrix[is.na(otu_matrix)] <- 0
  
  # Define all unique sectors for the plot after standardization
  all_sectors <- unique(c(rownames(otu_matrix), colnames(otu_matrix)))
  sectors_length = round(length(all_sectors)/8)
  
  preAllocateTracks = replicate(sectors_length, list(track.height = 0.02), simplify = FALSE)
  #cat("Sectors:\n")
  #print(sectors_length)
  
  # Print all unique sectors to confirm names in the matrix
  #cat("All unique sectors (categories) in the OTU matrix:\n")
  #print(all_sectors)
  
  # Print row names and column names separately for verification
  #cat("\nRow names (Bee species or other row categories):\n")
  #print(rownames(otu_matrix))
  
  #cat("\nColumn names (Plant species or other column categories):\n")
  #print(colnames(otu_matrix))
  
  # Identify and remove any empty rows and columns (those containing only zeros)
  non_zero_rows <- rowSums(otu_matrix) > 0
  non_zero_cols <- colSums(otu_matrix) > 0
  
  # Check for empty rows and columns and remove them
  if (any(!non_zero_rows, na.rm = TRUE)) {
    cat("\nRemoving empty rows (species with no links):\n")
    print(rownames(otu_matrix)[!non_zero_rows])
  }
  if (any(!non_zero_cols, na.rm = TRUE)) {
    cat("\nRemoving empty columns (species with no links):\n")
    print(colnames(otu_matrix)[!non_zero_cols])
  }
  
  # Filter the matrix to keep only non-empty rows and columns
  otu_matrix <- otu_matrix[non_zero_rows, non_zero_cols]
  
  # Define colors for each sector
  grid.col.plants <- shuffle(viridis(length(colnames(otu_matrix))))  # Colors for plant species (columns)
  grid.col.bees <- shuffle(rainbow(length(rownames(otu_matrix))))    # Colors for bee species (rows)
  grid.col <- c(grid.col.bees, grid.col.plants)
  names(grid.col) <- unique(c(rownames(otu_matrix), colnames(otu_matrix)))
  
  # Clear any previous circos settings
  circos.clear()
  
  # Plot the main chord diagram with pre-allocated tracks
  chordDiagram(
    x = otu_matrix,  # Use directly without additional transpose, as we transposed already
    annotationTrack = "grid",
    grid.col = grid.col,
    transparency = 0.4,
    big.gap = 5,
    preAllocateTracks = preAllocateTracks
  )
  
  # Plotting names in sectors for bees (columns)
  trackid <- sectors_length
  for (i in seq_len(ncol(otu_matrix))) {
    trackid <- trackid - 1
    if (trackid < 4) { trackid <- sectors_length }
    myFactor <- colnames(otu_matrix)[i]
    #print(paste("Processing column:", myFactor))
    if (myFactor %in% all_sectors) {
      circos.trackPlotRegion(
        track.index = trackid, factor = myFactor,
        panel.fun = function(x, y) {
          xlim <- get.cell.meta.data("xlim")
          ylim <- get.cell.meta.data("ylim")
          circos.text(
            mean(xlim), ylim[1] + .1, myFactor,
            cex = 1, font = 1, col =  grid.col[myFactor],
            facing = "bending", niceFacing = FALSE
          )
        },
        bg.border = NA
      )
    } else {
      message("Skipping missing sector (column): ", myFactor)
    }
  }
  
  # Plotting names in sectors for plant species (rows)
  trackid <- sectors_length
  for (i in seq_len(nrow(otu_matrix))) {
    myFactor <- rownames(otu_matrix)[i]
    myCol <- "black"
    trackid <- trackid - 1
    if (trackid < sectors_length-5) { trackid <- sectors_length }
    
    #print(paste("Processing row:", myFactor))
    if (myFactor %in% all_sectors) {
      circos.trackPlotRegion(
        track.index = trackid, factor = myFactor,
        panel.fun = function(x, y) {
          xlim <- get.cell.meta.data("xlim")
          ylim <- get.cell.meta.data("ylim")
          circos.text(
            mean(xlim), ylim[1] + .1, myFactor,
            cex = 1.5, col =  grid.col[myFactor], font = 2,
            facing = "bending", niceFacing = TRUE
          )
        },
        bg.border = NA
      )
    } else {
      message("Skipping missing sector (row): ", myFactor)
    }
  }
}
