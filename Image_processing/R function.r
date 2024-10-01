library(EBImage)
library(imager)
library(dplyr)
## adaptive thresholding in R
matsplitter<-function(M, r, c) {
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  cv
} 
megerdata = function(matlist, along = 5, orignal_x_dimension = 1024) {
  num_blocks  = length(matlist)
  rep_numb = num_blocks/along
  split_list = base::split(1:num_blocks, rep(1:rep_numb, each = along))
  new_mat = matrix(ncol = orignal_x_dimension) 
  for( i in 1:length(split_list)){
    dimension2_extract = split_list[[i]]
    temp_mat = do.call(cbind,matlist[dimension2_extract])
    new_mat = rbind(new_mat, temp_mat)
  }
  return(new_mat[-1,])
}
adaptiveThreshold_custome = function(matrix_image, block_size = c(64, 64), orignal_x_dimension = 2048, offset = 0, offset_x = 0.75) {
  gloabla_thresh = otsu(matrix_image)
  imag_blocks = matsplitter(matrix_image, block_size[1], block_size[2])
  imag_blocks_list <- lapply(seq(dim(imag_blocks)[3]), function(x) imag_blocks[, , x])
  imag_blocks_thresh_list = lapply(imag_blocks_list, function(x) x > min(otsu(x)*offset_x + offset, gloabla_thresh))
  #imag_blocks_thresh_list = lapply(imag_blocks_list, function(x) x > quantile(x,0.4) )
  along_num = orignal_x_dimension / block_size[1]
  adaptive_imag = megerdata(imag_blocks_thresh_list, along = along_num, orignal_x_dimension = orignal_x_dimension)
  return(adaptive_imag)
}
## beads capture rate
beads_rate <- function(dir, SD = 3, test_range = 1, threshold = 0.9) {
  # List .czi files in the specified directory
  figures <- list.files(dir, pattern = "\\.czi$", full.names = TRUE)
  sample_ids <- gsub(".czi", "", basename(figures))
  
  # Initialize the results data frame
  beads_rate <- data.frame(
    row.names = sample_ids,
    Sample_ID = sample_ids,
    Assay = remove_pattern_after_hyphen(figures),
    dropnum = NA,
    beads_num = NA,
    cap_rate = NA
  )
  
  zis <- reticulate::import("czifile")
  
  for (i in test_range:length(figures)) {
    # Load the current image
    sample_a <- figures[i]
    image_loaded <- zis$imread(sample_a)
    
    # Process the second image
    image_fam <- t(as.matrix(image_loaded[2, 1:2048, 1:2048, 1]))
    image_fam <- EBImage::normalize(image_fam)
    EBImage::display(image_fam, method = "raster", all = TRUE)
    
    # Adaptive thresholding
    adapted_image <- adaptiveThreshold_custome(image_fam, block_size = c(512, 512), orignal_x_dimension  = 2048, offset_x = 0.5, offset = 0)
    #EBImage::display(adapted_image, method = "raster", all = TRUE)
    
    # Edge detection
    edges <- cannyEdges(as.cimg(1 - adapted_image), alpha = 1, sigma = 1)
    EBImage::display(as.matrix(edges), method = "raster", all = TRUE)
    
    # Morphological closing
    closed_img <- closing(as.matrix(edges), makeBrush(100, shape = 'disc'))
    EBImage::display(as.matrix(closed_img), method = "raster", all = TRUE)
    
    closed_img[closed_img == 1] <- 0.5
    closed_img[closed_img == 0] <- 1
    closed_img[closed_img == 0.5] <- 0
    EBImage::display(as.matrix(closed_img), method = "raster", all = TRUE)
    
    # Watershed segmentation
    nmask <- EBImage::watershed(distmap(closed_img), makeBrush(1, shape = 'disc'))
    
    # Compute droplet sizes and filter
    size_droplets <- computeFeatures.shape(nmask)
    pass_filter <- rownames(size_droplets)[size_droplets[, "s.radius.max"] < 40 | size_droplets[, "s.radius.sd"] > SD]
    sG <- rmObjects(nmask, index = pass_filter)
    EBImage::display(colorLabels(sG), all = TRUE, method = "raster")
    
    # Compute features of segmented objects
    size_droplets_2 <- computeFeatures.shape(sG)
    # Create a list of thresholded images using lapply
    thresholded_images <- lapply(rownames(size_droplets_2), function(obj) {
      # Create a mask for the object
      obj_mask <- rmObjects(sG, index = (1:nrow(size_droplets_2))[!(rownames(size_droplets_2) %in% obj)])
      
      # Extract and threshold the object image
      obj_image <- image_fam * obj_mask
      thresholded <- thresh(obj_image, w = 30, h = 30)
      
      return(thresholded)  # Return the thresholded image
    })
    
    # Combine thresholded images
    combined_matrix <- Reduce(`+`, thresholded_images)
    EBImage::display(combined_matrix, all = TRUE, method = "raster")
    
    # Compute intensity features
    intensity_2nd <- computeFeatures.basic(sG, combined_matrix)
    detected_2nd <- intensity_2nd[, 1] < threshold 
    pass_filter <- rownames(intensity_2nd)[intensity_2nd[, 1] > threshold]
    sG2 <- rmObjects(sG, index = pass_filter)
    EBImage::display(sG2, all = TRUE, method = "raster")
    
    # Segment and display the results
    normalized_image <- toRGB(image_fam)
    segmented = paintObjects(sG2, normalized_image,col=c(NA, 'blue'),thick = TRUE, opac=c(1, 0.2))
    segmented2 <- paintObjects(sG, segmented, col = c('red'), thick = TRUE, opac = c(1))
    EBImage::display(segmented2, all = TRUE, method = "raster")
    
    # Save the labeled image
    filename <- paste0(gsub("\\.czi", "", sample_a), "_labeled.png")
    dev.print(png, filename = filename, width = dim(segmented)[1], height = dim(segmented)[2])
    
    # Update beads_rate data frame
    Large_number <- nrow(intensity_2nd)
    first_rate <- sum(detected_2nd)
    cap_rate <- first_rate / Large_number
    beads_rate[i, c("dropnum", "beads_num", "cap_rate")] <- c(Large_number, first_rate, cap_rate)
    
    svMisc::progress(i, length(figures))
    Sys.sleep(0.01)
    
    if (i == length(figures)) message("Done!")
  }
  
  return(beads_rate)
}
## 6. Gray loading rate 
gray_loading = function(dir) {
  figures <- list.files(dir)[grep("*.jpg", list.files(dir))]
  input_file <- paste0(dir, figures)
  load_1st_rate =data.frame(sample = gsub(x = figures, pattern = "*.jpg", replacement = ""),
                            droplets_rate = NA)
  for (i in c(1:length(input_file))) {
    test_loading = imager::load.image(input_file[i])
    test_run = test_loading[1:1280,1:720,1,2]
    normalized_image <- EBImage::normalize(test_run, ft = c(0, 1))
    EBImage::display(normalized_image, all=TRUE,method = "raster")
    
    nmask1 = adaptiveThreshold_custome(normalized_image, block_size = c(40, 40), orignal_x_dimension = 720,offset_x = 0.6)
    nmask1 <- EBImage::watershed( distmap(nmask1))
    EBImage::display(nmask1, all=TRUE,method = "raster")
    
    size.droplets <- cbind(computeFeatures.shape(nmask1) )
    pass_filter = rownames(size.droplets)[size.droplets[,"s.radius.max"] > 60 | size.droplets[,"s.radius.max"] < 30 | size.droplets[,"s.radius.min"] < 20|  size.droplets[,"s.radius.sd"] > 4]
    sG1 = rmObjects(nmask1,index =  pass_filter)
    EBImage::display(sG1, all=TRUE,method = "raster")
    
    size.droplets <- computeFeatures.shape(sG1)
    
    load_1st_rate[i,][c(2)] = dim(size.droplets)[1]
    
    normalized_image_over = toRGB(normalized_image)
    segmented = paintObjects(sG1, normalized_image_over, col=c('red', 'blue'),thick = TRUE, opac=c(1, 0.3))
    EBImage::display((segmented),all = TRUE, method = "raster")
    filename = paste0(gsub(x = input_file[i], pattern = "*.jpg", replacement = ""), " labled.png")
    dev.print(png, filename = filename , width = dim(segmented)[1]/2, height = dim(segmented)[2]/2)
    
    svMisc::progress(i, length(input_file))
    Sys.sleep(0.001)
    if (i == length(input_file)) message("Done!")
  }
  return(load_1st_rate)
} 
