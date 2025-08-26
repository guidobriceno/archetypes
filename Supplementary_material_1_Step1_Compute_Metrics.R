##%######################################################%##
#                                                          #
####                     Parameters                     ####
#                                                          #
##%######################################################%##

SITES = c("Paragominas", "Guaviare", "Cotriguacu")

ROI_BUFFER <- 1000 # meters if ROI CRS is lat/lon, projection units otherwise
MASK_ROI <- TRUE # if TRUE, will mask using ROI after cropping
WRITE_INTERMEDIATE <- TRUE

YEARS <- 1990:2021
AGG_LEVEL <- 30
YEAR_ACTIVENESS <- 2016


##%######################################################%##
#                                                          #
####                    Prepare data                    ####
#                                                          #
##%######################################################%##

library(terra)
library(data.table)
library(stringr)
terraOptions(tempdir="D:/temp")
terraOptions()

for ( i in 1:length(SITES)){
  
  # Site
  Site = SITES[i]
  
  cat("\n Prepare data for", Site, "\n")
  
  # study site
  ROI <- list.files(stringr::str_glue("data/{Site}"), pattern = "shp", full.names = TRUE) # filename or SpatVector
  # Open annual change collection
  ANNUAL_CHANGE <- list.dirs(stringr::str_glue("data/{Site}"), recursive = F, full.names = TRUE) # filename(s) or SpatRaster
  OUT_FOLDER <- stringr::str_glue("results/{Site}")
  
  dir.create(OUT_FOLDER, showWarnings = FALSE, recursive = TRUE)
  
  # ROI
  # ---------------------------------------------------------------- - - -
  if(is.character(ROI)) {
    roi <- vect(ROI)
  } else {
    roi <- ROI
  }
  
  # Crop annual change using ROI
  # ---------------------------------------------------------------- - - -
  
  if(is.character(ANNUAL_CHANGE)) {
    if(length(ANNUAL_CHANGE)==1) {
      terra::vrt(list.files(path= ANNUAL_CHANGE, pattern=".tif$", full.names = T, recursive = T),  
                 paste0(OUT_FOLDER, "/T1.vrt"),options="-separate", overwrite=T)
      
    } else {
      terra::vrt(list.files(path= ANNUAL_CHANGE[1], pattern=".tif$", full.names = T, recursive = T),  
                 paste0(ANNUAL_CHANGE[1], "/T1.vrt"),options="-separate", overwrite=T)
      terra::vrt(list.files(path= ANNUAL_CHANGE[2], pattern=".tif$", full.names = T, recursive = T), 
                 paste0(ANNUAL_CHANGE[2], "/T2.vrt"), options="-separate", overwrite=T)
      terra::vrt(c(paste0(ANNUAL_CHANGE[1], "/T1.vrt"), paste0(ANNUAL_CHANGE[2], "/T2.vrt")), paste0(OUT_FOLDER, "/Annual_change.vrt"), overwrite=T)
    }
  } else {
    annual_change <- ANNUAL_CHANGE
  }
  
  ANNUAL_CHANGE <- list.files(OUT_FOLDER, pattern = "vrt", full.names = T) # filename(s) or SpatRaster
  annual_change <- vrt(ANNUAL_CHANGE)
  # crop annual change using ROI
  roi <- project(roi, annual_change)
  roi_buf <- buffer(roi, ROI_BUFFER)
  set.crs(roi_buf, crs(roi))
  annual_change <- crop(annual_change, roi_buf, mask=MASK_ROI)
  names(annual_change) <- YEARS
  
  # initial forest cover
  # ---------------------------------------------------------------- - - -
  fc_initial <- annual_change[[1]]
  fc_initial <- subst(fc_initial, c(NA, 1, 2), c(NA, 1, 1), others=0)
  names(fc_initial) <- first(YEARS)
  
  # event (deforestation, regrowth, degradation) to year and aggregate 
  # ---------------------------------------------------------------- - - -
  event_to_agg_years <- function(init, stack, code, years) {
    stack <- subst(stack, c(NA, code), c(NA, 1), others=0)
    mask_stack <- sum(stack)
    mask_stack[mask_stack != 0] <- 1
    stack <- app(stack, "which.max")
    stack <- stack + first(years)-1
    stack <- stack * mask_stack
    #stack <-stack *years[33]
    stack <- segregate(stack, years[-1])
    names(stack) <- years[-1]
    aggregate(c(init, stack), AGG_LEVEL, fun="sum", na.rm=FALSE)
  }
  
  agg_def <- event_to_agg_years(fc_initial, annual_change, 3, YEARS)
  agg_deg <- event_to_agg_years(fc_initial, annual_change, 2, YEARS)
  agg_reg <- event_to_agg_years(fc_initial, annual_change, 4, YEARS)
  
  # build raw grid from initial forest cover, deforestation,degradation and regrowth events
  # ---------------------------------------------------------------- - - -
  if(WRITE_INTERMEDIATE) {
    writeRaster(agg_def, file.path(OUT_FOLDER, "agg_def.tif"), overwrite=TRUE)
    writeRaster(agg_deg, file.path(OUT_FOLDER, "agg_deg.tif"), overwrite=TRUE)
    writeRaster(agg_reg, file.path(OUT_FOLDER, "agg_reg.tif"), overwrite=TRUE)
  }
  
  # grid parameters
  # ---------------------------------------------------------------- - - -
  pixel_count <- AGG_LEVEL * AGG_LEVEL
  
  # empty raster to use as template for outputs
  empty_agg <- rast(agg_def, nlyrs=1, vals=NA)
  
  
  
  ##%######################################################%##
  #                                                          #
  ####                  compute metrics                   ####
  #                                                          #
  ##%######################################################%##
  
  cat("\n Compute metrics for", Site, "\n")
  
  compute_metrics <- function(grid_agg, metrics = c("baseline", "forestleft", "forestloss", "percreg", "percdeg", "activeness", "speed")) {
    
    # transform grid into a data.table for further processing
    # one line per pixel/year (wide=FALSE)
    # first year's value is initial forest cover, remaining values are deforestation events
    # values are in number of pixels with a maximum value of pixel_count
    # ---------------------------------------------------------------- - - -
    df_agg <- setDT(as.data.frame(grid_agg, cells=TRUE))
    df_agg <- melt(df_agg, id.vars = "cell", variable.factor = FALSE)
    
    setnames(df_agg, 1:3, c("cell", "year", "values"))
    df_agg[, year := as.integer(year)]
    
    setorder(df_agg, cell, year)
    
    area_agg <- cellSize(grid_agg, unit="km")/pixel_count
    df_agg[, area:=area_agg[cell]]
    
    # precompute some values
    # ---------------------------------------------------------------- - - -
    
    # compute cumulated deforestation
    df_agg[, by=cell, cum_def :=  c(0, cumsum(values[-1]))]
    
    # compute remaining forest
    df_agg[, by=cell, rem_for := first(values) - cum_def]
    
    # compute percentage (ratio ?) of woodland loss
    df_agg[, by=cell, pct_forlost := c(NA, values[-1]/rem_for[-.N])] # x[-.N] ~ x but it's last value
    df_agg[is.nan(pct_forlost), pct_forlost := 0]
    
    # we can remove first "values" now
    df_agg[year==first(YEARS), values := NA]
    
    # and rename it to yeardef
    setnames(df_agg, "values", "year_def")
    
    # Pixel surface
    # area_pix <- setDT(as.data.frame(area_agg, cells=TRUE))
    
    # compute metrics
    # ---------------------------------------------------------------- - - -
    
    results <- df_agg[,by=cell, .N]
    # results[, area_pix := terra::extract(area_agg, results$cell)]
    
    # 1 - baseline woodland
    # ---------------------------------------------------------------- - - -
    
    if ("baseline" %in% metrics) {
      df_baseline <- df_agg[year==first(YEARS)]
      results <- results[df_baseline, on="cell", baseline:=i.rem_for]
      results[, baseline:=100*baseline/pixel_count]
    }
    
    # 2 - percentage of woodland loss
    # ---------------------------------------------------------------- - - -
    if ("forestloss" %in% metrics || "percdeg" %in% metrics) {
      df_pct_loss <- df_agg[, by=cell, .(pct_loss = 100*last(cum_def)/first(rem_for))]
      results <- results[df_pct_loss, on="cell", forest_loss:=i.pct_loss]
    }
    
    # 2 - percentage of reg
    # ---------------------------------------------------------------- - - -
    if ("percreg" %in% metrics) {
      df_pct_gain <- df_agg[, by=cell, .(pct_gain = 100*last(cum_def)/pixel_count)]
      results <- results[df_pct_gain, on="cell", forest_gain:=i.pct_gain]
    }
    
    # 3 - woodland left
    # ---------------------------------------------------------------- - - -
    
    if ("forestleft" %in% metrics){
      df_wl_left <- df_agg[, by=cell, .(wl_left = last(rem_for))]
      results <- results[df_wl_left, on="cell", forest_left:=i.wl_left]
      results[, forest_left:=100*results$forest_left/pixel_count]
      
    }
    
    # 4 - Activeness
    # ---------------------------------------------------------------- - - -
    
    if ("activeness" %in% metrics){
      
      # compute rolling means using a 5 year window
      df_activeness <- df_agg[, by=cell, .(
        year_from = year[-1],  # start and end of current window
        year_to = year[-1] + 4,
        T = 100*frollmean(pct_forlost[-1], 5, align="left")
      )][!is.na(T)] # to remove incomplete window at start/end
      
      # detect
      df_activeness <- df_activeness[, by=cell, .(
        old = any(T[year_to < YEAR_ACTIVENESS] > 0.5),
        active = any(T[year_from <= YEAR_ACTIVENESS & YEAR_ACTIVENESS <= year_to] > 0.5),
        new = any(T[year_from > YEAR_ACTIVENESS] > 0.5)
      )]
      
      # recode
      if(TRUE) { 
        df_activeness[(old & !active & !new), activeness := 1] # Old
        df_activeness[(!old & !active & new), activeness := 3] # Emergent
        df_activeness[(!old & active & new), activeness := 3] # Emergent
        df_activeness[(old & active & new), activeness := 3] # Emergent
        df_activeness[(!old & active & !new), activeness := 2] # Active
        df_activeness[(old & active & !new), activeness := 2] # Active
        df_activeness[(!old & !active & !new), activeness := 0] # Forest
        df_activeness[(old & !active & new), activeness := 3] # Active
      }
     

      
      results <- results[df_activeness,on="cell", activeness:=i.activeness]
    }
    
    # 5 - Speed
    # ---------------------------------------------------------------- - - -
    if ("speed" %in% metrics) {
      df_speed <- df_agg[year > first(YEARS), by=cell, .(
        speed = max(diff(loess(year_def ~ year, span=0.3)$fitted))*area
      )]
      results <- results[df_speed,on="cell", speed:=i.speed]
    }
    
    # results[ ,area_pix:=NULL]
    results[ ,N:=NULL]
    melt(results, id.vars = "cell", variable.factor = FALSE)
  }
  
  def_metrics <- compute_metrics(agg_def, metrics = c("baseline", "forestleft", "forestloss", "activeness", "speed"))
  def_metrics[, source:="def"]
  
  deg_metrics <- compute_metrics(agg_deg, metrics = c("percdeg", "activeness", "speed"))
  deg_metrics[, source:="deg"]
  
  reg_metrics <- compute_metrics(agg_reg, metrics = c("percreg"))
  reg_metrics[, source:="reg"]
  
  metrics_all <- rbindlist(list(def_metrics, deg_metrics, reg_metrics))
  img_list <- split(metrics_all, by=c("source", "variable"))
  
  for(i in seq_along(img_list)) {
    dat <- img_list[[i]]
    img_name <- file.path(OUT_FOLDER, stringr::str_glue("{first(dat$variable)}_{first(dat$source)}.tif"))
    img_r <- rast(empty_agg)
    set.values(img_r, dat$cell, dat$value)
    writeRaster(img_r, img_name, overwrite=TRUE)
  }
  
  
  metrics_all <- dcast(metrics_all, cell ~ source+variable, fill=0)
  metrics_all[, site:= Site]
  saveRDS(metrics_all, file.path(OUT_FOLDER, str_glue("{Site}_metrics.rds")))
}
