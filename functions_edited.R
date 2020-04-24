# function that can be run for each folder "time1" "time2" "time3"
# 
# there's a time-consuming file parsing step that only has to be run once
# if you want to reparse the files for some reason, you can set force_make to TRUE
#
#
extract_data_from_files <- function(folder,force_make=FALSE)
{
  if(!file.exists(paste("results/",folder,".RData",sep="")) | force_make)
  {
    cat("Preparing",folder,"\n")
    
    # find all txt files in the data/ directory
    file_paths <- list.files(path = paste("data/",folder,sep=""),pattern = "txt",full.names = TRUE)
    files <- list.files(path = paste("data/",folder,sep=""),pattern = "txt",full.names = FALSE)
    
    # read in an arbitrary file to see what the wavelengths are (they're the same across all files)
    templine <- readLines(file_paths[1])
    templine <- read.table(text = templine[-(1:grep("<<<",templine))],sep="\t",header = FALSE)
    wavelengths <- templine[,1]
    rm(templine)
    
    
    # read in data from each file, store it in a data fram called dat
    
    dat <- data.frame(matrix(data = character(0),nrow = length(files),ncol = 1),
                      matrix(data = numeric(0),nrow = length(files),ncol = length(wavelengths)),
                      stringsAsFactors = FALSE)
    colnames(dat) <- c("plant_id",paste("X",wavelengths,sep=""))
    
    for(i in 1:length(files))
    {
      file_path <- file_paths[i]
      file <- files[i]
      line <- readLines(file_path)
      line <- read.table(text = line[-(1:grep("<<<",line))],sep="\t",header = FALSE)
      dat[i,1] <- strsplit(file,"_Splice")[[1]][1]
      dat[i,-1] <- line[,2]
    }
    
    save(dat,file = paste("results/",folder,".RData",sep=""))
  } else 
  {
    # if the above step has already been completed, just load the save file instead
    load(paste("results/",folder,".RData",sep=""))
  }
  
  # return hyperspec data frame
  dat
}

get_wavelengths <- function(data)
{
  as.numeric(sapply(strsplit(colnames(data)[-1],"X"),function(X) X[2]))
}

# visualize data
# optional: xlim is a numeric vector of length 2: the minimum and maximum wavelength you wish to visualize
# optional: ylim is a numeric vector of length 2; the minimum and maximum reflectance you wish to visualize
# 
# optional: to_plot is a character vector of length 1 or more; options include "random" "mean" and "all"
# if you select multiple, they will be included
# 
# if you wish to plot an additional random 100, set to_plot="random" and add=TRUE
#
visualize_data <- function(data,to_plot = c("random","mean","all")[1:2],add = FALSE,xlim = NULL,ylim=NULL,nrand_to_plant=200)
{
  wavelengths <- get_wavelengths(data)
  # plot all spectra
  if(!add) plot(expand.grid(range(wavelengths),range(data[,-1],na.rm=TRUE)),type="n",ylim=ylim,xlim=xlim)
  if(any(grepl("random",to_plot,ignore.case = TRUE)))
  {
    cap <- sapply(sample(x = 2:nrow(data[,-1]),size = min(c(nrow(data)-1,nrand_to_plant))),function(X) points(wavelengths,data[,-1][X,],col=rainbow(min(c(nrow(data)-1,nrand_to_plant)))[X],type="l"))
  }
  if(any(grepl("all",to_plot,ignore.case = TRUE)))
  {
    if(nrow(data)>500)
    {
      readline(prompt = "\n\nWarning: This migth take a while and might use a lot of RAM.\nYou might want to_plot='random' instead.\n\nContinue anyway? (press Esc to quit or Enter to continue)")
    }
    cap <- sapply(2:nrow(data[,-1]),function(X) points(wavelengths,data[,-1][X,],col=rainbow(nrow(data[,-1]))[X],type="l"))
  }
  if(any(grepl("mean",to_plot,ignore.case = TRUE)))
  {
    points(wavelengths,colMeans(data[,-1]),type="l",lwd=2,col="black")
  }
}

# removes overlap (severely conflicts) in UV_Vis and NIR scans
eliminate_overlap <- function(data,dominant_spectra_method = c("middle","NIR","UV_Vis")[3],avg_by_id_first=TRUE,resolve = TRUE)
{
  wavelengths <- get_wavelengths(data)
  UV_Vis_inds <- which(diff(wavelengths)<2)
  NIR_inds <- c(which(diff(wavelengths)>2),length(wavelengths))
  
  UV_Vis <- wavelengths[UV_Vis_inds]
  NIR <- wavelengths[NIR_inds]
  
  
  if(grepl("UV_Vis",dominant_spectra_method,ignore.case = TRUE)) 
  {
    NIR <- c(NIR[which(NIR>max(UV_Vis))])
  } else if(grepl("NIR",dominant_spectra_method,ignore.case = TRUE)) 
  {
    UV_Vis <- c(UV_Vis[which(UV_Vis<min(NIR))])
  } else if(grepl("middle",dominant_spectra_method,ignore.case = TRUE))
  {
    UV_Vis_overlap_inds <- UV_Vis_inds[which(UV_Vis>=min(NIR))]
    NIR_overlap_inds <- NIR_inds[NIR<=max(UV_Vis)]
    UV_Vis <- wavelengths[1:round(mean(UV_Vis_overlap_inds)-5)]
    NIR <- wavelengths[round(mean(NIR_overlap_inds)):length(wavelengths)]
  }
  
  data <- data[,c(colnames(data)[1],paste("X",c(UV_Vis,NIR),sep=""))]
  
  if(avg_by_id_first)
  {
    data <- average_dat_by_id(data)
  }
  if(resolve)
  {
    data[,paste("X",NIR,sep="")] <-   data[,paste("X",NIR,sep="")] /( (cbind(data[,paste("X",NIR[1],sep="")] / data[,paste("X",UV_Vis[length(UV_Vis)],sep="")]) %*% rep(1,length(NIR))))
  }
  
  data
}

trim_bad_spectra <- function(data,bad_wavelengths = 180:199)
{
  data[,-unique(unlist(lapply(paste("X",bad_wavelengths,sep=""),function(X) grep(paste(X,"\\.",sep=""),colnames(data)))))]
}  

average_dat_by_id <- function(dat)
{
  # average based on plant_id
  dat <- data.frame(plant_id = unique(dat$plant_id),apply(X = dat[,-(1)],MARGIN = 2,FUN = function(X) tapply(X = X,INDEX = dat$plant_id,FUN = function(X) mean(x = X,na.rm = TRUE))),stringsAsFactors = FALSE)
  dat
}

rescale_between_0_1 <- function(dat,convert_to_absorbance=FALSE)
{
  mindat <- min(dat[,-1],na.rm=TRUE)
  maxdat <- max(dat[,-1],na.rm=TRUE)
  dat <- data.frame(plant_id=dat[,1],apply(dat[,-1,drop=FALSE],2, 
    function(x) ((x-mindat)/(maxdat-mindat))*.999+.001),stringsAsFactors = FALSE)
  if(convert_to_absorbance)
  {
    dat[,-1] <- log(1/dat[,-1])
  }
  dat
}

remove_scans_with_bad_y_values <- function(data,good_y_range = c(0,Inf),drop_scans = FALSE)
{
  drop_inds <- which(apply(data[,-1],1,function(X) any(X<min(good_y_range))))
  drop_inds <- c(drop_inds,which(apply(data[,-1],1,function(X) any(X>max(good_y_range)))))
  drop_inds <- unique(drop_inds)
  if(!drop_scans & length(drop_inds)>0)
  {
    temp_dat <- as.matrix(data[,-1])
    temp_dat[temp_dat<min(good_y_range)] <- NA
    temp_dat[temp_dat>max(good_y_range)] <- NA
    data[,-1] <- temp_dat
    data
  } else if(length(drop_inds)>0) 
  {
    data[-drop_inds,]
  } else
  {
    data
  }
}

smooth_scans <- function(data)
{
  wavelengths <- get_wavelengths(data)
  for(i in 1:nrow(data[,-1]))
  {
    temp_dat <- data[i,-1]
    if(any(is.na(temp_dat)))
    {
      run_dat <- temp_dat[-which(is.na(temp_dat))]
      data[i,-1] <- predict(smooth.spline(x = wavelengths[-which(is.na(temp_dat))],y = run_dat),x = wavelengths)$y
    } else
    {
      run_dat <- temp_dat
      data[i,-1] <- predict(smooth.spline(x = wavelengths,y = run_dat),x = wavelengths)$y
    }
    
  }
  data
}

auto_process <- function(folder)
{
  data <- extract_data_from_files(folder)
  data <- trim_bad_spectra(data = data,bad_wavelengths = 180:229)
  data <- remove_scans_with_bad_y_values(data = data,good_y_range = c(0,Inf))
  data <- trim_bad_spectra(data = data,bad_wavelengths = 2300:3000)
  data <- eliminate_overlap(data = data,dominant_spectra_method = "UV_Vis",avg_by_id_first = TRUE)
  data <- rescale_between_0_1(data)
  data <- smooth_scans(data)
  data
}

write_data_in_raw_read_format <- function(data,filename,sep = ",")
{
  if((sep!=",") & (sep!="\t")) warning("You probably want you sep variable to be either , or \\t")
  dat <- t(data)
  colnames(dat) <- dat[1,]
  dat <- dat[-1,]
  dat <- cbind(as.numeric(sapply(strsplit(x = rownames(dat),split = "X"),function(X) X[2])),dat)
  colnames(dat)[1] <- "wavelength"
  write.table(x = dat,file = paste("results/",filename,sep=""),sep = sep,row.names = FALSE)
}

write_data_in_transposed_format <- function(data,filename,sep = ",")
{
  if((sep!=",") & (sep!="\t")) warning("You probably want you sep variable to be either , or \\t")
  colnames(data) <- c(as.numeric(sapply(strsplit(x = colnames(data),split = "X"),function(X) X[2])))
  colnames(data)[1] <- "plant_id"
  write.table(x = data,file = paste("results/",filename,sep=""),sep = sep,row.names = FALSE)
}
