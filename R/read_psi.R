#' @title Function to import and parse PSI text file
#' @description  the result from this function should be an object that will be used by the fitting functions
#' @param file_name name of the file to process, results will be written in this folder
#' @param out_name name for results files (default "results")
#' @param samples replicate names, if present should be in the form c("name1","name2")
#' otherwise leave at FALSE (default)
#' @param num_steps Manually input the number of light steps in each light curve (needs to be the same for all light curves in a data file). Defaults to 11, which is the maximum
#' @param calib_file select the calibration file
#' @return list all the objects the function returns plus the output files
#' @export
read_psi<-function(file_name,
                   out_name,
                   samples,
                   num_steps,
                   calib_file
){  
  
  
  #for the moment some of fitting settings are here but they should be moved elsewhere
  #this function should only parse the text and prepare the data.frames.  
  
  #name of the input file
  filename = file_name
  
  #input replicate names, e.g. 1 name per complete light curve in the input file - input them as a character vector... OR - set to false if you dont want to name your replicates....
  #e.g. rep_names<-c("LL_EP","A","B")
  rep_names<-samples
  
  
  #Manually input the number of light steps in each light curve (needs to be the same for all light curves in a data file)
  #for the moment I set it up to the maximum, i.e. 11, but it would be nice if it was calculated automatically by counting the number of samples and the number of replicates
  no.light.steps<-num_steps
  
  #input a name for the output files, e.g. Experiment_1 will be used to make e.g. Experiment_1_sti_str_data.csv...
  nam<-out_name
  
  
  
  
  
  # List of constants
  # BJ: are these constants specific of ech machine, i.e. does the user needs to change them?
  # BJ: I added them to the function parameters but they be removed if they are real constants 
  #BJ: moved them to a function outside (seems to work!)
  
  
  
  #install.packages("minpack.lm", dependencies=T) #delete the first # if the R package minpack.lm is not installed
  library("minpack.lm") #turns on the minpack.lm package
  
  
  
  ##Run through the psiworx_read commands to parse the data from the raw Flurowin input.....
  ##### Determine number of data sets there are in the file
  # read FluorWin data file specified above
  skipme<-9;
  dat.num <- read.delim(filename, skip=skipme, header=F) #must be text file for this to work!! not csv....
  
  # the number of data sets is the row number before the first "-------"
  #BJ:this could come as a parameter in the function
  #BJ:I could make a smaller function just to determine the number of datasets
  #that could be called separately or inside the bigger function
  data.sets <- which(dat.num[,1]=="-----------------------------------")-1
  
  # Observation data is stored after 12 standard heading rows plus 5 sub-heading rows for each data set
  #BJ: this could be replaced by adding the number of RLC step in the function parameters
  skip.row <- 12+(data.sets[1]*5)
  
  # Extract observation data from the focus file starting at row determined above
  # all.data <- read.table(filename , sep="", skip=skip.row, header=F);
  all.data <- read.table(filename , sep="\t", skip=skip.row, header=F);
  
  #### Extract Meta-Data From FluorWin Output File ####
  # extract time data from FluorWin output file
  #BJ:this produces factors, any reason why we should extractly immediatly the
  #time variable?
  
  time.data <- read.delim(filename , skip=skipme, header=F, nrows=data.sets[1])
  #print(time.data)
  
  # extract meta.data from FluorWin output file, puts all the info about the run in meta.data vector
  #BJ:this is a bit confusing here, why not create a section in the begginning
  #dedicated to extract the metadata, or even better, a function that does that
  n.runs <- ceiling(data.sets[1]/no.light.steps)
  meta.skip <- 11 + data.sets[1] # this number is where the info about the settings on the PSI start
  meta.data.names <- read.delim(filename , skip=meta.skip, header=T, nrows=1)
  meta.data <- as.data.frame(matrix(nrow=1,ncol=ncol(meta.data.names)))
  names(meta.data) <- names(meta.data.names)
  
  #populate meta data frame
  for(n in 1:n.runs){
    meta.n.skip <- 11+data.sets[1]+(40*(n-1))#for each n, 40 lines of code are before the data needed
    meta.data[n,] <- read.delim(filename , skip=meta.n.skip, header=T, nrows=1)
  }
  
  #********************************modified the script a bit here to check if background measurement taken, if not skip this part, if so, annotate and subtract background*************************#
  
  #check to see if PreFlash is > 0, i.e. a background measurement was made per flashlet (if so, subtract from flashlet, if not, skip and continue):
  
  if(meta.data$PreFlash[1]>0){
    
    #Add column to separate sample data ("measurement") and back-ground measurements ("background")
    all.data$Cond<-rep(c("background", "measurement"));
    
    # Subtract background from flashes (but not V1 which is time column and not Cond column which shows light or     background)
    cor.data<-all.data[which(all.data$Cond=="measurement"),!names(all.data)=="V1" & !names(all.data)=="Cond"]-all.data[which(all.data$Cond=="background"),!names(all.data)=="V1" & !names(all.data)=="Cond"];
    names(cor.data)<-paste("V",ncol(cor.data):1,sep="");
    
  } else {
    all.data$Cond<-rep("measurement")
    cor.data<-all.data[which(all.data$Cond=="measurement"), !names(all.data)=="V1" & !names(all.data)=="Cond"];
    names(cor.data)<-paste("V",ncol(cor.data):1,sep="");
    
  }#end of check for PreFlash > 0 i.e. was there background measurements before each flash?
  
  # Create time column (in micro-seconds)
  cor.data$time<-all.data$V1[which(all.data$Cond=="measurement")] * 1000 * 1000;
  # reorder columns so they are in the proper order, first flash is in last column
  cor.data<-cor.data[ncol(cor.data):1];
  
  
  #make the names of cor.data reflect the replicate that they belong to (e.g. did we do 3 light curves, with 10 steps? therefore 30 columns total)
  #BJ: the number of collumns can came from the user input
  #BJ: I think that we need to force the user to save each replicate separately
  no_curves<-length(which(grepl("V",names(cor.data))))
  no_light_steps<-meta.data$RLC_Levels[1]
  no_reps<-no_curves/no_light_steps
  names<-paste('lc', rep(1:no_reps, each=no_light_steps), '_step',rep(1:no_light_steps), sep='')
  names(cor.data)[2:ncol(cor.data)]<-names
  
  
  #Determine number of sti and str flashlets
  #BJ:sti -induction flashlets
  #str - recovery flashlets
  sti_flashes<-meta.data$ST1_num[1]   #number of st1 flashes
  str_flashes<-nrow(cor.data)-sti_flashes  #number of str flashes
  cor.data$phase<-rep(c('sti','str'),times=c(sti_flashes,str_flashes)) #add in phase (sti or str) identifier column
  
  #add in a colour column based on these identifiers, so we can plot our whole light (s
  cor.data$curve<-rep('light')
  #add an STI or STR phase identifierti+str) and dark (sti+str) traces with phases emphaised
  light_sti_col<-'blue'
  light_str_col<-'red'
  for(i in 1:nrow(cor.data)){
    if(cor.data$curve[i]=='light' & cor.data$phase[i]=='sti'){cor.data$col[i]<-light_sti_col} else {cor.data$col[i]<-light_str_col}
  }
  
  #split into different parts for models below:
  
  #data.sti
  data.light.sti<-cor.data[cor.data$phase=='sti' , which(grepl('step',names(cor.data)))]
  
  data.sti<-cbind(cor.data$time[cor.data$phase=='sti'],data.light.sti)
  names(data.sti)[1]<-'time'
  
  #data.str
  data.light.str<-cor.data[cor.data$phase=='str', which(grepl('step',names(cor.data)))]
  data.str<-cbind(cor.data$time[cor.data$phase=='str']-meta.data[1,"ST1_end"]*1000*1000, #corrected dark time in us
                  data.light.str)
  names(data.str)[1]<-'time'
  
  #BJ: It might be useful to find a way to make go from the reading function to the 
  #fit function without being exported to the workspace
  
  time.data<<-time.data
  all.data<<-all.data
  meta.data<<- meta.data
  cor.data<<-cor.data
  data.sti<<- data.sti
  data.str<<-data.str
  #calib_file<<-calib_file

  ##BJ assigning these object to a separate environment so that they don't conflit with the current environment while simultaneously being acesssible to the other packages.
  e2 <<- new.env()
  
  #rep_names<-assign("rep_names",rep_names, envir = e2)
  no_light_steps<-assign("no_light_steps",no_light_steps, envir = e2)
  
  #nam<-assign("nam",nam, envir = e2)
  #calib_file<-assign("calib_file",calib_file, envir = e2)
 
  
  }


