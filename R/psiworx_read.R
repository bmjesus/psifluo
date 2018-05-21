# psiworx_read
#
# - Read in data file
# - Parse metadata
# - Subtract darks from light flashlets (assumed alternating)
# - Split columns into "light" and "dark" curves
# - Split data by sti and str
# 
# Return: 
# all.data  = raw data input, dataframe
# meta.data = acquisiton settings, data frame
# data.sti  = level 1 data, induction only, data frame
# data.str  = level 1 data, relaxation only, data frame
#
# Columns of data.sti and data.str match. 

## License information
# 
# PSIworxR - Data processing from the PSI fluorometer in R
#     Copyright (C) 2017 by Audrey B. Ciochetto, M. Roodvoets, E. J. Austen, and 
#     Douglas A. Campbell. 
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Lesser General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Lesser General Public License for more details.
# 
#     You should have received a copy of the GNU Lesser General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# Contact Audrey B. Ciochetto at audrey.ciochetto@gmail.com

psiworx.read<-function(filename){
  ####  Read in data #### 
  # Note: not all data files are of the same format, I have tried
  # to make this as generic as possible by reading in with scan()
  # and then parsing the data. AB 2014-11-02  
  data.raw<-scan(filename,character(0),sep='\n',quiet=TRUE);
  
  # There are 3 dashed lines:
  # 1) The start of the initial header
  # 2) The end of the inital header, start of acquisiton settings
  # 3) The end of acquisiton settings, start of the data matrix
  # Data start 2 lines after the last dashed separator
  data.sep <- which(data.raw=="-----------------------------------");
  data.start <- grep(pattern="^Selected Data Sets:",x=data.raw);
  data.acq <- grep(pattern="^Data",x=data.raw);
  n.data.sets <- data.sep[2]-data.sep[1]-8; # number of data sets
  
  # Convert scans into data frame
  # Make text data into a list of numeric values
  tst2<-lapply(strsplit(data.raw[(data.sep[length(data.sep)]+2):length(data.raw)],split='\t'),as.numeric);
  # Take list of numeric values and make a data frame, should have n.data.sets+1 b/c one col is time
  all.data<-data.frame(do.call(rbind,lapply(tst2,matrix,ncol=length(tst2[[1]]),byrow=TRUE)));
  
  # Get sample times
  tst2<-lapply(strsplit(data.raw[(data.start+1):(data.sep[length(data.sep)-1]-1)],split='\t'),as.character);
  data.samples<-data.frame(do.call(rbind,tst2));
  
  # Extract Meta-Data
  meta.data.names<-unlist(strsplit(data.raw[data.sep[2]+2],split='\t'),use.names=FALSE);
  #tst<-lapply(strsplit(data.raw[seq(data.sep[2]+3,data.sep[3]-2,3)],split='\t'),as.numeric);
  #meta.data<-data.frame(do.call(rbind,lapply(tst,matrix,ncol=length(tst[[1]]),byrow=TRUE)));
  #names(meta.data)<-meta.data.names[1:length(meta.data.names)];
  # The old way, keeping here for now incase
  tst<-unlist(strsplit(data.raw[data.sep[2]+3],split='\t'));
  meta.data<-as.data.frame(as.list(as.numeric(tst[2:length(tst)])));
  names(meta.data)<-meta.data.names[2:length(meta.data.names)];
  
  #Add column to separate sample data ("light") and blanks ("dark")
  all.data$Cond<-rep(c("Dark","Light"));
  
  # Subtract darks from flashes
  cor.data<-all.data[which(all.data$Cond=="Light"),!names(all.data)=="X1" & !names(all.data)=="Cond"]-all.data[which(all.data$Cond=="Dark"),!names(all.data)=="X1" & !names(all.data)=="Cond"];
  names(cor.data)<-paste("V",ncol(cor.data):1,sep="");
  # Create time column
  cor.data$time<-all.data$X1[which(all.data$Cond=="Light")] * 1000 * 1000;
  # reorder columns so they are in the proper order, first flash is in last column
  cor.data<-cor.data[ncol(cor.data):1];
  
  # find duration of dark period between the two induction curves
  dark.period <- meta.data[1,"Dark_period"] *1000 * 1000;
  
  if(is.na(dark.period)){
    stop('Read function can only handle mutliple curves.');
  } else {
    # split data curves
    # For now, assuming only two curves
    halfway<-nrow(cor.data)/2;
    data.light<-cor.data[1:halfway,];  #the first set of rows (with time < 2) describe the first 'read' in each column. put this in a dataframe of their own
    names(data.light)<-gsub("V","light_",names(data.light));
    data.dark<-cor.data[(halfway+1):nrow(cor.data),];  #the second set (with time > 2) describes the second 'read'. Put these in their own dataframe, too.  The second read is usually better (according to Doug, Mitchell), so work with this one
    names(data.dark)<-gsub("V","dark_",names(data.dark));
    data.dark<-data.dark[,!(names(data.dark)=="time")];
    
    #data.light<-cor.data[which(cor.data$time<dark.period),];  #the first set of rows (with time < 2) describe the first 'read' in each column. put this in a dataframe of their own
    #names(data.light)<-gsub("V","light",names(data.light));
    #data.dark<-cor.data[which(cor.data$time>dark.period),];  #the second set (with time > 2) describes the second 'read'. Put these in their own dataframe, too.  The second read is usually better (according to Doug, Mitchell), so work with this one
    #names(data.dark)<-gsub("V","dark",names(data.dark));
    #data.dark<-data.dark[,!(names(data.dark)=="time")];
    
    # Combine data frames
    merge.data<-cbind(data.light,data.dark);

    # split data set into STI and STR
    sti.length<-meta.data[1,"ST1_duration"]*1000*1000;
    data.sti<-merge.data[which(merge.data$time<=sti.length),];
    data.str<-merge.data[which(merge.data$time>sti.length),];
    data.str$time<-data.str$time-sti.length;
  }
  
  ##Bruno's code
  #to export all the above objects to the workspace at the end of the import function
  all.data<<-all.data
  meta.data<<-meta.data
  data.sti<<-data.sti
  data.str<<-data.str
  data.samples<<-data.samples
  
  e2 <<- new.env()
  
  #rep_names<-assign("rep_names",rep_names, envir = e2)
  
  #this is just a temporary fix to get Doug's code working.
  #the cor.data needs to be created so that it reads only the corrected data instead
  #having twice the columns, i.e.light and dark columns
  no_light_steps<-length(data.sti[1,])-1
  no_light_steps<-assign("no_light_steps",no_light_steps, envir = e2)
  
  #creating the time.data object to imitate Chris code. Not sure why we need this
  #also, I multiply it by 2 because inthis method we have twice the number
  #of columns. The corrected values (dark) are attached to the non-corrected (light)
  dat.num <- read.delim(filename, skip=9, header=F) #must be text file for this to work!! not csv....
  #print(dat.num)
  data.sets <- which(dat.num[,1]=="-----------------------------------")-1
  time.data <- read.delim(filename , skip=9, header=F, nrows=data.sets[1])
  #print(time.data)
  time.data<<-time.data
  
  #nam<-assign("nam",nam, envir = e2)
  #calib_file<-assign("calib_file",calib_file, envir = e2)
  
  
  #return(list(all.data=all.data,meta.data=meta.data,data.sti=data.sti,data.str=data.str,data.samples=data.samples));
  
 
  
  
} # end of psiworx.read