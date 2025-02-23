specify_decimal <- function(x, k) format(round(x, k), nsmall=k)
inf2NA <- function(x) { x[is.infinite(x)] <- NA; x }
is.integer0 <- function(x) {is.integer(x) && length(x) == 0L}
numtoal <- function(x) {chartr("123456789", "abcdefghij", x)}
myMean03 <- function (a1,k) {
  a2 <- mean(a1,na.rm=TRUE)
  a3 <- sd(a1,na.rm=TRUE)
  a4 <- paste(as.character(specify_decimal(a2,k))," \u00b1 ",as.character(specify_decimal(a3,k)),sep="")
  return(c(a2,a3,a4))
}
GGpltRemLeg <- function (bp) {
  bp <- bp + theme(legend.position="none")
  return (bp)
}
cor.mtest <- function(mat) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  r.mat <- p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(r.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      dat <- cbind(mat[, i], mat[, j])
      is.na(dat) <- !dat
      idx1 <- which(!is.na(dat[,1]) & !is.na(dat[,2]))
      if (length(idx1) > 3) {
        dat <- dat[idx1,]
        tmr <- cor(dat[,1],dat[,2])
        tmp <- cor.test(dat[,1],dat[,2])
        r.mat[i, j] <- r.mat[j, i] <- tmr
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
  }
  colnames(r.mat) <- rownames(r.mat) <- 
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(list(r=r.mat,p=p.mat))
}
DBE_AI_mod <- function(d1) {
  d1$DBE <- (2 + 2*d1$C - d1$H + d1$N + d1$P)/2
  d1$`AI_mod` <- (1 + d1$C - 0.5*d1$O - d1$S - 0.5*d1$H)/(d1$C - 0.5*d1$O - d1$S - d1$N - d1$P)
  return(d1)
}




#####################
# MS_02_Process.r
#####################

PeakIdentify <- function(a2) {
  a3 <- as.data.frame(matrix(nrow=nrow(a2)-1,ncol=2))
  for (j in (2:nrow(a2))) {
    a3[j-1,1] <- as.numeric(a2[j,2])
    a3[j-1,2] <- sum(a2[j,4:ncol(a2)])
  }
  colnames(a3) <- c('time','abundance')
  
  
  # {baseline}
  array1 <- t(a3[,2])
  colnames(array1) <- as.character(t(a3[,1]))
  bs.ar <- baseline(array1)
  bsl <- getBaseline(bs.ar)
  spt <- getSpectra(bs.ar)
  cor <- getCorrected(bs.ar)
  
  
  # peak detection
  p1 <- as.data.frame(matrix(nrow=ncol(cor),ncol=7))
  p1[,1] <- colnames(cor)
  p1[,2] <- t(bsl)
  p1[,3] <- as.numeric(spt)
  p1[,4] <- as.numeric(cor)
  colnames(p1) <- c('time','baseline','abundance','corrected','Lpeak','Hpeak','nvalue')
  idx <- 4
  for (j in 1:nrow(p1)) {
    if(p1[j,idx]<0) {p1[j,]$nvalue=1} else {p1[j,]$nvalue=0}
    if (j==1 | j==nrow(p1)) {
      p1[j,]$Lpeak <- p1[j,]$Hpeak <- 0
    } else {
      if(p1[j,idx] < p1[j-1,idx] & p1[j,idx] < p1[j+1,idx]) {p1[j,]$Lpeak=1} else {p1[j,]$Lpeak=0}
      if(p1[j,idx] > 0 & p1[j,idx] > p1[j-1,idx] & p1[j,idx] > p1[j+1,idx]) {p1[j,]$Hpeak=1} else {p1[j,]$Hpeak=0}
    }
  }
  #plot(p1[121:155,c(1,4)],type='p')
  p2 <- p1[which(p1$Hpeak==1),]
  #addcol <- c('chem','formula','MF','RMF','Prob','CAS','MW','Lib','ID','MajorP','area1','area2')
  #p2[,addcol] <- NA
  
  
  # signal-to-noise
  p3 <- p2[which(p2$corrected>mean(p1$corrected)*3),]
  
  return(p=list("p1"=p1,"p2"=p2,"p3"=p3))
}




# revised from SearchNIST
SearchNIST2 <- function(path02){
  
  nistpath <- "C:\\NISTDEMO\\MSSEARCH\\"
  mspfile <- paste(nistpath,"library.msp",sep='')
  
  #create AUTIMP.MSD file
  secondlocatorpath<-paste(nistpath,"FILESPEC.FIL", sep="")
  #secondlocatorpath<-"FILESPEC.FIL"
  zz<-file(file.path(nistpath,"AUTOIMP.MSD"))
  cat(secondlocatorpath, file = zz, sep = "\n")
  close(zz)
  firstlocatorpath<-file.path(nistpath,"AUTOIMP.MSD")
  
  if(file.exists(file.path(nistpath,"SRCREADY.TXT"))==TRUE){unlink(file.path(nistpath,"SRCREADY.TXT"))}
  
  #create FILSPEC.FIL
  #zz<-file(paste(nistpath,"FILESPEC.FIL", sep=""), "w")
  zz<-file(paste(nistpath,"FILESPEC.FIL",sep=''), "w")
  cat(paste(file.path(mspfile),"Overwrite",sep=" "), file = zz, sep = "\n")
  cat(paste(23,62789), file = zz, sep = "\n")
  close(zz)
  
  
  code<-paste(file.path(nistpath,"nistms$.exe"), " /INSTRUMENT /PAR=2", sep="")
  system(code)
  stopflag=FALSE
  Sys.sleep(0.5)
  while(stopflag==FALSE){
    if(file.exists(file.path(nistpath,"SRCREADY.TXT"))==TRUE){
      stopflag=TRUE
    } else {
      system(code)
      Sys.sleep(0.5)
    }
  }
  
  unlink(file.path(nistpath,"*.HLM"), recursive = FALSE, force = FALSE)
  
  
  r1 <- readLines(file.path(nistpath,"SRCRESLT.TXT"))
  for (j in 2:length(r1)) {
    z1 <- r1[j]
    z2 <- substr(z1,1,gregexpr(';<<',z1)[[1]][1])
    z3 <- gsub('#','',gsub(';',',',z2))
    r1[j] <- paste(z3,substr(z1,gregexpr(';<<',z1)[[1]][1],nchar(z1)),sep='')
  }
  fileConn<-file(file.path(nistpath,"SRCRESLT.TXT"))
  writeLines(r1,fileConn)
  close(fileConn)
  
  
  output <- read.table(file.path(nistpath,"SRCRESLT.TXT"),header=F,skip=1,sep=';',quote="\"",stringsAsFactors=F)
  output2 <- as.data.frame(matrix(nrow=nrow(output),ncol=ncol(output)+1))
  if (nrow(output)!=0) {
    for (j in 1:nrow(output)) {
      output2[j,1] <- trimws(substr(output[j,1],1,gregexpr(':',output[j,1])[[1]][1]-1))
      output2[j,2] <- rm_between(output[j,1],"<<",">>",extract=T)[[1]][1]
      output2[j,3] <- rm_between(output[j,2],"<<",">>",extract=T)[[1]][1]
      output2[j,4] <- gsub('MF: ','',output[j,3])
      output2[j,5] <- gsub('RMF: ','',output[j,4])
      output2[j,6] <- gsub('Prob: ','',output[j,5])
      output2[j,7] <- gsub('CAS:','',output[j,6])
      output2[j,8] <- gsub('Mw: ','',output[j,7])
      output2[j,9] <- gsub('Lib: ','',output[j,8])
      output2[j,10] <- gsub('Id: ','',output[j,9])
    }
    output2[,4] <- as.numeric(output2[,4])
    output2[,5] <- as.numeric(output2[,5])
    output2[,6] <- as.numeric(output2[,6])
    output2[,8] <- as.numeric(output2[,8])
    output2[,10] <- as.numeric(output2[,10])
    
    #idx <- which(output[,5]==max(output[,5]))
    #output2 <- output[idx[1],]
  } else {
    #output2 <- matrix(nrow=9,ncol=1)
  }
  
  return(output2)
}




CheckCompList <- function(complist) {
  check1 <- as.data.frame(matrix(0,nrow=nrow(complist),ncol=2))
  for (i in 1:nrow(complist)) {
    if (complist[i,]$rt1>complist[i,]$rt2) {check1[i,1] <- 1}
    if (i > 1) {if (complist[i,]$rt1 <= complist[i-1,]$rt2) {check1[i,2] <- 1}}
  }
  return(check1)
}




RTSummary <- function(local2,path04,path06,label) {
  fdpos <- gregexpr('/',local2)[[1]]
  foldername <- substr(local2,fdpos[length(fdpos)-1]+1,fdpos[length(fdpos)]-1)
  
  filenames <- list.files(path=path04,pattern='.csv')
  rt_round <- 3
  
  rt1 <- as.data.frame(matrix(nrow=0,ncol=1))
  for (i in 1:length(filenames)) {
    a1 <- filenames[i]
    a2 <- read.csv(paste(path04,a1,sep=''),header=T,stringsAsFactors = F)
    a2 <- a2[!is.na(a2$chem),]
    a3 <- as.data.frame(trimws(as.character(specify_decimal(a2$time,rt_round)),which='both'),stringsAsFactors = F)
    if (nrow(a3)!=nrow(unique(a3))) {print('False')}
    colnames(rt1) <- colnames(a3) <- 'RetentionTime'
    rt1 <- rbind(rt1,a3)
  }
  rt2 <- as.data.frame(rt1[order(rt1[,1]),],stringsAsFactors = F)
  rt3 <- as.data.frame(unique(rt2[,1]),stringsAsFactors = F)
  colnames(rt3) <- colnames(rt2) <- 'RetentionTime'
  
  
  ot1 <- as.data.frame(matrix(nrow=nrow(rt3),ncol=1))
  ot1[,1] <- rt3
  colnames(ot1) <- 'time'
  for (i in 1:length(filenames)) {
    a1 <- filenames[i]
    a2 <- read.csv(paste(path04,a1,sep=''),header=T,stringsAsFactors = F)
    a2 <- a2[!is.na(a2$chem),]
    a2$time <- trimws(as.character(specify_decimal(a2$time,rt_round)),which='both')
    ot1 <- merge(ot1,a2[,c(2,10)],by='time',all.x=T)
    colnames(ot1)[i+1] <- a1
  }
  
  
  ot2 <- as.data.frame(matrix(nrow=nrow(rt3),ncol=12))
  ot2[,1] <- rt3
  colnames(ot2) <- 'time'
  for (i in 1:nrow(ot1)) {
    a1 <- ot1[i,2:ncol(ot1)]
    a2 <- a1[1,!is.na(a1)]
    if (sum(!is.na(a1))!=0) {
      a3 <- as.data.frame(table(t(a2)),stringsAsFactors = F)
      a4 <- a3[order(-a3[,2]),]
      a5 <- as.data.frame(matrix(nrow=1,ncol=nrow(a4)*2))
      for (j in 1:nrow(a3)) {a5[1,(j*2-1):(j*2)] <- a4[j,1:2]}
      
      ot2[i,2] <- nrow(a4)
      if (nrow(a4)<=5) {ot2[i,3:(ncol(a5)+2)] <- a5} else {ot2[i,3:12] <- a5[1,1:10]}
      
    } else {
      ot2[i,2] <- 0
    }
  }
  colnames(ot2) <- c('rt','NumofChem','Chem1','Freq1','Chem2','Freq2','Chem3','Freq3','Chem4','Freq4','Chem5','Freq5')
  aa <- ot2[which(ot2[,2]==0),]
  
  
  ot1[,1] <- as.numeric(ot1[,1])
  ot2[,1] <- as.numeric(ot2[,1])
  ot1 <- ot1[order(ot1[,1]),]
  ot2 <- ot2[order(ot2[,1]),]
  
  # {openxlsx}
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "rt1")
  openxlsx::addWorksheet(wb, "rt2")
  openxlsx::writeData(wb, sheet = "rt1", ot1, colNames = TRUE)
  openxlsx::writeData(wb, sheet = "rt2", ot2, colNames = TRUE)
  openxlsx::saveWorkbook(wb, paste(path06,foldername,'_rt_',label,'.xlsx',sep=''), overwrite = TRUE)
}




majorpeaks <-function(con,p3,j) {
  s0 <- trimws(gsub('>>','',gsub('<<','',p3[j,]$Lib)))
  #if (trimws(p2[j,]$Lib,which=c('both')) == '<<mainlib>>') {
  sqlquery <- paste("select * from raw01NIST where DB=",p3[j,]$ID," and library='",s0,"';",sep='')
  sq1 <- dbSendQuery(con, sqlquery)
  s1 <- fetch(sq1)
  #s1 <- lib[which(lib$DB==p3[j,]$ID & lib$library==s0),]
  
  ind1 <- which(colnames(s1)=='p1')
  
  s2 <- as.numeric(s1[,ind1:(ind1+9)][,which(s1[,(ind1+10):(ind1+19)]>=300)])
  return(s2)
}












