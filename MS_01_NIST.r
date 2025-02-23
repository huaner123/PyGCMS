rm(list = ls())

wants <- c("xlsx",'XML','flux','qdapRegex','openxlsx','RMySQL','stringdist')
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
for (i in (1:length(wants))){library(wants[i],character.only = TRUE)}

source("https://bioconductor.org/biocLite.R")
biocLite(pkgs=c('ChemmineR','ChemmineOB'),suppressUpdates=T,suppressAutoUpdate=T)
library(ChemmineR)
library(ChemmineOB)

if (Sys.info()[1]!="Linux") {
  local <- 'D:/Career_13_Instrument/004PyGCMS'
} else {
  local <- '/media/huan/HC/Job_7_CHOW/20170422_MS' # linux system
}

path01 <- paste(local,'/01Others/01MSP_NIST_MS_Library/',sep="")

except <- c('except','has','local','wants','path01')

# total 191436 in mainlib
# total  28307 in replib
#-------------------------
#       219743

#################################
# Major Peak
#################################

rm(list= ls()[!(ls() %in% except)])
source(paste(local,"/03Rscript/functions.R",sep=''))


filenames1 <- list.files(path=paste(path01,'msp/',sep=''),pattern='.MSP')
for (i in 1:length(filenames1)) {
  a1 <- filenames1[i]
  a2 <- as.data.frame(readLines(paste(path01,'msp/',a1,sep='')),stringsAsFactors=F)
  colnames(a2) <- NULL
  
  pos1 <- gregexpr('_', a1)[[1]][1]
  pos2 <- gregexpr('-', a1)[[1]][1]
  lib <- as.character(substr(a1, 1, pos1-1))
  stt <- as.numeric(substr(a1, pos1+1, pos2-1))
  end <- as.numeric(substr(a1, pos2+1,nchar(gsub('.MSP','',a1))))
  rge <- paste(stt,end,sep='-')
  
  y <- matrix(0,nrow(a2),1)
  for (j in 1:nrow(a2)) { if (as.character(a2[j,1])=="") {y[j,1] <- j}}
  y <- y[which(y!=0),1]
  
  ot1 <- as.data.frame(matrix(nrow=length(y),ncol=26))
  ot1[,1] <- a1
  ot1[,2] <- lib
  ot1[,3] <- rge
  idx <- 4
  for (j in 1:length(y)) {
    if (j==1) {st <- 1} else {st <- y[j-1]+1}
    if (j==length(y)) {en <- nrow(a2)} else {en <- y[j]-1}
    
    b1 <- a2[st:en,1]
    ot1[j,idx] <- gsub('DB#: ','',b1[which(grepl('DB#: ',b1[,1])),1])
    ot1[j,idx+1] <- gsub('Name: ','',b1[which(grepl('Name: ',b1[,1])),1])
    b2 <- gsub('CAS#: ','',b1[which(grepl('CAS#: ',b1[,1])),1])
    ot1[j,idx+2] <- substr(b2,1,as.numeric(gregexpr(';',b2)[[1]])-1)
    
    b3 <- ' '
    ind <- which(grepl('Num Peaks:',b1[,1])) + 1 + (st - 1)
    while (ind <= en) {b3 <- paste(b3,a2[ind,1],sep=' '); ind <- ind + 1}
    pos1 <- as.numeric(gregexpr(';',b3)[[1]])
    b4 <- as.data.frame(matrix(nrow=length(pos1),ncol=2))
    for (l in 1:length(pos1)) {
      if (l==1) {st2 <- 1} else {st2 <- pos1[l-1]+1}
      en2 <- pos1[l]-1
      l1 <- gsub(' ',';',substr(b3,st2,en2))
      while (grepl(';;',l1)) {l1 <- gsub(';;',';',l1)}
      pos2 <- as.numeric(gregexpr(';',l1)[[1]])
      b4[l,1] <- substr(l1,pos2[1]+1,pos2[2]-1)
      b4[l,2] <- substr(l1,pos2[2]+1,nchar(l1))
    }
    b4[,1] <- as.numeric(b4[,1])
    b4[,2] <- as.numeric(b4[,2])
    b5 <- b4[order(-b4[,2]),]
    if (nrow(b5)>=10) {
      ot1[j,(idx+3):(idx+12)] <- t(b5[1:10,1])
      ot1[j,(idx+13):(idx+22)] <- t(b5[1:10,2])
    } else {
      ot1[j,(idx+3):(idx+12)] <- cbind(t(b5[,1]),matrix(0,1,10-nrow(b5)))
      ot1[j,(idx+13):(idx+22)] <- cbind(t(b5[,2]),matrix(0,1,10-nrow(b5)))
    }
  }
  colnames(ot1) <- c('file','library','DBrange','DB','chem','cas',
                    'p1','p2','p3','p4','p5','p6','p7','p8','p9','p10',
                    'i1','i2','i3','i4','i5','i6','i7','i8','i9','i10')
  write.table(ot1,paste(path01,'msp/',gsub('.MSP','.txt',a1),sep=''),sep="\t",row.names=FALSE,col.names = TRUE)
  print(paste('Finished: ',i,'/',length(filenames1),sep=''))
}


filenames2 <- list.files(path=paste(path01,'msp/',sep=''),pattern='.txt')
ot <- as.data.frame(matrix(nrow=0,ncol=26))
for (i in 1:length(filenames2)) {
  a1 <- read.table(paste(path01,'msp/',filenames2[i],sep=''),header=T,stringsAsFactors=F)
  ot <- rbind(ot,a1)
  print(paste('Finished: ',i,'/',length(filenames2),'/',nrow(ot),'/',filenames2[i],sep=''))
}

ot2 <- ot
ot2$DB <- as.numeric(ot2$DB)
ot2 <- ot2[order(ot2$library,ot2$DB),]

aa <- ot2[duplicated(ot2[,c(2,5)]),]

write.table(ot2,paste(path01,'NIST_peaks.txt',sep=''),sep="\t",row.names=FALSE)




#################################
# Chemical structure
#################################

rm(list= ls()[!(ls() %in% except)])
source(paste(local,"/03Rscript/functions.R",sep=''))


filenames1 <- list.files(path=paste(path01,'sdf/',sep=''),pattern='.SDF')
for (i in 1:length(filenames1)) {
  a1 <- filenames1[i]
  pos1 <- gregexpr('_', a1)[[1]][1]
  pos2 <- gregexpr('-', a1)[[1]][1]
  lib <- as.character(substr(a1, 1, pos1-1))
  stt <- as.numeric(substr(a1, pos1+1, pos2-1))
  end <- as.numeric(substr(a1, pos2+1,nchar(gsub('.SDF','',a1))))
  rge <- paste(stt,end,sep='-')
  print(paste('Processing: ',a1,'(',i,'/',length(filenames1),')',sep=''))
  
  
  sdfset <- ChemmineR::read.SDFset(paste(path01,'sdf/',a1,sep=''))
  
  #sm <- as.data.frame(matrix(nrow=length(sdfset),ncol=1))
  #colnames(sm) <- 'smiles'
  #for (j in 1:nrow(sm)) {sm[j,1] <- as.character(sdf2smiles(sdfset[j])[[1]])}
  
  # {ChemmineR}
  propma <- data.frame(file=rep(a1,times=length(sdfset)),
                       library=rep(lib,times=length(sdfset)),
                       DBrange=rep(rge,times=length(sdfset)),
                       DB=stt:end,
                       chem=sdfid(sdfset),
                       #smiles=sm,
                       MF=MF(sdfset, addH=FALSE),
                       MW=MW(sdfset, addH=FALSE),
                       Ncharges=sapply(bonds(sdfset, type="charge"), length),
                       atomcountMA(sdfset, addH=TRUE), 
                       groups(sdfset, type="countMA"), 
                       rings(sdfset, upper=6, type="count", arom=TRUE))
  
  write.table(propma,paste(path01,'sdf/',gsub('.SDF','.txt',a1),sep=''),sep="\t",row.names=FALSE,col.names = TRUE)
}

#sdfset[[0]]
#header(sdfset[[1]])
#atomblock(sdfset[1])
#sdfid(sdfset)[1]
#cid(sdfset)[1]
#unique_ids <- makeUnique(sdfid(sdfset))
#plot(sdfset[1:10], regenCoords=TRUE,print=FALSE)

#sdf.visualize(sdfset[1:100])



filenames2 <- list.files(path=paste(path01,'sdf/',sep=''),pattern='.txt')
col <- as.data.frame(matrix(nrow=0,ncol=0))
for (i in 1:length(filenames2)) {
  a1 <- read.table(paste(path01,'sdf/',filenames2[i],sep=''),header=T,stringsAsFactors=F)
  a2 <- as.data.frame(as.character(colnames(a1)),stringsAsFactors=F)
  colnames(a2) <- 'columnnames'
  col <- rbind(col,a2,stringsAsFactors=F)
  col <- as.data.frame(col[!duplicated(col[,1]),],stringsAsFactors=F)
  colnames(col) <- 'columnnames'
}


ot <- as.data.frame(matrix(0,nrow=0,ncol=nrow(col)))
colnames(ot) <- col[,1]
for (i in 1:length(filenames2)) {
  a1 <- read.table(paste(path01,'sdf/',filenames2[i],sep=''),header=T,stringsAsFactors=F)
  a2 <- as.data.frame(matrix(0,nrow=nrow(a1),ncol=nrow(col)))
  colnames(a2) <- col[,1]
  for (j in 1:ncol(a1)) {a2[,which(col==colnames(a1)[j])] <- a1[,j]}
  ot <- rbind(ot,a2)
  print(paste('Finished: ',i,'/',length(filenames2),sep=''))
}

write.table(ot,paste(path01,'NIST_structure.txt',sep=''),sep="\t",row.names=FALSE)



#################################
# Combine
#################################

rm(list= ls()[!(ls() %in% except)])
source(paste(local,"/03Rscript/functions.R",sep=''))


raw01 <- read.table(paste(path01,'NIST_peaks.txt',sep=''),header=T,stringsAsFactors=F)
raw02 <- read.table(paste(path01,'NIST_structure.txt',sep=''),header=T,stringsAsFactors=F)

raw01$library <- trimws(raw01$library, which = c("both"))
raw01$DBrange <- trimws(raw01$DBrange, which = c("both"))
raw01$chem <- trimws(raw01$chem, which = c("both"))

raw02$library <- trimws(raw02$library, which = c("both"))
raw02$DBrange <- trimws(raw02$DBrange, which = c("both"))
raw02$chem <- trimws(raw02$chem, which = c("both"))

ot01 <- merge(raw01,raw02,by=c('library','DB'),all.x=T)

y <- ot01[ot01$chem.x!=ot01$chem.y,c(2,3,5,27,29,30)]
y <- y[which(!is.na(y[,1])),]
y[,'dis'] <- NA
for (i in 1:nrow(y)) {
  if (nchar(y[i,]$chem.x)>80) {sch <- substr(y[i,]$chem.x,1,80)} else {sch <- y[i,]$chem.x}
  y[i,]$dis <- stringdist(sch,y[i,]$chem.y,method = 'jw')
}
y <- y[order(-y$dis),]


wb <- createWorkbook()
addWorksheet(wb, 'NIST')
addWorksheet(wb, 'Namenoequal')
writeData(wb, sheet = "NIST", ot01, colNames = TRUE)
writeData(wb, sheet = "Namenoequal", y, colNames = TRUE)
saveWorkbook(wb, paste(path01,'NIST.xlsx',sep=''), overwrite = TRUE)


con <- dbConnect(MySQL(), user="root", password="205090", dbname='dom', host="localhost")
dbSendQuery(con, "DROP TABLE IF EXISTS `raw01NIST`;")
#sqlquery <- 'CREATE TABLE `raw01NIST` ('
#for (i in 1:ncol(ot01)) {
#  sqlquery <- paste(sqlquery,'`',colnames(ot01)[i],'` varchar(128)',sep='')
#  if (i!=ncol(ot01)) {sqlquery <- paste(sqlquery,',',sep='')} else {sqlquery <- paste(sqlquery,')',sep='')}
#}
#dbSendQuery(con, sqlquery)
dbWriteTable(con, value = ot01, name = "raw01NIST", overwrite = TRUE )

