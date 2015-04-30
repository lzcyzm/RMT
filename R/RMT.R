RMT <-
function(
  INPUT_BAM,IP_BAM,
  INPUT2IP = NA,
  EXOME_OUTPUT_DIR = NA,
  GO_OUTPUT_DIR = NA,
  GO_EXPERIMENT_NAME = "GO_Result"){
  # Wrap parameters ##################################################
  PARAMETERS = list();
  PARAMETERS$INPUT_BAM = INPUT_BAM
  PARAMETERS$IP_BAM = IP_BAM
  PARAMETERS$INPUT2IP = INPUT2IP
  PARAMETERS$EXOME_OUTPUT_DIR = EXOME_OUTPUT_DIR
  PARAMETERS$GO_OUTPUT_DIR = GO_OUTPUT_DIR
  PARAMETERS$GO_EXPERIMENT_NAME = GO_EXPERIMENT_NAME
  
  # dependent variables
  if (is.na(PARAMETERS$EXOME_OUTPUT_DIR)) {
    PARAMETERS$EXOME_OUTPUT_DIR <- getwd()
  }
  if (is.na(PARAMETERS$GO_OUTPUT_DIR)) {
    PARAMETERS$GO_OUTPUT_DIR <- getwd()
  }
  
  # algrithm ##################################################
  #peak calling
  input_length <- length(PARAMETERS$INPUT_BAM)
  ip_length <- length(PARAMETERS$IP_BAM)
  all_filepath <- vector()
  if(is.na(PARAMETERS$INPUT2IP)){
    for(i in 1:ip_length){
      experiment_name = PARAMETERS$IP_BAM[i]
      all_filepath[i] = PARAMETERS$IP_BAM[i]
      exomepeak(GENOME = "hg19", IP_BAM = c(PARAMETERS$IP_BAM[i]), 
                INPUT_BAM = c(PARAMETERS$INPUT_BAM[i]), 
                OUTPUT_DIR = PARAMETERS$EXOME_OUTPUT_DIR, 
                POISSON_MEAN_RATIO = 1,
                EXPERIMENT_NAME = experiment_name)
      }
  } else {
    for(i in 1:ip_length){
      experiment_name = PARAMETERS$IP_BAM[i]
      input_bam = PARAMETERS$INPUT_BAM[PARAMETERS$INPUT2IP[[i]][1]]
      if(length(PARAMETERS$INPUT2IP[[i]])>1){
        for(j in 2:length(PARAMETERS$INPUT2IP[[i]])){
          input_bam=c(input_bam,c(PARAMETERS$INPUT_BAM[PARAMETERS$INPUT2IP[[i]][j]]))
        }
      }      
        all_filepath[i] <- paste(PARAMETERS$EXOME_OUTPUT_DIR,experiment_name,sep = '/')
        exomepeak(GENOME = "hg19", IP_BAM = c(PARAMETERS$IP_BAM[i]), 
                  INPUT_BAM = input_bam, 
                  OUTPUT_DIR = PARAMETERS$EXOME_OUTPUT_DIR, 
                  POISSON_MEAN_RATIO = 1, 
                  EXPERIMENT_NAME = experiment_name)   
    }
  }
  
  #merge all peak
  for(i in 1:length(all_filepath)){
    path <- paste(all_filepath[i],"peak.xls",sep = '/')
    bam_file <- read.table(path,sep = "\t",header = FALSE)
    if(i == 1){
      all_peak <- bam_file
    } else {
      all_peak = rbind(all_peak,bam_file[2:nrow(bam_file), ])
    } 
  }
  dir.create(paste(PARAMETERS$GO_OUTPUT_DIR,
                   PARAMETERS$GO_EXPERIMENT_NAME,sep = '/'),
             recursive = TRUE,showWarnings = FALSE)
  dir = paste(PARAMETERS$GO_OUTPUT_DIR,
              PARAMETERS$GO_EXPERIMENT_NAME,sep = '/')
  write.table(all_peak,paste(dir,"all_peak.xls",sep = '/'),
              sep = "\t",quote = FALSE,col.names = FALSE,
              row.names = FALSE)
  
  #read all peaks and transform GRangesList
  filepath = paste(dir,"all_peak.xls",sep = '/')
  matrix_peak_read = read.table(filepath,header = TRUE,stringsAsFactors = FALSE)
 
  ori_peak <- xls2Grangeslist(filepath)
  
  #remove overlaps
  overlaps <- findOverlaps(ori_peak,ori_peak)
  matrix_overlap <- matrix(0,length(overlaps),2)
  matrix_overlap[,1] <- overlaps@queryHits
  matrix_overlap[,2] <- overlaps@subjectHits
 
  num <- c(rep(0,length(ori_peak)))
  
  for(i in 1:nrow(matrix_overlap)){
    if(matrix_overlap[i, 1][1][[1]] != matrix_overlap[i, 2][1][[1]]){
      x <- ranges(ori_peak[matrix_overlap[i, 1]])[[1]]@width
      y <- ranges(ori_peak[matrix_overlap[i, 2]])[[1]]@width
      #num[i]=if(x<y) matrix_overlap[i,1][1][[1]] else matrix_overlap[i,2][1][[1]]
      if(x < y){
        num[matrix_overlap[i, 1][1][[1]]] <- matrix_overlap[i, 1][1][[1]] 
      } else if(x > y){
        num[matrix_overlap[i, 2][1][[1]]] <- matrix_overlap[i, 2][1][[1]]
      } else {
        num[matrix_overlap[i, 1][1][[1]]] <- matrix_overlap[i, 1][1][[1]]
        num[matrix_overlap[i, 2][1][[1]]] <- matrix_overlap[i, 2][1][[1]]
      }
    }
  }
  peak <- ori_peak[num == 0]
  tab <- which(num == 0)
  write.table(matrix_peak_read[tab, ],paste(dir,"rm_peak.xls",sep = '/'),sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  
  # read gene annotation from internet
  txdb <- makeTranscriptDbFromUCSC(genome="hg19")
  exonRanges <- exonsBy(txdb, "gene")
  
  #rpkm in all inputbams and ipbams
  rpkm_input <- matrix(0,length(peak),length(PARAMETERS$INPUT_BAM))
  rpkm_ip <- matrix(0,length(peak),length(PARAMETERS$IP_BAM))
  count_input <- matrix(0,length(peak),length(PARAMETERS$INPUT_BAM))
  count_ip <- matrix(0,length(peak),length(PARAMETERS$IP_BAM))
  #rpkm in inutbams
  for(i in 1:length(PARAMETERS$INPUT_BAM)){
    filepath <- PARAMETERS$INPUT_BAM[i]
    aligns <- readBamGappedAlignments(filepath)
    para <- ScanBamParam(what="mapq")
    mapq <- scanBam(filepath, param=para)[[1]][[1]]
    # filter reads with mapq smaller than 30. 
    mapq[is.na(mapq)] <- 255  # Note: mapq "NA" means mapq = 255
    ID_keep <- (mapq >30)
    filtered <- aligns[ID_keep]
    id <- countOverlaps(filtered,exonRanges)
    transcriptome_filtered_aligns <- filtered[id>0]
    counts <- countOverlaps(peak, transcriptome_filtered_aligns)
    count_input[,i] <- counts
    numBases <- sum(width(peak))
    geneLengthsInKB <- numBases / 1000
    millionsMapped <- length(transcriptome_filtered_aligns) / 1000000
    rpm <- counts / millionsMapped
    rpkm_input[,i] <- rpm / geneLengthsInKB
    
  }
  write.table(count_input,paste(dir,"count_input.xls",sep = '/'),sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  write.table(count_input,paste(dir,"rpkm_input.xls",sep = '/'),sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  #rpkm in ipbams
  for(i in 1:length(PARAMETERS$IP_BAM)){
    filepath <- PARAMETERS$IP_BAM[i]
    aligns <- readBamGappedAlignments(filepath)
    para <- ScanBamParam(what="mapq")
    mapq <- scanBam(filepath, param=para)[[1]][[1]]
    # filter reads with mapq smaller than 30. 
    mapq[is.na(mapq)] <- 255  # Note: mapq "NA" means mapq = 255
    ID_keep <- (mapq >30)
    filtered <- aligns[ID_keep]
    id <- countOverlaps(filtered,exonRanges)
    transcriptome_filtered_aligns <- filtered[id>0]
    counts <- countOverlaps(peak, transcriptome_filtered_aligns)
    count_ip[,i] <- counts
    numBases <- sum(width(peak))
    geneLengthsInKB <- numBases / 1000
    millionsMapped <- length(transcriptome_filtered_aligns) / 1000000
    rpm <- counts / millionsMapped
    rpkm_ip[,i] <- rpm / geneLengthsInKB 
  }
  write.table(count_ip,paste(dir,"count_ip.xls",sep = '/'),sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  write.table(count_ip,paste(dir,"rpkm_ip.xls",sep = '/'),sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  
  #return to all conditions
  matrix_log2_fc = matrix(0,length(peak),ip_length)
  if(is.na(PARAMETERS$INPUT2IP)){
    for(i in 1:ip_length){
      matrix_log2_fc[,i] <- log2(rpkm_ip[,i]+0.01)-log2(rpkm_input[,i]+0.01)
    }
  } else {
    for(i in 1:ip_length){
      if(length(PARAMETERS$INPUT2IP[[i]])==1){
        matrix_log2_fc[,i] <- log2(rpkm_ip[,i]+0.01)-log2(rpkm_input[,PARAMETERS$INPUT2IP[[i]]]+0.01)
      } else {
        sum = matrix(0,length(peak),1)
        for(j in 1:length(PARAMETERS$INPUT2IP[i])){
          sum = sum + rpkm_input[,PARAMETERS$INPUT2IP[[i]][j]]
        }
        matrix_log2_fc[,i] <- log2(rpkm_ip[,i]+0.01)-log2(sum + 0.01)
      }
    }
  }
  
  #feature selection
  v <- matrix_log2_fc[,1]
  for(i in 1:nrow(matrix_log2_fc)) {v[i] <- var(matrix_log2_fc[i,])}
  hist(v)
  matrix_log2_fc <- cbind(tab,matrix_log2_fc[1:nrow(matrix_log2_fc),],matrix_peak_read[tab,])
  write.table(matrix_log2_fc,paste(dir,"all_information.xls",sep = '/'),sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
}
