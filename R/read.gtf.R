
# read gtf file
.read.gtf <- function(PARAMETERS){
  # download the annotation
  if ((!is.na(PARAMETERS$GENOME)) & (!is.na(PARAMETERS$UCSC_TABLE_NAME)) & is.na(PARAMETERS$GENE_ANNO_GTF) & is.na(PARAMETERS$TXDB)) {
    op <- options(warn = (-1))
    txdb =makeTxDbFromUCSC(genome=PARAMETERS$GENOME,
                                   tablename=PARAMETERS$UCSC_TABLE_NAME)
    options(op)
  }
  
  # use provided annotation data file
  if (!is.na(PARAMETERS$GENE_ANNO_GTF) & is.na(PARAMETERS$TXDB) ) {
    op <- options(warn = (-1))
    txdb=makeTxDbFromGFF(PARAMETERS$GENE_ANNO_GTF,format="gtf")
    options(op)
  }
  
  # use provided annotation data file
  if (!is.na(PARAMETERS$TXDB) ) {
    txdb=PARAMETERS$TXDB
  }

  print("Available keytypes:")
  print(keytypes(txdb))
  
  print("Available columns:")
  print(columns(txdb))
  
  # try the internet method
  op <- options(warn = (-1))
  ID <- keys(txdb, keytype = "GENEID")
  temp <- select(txdb,
                 keys = ID,
                 columns = c("GENEID", "TXID", "TXCHROM", "TXSTART", "TXEND", "TXSTRAND"),
                 keytype = "GENEID")

  temp$feature <- "exon"
  options(op)

  
  
  gtf <- data.frame(
    chr        = temp$TXCHROM,
    feature    = temp$feature,
    start      = as.numeric(temp$TXSTART),
    stop       = as.numeric(temp$TXEND),
    strand     = temp$TXSTRAND,
    gene       = temp$GENEID,
    transcript = temp$TXID,
    stringsAsFactors = FALSE
  )
  gtf <- gtf[complete.cases(gtf[, c("chr", "start", "stop", "strand", "gene")]), ]
  
  # fix duplication issues
  gtf <- gtf[!duplicated(gtf[, c("gene", "start", "stop", "strand")]), ]

  print("Using modified .read.gtf()")
  print(head(gtf))
  # return data
  return(gtf)
}
