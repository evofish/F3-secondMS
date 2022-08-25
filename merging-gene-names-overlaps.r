##### Merging the gene annotation with the contig name 
setwd("~/Desktop/F3-all/Pairwise-secondMS-comparisons/")
w<-read.table("~/Desktop/F3-all/annotation/apoly-annot-columns-jul18.txt", header=TRUE, sep='\t')

files <- list.files(pattern="*-sig.csv", full.names=TRUE, recursive=FALSE)

#make sure that the headers match

lapply(files, function(x) {
  y <- read.csv(x, header=TRUE, sep=',') # load file
  # apply function
  out <- merge(w, y, by.y="x" )
  # write to file
   write.csv(out, file=paste(x,"_merged.csv"), quote=TRUE, row.names=FALSE, col.names=TRUE)
}) #this generates the merged files, the name of the files needs to be edited further

##### Example of how to evaluate overlap in DEGs between two categories
#Overlap between A-LCAT and A-LCATaC
x<-read.table('A-LCAt-sig.csv', header=TRUE, sep=',')
y<-read.table('A-LCAtaC-sig.csv', header=TRUE, sep=',')
xy<-merge(x,y, by='x')
nrow(xy) #10 overalp 
write.table(xy, file="overlaps/overlap-A-LCAt-A-LCAtaC.txt", quote=FALSE, sep='\t')
