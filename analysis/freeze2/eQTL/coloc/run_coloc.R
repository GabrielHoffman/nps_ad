
my_colc<-function(input,pair)
 {

 library(coloc)
 
 data_ori<-read.table(input)
 pair_list<-read.table(pair)
# gene_list<-read.table(gene) 
 data_ori[,4]<-data_ori[,3]*data_ori[,4]
 data_ori[,10]<-data_ori[,9]*data_ori[,10]
 for(x in 1:dim(pair_list)[1])
#   for(j in gene_list[,1])
     {
   i<-pair_list[x,1]
   j<-pair_list[x,2]
 cat("explored probe is: ",i, " and explored gene is: ",j,"\n")
 data<-data_ori[which(data_ori[,1]==i & data_ori[,7]==j),]

 if(dim(data)[1]==0)
   {
     next
   }
  if(sum(which(duplicated(data[,2])))>0)
   {
  
     data<-data[-which(duplicated(data[,2])),]
   }
# cat("I am here\n")
# data<-read.table(paste0("./summary_for_",probe,"_",gene))
 data[,3]<-data[,3]*data[,3]
 data[,9]<-data[,9]*data[,9]

# cat("I am here1\n")
 Current.eqtl<-data[,c(2,4,3)]
 Current.eqtl[,4]<-seq(1:dim(data)[1])
# cat("I am here1_1\n")
 Current.eqtl[,5]<-NA
 Current.eqtl[,6]<-data[,5]
 Current.eqtl[,7]<-data[,6]
 Current.eqtl[,8]<-1
# cat("I am here1_2\n")
 colnames(Current.eqtl)<-c("snp","beta","varbeta","position","p","N","MAF","sdY")
# cat("I am here2\n")
 Current.gwas<-data[,c(2,10,9)]
 Current.gwas[,4]<-seq(1:dim(data)[1]) 
 Current.gwas[,5]<-NA
 Current.gwas[,6]<-data[,11]
 Current.gwas[,7]<-data[,12]

 cat("I am here3\n")
 colnames(Current.gwas)<-c("snp","beta","varbeta","position","p","N","s")
Current.gwas$snp <- as.character(Current.gwas$snp)
Current.eqtl$snp <- as.character(Current.eqtl$snp)
                        Current.gwas <- as.list(as.data.frame(Current.gwas)); Current.gwas[["type"]] <- "cc"; str(Current.gwas)
                        Current.eqtl <- as.list(as.data.frame(Current.eqtl)); Current.eqtl[["type"]] <- "quant"; str(Current.eqtl)
                        coloc <- coloc.abf(dataset1 = Current.gwas, dataset2 = Current.eqtl, p1=1e-4, p2=1e-4, p12=1e-6)     ### used for further analysis (default parameters)



cat("Result for: ","probe ",i," and gene ",j," is: ",coloc$summary,"\n")
}

}

Args <- commandArgs()

cat("Useage: R --vanilla --slave --input inpute_file --output methylation_profile_file  < /nv/hp10/jzhao48/scratch/zengbiao/script/BH_FDR_adjust.R\n")

for (i in 1:length(Args)) {
    if (Args[i] == "--input") file_name = Args[i+1]
    if (Args[i] == "--pair") probe = Args[i+1]
}


my_colc(file_name, probe)
