#You need these libraries:
library(ggplot2)
library(reshape2)
library(topGO)
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
# Modified version of tAI  version 0.21 (Dec 2016) Mario dos Reis
#Reis, M. D., Savva, R., & Wernisch, L. (2004). Solving the riddle of codon usage preferences: a test for translational selection. 
#Nucleic acids research, 32(17), 5036-5044.
# This file forms part of tAI, a package for the analysis of codon
# usage in DNA coding sequences.
# Function to calculate relative adaptiveness values
get.ws <- function(tRNA,      # tRNA gene copy number
                   s = NULL,  # selective constraints
                   sking)     # super kingdom: 0-eukaryota, 1-prokaryota
{
  # optimised s-values:
  if(is.null(s)) s <- c(0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68, 0.89)
  
  p = 1 - s
  
  # initialise w vector
  W = NULL  # don't confuse w (lowercase) and W (uppercase)
  
  # obtain absolute adaptiveness values (Ws)
  for (i in seq(1, 61, by=4))
    W = c(W,
          p[1]*tRNA[i]   + p[5]*tRNA[i+1],     # INN -> NNT, NNC, NNA
          p[2]*tRNA[i+1] + p[6]*tRNA[i],       # GNN -> NNT, NNC
          p[3]*tRNA[i+2] + p[7]*tRNA[i],       # TNN -> NNA, NNG
          p[4]*tRNA[i+3] + p[8]*tRNA[i+2])     # CNN -> NNG
  
  # check methionine
  W[36] = p[4]*tRNA[36]
  
  # if bacteria, modify isoleucine ATA codon
  if(sking == 1) W[35] = p[9]
  
  # get rid of stop codons (11, 12, 15) and methionine (36)
  W = W[-c(11,12,15,36)]
  
  # get ws
  w = W/max(W)
  
  if(sum(w == 0) > 0) {
    ws <- w[w != 0] # zero-less ws
    gm <- exp(sum(log(ws))/length(ws)) # geometric mean
    w[w == 0] = gm # substitute 0-ws by gm
  }
  
  return(w)
}
# After calculating all w's, the tRNA adaptation index (tAI) can then be
# calculated. x is a matrix with all the codon frequencies per ORF in any
# analysed genome
#' The tRNA adaptation index
#' 
#' Calculates the tRNA adaptation index (tAI) of dos Reis et al. (2003, 2004).

get.tai <- function(x, #an n by 60 matrix of codon frequencies for n open reading frames.
                    w) #a vector of length 60 of relative adaptiveness values for codons. It's the output of get.ws 
  
{
  
  w = log(w)              #calculate log of w
  n = apply(x,1,'*',w)    #multiply each row of by the weights
  n = t(n)                #transpose
  n = apply(n,1,sum)      #sum rows
  L = apply(x,1,sum)      #get each ORF length
  tAI = exp(n/L)          #get tai
  return(tAI)
}

# After calculating tAI, we can then calculate R_ENC', i.e. the correlation
# between tAI and Nc' John Novembre (Feb 2006)
#Novembre, J. A. (2002). Accounting for background nucleotide composition when measuring codon usage bias. 
#Molecular biology and evolution, 19(8), 1390-1394.
#' Correlation between tAI and Nc adjusted (Nc')
get_renc <- function(tAI, #a vector of length n with tAI values for genes
                     ncp)  #a vector of length n with Nc' values for genes
{
  ts <- cor.test(tAI, ncp, use='p')
  return(ts)
}

# Statistical test for tAI
#' Monte Carlo test of correlation between tAI and Nc adjusted (Nc')
#' 
#' Calculates the p-value (using a Monte Carlo or randomisation test) that the 
#' correlation (the R_ENC' value) between tAI and the adjusted Nc (Nc') for a set of 
#' genes is different from zero.
#' The Monte Carlo test is described in dos Reis et al. (2004). When
#' working with complete genomes, matrix \code{m} can have a very large number 
#' of rows (large k). In this case it may be advisable to choose \code{samp.size}
#' < k to speed up the computation.
#' 
#' @return A list with elements \code{p.value}, the p-value for the test, and 
#' \code{ts.simulated}, a vector of length \code{n} with the simulated 
#' correlations between tAI and adjusted Nc.
ts.test <- function(m, #a k by 60 matrix of codon frequencies for k genes
                    ws, #ws vector of length 60 of relative adaptiveness values of codons
                    ncp, #ncp vector of length k of Nc values for genes
                    ts.obs, #vector of length 1 with observed correlation between tAI and Nc adjusted for the k genes
                    samp.size, #a vector of length 1 with the number of genes to be sampled from m (see details)
                    n=1000) #n the number of permutations of ws in the randomisation test
{
  
  # create a matrix of randomly permuted w-values:
  ws.permuted <- matrix(rep(ws, n), ncol = 60, byrow = T)
  ws.permuted <- t(apply(ws.permuted, 1, sample))
  
  # initialise 'translational selection' (ts) vector:
  ts = numeric(length(ws.permuted[,1]))
  # work with a smaller m matrix to reduce computation time:
  samp <- sample(nrow(m), samp.size)
  m.samp <- m[samp,]
  ncp <- ncp[samp]
  
  # calculate simulated ts values:
  for(i in 1:length(ws.permuted[,1])) {
    tai <- get.tai(m.samp, ws.permuted[i,])
    co <- cor(ncp, tai, use = 'p')
    ts[i] <- co
  }
  
  # p-values is:
  p.value <- sum(ts < ts.obs)/n

  return(p.value)
}
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
#Function to calculate R_ENC' and other metrics from a list of genomes
#
#genomes<-read.table("genomes.txt", header = T, stringsAsFactors=FALSE)
#The table needs to have a column containing ids named "FTP"

RENC<-function(genomes, #a table containing ids of the genomes
               type,  #if prokaryote 1, eukaryote 0
               output) #name of the output table
{
  R_ENC<-NULL #vector containing R_ENC' values
  ENCp<-NULL #vector containing R_ENC' p.values
  mean_tAI<-NULL #vector containing the mean of tAI
  sd_tAI<-NULL #vector containing the standard deviation of tAI values
  mean_ENc<-NULL #vector containing the mean of Nc'
  mean_Nc<-NULL #vector containing the mean of Nc
  sd_ENc<-NULL #vector containing the standard deviation of Nc' values
  sd_Nc<-NULL #vector containing the standard deviation of Nc values
  trnacn<-NULL #vector containing the number of tRNAs per genome
  
  
  for (i in 1:length(genomes$FTP)){
    trna <- scan(paste("trnas_tables/", genomes$FTP[i],"_genomic.fna.out.trna", sep = ''))
    #reading the output of tRNA-scan
    trnacn[i]<-sum(trna)
    #counting the number of tRNAs
    ws <- get.ws(tRNA=trna, sking=type)
    #calculating weights of each tRNA
    m <- matrix(scan(paste("triplets/", genomes$FTP[i],"_cds_from_genomic.fna.table_triplets.txt", sep = '')), ncol=61, byrow=TRUE)[,-33]
    #reading codon frequencies per each tRNA
    codon_bias<-read.csv(paste("CU_cvs/", genomes$FTP[i],"_cds_from_genomic.fna.table_codon_usage.cvs", sep = ''), header = T, sep="\t", 
                         fill = T, stringsAsFactors = FALSE)
    #reading the table with protein information
    tENC<-read.csv(file = paste("ENC_table/ENC_",genomes$FTP[i],"_table", sep= ''), header = F, sep='')
    colnames(tENC)<-c("Protein_ID","Nc_n","Ncp", "caledChi", "SumChi", "df", "p", "B_KM", "n_codons")
    #reading the table output of ENCprime with the Nc' values
    tENC<-tENC[which(tENC$n_codons>99),]
    #filtering proteins with more than 100 aminoacids
    
    <-codon_bias[!is.na(codon_bias$Nc)&&!is.na(codon_bias$GC3s),]
    #deleting NAs
    cods$Nce<-(-6) + cods$GC3 + 34/(cods$GC3^2 + (1.025 - cods$GC3)^2)
    #calculating the codon usage given the GC3 content
    cods$tai <- get.tai(m, ws)
    #calculating tAI
    cods<-merge(cods, tENC,by="Protein_ID")
    #building a table with all protein information
    write.csv(x=cods, file=paste("CU_table/protein",genomes$FTP[i],"_table", sep= ''))
    #output table
    
    #plotting the nucleotic background and Nc
    png(filename = paste("CU_plots/Nc_vs_GC3_",genomes$FTP[i],".png", sep = ''))
    plot(cods$GC3s,cods$Nc, xlab = "GC3s", ylab="Nc", ylim = c(20,61), xlim=c(0,1), main = genomes$Name[i])
    curve(2+x+(29/(x^2+(1-x)^2)), 0, 1, add = T, col="blue")
    dev.off()
    
    #plotting the nucleotic background and Nc'
    png(filename = paste("CU_plots/Ncp_vs_GC3_",genomes$FTP[i],".png", sep = ''))
    plot(cods$GC3s,cods$Ncp, xlab = "GC3s", ylab="Nc'", ylim = c(20,61), xlim=c(0,1), main = genomes$Name[i])
    curve(2+x+(29/(x^2+(1-x)^2)), 0, 1, add = T, col="blue")
    dev.off()
    
    #R_ENC
    R_ENC[i]<-get_renc(cods$tai,cods$Ncp)
    ENCp[i]<-ts.test(m,ws,cods$Ncp, R_ENC[i],length(cods$Ncp))
    
    mean_tAI[i]<-mean(cods$tai, na.rm = T)
    mean_Nc[i]<-mean(cods$Nc_n, na.rm = T)
    mean_ENc[i]<-mean(cods$Ncp, na.rm = T)
    
    sd_tAI[i]<-sd(cods$tai, na.rm = T)
    sd_Nc[i]<-sd(cods$Nc_n, na.rm = T)
    sd_ENc[i]<-sd(cods$Ncp, na.rm = T)
    
    
    #Plotting the correlation
    x <-cods$Ncp
    y <- cods$tai
    mydf <- data.frame(x = x, y = y)
    
    png(filename = paste("CU_plots/Ncp_vs_tAI_",genomes$FTP[i],".png", sep = ''))
    plot(y ~ x, data = mydf,ylab = "tAI", xlab="Nc'",cex.main=.9 ,main = paste(genomes$Name[i],"\t","r=",round(R_ENC,4),"\t","p.value=",ENCp))
    model <- lm(y ~ x, data = mydf)
    abline(model, col = "red")
    dev.off()
    
    
  }
  
  genomes$R_ENC<-R_ENC
  genomes$ENCp<-ENCp
  genomes$
  genomes$mean_tAI<-mean_tAI
  genomes$mean_Nc<-mean_Nc
  genomes$mean_ENc<-mean_ENc
  genomes$sd_tAI<-sd_tAI
  genomes$sd_dN<-sd_Nc
  genomes$sd_ENc<-sd_ENc
  genomes$trnacn<-trnacn
  
  write.csv(x = genomes, file = output)
  return(genomes)
}

###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
#Function to do a Gene Set Enrichment Analysis and
#Overlap statistic (hypergeometric) analysis
library(topGO) #This function uses topGO package which provides tools for testing GO terms while accounting for the topology of
#the graph.
#Alexa A, Rahnenfuhrer J (2018). topGO: Enrichment Analysis for Gene Ontology. R package version 2.34.0.
CUB_gsea<-function(genomes, #a table containing ids of the genomes
                   test,  #"gsea" for Kolmogorov-Smirnov, 
                   #"quart" for fisher's exact test using ~25% of the genes as threshold  
                   #"ec" for fisher's exact test usign E. coli's ribosomal proteins as threshold
                   out_put_cub) #name of the output table
{
  genome_anot<-data.frame(0,0,0,0,0)
  j<-1
  for (i in 1:length(genomes$FTP)){
    tryCatch({
      geneID2GO<-NULL
      #reading the GO annotations
      geneID2GO <- readMappings(file = paste("GO_table/GO_",genomes$FTP[i],"_table_map.txt", sep = ""))
      info_gene<-NULL
      #reading the protein information (ids, Nc' and tAI)
      info_gene<-read.table(file = paste("CU_table/protein_",genomes$FTP[i],"_table", sep = ""), sep =",", header = T)
      #counting hypothetical proteins
      nhp<-nrow(info_gene[grep("hypothetical",info_gene$NAME),])
      
      gonames<-names(geneID2GO)
      #parsing the protein ids
      fonames<-sub("\\..*","",info_gene$Protein_ID)
      #counting the annotated genes 
      num_gene_anot<-length(intersect(gonames,fonames))
      #counting the total of genes
      num_gene<-length(fonames)
      #estimating how well the genome was annotated
      porncent_anot<-num_gene_anot/num_gene
      genome_anot[j,]<-c(genomes$FTP[i],num_gene_anot,porncent_anot,nhp,num_gene)
      j<-j+1
      
      genesOfInterest<-data.frame(names=info_gene$Protein_ID,tai=info_gene$tai,ncp=info_gene$Ncp, ant=info_gene$NAME)
      
      #fitting the cub values to a Z distribution with mean 0 and sd 1
      #excluding outliers
      outliers <- boxplot(genesOfInterest$tai, plot=FALSE)$out
      tais_o<-genesOfInterest$tai
      if (length(outliers)>0){
        tais_o<-genesOfInterest$tai[-which(genesOfInterest$tai %in% outliers)]
      }
      outliers <- boxplot(genesOfInterest$ncp, plot=FALSE)$out
      ncps_o<-genesOfInterest$ncp
      if (length(outliers)>0){
        ncps_o<-genesOfInterest$ncp[-which(genesOfInterest$ncp %in% outliers)]
      }
      genesOfInterest$tai<-(genesOfInterest$tai-mean(tais_o))/sd(tais_o)
      genesOfInterest$ncp<-(genesOfInterest$ncp-mean(ncps_o))/sd(ncps_o)
      
      
      ##plot
      makeTransparent<-function(someColor, alpha=100)
      {
        newColor<-col2rgb(someColor)
        apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                                    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
      }
      
      ###Kolmogorov-Smirnov   
      #####
      if(test=="gsea"){
        #KS_ncp
        
        geneList <- genesOfInterest$ncp
        
        names(geneList) <- sub("\\..*","",genesOfInterest$names)
        
        myGOdata <- new("topGOdata",
                        description="My project", 
                        ontology="BP",
                        allGenes=geneList,
                        nodeSize = 1,
                        geneSelectionFun= function(x)x,
                        annot = annFUN.gene2GO,
                        gene2GO = geneID2GO)
        
        
        resultks_c_u <- runTest(myGOdata, algorithm="classic", statistic="ks")
        
        
        s<-(score(resultks_c_u))
        
        p_adjust_FDR<-p.adjust(s, method = 'fdr', n = length(s))
        p_adjust_BNF<-p.adjust(s, method = 'bonferroni', n = length(s))
        
        p_val_g<-s[which(s<=0.05)]
        
        padj_FDR<-data.frame(FDR=p_adjust_FDR,GO.ID=names(p_adjust_FDR))
        padj_BNF<-data.frame(Bonferroni=p_adjust_BNF,GO.ID=names(p_adjust_BNF))
        
        allRes <- GenTable(myGOdata,
                           classicFisherU = resultks_c_u,
                           orderBy = "classicFisherU", 
                           topNodes=length(p_val_g)+10)
        
        allRes<-merge(allRes, padj_FDR, by="GO.ID")
        allRes<-merge(allRes, padj_BNF, by="GO.ID")
        allRes_ncp<-allRes
        #ks_tai    
        
        geneList <- (-1)*(genesOfInterest$tai)
        
        names(geneList) <- sub("\\..*","",genesOfInterest$names)
        
        myGOdata <- new("topGOdata",
                        description="My project", 
                        ontology="BP",
                        allGenes=geneList,
                        nodeSize = 1,
                        geneSelectionFun= function(x)x,
                        annot = annFUN.gene2GO,
                        gene2GO = geneID2GO)
        
        
        resultks_c_u <- runTest(myGOdata, algorithm="classic", statistic="ks")
        
        
        s<-(score(resultks_c_u))
        
        p_adjust_FDR<-p.adjust(s, method = 'fdr', n = length(s))
        p_adjust_BNF<-p.adjust(s, method = 'bonferroni', n = length(s))
        
        p_val_g<-s[which(s<=0.05)]
        
        padj_FDR<-data.frame(FDR=p_adjust_FDR,GO.ID=names(p_adjust_FDR))
        padj_BNF<-data.frame(Bonferroni=p_adjust_BNF,GO.ID=names(p_adjust_BNF))
        
        allRes <- GenTable(myGOdata,
                           classicFisherU = resultks_c_u,
                           orderBy = "classicFisherU", 
                           topNodes=length(p_val_g)+10)
        
        allRes<-merge(allRes, padj_FDR, by="GO.ID")
        allRes<-merge(allRes, padj_BNF, by="GO.ID")
        allRes_tai<-allRes
        
        allRes<-merge(allRes_ncp, allRes_tai, by="GO.ID")
        
        write.table(file=paste("topGO/topGO_",genomes$FTP[i],"_ks_txt",sep = ""),sep="\t",allRes,col.names=T, row.names = F)
        
        genesOfInterest$Protein<-rep("other", nrow(genesOfInterest))
        groel<-grep("GroEL",genesOfInterest$ant)  
        ribo<-grep("(riboso|Riboso",genesOfInterest$ant)
        ribo<-ribo[grep("(factor|matur|phosp|methyl|protease|trans|ase|modulation|involved)", genesOfInterest$ant[ribo], invert = T)]
        genesOfInterest$Protein[ribo]<-"RP"
        genesOfInterest$Protein[groel]<-"GroEL"
        
        
        
        ribs<-ggplot(genesOfInterest, aes(ncp,tai)) + 
          geom_point(aes(colour=Protein)) + 
          scale_color_manual(values = c("RP" = "red", "GroEL" = "blue", "other" = makeTransparent("grey")))
        
        svg(file=paste("topGO/Gene_space_ks_",genomes$FTP[i],".svg",sep=""))
        print(ribs)
        dev.off()
        
        pdf(file=paste("topGO/Gene_space_ks_",genomes$FTP[i],".pdf",sep=""))
        print(ribs)
        dev.off()
        
      }
      ###Fisher's exact test using ~25% of the genes as threshold  
      #####
      if(test=="quart"){
        
        porc_gene<-0
        a<-quantile(genesOfInterest$tai, probs=.75)
        b<-quantile(genesOfInterest$ncp, probs=.25)
        while(porc_gene<.25){
          a<-a-.1
          b<-b+.1
          BFs<-genesOfInterest[which(genesOfInterest$tai>a&genesOfInterest$ncp<b),]$names
          porc_gene<-length(BFs)/nrow(genesOfInterest)
        }
        
        BFs<-sub("\\..*","",genesOfInterest[which(genesOfInterest$tai>a &
                                                    genesOfInterest$ncp<b),]$names)
        
        
        geneUniverse <- sub("\\..*","",genesOfInterest$names)
        
        geneList <- factor(as.integer(geneUniverse %in% BFs))
        
        names(geneList) <- geneUniverse
        
        myGOdata <- new("topGOdata",
                        description="My project", 
                        ontology="BP",
                        allGenes=geneList,
                        nodeSize = 1,
                        annot = annFUN.gene2GO,
                        gene2GO = geneID2GO)
        
        #Fisher_BC
        resultFisher_c_u <- runTest(myGOdata, algorithm="classic", statistic="fisher")
        s<-(score(resultFisher_c_u))
        
        p_adjust_FDR<-p.adjust(s, method = 'fdr', n = length(s))
        p_adjust_BNF<-p.adjust(s, method = 'bonferroni', n = length(s))
        
        p_val_g<-s[which(s<=0.05)]
        
        padj_FDR<-data.frame(FDR=p_adjust_FDR,GO.ID=names(p_adjust_FDR))
        padj_BNF<-data.frame(Bonferroni=p_adjust_BNF,GO.ID=names(p_adjust_BNF))
        
        allRes <- GenTable(myGOdata,
                           classicFisherU = resultFisher_c_u,
                           orderBy = "classicFisherU", 
                           topNodes=length(s))
        
        allRes<-merge(allRes, padj_FDR, by="GO.ID")
        allRes<-merge(allRes, padj_BNF, by="GO.ID")
        
        write.table(file=paste("topGO/topGO_",genomes$FTP[i],"_25_txt",sep = ""),sep="\t",allRes,col.names=T, row.names = F)
        
        
        
        genesOfInterest$Protein<-rep("other", nrow(genesOfInterest))
        groel<-grep("GroEL",genesOfInterest$ant)
        ribo<-grep("(riboso|Riboso",genesOfInterest$ant)
        ribo<-ribo[grep("(factor|matur|phosp|methyl|protease|trans|ase|modulation|involved)", genesOfInterest$ant[ribo], invert = T)]
        genesOfInterest$Protein[ribo]<-"RP"
        genesOfInterest$Protein[groel]<-"GroEL"
        
        
        
        ribs<-ggplot(genesOfInterest, aes(ncp,tai)) + 
          geom_point(aes(colour=Protein)) + 
          geom_segment(aes(x =b, y = a, xend = b, yend = max(genesOfInterest$tai)), linetype="dashed") +
          geom_segment(aes(x = b, y = a, xend = min(genesOfInterest$ncp), yend = a), linetype="dashed") +
          stat_smooth(method = "lm", se = FALSE) + 
          scale_color_manual(values = c("RP" = "red", "GroEL" = "blue", "other" = makeTransparent("grey")))
        
        svg(file=paste("topGO/Gene_space_25_",genomes$FTP[i],".svg",sep=""))
        print(ribs)
        dev.off()
        
        pdf(file=paste("topGO/Gene_space_25_",genomes$FTP[i],".pdf",sep=""))
        print(ribs)
        dev.off()
        
        
      }
      
      ###Fisher's exact test using E. coli's ribosomal proteins as threshold 
      #####
      if (test=="ec"){
        
        BFs<-sub("\\..*","",genesOfInterest[which(genesOfInterest$tai>1.1 &
                                                    genesOfInterest$ncp<(-.44)),]$names)
        porc_gene<-length(BFs)/nrow(genesOfInterest)
        
        geneUniverse <- sub("\\..*","",genesOfInterest$names)
        
        geneList <- factor(as.integer(geneUniverse %in% BFs))
        
        names(geneList) <- geneUniverse
        
        myGOdata <- new("topGOdata",
                        description="My project", 
                        ontology="BP",
                        allGenes=geneList,
                        nodeSize = 1,
                        annot = annFUN.gene2GO,
                        gene2GO = geneID2GO)
        
        #Fisher_BC
        resultFisher_c_u <- runTest(myGOdata, algorithm="classic", statistic="fisher")
        s<-(score(resultFisher_c_u))
        
        p_adjust_FDR<-p.adjust(s, method = 'fdr', n = length(s))
        p_adjust_BNF<-p.adjust(s, method = 'bonferroni', n = length(s))
        
        p_val_g<-s[which(s<=0.05)]
        
        padj_FDR<-data.frame(FDR=p_adjust_FDR,GO.ID=names(p_adjust_FDR))
        padj_BNF<-data.frame(Bonferroni=p_adjust_BNF,GO.ID=names(p_adjust_BNF))
        
        allRes <- GenTable(myGOdata,
                           classicFisherU = resultFisher_c_u,
                           orderBy = "classicFisherU", 
                           #topNodes=length(p_val_g))
                           topNodes=length(s))
        
        allRes<-merge(allRes, padj_FDR, by="GO.ID")
        allRes<-merge(allRes, padj_BNF, by="GO.ID")
        
        write.table(file=paste("topGO/topGO_",genomes$FTP[i],"_rib_ec_txt",sep = ""),sep="\t",allRes,col.names=T, row.names = F)
        
        
        genesOfInterest$Protein<-rep("other", nrow(genesOfInterest))
        groel<-grep("GroEL",genesOfInterest$ant)
        ribo<-grep("(riboso|Riboso",genesOfInterest$ant)
        ribo<-ribo[grep("(factor|matur|phosp|methyl|protease|trans|ase|modulation|involved)", genesOfInterest$ant[ribo], invert = T)]
        genesOfInterest$Protein[ribo]<-"RP"
        genesOfInterest$Protein[groel]<-"GroEL"
        
        
        
        ribs<-ggplot(genesOfInterest, aes(ncp,tai)) + 
          geom_point(aes(colour=Protein)) + 
          geom_segment(aes(x =-.44, y = 1.1, xend = -.44, yend = max(genesOfInterest$tai)), linetype="dashed") +
          geom_segment(aes(x = -.44, y = 1.1, xend = min(genesOfInterest$ncp), yend = 1.1), linetype="dashed") +
          stat_smooth(method = "lm", se = FALSE) + 
          scale_color_manual(values = c("RP" = "red", "GroEL" = "blue", "other" = makeTransparent("grey")))
        
        svg(file=paste("topGO/Gene_space_rib_ec_",genomes$FTP[i],".svg",sep=""))
        print(ribs)
        dev.off()
        
        pdf(file=paste("topGO/Gene_space_rib_ec_",genomes$FTP[i],".pdf",sep=""))
        print(ribs)
        dev.off()
      }
      
      #
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
  }
  
  colnames(genome_anot)<-c("FTP","Num_gene_anot","Porcentage_anot","Hypothetical_P", "Num_tot_gene")
  write.csv(x=genome_anot,file = out_put_cub)
  return(genome_anot)
}