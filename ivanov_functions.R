#source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

# DMP validation function
validate_dmps = function(GRset.funnorm,studyID) {

GRset.funnorm_subset=GRset.funnorm[,which(pData(GRset.funnorm)$studyID==studyID)]
beta=getBeta(GRset.funnorm_subset)

pd=pData(GRset.funnorm_subset)

mod=model.matrix(~as.factor(pd$sex))
svaobj = sva(beta, mod)

mod = model.matrix(~as.factor(pd$sex))
mod=cbind(mod,svaobj$sv)
probe_fit=lmFit(object=beta,design=mod)
eb=eBayes(probe_fit)

slope=probe_fit$coefficients[,2]
intercept=probe_fit$coefficients[,1]
p.value=eb$p.value[,2]
t=eb$t[,2]
fdr=p.adjust(as.vector(eb$p.value[,2]),method='fdr')

probeFit=data.frame(probes=rownames(eb$p.value), slope=slope, intercept=intercept, 
  p.value=p.value, fdr=fdr, t=t)
rownames(probeFit)=NULL

file=paste('/athena/masonlab/scratch/users/nai2008/DNAm/rdas/probeFit_',studyID,'.rda',sep='')
save(probeFit,file=file)

dmps=which(probeFit$fdr<=.05)

return(length(dmps))

}


# Gene Ontology functions

dogo_entrez <- function(names,universe,species="human", goP = 0.01, #input genes (user set and universe) as entrez
  cond=FALSE, ontology = "BP"){
    if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Hs.egREFSEQ2EG
  } else  if (species == "mouse") {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egREFSEQ2EG
  } else if (species == "rat") {
    golib="org.Rn.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Rn.egREFSEQ2EG
  }
  require(GOstats)
  x=names
  x=x[!is.na(x)]

  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = universe,
                annotation = golib,
                ontology = ontology, pvalueCutoff = goP, conditional = cond,
                testDirection="over")
  
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)

}

dogo_RefSeq <- function(names,universe,species="human", goP = 0.01, #input genes (user set and universe) as RefSeq
	cond=FALSE, ontology = "BP"){
    if(species=="human"){
		golib="org.Hs.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Hs.egREFSEQ2EG
  } else  if (species == "mouse") {
		golib="org.Mm.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Mm.egREFSEQ2EG
  } else if (species == "rat") {
		golib="org.Rn.eg.db"
		library(golib,character.only=TRUE)
		gomap= org.Rn.egREFSEQ2EG
	}
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA)) # convert inputted RefSeq genes (user set) to entrez 
  x=x[!is.na(x)]
  Universe=unlist(mget(as.character(universe),gomap,ifnotfound = NA)) # convert inputted RefSeq genes (universe) to entrez  
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))

  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = ontology, pvalueCutoff = goP, conditional = cond,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)

}

dogo_GeneSymbols <- function(names,universe,species="human", goP = 0.01, #input genes (user set and universe) as Gene Symbol
  cond=FALSE, ontology = "BP"){
    if(species=="human"){
    golib="org.Hs.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Hs.egSYMBOL2EG
  } else  if (species == "mouse") {
    golib="org.Mm.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Mm.egSYMBOL2EG
  } else if (species == "rat") {
    golib="org.Rn.eg.db"
    library(golib,character.only=TRUE)
    gomap= org.Rn.egSYMBOL2EG
  }
  require(GOstats)
  x=unlist(mget(as.character(names), gomap,ifnotfound = NA)) # convert inputted Gene Symbols (user set) to entrez 
  x=x[!is.na(x)]
  Universe=unlist(mget(as.character(universe),gomap,ifnotfound = NA)) # convert inputted Gene Symbols (universe) to entrez  
  Universe=unique(c(Universe[!is.na(Universe)],unique(x)))

  params <- new("GOHyperGParams", geneIds = unique(x),
                universeGeneIds = Universe,
                annotation = golib,
                ontology = ontology, pvalueCutoff = goP, conditional = cond,
                testDirection="over")
  ht=hyperGTest(params)
  tab=summary(ht)
  tmp1=geneIdsByCategory(ht)
  tmp1=tmp1[tab[,1]]
  tab$IDs=sapply(tmp1,function(y) paste(names(x)[x%in%y],collapse=";"))
  return(tab)

}

# wrapper for string split and sapply
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

# command that will return % of variance explained by each PC
getPcaVars = function(pca)  signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100

# function that takes 2 arguments x1 and x2, and matches all elements of x1 to ALL elements of x2, and returns an index vector corresponding to x2

match_all = function (x1, x2){

indexes=list()

  for (i in 1:length(x1)){ 

    aa = which( x2 %in% x1[i]) 

    if(length(aa)>0){indexes[[i]] = aa 
      } else if (length(aa)==0) { indexes[[i]] = NA }
  }

oo=unlist(indexes)

return(oo)

}

# function to generate a summary file a from fastqc files

fastqcSummary <- function(files){ 
#files=paths to unzipped fastqc folders 

ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), "[", slot)

txt_files=paste0(files,'/fastqc_data.txt')

  for (i in 1:length(txt_files)){
  
    scan=as.vector(sapply(txt_files[i], function(x) scan(x,"",sep="\n")))

    ## Total Sequences
    ts=scan[grep('Total Sequences',scan)]
    Total_Sequences=as.numeric(ss(ts,'Sequences\t',2))

    ## Per base sequence quality 
    Per_base_sequence_quality=read.table(txt_files[i],sep='\t',comment.char="",skip=grep('Per base sequence quality',scan)+1, nrows=(grep('Per tile sequence quality',scan)-grep('Per base sequence quality',scan)-3))
    if (ncol(Per_base_sequence_quality)==7){
      colnames(Per_base_sequence_quality)=c('Base', 'Mean', 'Median', 'Lower_Quartile', 'Upper_Quartile', '10th_Percentile', '90th_Percentile') 
    }

    ## Adapter Content
    Adapter_Content=read.table(txt_files[i],sep='\t',comment.char="",skip=grep('Adapter Content',scan)+1, nrows=(grep('Kmer Content',scan)-grep('Adapter Content',scan)-3)) 
    if(ncol(Adapter_Content)==4){
    colnames(Adapter_Content)=c('Position', 'Illumina_Universal_Adapter', 'Illumina_Small_RNA_Adapter', 'Nextera_Transposase_Sequence')
    } else if (ncol(Adapter_Content)==5){
    colnames(Adapter_Content)=c('Position', 'Illumina_Universal_Adapter', 'Illumina_Small_RNA_Adapter', 'Nextera_Transposase_Sequence', 'SOLID_Small_RNA_Adapter')
    }

    fastqcSummary=list(Total_Sequences,Per_base_sequence_quality,Adapter_Content)
    names(fastqcSummary)=c('total_sequences','Per_base_sequence_quality','Adapter_Content')

    output_dir=files[i]
    o=paste0(output_dir,'/fastqcSummary.rda')
    save(fastqcSummary,file=o)

  } 

  #if (length(files)==1) { return(fastqcSummary) }

} #end of function


## from Andrew
## cset: result of cpgCollapse()$cobj
## block450: results of blockFinder
## coi = covariate of interest
## N: the number of blocks to plot, default=10
## blockname = name of block track, default='coi'
## filename = where to save plots
## scale, in kb. default = 100
blockPlot = function(cset, blocks450, coi, N=10,
  blockname = "coi", filename=paste0(blockname,"_blocks.pdf"),scale=100,
  showMethPanel = TRUE, showGenePanel=TRUE, showDiffPanel=TRUE, 
  showCancerPanel = TRUE, bty= "o") {
  panels = c(showMethPanel,showGenePanel, showDiffPanel, showCancerPanel)
  
  require(GenomicRanges)

  blocksTable=with(blocks450$table, GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
  colIds=match(c("chr","start","end"),names(blocks450$table))
  mcols(blocksTable)=blocks450$table[-colIds]

  plotRegion = blocksTable[1:N]

  ## annotation based on ensembl
  cat("Loading Annotation.\n")
  load("/home/epi/ajaffe/GenomicStates/GenomicState.Hsapiens.ensembl.GRCh37.p12.rda")
  gs = GenomicState.Hsapiens.ensembl.GRCh37.p12$fullGenome
  oo = findOverlaps(blocksTable, gs)
  anno = split(gs[subjectHits(oo)], queryHits(oo))

  # cancer blocks
  if(showCancerPanel) {
    load("/home/epi/ajaffe/Lieber/Projects/450k/devPaper/cancer_blocks_hansen.rda")
    genome(blocks)="hg19"
    cancerBlocks = blocks
  }

  cat("Ploting.")
  pdf(filename,height=5,width=10)
  par(bty=bty)

  for(i in seq(along=plotRegion)) {
    cat(".")
    r = plotRegion[i]
    tmp=subsetByOverlaps(cset,r)
    tmp450=sort(subsetByOverlaps(blocksTable,r))
    if(showCancerPanel) tmpBsmooth=subsetByOverlaps(cancerBlocks,r)
                      
    beta=getBeta(tmp)
    x=start(tmp)

    ii=cset %over% r
    d=blocks450$coef[ii]
    sd=blocks450$fitted[ii]
    ## which rows
    Index=split(seq_along(coi),coi)
    mns=sapply(Index,function(ind) rowMeans(beta[,ind]))  
    smns=apply(mns,2,function(y) limma::loessFit(y,x,span=.2)$fitted)

    ## paneling, from hector
    mypar(1,1, brewer.name = "Set1")
    par(mar=par()$mar+c(0,3,0,0))
    omar=par()$mar
    cmar=omar
    cmar[1]=.5
    par(mar=cmar)
    
    layout(cbind(1:sum(panels)),height=c(1+2/3,1, 0.75,0.75)[1:sum(panels)])
    if(showMethPanel) {
      matplot(x,beta,col=as.numeric(factor(coi)),type="p",pch=".",
      xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,1))
      matplot(x,smns,col=1:2,type="l",lwd=2.5,add=TRUE,lty=1)
      legend("bottomright",col=seq(along=levels(factor(coi))),
        lty=1,lwd=2,legend=levels(factor(coi)),cex=.8, bty=bty)
        
      segments(min(x), .2, min(x)+scale*1000, .2)
      text(min(x),.05,labels=sprintf("%dkb",scale),pos=4,offset=0)
      axis(side=2,at=c(.2,.5,.8), cex.axis=1.8)
      mtext(sprintf("%s:%d-%d", seqnames(r), start(r), end(r)), side=3)
      mtext("Methylation",side=2, line = 2.5,cex=1.5)

      legend("topright", paste0("fwer = ",
        signif(blocks450$tab$fwer[i],3)),bty=bty)
    }   
    # annotation
    if(showGenePanel) {
      cmar=omar
      if(!is.na(showDiffPanel)) {
        cmar[3]=0.5
      } else {
        cmar[c(1,3)]=c(0.5,0.5)
      }
      par(mar=cmar)

      plot(x,rep(0,length(x)), type="n",ylim=c(-1.5,1.5),yaxt="n",ylab="",
         xlab="",cex.axis = 1.5, cex.lab =1.5,xaxt="n")
      a = as.data.frame(anno[[i]])
      Strand = ifelse(a$strand == "+", 1, ifelse(a$strand=="-", -1, 0))
      Col = ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue","white"))
      Lwd = ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron",0.5,0))
      axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1,cex.axis = 3)
      abline(h=0,lty=3)
      for(k in 1:nrow(a)) {
        polygon(c(a$start[k],a$end[k],a$end[k],a$start[k]),
          Strand[k]/2+c(-0.3,-0.3,0.3,0.3)*Lwd[k],col=Col[k])
      }

      ## by gene
      g = split(a, sapply(a$symbol,"[", 1))
      # g = split(a, a$Symbol)
      s2 = ifelse(sapply(g, function(x) unique(x$strand))=="+",1,-1)
      g = sapply(g, function(x) (max(x$end) - min(x$start))/2 + min(x$start) )
      
      if(length(g) > 0) text(g, y=s2, names(g),font=1,pos=s2+2,cex=0.8)
          
      mtext("Genes",side=2, line = 2.5,cex=1.5)
      if(!showDiffPanel) {
        xtick=pretty(x)
        axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
      }
    }

    ## mean diff
    if(showDiffPanel) {
      cmar=omar
      cmar[c(1,3)]=c(.5,.5)
      par(mar=cmar)

      zz=granges(tmp)

      matplot(x,sd,xaxt="n",ylab="",xlab="",type="n",lty=1,ylim=c(-.6,.6),yaxt="n",pch=21)
      axis(side=2,at=c(-.3,0,.3),labels=c("-.3","0",".3"))
      xtick=pretty(x)
      axis(side=1,at=xtick,labels=sprintf("%.1fMb",xtick / 1e6))
      mtext("Diff",side=2, line = 2.5,cex=1.5)

      ii=which(zz$type=="OpenSea")
      blockgroup=zz$blockgroup[ii]

      blockIndexes=split(seq(along=blockgroup),blockgroup)
      for (ind in blockIndexes) {
        ind=ii[ind]
        lines(x[ind], sd[ind], lwd=2.5,col="black")
      }


      points(x[ii],d[ii],pch=21,cex=1.4,bg="black")
      axis(side=2,at=c(-2,0,2))
      abline(h=0,lty=2,col="black")

      cmar=omar
      cmar[3]=.5
      par(mar=cmar)
      matplot(x,beta,type="n",xaxt="n",yaxt="n",xlab="",
        ylab="",ylim=c(0,2),bty="n")
    }
    
    #  browser()
    if(showCancerPanel) {
      col=ifelse(tmp450$value<0 & tmp450$p.value<.05,"blue",ifelse(tmp450$value>0 & tmp450$p.value<.05,"red","black"))
      rect(start(tmp450),1+1/3,end(tmp450),1+2/3,col=col)
      if(length(tmpBsmooth) > 0)  rect(start(tmpBsmooth),1/3,end(tmpBsmooth),2/3,col=ifelse(tmpBsmooth$direction=="hypo","blue","red"))
      axis(side=2,at=c(.5,1.5),labels=c("Hansen et al.",blockname),las=1,lwd=0)
      legend("bottomleft",pt.bg=c("blue","red"),legend=c("hypo","hyper"),pch=22,cex=.8)
    }
  }
  dev.off()
}


## from Andrew
# "http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg18.txt"
# "http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg19.txt"
dmrPlot = function(regions, p, chr, pos, cluster, genes, coi, gs=gs,
  genomicState, build="hg19", species = "human", Jitter = FALSE,
  number=100,cols=NULL, lines = FALSE, linesSmooth = TRUE,
  title = TRUE, Legend = TRUE, colorRamp = FALSE, 
  meanSmooth=TRUE, plotCpG = TRUE, geneAnno = "gene") {
  
  require(bumphunter)
  require(derfinder)
  gr = GRanges(regions$chr, IRanges(regions$start, regions$end))
  
  if(build == "hg18") {
    cpg.cur = read.table("http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-hg18.txt",
      header = TRUE, as.is=TRUE)
    library("BSgenome.Hsapiens.UCSC.hg18")
  }
  
  if(build == "hg19") {
    cpg.cur = read.table("http://web.stanford.edu/class/bios221/data/model-based-cpg-islands-hg19.txt",
      header = TRUE, as.is=TRUE)
    library("BSgenome.Hsapiens.UCSC.hg19")
  }
  
  if(build == "mm9") {
    cpg.cur = read.table("http://rafalab.jhsph.edu/CGI/model-based-cpg-islands-mm9.txt",
      header = TRUE, as.is=TRUE)
    library("BSgenome.Mmusculus.UCSC.mm9")
  }
    
  ocpgi=data.frame(chr=I(cpg.cur[,1]), 
    start=as.numeric(cpg.cur[,2]), 
    end=as.numeric(cpg.cur[,3]))
  ADD1 = 1500; PAD = 10
  
  gr2 = GRanges(regions$chr, IRanges(regions$start - ADD1, regions$end + ADD1))
  anno = annotateRegions(gr2, gs)$annotationList
  
  # check regions
  if(is.numeric(coi)) {
    groups=coi
    gNames= sort(unique(coi))
  }

  if(is.character(coi) | is.factor(coi)) {
    groups = factor(coi)
    gNames= levels(groups)
  }
  gIndexes=split(1:length(groups),groups)
  
  brewer.n=max(3,min(length(unique(coi)),11))

  
  if(is.null(cols)) {
    mypar(brewer.n=brewer.n)
  } else if(length(cols) == 1 & cols[1] %in% rownames(brewer.pal.info)) {
    mypar(brewer.name=cols, brewer.n=brewer.n)
  } else {
    mypar()
    palette(cols)
  } 
    
  if(colorRamp) {
    pal = colorRampPalette(palette()[3:brewer.n])
    palette(pal(brewer.n))
  }
  
  cat("Plotting.\n")
  for(j in 1:(min(nrow(regions), number))) {
  
    layout(matrix(1:2,ncol=1),heights=c(0.7,0.3))

    # first plot, region
    par(mar=c(0,4.5,0.25,1.1),oma=c(0,0,2,0))
    
    Index=(regions[j,7]-PAD):(regions[j,8]+PAD)
    Index = Index[Index %in% seq(along=cluster)]
    Index=Index[cluster[Index]==regions[j,6]]
    
    # make DMR plots
    if(Jitter) {
      posx = matrix(pos[Index], nc = ncol(p), nr = length(Index),
        byrow=  FALSE)
      posx = t(apply(posx, 1, function(x) jitter(x,amount = 12)))
      
    } else posx = pos[Index]
          
    matplot(posx, p[Index,], ylim = c(0,1),
      ylab = "", xlab = "",xaxt = "n",cex=0.7,
      cex.axis = 1.7, cex.lab = 1.7, pch=21,
      bg = as.numeric(factor(groups)),col="black",
      xlim = range(pos[Index]), yaxt="n")
    axis(2, at = c(0.2, 0.5, 0.8), cex.axis=1.7)

    xx=pos[Index]
    for(k in seq(along=gIndexes)){
      if(length(gIndexes[[k]]) == 1) yy=p[Index,gIndexes[[k]]]
      if(length(gIndexes[[k]]) > 1)   yy=rowMeans(p[Index,gIndexes[[k]]])
      if(meanSmooth) { 
        fit1=loess(yy~xx,degree=1,span=300/(36*length(xx)),
          family="symmetric")
        lines(xx,fit1$fitted,col=k,lwd=2)
      } else  lines(xx,yy,col=k,lwd=2)

    }
    
    mtext("Methylation",side=2, line = 2.5,cex=1.8)

    if(Legend) {
      if(length(unique(coi)) < 4) {
        legend("topleft",legend=gNames,col=1:length(gNames),
        lty=1, lwd = 4,cex=1,bty="n")
      } else {
        legend("topleft",legend=gNames,col=1:length(gNames),
          pch=19, pt.cex = 2,cex=1.1, nc = 6,bty="n")
      } 
    }
        
    abline(v=(regions$start[j]-15),lty=1)
    abline(v=(regions$end[j]+15),lty=1)

    if(title) mtext(paste0(genes$name[j],"; ", 
      genes$distance[j],"bp from tss:",genes$description[j]), outer=T,cex=1.3)
      
    # add cpgs
    if(plotCpG) {
      thechr=as.character(regions$chr[j])
      chrName = strsplit(thechr, "r")[[1]][2]
      chrName = paste("Chromosome",chrName)
      
      start = pos[Index[1]]
      end = pos[Index[length(Index)]]
      ocpgi2=ocpgi[ocpgi$chr%in%unique(as.character(thechr)),]
      
      ##PLOT CPG ISLANDS
      if(species=="human") seq<-Hsapiens[[as.character(thechr) ]]
      if(species=="mouse") seq<-Mmusculus[[as.character(thechr) ]]
      
      subseq<-subseq(seq,start=start,end=end)
      cpgs=start(matchPattern("CG",subseq))+start-1

      if(plotCpG) Rug(cpgs,col="black")

      Index1 = which(ocpgi2[,1]==as.character(thechr) &
           ((ocpgi2[,2] > start & ocpgi2[,2]< end) |
            (ocpgi2[,3] > start & ocpgi2[,3]< end)))
      if(length(Index1)>0) sapply(Index1,function(j) Rug(unlist(ocpgi2[j,2:3]),
             col="darkgreen",lwd=3,side=1))
    }
    
    # plot 3
    ##PLOT GENES
    par(mar=c(3.5,4.5,0.25,1.1))

    plot(0,0, type="n", xlim=range(xx),ylim=c(-1.5,1.5),yaxt="n",ylab="",
       xlab="",cex.axis = 1.5, cex.lab =1.5)
    a = as.data.frame(anno[[j]])
    Strand = ifelse(a$strand == "+", 1, ifelse(a$strand=="-", -1, 0))
    Col = ifelse(a$theRegion == "exon", "blue", ifelse(a$theRegion == "intron", "lightblue","white"))
    Lwd = ifelse(a$theRegion == "exon", 1, ifelse(a$theRegion == "intron",0.5,0))
    axis(2,c(-1,1),c("-","+"),tick=FALSE,las=1,cex.axis = 3)
    abline(h=0,lty=3)
    for(k in 1:nrow(a)) {
      polygon(c(a$start[k],a$end[k],a$end[k],a$start[k]),
        Strand[k]/2+c(-0.3,-0.3,0.3,0.3)*Lwd[k],col=Col[k])
    }
    
    if(sum(a$theRegion=="exon") > 0) {
      e = a[a$theRegion=="exon",]
      s2 = Strand[a$theRegion=="exon"]
      g = unlist(e$symbol)
      g[is.na(g)] = ""
      if(length(g) > 0) text(x=e$start + e$width/2,
        y=s2*0.75, g,font=1,pos=s2+2,cex=1.2)
    }
        
    mtext("Genes",side=2, line = 2.5,cex=1.5)
    mtext(chrName,side=1, line = 2,cex=1.4)

    abline(v=(pos[regions[j,7]]-15),lty=1)
    abline(v=(pos[regions[j,8]]+15),lty=1)
    
  }
}


## finding DMRs

DMRs_v2 = function(input_GRset, DMRs_filename, B) {

# Enable parallelization
require(doParallel)
registerDoParallel(cores = 3)

# Find bumps
library(bumphunter)

pd=pData(input_GRset)
neg_control_PCs=pd[,10:16]
sex=factor(as.vector(pd$sex), levels=c('Female','Male')) # Female = 0; Male=1;

#Arguments for bumphunter
mod=model.matrix(~sex)
mod=cbind(mod, neg_control_PCs)
p=getBeta(input_GRset)

anno=getAnnotation(input_GRset)
chr=as.vector(anno$chr)
pos=as.vector(anno$pos)

bumps = bumphunterEngine(p, mod, chr = chr, 
pos = pos, cutoff= 0.1, nullMethod = "bootstrap",
smooth=TRUE, B=B)

dat=data.frame(chr=bumps$tab$chr,start=bumps$tab$start,end=bumps$tab$end,p.value=bumps$tab$p.value, FWER=bumps$tab$fwer)

save(bumps, dat, file=DMRs_filename)

num_dmrs=length(which(dat$FWER<0.05))

return(num_dmrs)

}

## finding DNAm blocks

blocks_v2 = function(input_GRset, blocks_filename, B) {

library(minfi)

# Enable parallelization
require(doParallel)
registerDoParallel(cores = 3)

# Find blocks
pd=pData(input_GRset)
neg_control_PCs=pd[,10:16]
sex=factor(as.vector(pd$sex), levels=c('Female','Male')) # Female = 0; Male=1;

mod=model.matrix(~sex)
mod=cbind(mod, neg_control_PCs)

cobj=cpgCollapse(input_GRset, what="Beta")

blocks=blockFinder(cobj$object, design=mod, coef = 2, what = 'Beta', cluster=NULL, nullMethod='bootstrap',
cutoff = 0.1, pickCutoff = FALSE, smooth = TRUE, smoothFunction = locfitByCluster,
B = B, verbose = TRUE, bpSpan = 2.5 * 10^5)

dat=out=data.frame(chr=blocks$tab$chr,start=blocks$tab$start,end=blocks$tab$end,p.value=blocks$tab$p.value, FWER=blocks$tab$fwer)

save(blocks, dat, file=blocks_filename)

num_blocks=length(which(dat$FWER<0.05))

return(num_blocks)

}







