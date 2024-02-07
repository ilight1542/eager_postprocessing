#!/usr/bin/env Rscript
library(getopt); # for args parsing
library(parallel); # sumstat calculation is parallized across libraries

options(error = quote({
  dump.frames(to.file=T, dumpto='last.dump')
  load('last.dump.rda')
  print(last.dump)
  q()
}))

## FUNCTIONs
extract.stats5 <- function(id,path,malt.mode,paried_end_mode){
    ind <- id
    out <- list()
    for (run in malt.mode){
        ed.dis <- read.table(paste(path,run,'/editDistance/',ind,'_editDistance.txt',sep=""),header=F,sep="\t",skip=1,row.names=1,comment.char='')[,1:6] # keep relevant columns only
        colnames(ed.dis) <- c('0','1','2','3','4','5') # R does not like reading numeric header's
        mp.dam <- read.table(paste(path,run,'/damageMismatch/',ind,'_damageMismatch.txt',sep="" ),header=T,sep="\t",row.names=1,check.names=F,comment.char='')
        rd.dis <- read.table(paste(path,'default','/readDist/',ind,'_alignmentDist.txt',sep="" ),header=T,sep="\t",row.names=1,check.names=F,comment.char='')
        reads <- read.table(paste(path,run,'/RunSummary.txt',sep=''),header=T,sep="\t",row.names=1,check.names=F,comment.char='')
        if( run == 'ancient' ){ rhoNM <- 2:6 ;keptDiff <- 1:2} else { rhoNM <- 1:6; keptDiff <- 1:3}
        for ( spec in unq.spec ){
            ## spec <- 'Yersinia_pestis'
            reads.spec <- as.numeric(reads[ spec , ind ])
            ed.dis.spec <- ed.dis[ spec , ]
            mp.dam.spec <- mp.dam[ spec , ]
            rd.dis.spec <- rd.dis[ spec , ]
            ## get RatioOfDifferences and N (editdistance0-4 and 0-6)
            ## NOTE: Could be shorter if indeed always only 1 Node presented
            res <- matrix(ncol=5,nrow=nrow(ed.dis.spec)); rownames(res) <- rownames(ed.dis.spec); colnames(res) <- paste(run,c('node','dr6','n6','dr4','n4'),sep=".")
            for (subset in rownames(ed.dis.spec)){
                a <- diff(as.numeric(ed.dis.spec[ subset , rhoNM ]))
                dr6 <- round(sum(abs(a[a<0]))/sum(abs(a)),3)
                b <- a[ keptDiff ] # only diffs 1:3 considered. When ancient only diffs 1:2
                dr4 <- round(sum(abs(b[b<0]))/sum(abs(b)),3) ## calculate sum of negative differences to neighboring col to right divided by sum of absolute values of all differences (1 if strictly declining edit distances)
                res[subset, ] <- c( subset, dr6 , sum(ed.dis.spec[ subset , 1:6 ]) , dr4 , sum(ed.dis.spec[ subset , 1:4 ]))
            }
            ## require minimum of 10 reads present dr4 analysis and pick the one with highest number of reads
            ## NOTE: Could be shorter if indeed always only 1 Node presented!
            rowMax <- which(as.numeric(res[,paste(run,'.n4',sep="")])==max(as.numeric(res[,paste(run,'.n4',sep="")])))[1]
            if( !is.na(rowMax) & as.numeric(res[ rowMax , paste(run,'.n4',sep="") ]) > 1 ){
                top.dr <- res[ rowMax , ]
            } else {
                top.dr <- rep(NA,5)
            }
            ## extract map.damage:sum(C>T or G>A pos 1 2 (-1 -2 respectively)) for TopScorer@EdDis
            # if paired end: value def.mapDam and anc.mapDam have propotion of reads with damage at first index of 3' end or last base of 5' end, depending on which is higher
            # if single end: value def.mapDam and anc.mapDam have propotion of reads with damage at first index of 3'
            if (paried_end_mode) {
                mp.dam.spec.max <- max(mp.dam.spec[ rowMax ,"C>T_1" ] , mp.dam.spec[ rowMax ,"G>A_20"])
            } else {
                mp.dam.spec.max <- mp.dam.spec[ rowMax ,"C>T_1" ]
            }
            ## extract max readDis:uniquePerReference for TopScorer@EdDis
            read.dis.uniq <- rd.dis.spec[ rowMax ,'uniquePerReference']
            if(length(read.dis.uniq) == 0){ read.dis.uniq <- NA }
            ## write results list
            if( paste(ind,spec,sep="_") %in% names(out) ){
                out[[ paste(ind,spec,sep="_") ]] <- c( out[[ paste(ind,spec,sep="_") ]] , top.dr , mp.dam.spec.max , read.dis.uniq , reads.spec )
            } else {
                out[[ paste(ind,spec,sep="_") ]] <- c( ind , spec , top.dr , mp.dam.spec.max , read.dis.uniq , reads.spec )
            }
        }
    }
    out2 <- do.call(rbind,out)
    if(length(malt.mode)==2){
        colnames(out2) <- c('id','spec','def.node','def.dr6','def.n6','def.dr4','def.n4','def.mapDam','def.rd', 'reads_all','anc.node','anc.dr6','anc.n6','anc.dr4','anc.n4','anc.mapDam','anc.rd', 'reads_ancient')
    } else {
        colnames(out2) <- c('id','spec','def.node','def.dr6','def.n6','def.dr4','def.n4','def.mapDam','def.rd', 'reads_all')
        }
    return(out2)
}

plot.editDis <- function(id,tax,folders){
    ## function writes barplot EditDistance: default (and ancient)
    switch <- FALSE
    for (i in 1:length(folders)){
        ## Edit Distance: read data and grep taxon(s) of interest
        res <- read.table(paste(folders[i],'editDistance/',id,'_editDistance.txt',sep=''),header=F,sep="\t",skip=1,comment.char='',row.names=1)[,1:12] # R does not get the numbers as headers
        rownames(res) <- chartr("><","..",rownames(res)) # Unusual character fix by JFY, based on plot_summary_rmaex_v05
        colnames(res) <- c('0','1','2','3','4','5','6','7','8','9','10','>10')
        res <- res[ tax , ]
        mc1 <- c(colorRampPalette(c("lightgreen","darkgreen"))(nrow(res))) # remnant from multiple spec
        barplot(as.matrix(res),beside=T,col=mc1,main=id,xlab=paste("edit distance:",basename(folders[i])),ylab="read count")
        legend("top",legend=paste('Node:',tax,sum(res)),fill=mc1, cex = 0.6, bty = "n")
    }
}

plot.mapDamage3 <- function(id,tax,folder){
    ## plot damage pattern
    dam <- read.table(paste(folder,'damageMismatch/',id,'_damageMismatch.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    dam <- dam[ tax , ]
    maxY <- (max(dam[, -dim(dam)[2] ]) * 1.1)
    ## check if mapDamage was calculated at all (apparenlty sometimes everything is 0)
    ## We plot the cumulative damage pattern over all species, normalized by the number of strains that have been detected. The numbers provided by RMAex are fullly corrected per strain!
    ## Plot C>T and G>A in one plot:
    plot("",col="grey",ylim=c(0,maxY), xlim=c(1,20),type="l",lwd=2,ylab="C2T/G2A rate",xlab="Read Position",xaxt='n',main=paste("Damage plot for",tax,"node"))
    lines(x=1:10,dam[tax,21:30],col='grey',lwd=1.5)
    lines(x=11:20,dam[tax,31:40],col='grey',lwd=1.5)
    lines(x=1:10,y=dam[tax,1:10],col='red',lwd=1.5)
    lines(x=11:20,dam[tax,11:20],col='blue',lwd=1.5)
    axis(1, at=seq(1,20,2), labels=c(seq(1,10,2),seq(-10,-1,2)))
    # legend("top",legend=paste(rownames(dam),dam[, dim(dam)[2] ] ), cex = 0.6, bty = "n")
}

plot.readDisTable5 <- function(id,tax,folder){
    ## Plot valuble info about sample-species pair
    library(gridBase)
    library(gridExtra)
    ## read stuff from readDis file
    rd <- read.table(paste(folder,'readDist/',id,'_alignmentDist.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    rd <- rd[ tax , ]
    topNode <- rownames(rd)
    ## read mapDamage stuff
    dam <- read.table(paste(folder,'damageMismatch/',id,'_damageMismatch.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    idxN <- dim(dam)[2]
    dam <- dam[ grep(tax,rownames(dam)) , ]
    ## length distribution info
    ld <- read.table(paste(folder,'readDist/',id,'_readLengthStat.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    ld <- paste( round(ld[ tax , 'Mean' ],0),' (',round(ld[ tax , 'StandardDev' ],3),')',sep="")
    ## destacking on/off?
    ds <- read.table(paste(folder,'filterInformation/',id,'_filterTable.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    ds <- ds[tax, "turnedOn?"]
    ## join output and plot table
    data <- c(topNode,rd[topNode,'Reference'],rd[topNode,'TotalAlignmentsOnReference'],rd[topNode,'nonDuplicatesonReference'],rd[topNode,'uniquePerReference'],rd[topNode,'nonStacked'],ds,round(dam[topNode,'C>T_1'],4),round(dam[topNode,'G>A_20'],4),ld)
    data <- cbind( c('Node','Top Reference','all reads','nonDup','readDis','nonStacked','destacking?','C>T_1','G>A_-1','mean length (sd)'), data)
    colnames(data)=NULL; rownames(data)=NULL
    plot.new()
    mytheme <- gridExtra::ttheme_default(base_size=8) # https://github.com/baptiste/gridextra/wiki/tableGrob
    grid.table(data, vp=baseViewports()$figure,theme=mytheme)
}

table.additionalNodeEntries1 <- function(id,tax,folder){
    library(gridBase)
    library(gridExtra)
    ## read file and format
    ar <- read.table(paste(folder,'readDist/',id,'_additionalNodeEntries.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='',sep="\t",quote="")
    ar <- ar[ tax , ]
    ar <- paste(sub(";_TOPREFPERCREADS"," ",ar) ,"%",sep="")
    ## plot table
    plot.new()
    mytheme <- gridExtra::ttheme_default(base_size=8) # https://github.com/baptiste/gridextra/wiki/tableGrob
    grid.table(ar, vp=baseViewports()$figure,theme=mytheme)
}

## CODE

## INFO
## This scripts gathers signatures of species presence at nodes interrogated by MALTextract.
## Evidence is plotted in a heatmap for all samples and solely the samples with species specific evidence.
## For all sample-species pairs with evidence a profile-signature-pdf is plotted.
## Please see profilePDF_explained.pdf for details!
## This script is not designed for 'scan' output and come with no warranty.
## Questions/comments >> key@shh.mpg.de

## USAGE
## Input requires outpath of MALTextract -f <def_anc,default,ancient> run and the taxon-of-interest list (node.list) [e.g. MALTextract taxon list ] .
## Rscript ./postprocessing.v5.r -h

## get options, using the spec as defined by the enclosed list.
## we read the options from the default: commandArgs(TRUE).

spec = matrix(c(
    "rmaex.out.fld",  "r" , 1, "character", "MALTextract output folder.",
    "maltex.filter",  "m" , 2, "character", "MALTextract filter mode: <default,def_anc>. This script is not designed for 'scan' output. Default: <def_anc>.",
    "threads",  "t" , 1, "integer", "Max number of cores used.",
    "help"    ,  "h" , 0, "logical", "Print this help.",
    "sequencestrategy", "s", 1, "character", "pe or se for needing damage on either forward or reverse overhang for flagging, default se",
    "node.list"   ,  "n" , 1, "character","List (\\n separated) of nodes to be reported on (aka input species/node list used for MALTextract).",
    "heatmap.json"   ,  "j", 2, "logical", "Optional exporting of heatmap data in json format.",
    "reads" ,   "b" , 2, "double", "Minimum number of reads on node to output",
    "dmgcutoff" ,   "d" ,   2,  "double",  "Cutoff threshold for 3 prime damage for outputting plot. Default: 0, no cutoff is used",
    "readdistcutoff","c", 2,  "double",  "Cutoff threshold for read distribution (stacking) for outputting plot. Default: 0, no cutoff is used",
    "defratio"  ,   "e", 2, "double", "Minimum ratio of value sums of edit distances between successive bars of default edit distance needed for candidate taxon on a given sample to pass threshold, lower value is more permissive. Default: 0.9",
    "ancratio"  ,   "a",    2,  "double", "Minimum ratio of value sums of edit distances between successive bars of ancient edit distance needed for candidate taxon on a given sample to pass threshold, lower value is more permissive. Default: 0.8",
    "firm"  , "f", 0, "logical", "Use firm cutoffs, only output if all thresholds met for outputting plots, eg dmg, read dist, def ratio and anc ratio",
    "outputall", "o", 0, "logical", "Output all pdfs of sample x node any with reads above threshold (default threshold is 1 read)"
), byrow=TRUE, ncol=5);
opt = getopt(spec);


## and exit with a non-zero error code
if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}

### ARG parsing and sanity checks
## assign args and modify node.vec (tr ' ' '_')
path <- opt$rmaex.out.fld
if ( substr(path,nchar(path),nchar(path)) != "/"){path <- paste(path ,"/",sep="")} # add trailing "/" if missing
## parsing maltex filter and return update about what is being used
if ( is.null(opt$maltex.filter) ) {maltex.mode <- c('default','ancient'); print("No filter type provided, using default malt filter mode <def_anc>")
} else if (opt$maltex.filter == 'def_anc') {maltex.mode <- c('default','ancient')
} else if (opt$maltex.filter == 'default') {maltex.mode <- 'default'
} else {maltex.mode <- c('default','ancient'); print('Non standard malt filter mode provided, defaulting to <def_anc>')}
if ( !is.null(opt$reads) ) {reads_threshold <- opt$reads} else {reads_threshold <- 1}
if ( !is.null(opt$dmgcutoff) ) {dmgcutoff <- opt$dmgcutoff} else {dmgcutoff <- 0}
if ( !is.null(opt$readdistcutoff) ) {readdistcutoff <- opt$readdistcutoff} else {readdistcutoff <- 0}
if ( !is.null(opt$defratio) ) {defratio <- opt$defratio} else {defratio <- 0.9}
if ( !is.null(opt$ancratio) ) {ancratio <- opt$ancratio} else {ancratio <- 0.8}
if ( is.null(opt$sequencestrategy) ) {paired_end_mode <- FALSE
} else if (opt$sequencestrategy =='pe') {paired_end_mode <- TRUE
} else { paired_end_mode <- FALSE }
if ( !is.null(opt$firm) ) {firm <- TRUE} else {firm <- FALSE}
if ( !is.null(opt$outputall) ) {outputall <- TRUE} else {outputall <- FALSE}

## check if custom filtering values are acceptable
if (dmgcutoff < 0 || dmgcutoff > 1) {stop("damage cutoff value should be within range of [0,1]")}
if (readdistcutoff < 0 || readdistcutoff > 1) {stop("coverage cutoff value should be within range of [0,1]")}
if (defratio < 0 || defratio > 1) {stop("default ratio value should be within range of [0,1] (but likely below 0.9)")}
if (ancratio < 0 || ancratio > 1) {stop("ancient ratio value should be within range of [0,1] (but likely below 0.8)")}

# Print all filtering variables:
print('Printing all filtering parameter values:')
print('Damage cutoff:' )
print( dmgcutoff)
print('Read distribution cutoff:' )
print( readdistcutoff)
print('Default edit distance ration:' )
print( defratio )
print('Ancient edit distance ratio:' )
print( ancratio )
print('Reads threshold:' )
print( reads_threshold)


unq.spec <- unique(gsub(" ","_",scan(file=opt$node,sep="\n",what='character'))) # scan nodes, kill ' ', unique is solely sanity control
## START DATA PROCESSING
all.inds <- colnames(as.matrix(read.table(paste(path,'/default/RunSummary.txt',sep=''),sep="\t",header=T,stringsAsFactors=F,row.names=1,check.names=FALSE,comment.char='')))

### Extract MetaData for all Sample-Species Pairs
out.lists <- mclapply(1:length(all.inds), function(j) extract.stats5( all.inds[j],path,maltex.mode,paired_end_mode ), mc.cores=opt$threads )
data <- do.call(rbind, out.lists)
data <- data.frame(data,stringsAsFactors=F)
if ( length(maltex.mode) == 1 ) {
    # default mode, only grab data associated with all reads
    data[, c(4:10) ] = apply(data[ , c(4:10)], 2, function(x) as.numeric(as.character(x)))
} else {
    # def_anc mode, also grab data for ancient-damage reads
    data[, c(4:10,12:18) ] = apply(data[ , c(4:10,12:18)], 2, function(x) as.numeric(as.character(x)))
}

print('data frame parsing OK')
#############
## Extract scores and build matrix
#############
## scores are based on edit distance ratios, and damage
## only nodes with one of these + a read distribution above the cutoff will be output to pdf format/in the heatmap
if ( outputall ) {
    trg1 <- data[ data[,'reads_all'] >= reads_threshold &  !is.na(data[,'reads_all']) , ]
    print(rownames(trg1))
    res <- matrix(1L,nrow=length(unq.spec),ncol=length(all.inds),dimnames=list(a=unq.spec,b=all.inds))
    for (p in rownames(trg1)){
        print(p)
        res[ trg1[p,'spec']  , trg1[p,'id'] ] <- 2
    }
} else {
    if(length(maltex.mode) == 2){
        ## Default-Ancient
        ## if firm only output step 3 (all thresholds met), default only needs to satisfy default ratio + read distribution
        if ( firm ) {
            positions <- data[ data[,'reads_all'] >= reads_threshold & data[,'def.dr4'] >= defratio & !is.na(data[,'def.dr4']) & data[,'def.rd'] >= readdistcutoff & data[,'def.mapDam'] >= dmgcutoff & !is.na(data[,'def.mapDam']) & data[,'anc.dr4'] >= ancratio & !is.na(data[,'anc.dr4']), ]
            trg1 <- positions
            trg2 <- positions
            trg3 <- positions
        } else {
            trg1 <- data[ data[,'reads_all'] >= reads_threshold & data[,'def.dr4'] >= defratio & !is.na(data[,'def.dr4']) & data[,'def.rd'] >= readdistcutoff, ] ## Step1: DiffRatio0-4: > defratio (default = 0.9) and read distribution > cutoff (default = 0)
            trg2 <- data[ data[,'def.mapDam'] >= dmgcutoff & !is.na(data[,'def.mapDam']) & data[,'def.rd'] >= readdistcutoff, ] ## Step2: Terminal Damage Present (default = 0)
            trg3 <- data[ data[,'anc.dr4'] >= ancratio & !is.na(data[,'anc.dr4']) & data[,'def.rd'] >= readdistcutoff, ] ## Step3: DiffRatio1-4: > ancratio (default = 0.8)
        }


        # Build Matrix for Heatmap
        res <- matrix(1L,nrow=length(unq.spec),ncol=length(all.inds),dimnames=list(a=unq.spec,b=all.inds))
        for (p in rownames(trg1)){
            if( !p %in% rownames(trg2) ){
                res[ trg1[p,'spec']  , trg1[p,'id'] ] <- 2
            } else if( !p %in% rownames(trg3) ){
                res[ trg2[p,'spec']  , trg2[p,'id'] ] <- 3
            } else {
                res[ trg3[p,'spec']  , trg3[p,'id'] ] <- 4
            }
        }
    } else {
        ## Default: Extract scores and build matrix
        trg1 <- data[ data[,'reads_all'] >= reads_threshold & data[,'def.dr4'] >= defratio & !is.na(data[,'def.dr4']) & data[,'def.rd'] >= readdistcutoff , ] ## Step1: DiffRatio0-4: > defratio (default = 0.9)
        trg2 <- data[ data[,'def.mapDam'] >= dmgcutoff & !is.na(data[,'def.mapDam']) & data[,'def.rd'] >= readdistcutoff , ] ## Step2: Terminal Damage Present (default = 0)

        # Build Matrix for Heatmap
        res <- matrix(1L,nrow=length(unq.spec),ncol=length(all.inds),dimnames=list(a=unq.spec,b=all.inds))
        for (p in rownames(trg1)){
            if( !p %in% rownames(trg2) ){
                res[ trg1[p,'spec']  , trg1[p,'id'] ] <- 2
            } else {
                res[ trg2[p,'spec']  , trg2[p,'id'] ] <- 3
            }
        }
    }
}

print('Identifying Candides OK')

##############
## Plot heatmap
##############
## plot reduced overview heatmap only for samples and tax with evidence
## Def-Anc shades of red and Default only shades of green
if(length(maltex.mode) == 2){
    mycol=c('lightgray','yellow','orange','red')
    leg.txt <- c('Edit distance','+Damage','+Dam. Edit Dist.')
} else {
    mycol=c('lightgray','lightgreen','darkgreen')
    leg.txt <- c('Edit distance','+Damage')
}
red.res <- res[, colSums(res) > dim(res)[1] , drop = FALSE ] ## drops columns with only 1 (of 1,2,3,4) --> lightgray on heatmap (so if all species for a given sample are lightgray then don't put them in heatmap)
red.res <- red.res[ rowSums(red.res) > dim(red.res)[2] , , drop = FALSE ]

if ( dim(red.res)[1] == 0 && dim(red.res)[2] == 0 ) {
    print("No samples identified with suffienct signal on candidate taxon!")
    print("postprocessing.AMPS.r script complete.")
} else {
    ## NT & organelle results hat <> in species name (equus), caused bug
    rownames(red.res)=chartr("><","..",rownames(red.res))

    pdf.height <- max(dim(red.res)[1]/2.5,20)
    pdf.width <- max(dim(red.res)[2]/10 , 20)
    pdf(paste(path,'heatmap_overview_Wevid.pdf',sep=""),height=pdf.height,width=pdf.width)
    par(mar=c(5.1,30.1,25.1,2.1))
    if(ncol(red.res)!=0){image(x=1:ncol(red.res),y=1:nrow(red.res),z=t(red.res),col=mycol,axes=F,ylab="",xlab="",zlim=c(1,4))
        axis(side=2,at=1:nrow(red.res),labels=rownames(red.res),las=1,cex.axis=2)
        axis(side=3,at=1:ncol(red.res),labels=colnames(red.res),las=2,cex.axis=2,tick=F)
        abline(h=1:length(rownames(red.res))+0.5,col='darkgrey') # add horizontal lines for improved vision
        abline(v=1:length(colnames(red.res))+0.5,col='darkgrey') # add vertical lines for improved vision
        xleg <- ncol(red.res)-(ncol(red.res)*1.35)
        yleg <- nrow(red.res)+5
        legend(x=xleg,y=yleg, legend=leg.txt, fill = mycol[-1],xpd=T,cex=3)
        dev.off()
    }
    ## Export table format
    red.res.tab <- cbind(rownames(red.res), data.frame(red.res, row.names=NULL))
    colnames(red.res.tab)[1] <- "node"
    write.table(red.res.tab, file = paste(path,'heatmap_overview_Wevid.tsv',sep = ""), sep = "\t", row.names = F)

    if (!is.null(opt$heatmap.json)) {
        library("jsonlite")

        ## Prepare for tab to list
        red.res.tab.json <- red.res.tab
        rownames(red.res.tab.json) <- red.res.tab.json$node
        red.res.tab.json <- subset(red.res.tab.json, select = -node)

        ## convert to list
        red.res.json <- Map(function(x) {
            value_list <- as.list(x)
            names(value_list) <- row.names(red.res.tab.json)
            value_list
        }, as.list(red.res.tab.json))

        ## convert and save list as json
        write_json(red.res.json, path = paste(path,'heatmap_overview_Wevid.json',sep = ""), pretty = T)
    }

    ## Export postprocessing parameters text file
    if (length(maltex.mode)==2) {mltexmd <- 'def_anc'} else {mltexmd <- 'default'}

    output_parameters_table <- data.frame(
        variable=c('malt_extract_mode','malt_extract_output_folder','paired_end_mode','node_list','damage_cutoff','read_distribution_cutoff','default_edit_distance_ratio','ancient_edit_distance_ratio','firm_cutoff_for_output'),
        values=c(mltexmd,opt$rmaex.out.fld,toString(paired_end_mode),opt$node,dmgcutoff,readdistcutoff,defratio,ancratio,firm)
    )
    write.table(output_parameters_table, file = paste(path,"post_processing_parameters.txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE,quote=FALSE)

    ########################
    ###### Candidate Profile PDFs
    ########################
    ## plot summary pdf's for candidates if outpath specified
    folder.names <- paste(path,maltex.mode,'/',sep="")

    for (spl in colnames(red.res)){
        for (tax in rownames(red.res)){
            if( red.res[tax,spl] > 1 ){
                system(paste('mkdir -p ',path,'pdf_candidate_profiles/',tax,sep='')) #mk pdf output folder
                pdf(paste(path,'pdf_candidate_profiles/',tax,'/stp',red.res[tax,spl]-1,'_',spl,'_',tax,'_summary.pdf',sep=''))
                par(mfrow=c(3,2))
                plot.editDis(spl,tax,folder.names) # plot default (and ancient) edit distance
                plot.mapDamage3(spl,tax,folder.names[1]) # only calculated for default ("ancient" mode bias damage)
                plot.readDisTable5(spl,tax,folder.names[1]) # table w/ detailled info on best ref
                table.additionalNodeEntries1(spl,tax,folder.names[1])
                dev.off()
            }
        }
    }
    print("Samples identified with suffienct signal on candidate taxon!")
    print("postprocessing.AMPS.r script complete.")
    
    #########
    ## Save RData
    #########
    save.image(paste(path,'analysis.RData',sep=''))
}