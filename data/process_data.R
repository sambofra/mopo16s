library(Biostrings)
library(xlsx)

### Parameters ###

# Choose cluster similarity level of the GreenGenes OTUs used as representative set
rep.set.level <- "85"

# For the initial primer pairs, choose:
# minimum and maximum desired amplicon length
max.alen <- 800
min.alen <- 700

# minimum and maximum primer length
max.plen <- 21
min.plen <- 17

# maximum total degeneracy
max.degen <- 20


### Data processing ###

# load the representative set and its updated taxonomy
rep.set <- readDNAStringSet(paste("gg_13_5_otus/rep_set/",rep.set.level,"_otus.fasta",sep=""))
rep.set.tax <- read.table(paste("gg_13_5_otus/taxonomy/",rep.set.level,"_upd_tax.txt",sep=""),
                          stringsAsFactors=F,header=F,sep="\t")
rep.set.tax <- rep.set.tax[match(names(rep.set),rep.set.tax[,1]),]

# select only bacterial sequences
rep.set.bact <- rep.set[grep("Bacteria",rep.set.tax[,2])]

# save file
writeXStringSet(rep.set.bact, "rep_set.fasta",format = "fasta")


# load initial primers and primer pairs
t.primers <- read.xlsx2("supplementary_table_01.xlsx",sheetIndex = 1, startRow = 2,
                        endRow=177, colIndex = c(2,5,8), stringsAsFactors=F)
init.primers <- t.primers[,2]
names(init.primers) <- t.primers[,1]

# function to remove non-letter characters
trim <- function (x) gsub("[^A-Z]", "",x)
init.primers <- trim(init.primers)

t.pairs <- read.xlsx2("supplementary_table_08.xlsx",sheetIndex = 1, startRow = 2,
                        endRow=514, colIndex = c(1,2,4), stringsAsFactors=F)
names(t.pairs) <- c("fwd","rev","alen")

# retrieve domain and primer length for forward and reverse primers
split.fwd <- simplify2array(strsplit(t.pairs$fwd,"-"))
split.rev <- simplify2array(strsplit(t.pairs$rev,"-"))

domain.fwd <- split.fwd[3,]
domain.rev <- split.rev[3,]
plen.fwd <- as.numeric(split.fwd[7,])
plen.rev <- as.numeric(split.rev[7,])
 
# select good pairs:
# - reference amplicon length between min.alen and max.alen
# - target domain != Arch
# - with primer length between min.plen and max.plen
# - with both primers in init.primers
# - with total degeneracy <= max.degen
good.pairs <- t.pairs[ t.pairs$alen <= max.alen & t.pairs$alen >= min.alen & 
                         domain.fwd != "Arch" & domain.rev != "Arch" &
                         plen.fwd >= min.plen & plen.rev >= min.plen &
                         plen.fwd <= max.plen & plen.rev <= max.plen, ]

# check that both primers are in init.primer
good.pairs <- good.pairs[!is.na(init.primers[good.pairs$fwd]) & 
                           !is.na(init.primers[good.pairs$rev]),]

# check maximum degeneracy
fwd.seq <- DNAStringSet( init.primers[good.pairs$fwd] )
rev.seq <- DNAStringSet( init.primers[good.pairs$rev] )

# function to compute non degen primer set from degen primer
degen.to.set <- function( degen.primer )
{
  DNAStringSet(
    apply( expand.grid ( 
      lapply(
        as.list(IUPAC_CODE_MAP[strsplit(toString( degen.primer ),"")[[1]]]),
        function(el)(strsplit(el,"")[[1]])), 
      stringsAsFactors=F
    ), 1, paste, collapse=""
    )
  )
}

good.pairs <- good.pairs[
  sapply(fwd.seq,function(el)(length(degen.to.set(el)))) + 
    sapply(rev.seq,function(el)(length(degen.to.set(el)))) <= 20,]

# save data in fasta format, alternating forward and reverse primers
writeXStringSet(
  DNAStringSet(init.primers[c(rbind(good.pairs$fwd,good.pairs$rev))]),
  "good_pairs.fasta",format="fasta")
