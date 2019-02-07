#Assignment 5 script. 
#Kai Ellis, 20019803
#Github link: https://github.com/kellisfm/Assignment-5-Protein-alignments

#library for this script:
library(BiocManager)
library(sangerseqR)
library(genbankr)
library(annotate)
library(ape)
#this script will:
#take the following sequences:
ProteinSeq="NPNQKITTIGSICMVIGIVSLMLQIGNIISIWVSHSIQTGNQHQ
AEPCNQSIITYENNTWVNQTYVNISNTNFLTEKAVASVTLAGNSSLCPISGWAVYSKD
NGIRIGSKGDVFVIREPFISCSHLECRTFFLTQGALLNDKHSNGTVKDRSPHRTLMSC
PVGEAPSPYNSRFESVAWSASACHDGTSWLTIGISGPDNGAVAVLKYNGIITDTIKSW
RNNIMRTQESECACVNGSCFTVMTDGPSNGQASYKIFRIEKGKVVKSAELNAPNYHYE
ECSCYPDAGEITCVCRDNWHGSNRPWVSFNQNLEYRIGYICSGVFGDNPRPNDGTGSC
GPVSPKGAYGIKGFSFKYGNGVWIGRTKSTNSRSGFEMIWDPNGWTGTDSNFSVKQDI
VAITDWSGYSGSFVQHPELTGLDCIRPCFWVELIRGRPKESTIWTSGSSISFCGVNSD
TVGWSWPDGAELPFTIDK"

DNASeq="AGCAAAAGCAGGAGTTCAAAATGAATCCAAATCAGAAGATAACAACCATTGGATCAATCTGTATGGTAAT
TGGAATAGTTAGCTTGATGTTACAAATTGGGAACATAATCTCAATATGGGTTAGTCATTCAATTCAAACA
GGGAATCAACACCAGGCTGAACCATGCAATCAAAGCATTATTACTTATGAAAACAACACCTGGGTAAACC
AGACATATGTCAACATCAGCAATACCAATTTTCTTACTGAGAAAGCTGTGGCTTCAGTAACATTAGCGGG
CAATTCATCTCTTTGCCCCATTAGTGGATGGGCTGTATACAGTAAGGACAACGGTATAAGAATCGGTTCC
AAGGGGGATGTGTTTGTTATAAGAGAGCCGTTCATCTCATGCTCCCACTTGGAATGCAGAACTTTCTTTT
TGACTCAGGGAGCCTTGCTGAATGACAAGCATTCTAATGGGACCGTCAAAGACAGAAGCCCTCACAGAAC
ATTAATGAGTTGTCCCGTGGGTGAGGCTCCTTCCCCATACAACTCGAGGTTTGAGTCTGTTGCTTGGTCG
GCAAGTGCTTGTCATGATGGCACTAGTTGGTTGACAATTGGAATTTCTGGCCCAGACAATGGGGCTGTGG
CTGTATTGAAATACAATGGCATAATAACAGACACTATCAAGAGTTGGAGGAACAACATAATGAGAACTCA
AGAGTCTGAATGTGCATGTGTAAATGGCTCTTGCTTTACTGTTATGACTGATGGACCAAGTAATGGGCAG
GCTTCATACAAAATCTTCAGAATAGAAAAAGGGAAAGTAGTTAAATCAGCCGAATTAAATGCCCCTAATT
ATCACTATGAGGAGTGCTCCTGTTATCCTGATGCTGGAGAAATCACATGTGTGTGCAGGGATAACTGGCA
TGGCTCAAATCGGCCATGGGTATCTTTCAATCAAAATTTGGAGTATCGAATAGGATATATATGCAGTGGA
GTTTTCGGAGACAATCCACGCCCCAATGATGGGACAGGCAGTTGTGGTCCGGTGTCCCCTAAAGGGGCAT
ATGGAATAAAAGGGTTCTCATTTAAATACGGCAATGGTGTTTGGATCGGGAGAACCAAAAGCACTAATTC
CAGGAGCGGCTTTGAAATGATTTGGGATCCAAATGGATGGACTGGTACGGACAGTAATTTTTCAGTAAAG
CAAGATATTGTAGCTATAACCGATTGGTCAGGATATAGCGGGAGTTTTGTCCAGCATCCAGAACTGACAG
GATTAGATTGCATAAGACCTTGTTTCTGGGTTGAGCTAATCAGAGGGCGGCCCAAAGAGAGCACAATTTG
GACTAGTGGGAGCAGCATATCCTTTTGTGGTGTAAATAGTGACACTGTGGGTTGGTCTTGGCCAGACGGT
GCTGAGTTGCCATTCACCATTGACAAGTAGTTTGTTCAAAAAACTCCTTGTTTCTACT"

#and analyze them with alignments and phylogentics
#to determine if there is anything to be concerned about with this seq

#use parameter database="nr",program="blastp" to find protein blast with blastSequences()

#we're gonna start by blasting, and seeing what happens
Pblast=blastSequences(ProteinSeq,database = "nr",program = "blastp", as ="data.frame", 
                      hitListSize = 10, timeout = 100 )
#Running a genbank check on the first three protein hits to see what comes up
PhitSeqs=read.GenBank(Pblast$Hit_accession[1:3])
attr(PhitSeqs,"species") #looks like we're dealing with influenza A (H5N1), but lets do some more analysis with the dna seq

#dna blast functions similarly to protein blast, but we can more easily create phylotrees from it
DNAblast=blastSequences(DNASeq,as ="data.frame", hitListSize = 50, timeout = 100 )
?blastSequences

#to make a phylo tree we're first gonna need to align the sequences with muscle, which requires a conversion to DNAbin format
DNAHitsDF=data.frame(ID=DNAblast$Hit_accession,seq=DNAblast$Hsp_hseq,stringsAsFactors=F)

#Since we only have around 40 hits, it doesnt take tooo long to run get the actual names of the 
#sequences we've found. Additionally setting these as the ids for our dna bin makes reading the data much easier.
DNAhitSeqs=read.GenBank(DNAblast$Hit_accession)
ID=attr(DNAhitSeqs,"species")
DNAHits=sapply(DNAHitsDF$seq,strsplit,split="")
#names(DNAHits)=paste(1:nrow(DNAHitsDF),DNAHitsDF$ID,sep="_")
names(DNAHits)=paste(ID)

#convert to DNAbin to allow muscle use
DNAHits=as.DNAbin(DNAHits)
DNAmuscle=muscle(DNAHits,quiet=F)

#Now we check the alignment to ensure theres no large gaps or other errors
checkAlignment(DNAmuscle,what=1)
checkAlignment(DNAmuscle,what=2)

#it seems the large majority of gaps are quite small, and are confined mostly to the start or end of the graph
#however I'm gonna see if i can eliminate the larger gaps anyway just to assist the algorithm
checkAlignment(DNAmuscle,what=3)
checkAlignment(DNAmuscle,what=4)

#this plots out the hit lengths so we can see if there are any outliers
SeqLen=as.numeric(lapply(DNAHits,length))
library(ggplot2)
qplot(SeqLen)+theme_bw()

#There are two large clusters of sequences, one at roughly 1430 and the other at roughly 1460
#We can cut out everything below these clusters to simplfy the alignment a bit
keepSeq=SeqLen>1430
DNASubset=DNAmuscle[keepSeq,]

checkAlignment(DNASubset,what=1)
checkAlignment(DNASubset,what=2)
#Now the largest gap is only 13, down from 42, and the edges are a bit more consistant.
#time to realign based on this subset to see if the algorythem makes any changes

DNAsubAlign=muscle(DNASubset,quiet=F)
checkAlignment(DNAsubAlign,what=1)
checkAlignment(DNAsubAlign,what=2)
#nothing seems to have changed, but it was worth a shot!
#Now that the alignment has been simplified, we can make a matrix
DNADM=dist.dna(DNAsubAlign,model="K80")
DNADMmat=as.matrix(DNADM)

#convert it from wide form to long form with the melt package
library(reshape2)
Ddat=melt(DNADMmat)

#make a heat plot with ggplot
ggplot(data=Ddat,aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_gradientn(colours=c("white","blue","green","red"))

#and a DNA tree with ggtree
library(ggtree)
DNAtree=nj(DNADM)
ggtree(DNAtree) +geom_tiplab()

