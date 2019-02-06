#Assignment 5 script. 
#Kai Ellis, 20019803
#Github link: https://github.com/kellisfm/Assignment-5-Protein-alignments

#library for this script:
library(BiocManager)
library(sangerseqR)
library(reshape2)
library(genbankr)
library(annotate)
library(ape)
library(ggplot2)
library(ggtree)
#this script will:
#take the following sequence:
ProteinSeq="NPNQKITTIGSICMVIGIVSLMLQIGNIISIWVSHSIQTGNQHQ
AEPCNQSIITYENNTWVNQTYVNISNTNFLTEKAVASVTLAGNSSLCPISGWAVYSKD
NGIRIGSKGDVFVIREPFISCSHLECRTFFLTQGALLNDKHSNGTVKDRSPHRTLMSC
PVGEAPSPYNSRFESVAWSASACHDGTSWLTIGISGPDNGAVAVLKYNGIITDTIKSW
RNNIMRTQESECACVNGSCFTVMTDGPSNGQASYKIFRIEKGKVVKSAELNAPNYHYE
ECSCYPDAGEITCVCRDNWHGSNRPWVSFNQNLEYRIGYICSGVFGDNPRPNDGTGSC
GPVSPKGAYGIKGFSFKYGNGVWIGRTKSTNSRSGFEMIWDPNGWTGTDSNFSVKQDI
VAITDWSGYSGSFVQHPELTGLDCIRPCFWVELIRGRPKESTIWTSGSSISFCGVNSD
TVGWSWPDGAELPFTIDK"
#and analyze it with alignments and phylogentics
#to determine if there is anything to be concerned about with this seq

#use parameter database="nr",program="blastp" to find protein blast with blastSequences()


