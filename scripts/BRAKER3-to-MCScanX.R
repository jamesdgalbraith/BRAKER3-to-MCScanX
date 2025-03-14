suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(plyranges))
suppressMessages(library(optparse))

# Script to reorganise BRAKER3 output for use with MCScanX
# Selects the longest transcript of each gene and adds prefix

option_list <- list(
  make_option(c("-a", "--in_aa"), default=NA, type = "character",
              help="Input braker.aa file"),
  make_option(c("-c", "--in_cds"), default=NA, type = "character",
              help="Input braker.cds file"),
  make_option(c("-g", "--in_gff"), default=NA, type = "character",
              help="Input braker.gff file"),
  make_option(c("-p", "--prefix"), default=NA, type = "character",
              help="Prefix to add to gene names"),
  make_option(c("-o", "--outdir"), default=NA, type = "character",
              help="Output directory")
)


# Parse and chek options
opt <- parse_args(OptionParser(option_list=option_list))
if(is.na(opt$outdir)){
  stop("Out directory must be specified")
}
if(!dir.exists(opt$outdir)){
  dir.create(opt$outdir, recursive = T)
}
if(is.na(opt$in_aa) | is.na(opt$in_cds) | is.na(opt$in_gff) | is.na(opt$prefix)){
  stop("One of required variables not specified")
}

if(!file.exists(opt$in_aa) | !file.exists(opt$in_cds) | !file.exists(opt$in_gff)){
  stop("One of required files does not exist")
}
if(nchar(opt$prefix) != 2){
  stop("Prefix must be two letters")
}

# Sort out output paths
out_aa_path <- paste0(opt$outdir, '/MCScanX_', tail(unlist(strsplit(opt$in_aa, "/")), 1))
out_gff_path <- paste0(opt$outdir, '/MCScanX_', sub("gff3", "gff", tail(unlist(strsplit(opt$in_gff, "/")), 1)))
out_gff3_path <- paste0(opt$outdir, '/MCScanX_', tail(unlist(strsplit(opt$in_gff, "/")), 1))
out_cds_path <- paste0(opt$outdir, '/MCScanX_', sub(".codingseq", ".cds", tail(unlist(strsplit(opt$in_cds, "/")), 1)))

# Read in coding seq
aa_in <- Biostrings::readAAStringSet(opt$in_aa)
# Find longest transcript
longest_transcripts <- tibble(seqnames = names(aa_in), width = width(aa_in)) %>%
  mutate(sep_names = seqnames) %>%
  tidyr::separate(sep_names, into = c("gene", "transcript"), sep = "\\.") %>%
  mutate(gene = as.numeric(sub("g", "", gene))) %>%
  dplyr::group_by(gene) %>%
  dplyr::arrange(gene, -width) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()
# Get longest transcript from gff
longest_transcripts_coords <- plyranges::read_gff(opt$in_gff) %>%
  filter(ID %in% longest_transcripts$seqnames) %>%
  as_tibble() %>%
  mutate(seqnames = as.character(seqnames),
         gene = sub("g", "G", sub("\\.t.*", "", ID))) %>%
  dplyr::select(seqnames, gene, start, end, ID, strand, source)
chr <- tibble(seqnames = base::unique(longest_transcripts_coords$seqnames)) %>%
  mutate(chr_no = row_number())
longest_transcripts_coords <- inner_join(longest_transcripts_coords, chr, by = "seqnames") %>%
  mutate(gene = paste0(opt$prefix, chr_no, gene))
longest_transcripts_names <- longest_transcripts_coords %>%
  dplyr::select(gene, ID) %>%
  dplyr::rename(seqnames = ID)
longest_transcripts_coords_gff <- longest_transcripts_coords[1:4]
longest_transcripts_coords_gff3 <- longest_transcripts_coords %>%
  mutate(X6 = ".", X8 = ".", type = "gene",
         X9 = paste0("ID=", gene, ";Name=", gene)) %>%
  dplyr::select(seqnames, source, type, start, end, X6, strand, X8, X9)

# Get longest transcript seq, add custom gene names
longest_transcripts_seq <- aa_in[names(aa_in) %in% longest_transcripts$seqnames]
longest_transcripts_seq_tibble <- tibble(seqnames = names(longest_transcripts_seq), seq = as.character(longest_transcripts_seq, use.names=TRUE)) %>%
  inner_join(longest_transcripts_names, by = "seqnames")
longest_transcripts_seq <- AAStringSet(longest_transcripts_seq_tibble$seq)
names(longest_transcripts_seq) <- longest_transcripts_seq_tibble$gene

# Repeat for cds
nt_in <- Biostrings::readAAStringSet(opt$in_cds)
longest_transcripts_nt_seq <- nt_in[names(nt_in) %in% longest_transcripts$seqnames]
longest_transcripts_nt_seq_tibble <- tibble(seqnames = names(longest_transcripts_nt_seq),
                                            seq = as.character(longest_transcripts_nt_seq, use.names=TRUE)) %>%
  inner_join(longest_transcripts_names, by = "seqnames")
longest_transcripts_nt_seq <- DNAStringSet(longest_transcripts_nt_seq_tibble$seq)
names(longest_transcripts_nt_seq) <- longest_transcripts_nt_seq_tibble$gene

# Output data in MCScanX require format
writeXStringSet(longest_transcripts_nt_seq, out_cds_path)
writeXStringSet(longest_transcripts_seq, out_aa_path)
write_tsv(longest_transcripts_coords_gff, out_gff_path, col_names = F)
write_tsv(tibble("##gff-version 3"), out_gff3_path, col_names = F)
write_tsv(longest_transcripts_coords_gff3, out_gff3_path, col_names = F, append = T)
