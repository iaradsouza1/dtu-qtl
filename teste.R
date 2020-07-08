library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(biomaRt)
library(Rsamtools)
library(VariantAnnotation)
library(tools)

# Annotation and metadata
gtf_file <- "data/gtf/teste.gtf"
gtf <- import(gtf_file, format = "gtf")

# Keep only protein coding genes
chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
keep <- mcols(gtf)$gene_biotype == "protein_coding" & mcols(gtf)$type == "gene" & seqnames(gtf) %in% chroms
gtf <- gtf[keep]
genes_bed <- data.frame(chr = seqnames(gtf), start = start(gtf),
                        end = end(gtf), geneId = mcols(gtf)$gene_id,
                        stringsAsFactors = FALSE)
# Save files
# for(i in as.character(1:22)){
#   genes_bed_sub <- genes_bed[genes_bed$chr == i, ]
#   write.table(genes_bed_sub, paste0("data/gtf/genes_chr", i ,".bed"), quote = FALSE,
#               sep = "\t", row.names = FALSE, col.names = FALSE)
# }

# Read metadata
ann <- read.table("data/meta/SraRunTable.txt", sep = ",", header = T)
ann <- ann %>% 
  filter(Organism == "Homo sapiens") %>% 
  mutate(region = case_when(
    tissue == "Orbitofrontal (OFC; BA11)" ~ "OFC",
    tissue == "Dorsolateral prefrontal cortex (dlPFC; BA8/9)" ~ "dlPFC",
    tissue == "Cingulate gyrus 25 (Cg25)" ~ "Cg25",
    tissue == "Anterior Insula (aINS)" ~ "aINS",
    tissue == "Nucleus Accumbens (Nac)" ~ "Nac",
    tissue == "Subiculum (Sub)" ~ "Sub"
  ),
  group = paste(region, phenotype, gender, sep = "_"),
  sample_id = Run) %>% 
  dplyr::select(sample_id, group)
  
# Transcript counts
load("data/expr/txi_tx.rda")
cts <- txi$counts
rownames(cts) <- gsub("\\.\\d+", "", rownames(cts))
cts <- as.data.frame(cts)
cts$feature_id <- rownames(cts)
rownames(cts) <- NULL
dict <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"), 
              mart = useMart("ensembl", "hsapiens_gene_ensembl"))
cts <- inner_join(cts, dict, by = c("feature_id" = "ensembl_transcript_id"))

# Genotype

#files <- paste0("data/vcf", ann$sample_id, ".ann.vcf")
files <- sapply(strsplit(list.files("data/vcf"), split = "\\."), "[[", 1)
files <- paste0("data/vcf/", files, ".ann.vcf")
for (i in 1:length(files)) {
  if(!file.exists(paste0(files[i], ".bgz.tbi"))) {
    zipped <- bgzip(files[i])
    idx <- indexTabix(zipped, format = "vcf")
  }
}

# Create tabix indexes
WINDOW <- 5000
gene_ranges <- resize(gtf, GenomicRanges::width(gtf) + 2 * WINDOW, fix = "center")
zipped <- paste0(files[1], ".bgz")
idx <- paste0(files[1], ".bgz.tbi")
tab <- TabixFile(zipped, idx)
### Explore the file header with scanVcfHeader
hdr <- scanVcfHeader(tab)
param <- GRanges(seqnames = Rle(values = c(paste0("chr", 1:22), "chrX", "chrY"), 
                                runLength(seqnames(gtf))), 
                 ranges = gtf@ranges)
vcf <- readVcf(file = tab, "hg38", param = param)

### Keep only the bi-allelic SNPs
# width of ref seq
rw <- width(ref(vcf))
# width of first alt seq
#aw <- sapply(alt(vcf), function(x) {width(x)})
# number of alternate genotypes
nalt <- elementNROWS(alt(vcf))
# select only bi-allelic SNPs (monomorphic OK, so aw can be 0 or 1) 
snp <- rw == 1 & nalt == 1
# subset vcf
vcfbi <- vcf[snp,]
rowdata <- rowData(vcfbi)

### Convert genotype into number of alleles different from reference
geno <- geno(vcfbi)$GT
geno01 <- geno
geno01[,] <- -1
geno01[geno %in% c("0/0", "0|0")] <- 0 # REF/REF
geno01[geno %in% c("0/1", "0|1", "1/0", "1|0")] <- 1 # REF/ALT
geno01[geno %in% c("1/1", "1|1")] <- 2 # ALT/ALT
mode(geno01) <- "integer"

genotypes <- unique(data.frame(chr = as.vector(seqnames(vcfbi), "character"),
                               start = start(vcfbi), end = end(vcfbi), snpId = rownames(geno01),
                               geno01, stringsAsFactors = FALSE))

geuv_genotypes <- GeuvadisTranscriptExpr::genotypes




