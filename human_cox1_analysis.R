# Project: Human Mitochondrial Genome Analysis
# Script: Analysis of the Homo sapiens COX1 Gene

# 1. SETUP & DATA LOADING
# Load the essential bioinformatics library
setwd() setwd("C:\\Users\\click\\OneDrive\\Pictures\\Desktop\\MiniProjectR")
list.files()

library(Biostrings)

# Load your Human COX1 FASTA file
# Note: Ensure "COX1_Human.fasta" is in your working folder
dna <- readDNAStringSet("COX1_Human.fasta")
dna_seq <- dna[[1]]

# Display Basic Info
print("Sequence Overview")
print(dna)
print(paste("Length:", width(dna_seq), "bp"))

# 2. GC CONTENT ANALYSIS
gc <- letterFrequency(dna_seq, letters = c("G","C"), as.prob = TRUE)
gc_percent <- sum(gc) * 100
print(paste("GC Content:", round(gc_percent, 2), "%"))

# Export GC Plot
png("human_gc_content.png", width=1200, height=900)
barplot(letterFrequency(dna_seq, c("G","C")), 
        main="Human COX1 GC Content", 
        col="coral", xlab="Base", ylab="Count")
dev.off()

# 3. NUCLEOTIDE FREQUENCY
freq <- letterFrequency(dna_seq, letters = c("A","T","G","C"))
png("human_nucleotide_frequency.png", width=1200, height=900)
barplot(freq, 
        main="Human COX1 Nucleotide Frequency", 
        col="skyblue", xlab="Nucleotide", ylab="Count")
dev.off()

# 4. MOTIF SEARCH (Regulatory Elements)
# Searching for Start codons, TATA promoters, PolyA signals, and CpG sites
motifs <- c("ATG", "TATA", "AATAAA", "CG")
motif_results <- lapply(motifs, function(m) matchPattern(m, dna_seq))
names(motif_results) <- motifs

print("--- Motif Counts ---")
print(sapply(motif_results, length))

# 5. ORF DETECTION & TRANSLATION
# Find potential Open Reading Frames (ORFs)
starts <- start(matchPattern("ATG", dna_seq))
stops <- sort(c(start(matchPattern("TAA", dna_seq)), 
                start(matchPattern("TAG", dna_seq)), 
                start(matchPattern("TGA", dna_seq))))

orf_list <- list()
for (st in starts) {
  possible_stops <- stops[stops > st]
  if (length(possible_stops) > 0) {
    end_pos <- min(possible_stops)
    orf_list[[length(orf_list) + 1]] <- dna_seq[st:(end_pos + 2)]
  }
}

# Translate the longest ORF into Protein
if (length(orf_list) > 0) {
  orf_lengths <- sapply(orf_list, length)
  longest_orf <- orf_list[[which.max(orf_lengths)]]
  protein <- translate(longest_orf)
  
  print("--- Protein Translation (Longest ORF) ---")
  print(protein)
  
  # Export Codon Usage Plot for the protein-coding region
  codon_freq <- oligonucleotideFrequency(longest_orf, width = 3)
  nonzero <- codon_freq[codon_freq > 0]
  png("human_codon_usage.png", width=1500, height=900)
  barplot(nonzero, las=2, main="Human Codon Usage (Non-zero)", 
          col="darkseagreen", ylab="Frequency")
  dev.off()
}

# 6. SEQUENCE MANIPULATION
# Reverse Complement of the sequence
rev_dna <- reverseComplement(dna_seq)

print("Analysis Complete. All plots have been saved as .png files.")