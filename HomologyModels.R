# HomologModels.R
# Date: 2019-11-27
# Author: Xiaowen Zhang
# Purpose: Produce two models (1BM8 and 4UX5) for MBP1_MYSPE, and use bio3d to superimpose the two models, calculate the
# RMSD between each residue pair, writes the RMSD to each residues' B Factor field for both of the two models, and save
# the resulting PDB files.
#
# ToDo:
# Note:
#
# Reference: The codes in the step 1 was retrieved from "BIN-SX-Homology_modelling" learning unit (Steipe, B. 2019)
# with modifications.
#
# ==============================================================================

# ====  PACKAGES  ==============================================================
# Check that required packages have been installed. Install if needed.
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (! requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

if (! requireNamespace("msa", quietly=TRUE)) {
  BiocManager::install("msa")
}

if (! requireNamespace("bio3d", quietly = TRUE)) {
  install.packages("bio3d")
}

# ==== Step 0: Initialize & Model1BM8 ==============================================================

source("makeProteinDB.R") # get myDB loaded

# From the previous learning unit, I've generated a model for MBP1_MYSPE based on 1BM8.
Model1BM8 <- bio3d::read.pdb("MBP1_MYSPE-APSESrenum.pdb") # read the pdb file of Model1BM8

# ====  Step 1:  Model4UX5  =============================================================
# Processes: Produce a second model for MBP1_MYSPE based on 4UX5.
# Objective: save a .pdb file containing the Model1UX5.

# ====  Step 1.1: Target Sequence =============================================================
# Process: Define the TARGET sequence by sub-selecting the "APSES fold" in MBP1_THETH sequence.

targetName <- sprintf("MBP1_%s", biCode(MYSPE)) #concatenate the targetName
# Get the protein IDs.
sel <- which(myDB$protein$name == targetName)
proID <- myDB$protein$ID[sel]
# Find the feature ID in the feature table
ftrID <- myDB$feature$ID[myDB$feature$name == "APSES fold"]
# Get the annotation ID.
fanID <- myDB$annotation$ID[myDB$annotation$proteinID == proID &
                               myDB$annotation$featureID == ftrID]

# Get the feature start and end:
start <- myDB$annotation$start[fanID]
end   <- myDB$annotation$end[fanID]

# Extract the feature from the sequence
targetSeq <- substring(myDB$protein$sequence[sel], first = start, last = end)

names(targetSeq) <- targetName # Name it

# ====  Step 1.2: MSA =============================================================
# Note: a similar step should has been done when generating the first Model1BM8

sel <- grep("^MBP1_", myDB$protein$name) # get all MBP1 Sequences
MBP1Set <- myDB$protein$sequence[sel] # extract all the MBP1 sequences
names(MBP1Set) <- myDB$protein$name[sel] # name the sequences
seq1BM8 <- dbSanitizeSequence(readLines("1BM8_A.fa")) # read the template sequence
names(seq1BM8) <- "1BM8_A"
seq4UX5 <- dbSanitizeSequence(readLines("4UX5_A.fa")) # read the template sequence
names(seq4UX5) <- "4UX5_A"

MBP1Set <- c(MBP1Set, seq1BM8, seq4UX5) # Add the template sequences to the MBP1set
MBP1Set <- Biostrings::AAStringSet(MBP1Set) # Turn it into an Biostrings::AAStringSet

MBP1msa <- msa::msaMuscle(MBP1Set) # Calculate an msa by msaMuscle

library(Biostrings)
writeMFA(fetchMSAmotif(MBP1msa, seq4UX5), myCon = "APSES-MBP1_new.fa") # Write the alignments to file

# ====  Step 1.3: Extract =============================================================
# Objective: extract the aligned target and 4UX5 template sequence, while masking gaps
# that are not needed for the aligned pair.

myT <- seq4UX5 #set template as 4UX5
targetSeq   <- as.character(fetchMSAmotif(MBP1msa, myT)[targetName])
templateSeq <- as.character(fetchMSAmotif(MBP1msa, myT)[names(myT)])

# Drop positions in which both sequences have hyphens.
targetSeq   <- unlist(strsplit(targetSeq,   ""))
templateSeq <- unlist(strsplit(templateSeq, ""))
gapMask <- ! ((targetSeq == "-") & (templateSeq == "-"))
targetSeq   <- paste0(targetSeq[gapMask], collapse = "")
templateSeq <- paste0(templateSeq[gapMask], collapse = "")

# Assemble sequences into a set
TTset <- character()
TTset[1] <- targetSeq
TTset[2] <- templateSeq
names(TTset) <- c(targetName, names(myT))

writeMFA(TTset)

# >MBP1_THETH
# ------MAADHPKA-GIYSATYSGIPVYEYQFGPDLKEHVMRRREDNWINATHILKAAGFDKPARTRILERDVQKDIHEK
# IQGGYGKYQGTWIPLEHGEALAQRNNVYERLRPIFEFQPGNESPPPAPRHASKPKVP-
#
# >4UX5_A
# MVKAAAAAASAPTGPGIYSATYSGIPVYEYQFG--LKEHVMRRRVDDWINATHILKAAGFDKPARTRILEREVQKDQHEK
# VQGGYGKYQGTWIPLEAGEALAHRNNIFDRLRPIFEFSPGPDSPPPAPRHTSKPKQPK

# ====  Step 1.4: SwissModel =============================================================
# The step is done on the online SwissModel server.

# Process: Paste the aligned sequences of the MYSPE target and the template into the form field. Target is THETH,
# and template is 4UX5. Build Model and save PDB file as "MBP1_MYSPE-APSES_new.pdb" in the directory.

# GMQE: 0.72	QMEAN: -1.08
# GMAE ranges from 0 to 1 and reflects the expected accuracy of a model built with that alignment and template
# and the coverage of the target. 0.72 is not bad.
# QMEAN Z-score indicates whether the QMEAN score of the model is comparable to what one would expect from
# experimental structures of similar size. It is below 0 but not far away from 0.
# Overall, the model is acceptable based on the two metrics.


# ==== Step 1.5: Renumbering the model ================================================
PDB_INFILE      <- "MBP1_MYSPE-APSES_new.pdb"
PDB_OUTFILE     <- "MBP1_MYSPE-APSESrenum_new.pdb"

iFirst <- 14  # residue number for the first residue if your template was 4UX5

# Read the MYSPE pdb file
MYSPEmodel <- bio3d::read.pdb(PDB_INFILE) # read the PDB file into a list

# Modify residue numbers for each atom
resNum <- as.numeric(MYSPEmodel$atom[,"resno"])
resNum <- resNum - resNum[1] + iFirst  # add offset
MYSPEmodel$atom[ , "resno"] <- resNum   # replace old numbers with new

# ==== Step 1.6: Write output to file ================================================
bio3d::write.pdb(pdb = MYSPEmodel, file="MBP1_MYSPE-APSESrenum_new.pdb")

# ==== Step 2: Compare two models ====================================================
# Process: superimpose the two models, calculate the RMSD between each residue pair,
# writes the RMSD to each residues' B Factor field for both models, and save the resulting PDB files.

library(bio3d)
Model4UX5 <- bio3d::read.pdb("MBP1_MYSPE-APSESrenum_new.pdb")
superimpose <- bio3d::pdbaln(list(Model1BM8, Model4UX5), fit = TRUE, web.args =
                                   list(email = "xiaow.zhang@mail.utoronto.ca")) # do sequence alignment from the list of two PDB objects; fit = TRUE means performing coordinate superposition on the input structures.
superimpose$xyz # two-row numeric matrix of aligned C-alpha coordinates
superimpose$resno # two-row matrix of aligned residue numbers
# Note: the first row corresponds to Model1BM8, and the second is Model4UX5.

RMSD <- numeric(ncol(superimpose$resno)) # create a vector to store RMSD for each superposed residue pair
for (i in seq_len(ncol(superimpose$resno))) {
  inds <- seq(from = 1, to = ncol(superimpose$xyz), by=3) # 3 numbers per xyz-coordinate; provide the first index numbers of xyz-coordinates
  index <- inds[i]
  RMSD[i] <- rmsd(superimpose$xyz[1, index:index + 2], superimpose$xyz[2, index:index + 2]) # calculate the RMSD between each residue pair

  # writes the RMSD to each residues' B Factor field for both of the two models according to the aligned residue numbers
  a.ids <- superimpose$resno[1, i] # the aligned residue number for Model1BM8
  b.ids <- superimpose$resno[2, i] # the aligned residue number for Model4UX5
  Model1BM8$atom$b[Model1BM8$atom$resno == a.ids] <- RMSD[i] # use the residue number to find the corresponding residue in Model1BM8 and replace its b-factor as RMSD
  Model4UX5$atom$b[Model4UX5$atom$resno == b.ids] <- RMSD[i] # do the same thing for Model4UX5

  #Since Model4UX5 is longer than Model1BM8, for those gaps, assign RMSD as the highest RMSD between the coordinates.
  Model4UX5$atom$b[is.na(Model4UX5$atom$b)] <- max(RMSD, na.rm = TRUE)
  }

# save the pdb files
write.pdb(pdb = Model1BM8, file = "MBP1_MYSPE-APSES_RMSD_1BM8.pdb")
write.pdb(pdb = Model4UX5, file = "MBP1_MYSPE-APSES_RMSD_4UX5.pdb")

# ==== Step 3: Chimera ====================================================
# Process:  Load the two models, and superimpose them. Color them by the RMSD values I have computed.
# Save the stereo image and submit.

# See the detailed process below.

# [END]


