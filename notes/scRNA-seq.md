# Single cell RNA sequencing (scRNA-seq) procedure

Estimates gene expression in a single cell at a point through extracting the mRNA transcripts in the cytosol and quantifying their association with certain genes. 

## Cell extraction

Collect a group of cells which are physically independent from one another (lack of cellular junctions). Can achieve this via in vitro cells (i.e. cell culture), harvesting from solid tissues, or harvesting from mobile cell populations (e.g.    hematopoietic cells). 

Harvesting from solid tissues is followed by enzymatic digestion and filtration of the living cellular components. In most cases, there are quality checks to ensure that a sufficient proportion of cells survived the digestion process otherwise the sample is excluded. 

## Cell capture

Microfluidic chip with three input streams: filtered cell mixture, gel beads containing biological machinery, and oil. The chip functions to route the inputs such that **Gel Bead-in-Emulsion (GEM) droplets** are formed. Ideally the droplets contain one gel bead and one cell within the oil solvent. 

The gel bead carries millions of copies of characteristic oligonucleotides structured as: 
- 5' [PCR handle]-[16bp cell barcode]-[12bp unique molecular identified (UMI)]-[poly-dThymine tail] 3'
    - The 16 bp cell barcode identifies the bead and therefore cell (interbead distinction)
    - The UMI identifies the specific olignonucleotide (intrabead distinction)


## mRNA capture (intradroplet) 

Within the droplet the gel bead dissolves releasing its oligonucleotides. The cell is also lysed and thus within a GEM droplet we have circulating gel bead oligonucleotides and cell contents (only mRNA molecules are relevant). **The poly-deoxythymidine tail on the oligonucleotide binds the poly-adenosine tail on the mRNA transcripts capturing it.**
 - 5' [PCR handle]-[16bp cell barcode]-[12bp (UMI)]-[poly-dT] 3'
 - 3' [poly-A]-[UTR]-[Coding sequence]-[UTR] 5' 

**Reverse transcription of captured mRNA.** Reverse transcriptase enzymes are also released from the gel bead which reverse transcibe the mRNA into complementary DNA (cDNA) which uses the poly-dT tail as a primer (this is why we must use deoxythymidine instead of thymidine). Note that only a portion of the coding sequence is reverse transcribed, not all of it. 

**The resulting structure: 5' [PCR handle]-[barcode]-[UMI]-[CDS] 3'**

## Amplification and processing of cDNA (extradroplet)

**Polymerase chain reaction (PCR)**: After reverse transcriptase generates the cDNA, the droplets are dissolved (recall the bead and cell existed within the droplet). We then use PCR to amplify the cDNA structure (including the barcode and UMI, etc.). 

The cDNA structures are then sequenced (generation of digital ordering of AGCT nucleotides). 

**Processing of sequencing**: the digital format of the sequences is matched to the human genes (producing a read for that gene upon every match). For each barcode and gene combination, reads with the same UMI are collapsed to a single read (these are duplicates from PCR). There are minor processing steps within this to ensure the count matrix is as accurate as possible. 

The result is a map of human genes and a positive integer representing how many unique cDNA structures were associated with that gene. The higher the count indicates that more mRNA transcripts were present in the cell and thus that cell was expressing that gene further (vice versa). 

## Limitations

scRNA-seq data is notoriously sparse, many genes have 0 mRNA transcipt reads, which is due to fundamental yield limitations of the technique, heterogeneity, and decreased expression.

The experimental procedure is a narrowing funnel which results in potential data loss at multiple steps. The largest loss occurs during mRNA capture where oligonucleotides may fail to bind mRNA transcripts (only 10-40% of mRNA transcripts are actually bound).

Biological heterogeneity also results in some genes having 0 reads despite being housekeeping genes (constitutively active). This is a result of the granularity of examining a single cell at a time. If you were to analyze thousands of cells housekeeping genes would certainly show high gene expression. 