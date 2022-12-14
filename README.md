# PhD Supplementary Files: *Tetramesa* phylogenetics

This GitHub repository serves as supplementary material to my PhD thesis, entitled:

**"A molecular investigation of stem-galling _Tetramesa_ Walker (Hymenoptera: Eurytomidae) on African grasses: an application to biological control"**

Department of Zoology and Entomology, Centre for Biological Control, Rhodes University, Grahamstown, South Africa   

Supervisors: Prof. Iain D. Paterson and Dr. Guy F. Sutton

<img src="https://github.com/clarkevansteenderen/PhD_files/blob/main/tetramesa.png" height = 120>

## Instructions :pen:

Download the repository, unzip the folder, and open the R project **"PhD_curated_code"**. The R script **"data_analyses"** contains all the code used for the analyses used in my thesis. 

The **PopArt** folders contain the input files for haplotype networks in the PopArt program, for both the COI and 28S data, and the folder containing other species delimitation results are the output files from the ABGD, ASAP, PTP, bPTP, and Haplowebs (COI only).

The **FIGURES** folder contains all the figures used in the thesis, organised by chapter.

The **input_files** folder contains the input files of aligned sequences (COI and 28S) ready for MrBayes and IQ-TREE.

The structure of the R code follows the order below:

#### CHAPTER 2 :page_with_curl:
--- 

**PART 1: PHYLOGENY PLOTTING**    
Plots of the Bayesian and Maximum Likelihood phylogenies for the COI and 28S gene regions, and the dated 28S BEAST tree.

**PART 2: GMYC SPECIES DELIMITATION TESTS**   
Species delimitation tests using the single-threshold GMYC model on the COI and 28S FASTA alignments.

**PART 3: P-DISTANCES**   
Genetic p-distances for the COI and 28S data, and the generation of heatmaps for pairwise comparisons between clades, and barcode gap plots. Various threshold values are used to find when groups become conspecific.

**PART 4: AMOVAS FOR THE GOLDEN-SHOULDERED TETRAMESA ON ERAGROSTIS CURVULA**      
Analysis of Molecular Variance (AMOVA) for the golden-shouldered *Tetramesa* on *Eragrostis curvula* across four provinces (Eastern Cape, Free State, Mpumalanga, and the Northern Cape), and within the Eastern Cape only.

**PART 5: MANTEL TESTS: GENETIC DISTANCE VS GEOGRAPHIC DISTANCE**   
Mantel tests checked whether there was a correlation between genetic distance and geographic separation for the golden-shouldered *Tetramesa* on *Eragrostis curvula*.

#### CHAPTER 3 :page_with_curl:
--- 

**PART 6: COPHYLOGENETIC ANALYSES**      
Global-fit cophylogenetic analyses using a wasp COI Maximum Likelihood phylogeny, and chloroplast (rpl32-trnL, rps16-trnK, and rps16) and nuclear (ITS) trees for the grasses. A PACO analysis is included with phylogenies containing representative sequences from clades. 

#### CHAPTER 4 :page_with_curl:
---

**PART 7: SPEDE-SAMPLER**          
Includes all the plotting for the results obtained from the running of the SPEDE-Sampler software on the COI and 28S data.

#### MISCELLANEOUS
--- 

**PART 8: MAPS** :globe_with_meridians:   

**(1) TETRAMESA SPECIES DESCRIPTIONS GLOBALLY**    
This world map shows the sampling effort for *Tetramesa* descriptions, where there is a strong Northern-Hemisphere bias.

**(2) GRASS DISTRIBUTIONS IN NATIVE AND INVADED RANGES**    
These maps show the distributions of a query grass species in its native and invasive range, using the GBIF database.

The scripts in (1) and (2) were modified from those kindly provided by Guy Sutton.

**PART 9: SUMMARY BARPLOT FOR THE SPECIMENS COLLECTED**         
This barplot shows how many wasp specimens were collected on each host grass of interest, and the province it was from.
