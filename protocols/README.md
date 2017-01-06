# Protocols for Zika sequencing

#### Allison Black<sup>1,2</sup>

<sup>1</sup>Department of Epidemiology, University of Washington, Seattle, WA, USA, <sup>2</sup>Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA

## Multiplex PCR for MinION sequencing of Zika virus

*This protocol is a bench-ready version of the protocol presented in the following paper:*

> Quick J, Grubaugh ND, Pullan ST, Claro IM, Smith AD, Gangavarapu K, Oliveira G, Robles-Sikisaka R, Rogers TF, Beutler NA, Burton DR, Lewis-Ximenez LL, Goes de Jesus J, Giovanetti M, Hill S, Black A, Bedford T, Carroll MW, Nunes M, Alcantara LCJ, Sabino EC, Baylis SA, Faria N, Loose M, Simpson JT, Pybus OG, Andersen KG, Loman NJ. 2017. Multiplex PCR method for MinION and Illumina sequencing of Zika and other virus genomes directly from clinical samples. In prep.

*Additional thanks to Josh Quick and Nathan Grubaugh for lots of useful advice, which I've tried to work in to the protocol as helpful hints.*

### Bench protocols

1. [Viral RNA extraction](viral-rna-extraction.md)
2. [Two-step RT-PCR for amplicon generation](amplicon-generation.md)
3. [Post-PCR clean-up and amplicon quantification](cleanup-and-amplicon-quantification.md)
4. [MinION library preparation](minion-library-preparation.md)
5. [MinION flowcell priming and library loading](minion-flowcell-loading.md)

### Equipment and Reagent Checklist

##### Tiling amplicon generation
* QIAamp Viral RNA Mini Kit (Qiagen, cat. no. 52906)
* Random Hexamers (50 µM) (Thermo Fisher, cat. no. N8080127)
* Protoscript II First Strand cDNA Synthesis Kit (NEB, cat. no. E6560)
* Deoxynucleotide (dNTP) Solution Mix (NEB, cat. no. N0447)
* Q5 Hot Start High-Fidelity DNA Polymerase (NEB, cat. no. M0493)
* Custom PCR primers [click here for primer sequences](zika-multiplex-primers.xls)
* Agencourt AMPure XP beads (Beckman Coulter, cat. no. A63881)
* Qubit dsDNA HS Assay Kit (Thermo Fisher, cat. no. Q32854)
* HyClone Molecular Biology Grade Water (GE Life Sciences, cat. no. SH30221.10)
* 100% Ethanol

##### MinION sequencing
* SpotON Flow Cell Mk I (R9.4) (Oxford Nanopore Technologies, cat. no. FLO-MIN106)
* Nanopore Sequencing Kit (R9) (Oxford Nanopore Technologies, cat. no. SQK-NSK007)
* Native Barcoding Kit (Oxford Nanopore Technologies, cat. no. EXP-NBD002)
* NEBNext Ultra II End-repair/dA-tailing Module (NEB, cat no. E7546)
* NEB Blunt/TA Ligase Master Mix (NEB, cat. no. M0367)
* MyOne C1 Streptavidin beads (Thermo Fisher, cat. no. 65001)
