# Protocols for Zika sequencing

#### Allison Black<sup>1,2</sup>

<sup>1</sup>Department of Epidemiology, University of Washington, Seattle, WA, USA, <sup>2</sup>Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, WA, USA

## Multiplex PCR for MinION sequencing of Zika virus

*This protocol is a bench-ready version of the protocol presented in the following paper:*

> Quick J, Grubaugh ND, Pullan ST, Claro IM, Smith AD, Gangavarapu K, Oliveira G, Robles-Sikisaka R, Rogers TF, Beutler NA, Burton DR, Lewis-Ximenez LL, Goes de Jesus J, Giovanetti M, Hill S, Black A, Bedford T, Carroll MW, Nunes M, Alcantara LCJ, Sabino EC, Baylis SA, Faria N, Loose M, Simpson JT, Pybus OG, Andersen KG, Loman NJ. 2017. Multiplex PCR method for MinION and Illumina sequencing of Zika and other virus genomes directly from clinical samples. In prep.

A pre-print of the paper can be found [here](http://biorxiv.org/content/early/2017/01/09/098913)

*Additional thanks to Josh Quick and Nate Grubaugh for lots of useful advice, which I've tried to work in to the protocol as helpful hints.*

### Bench protocols

Regardless of which platform you'll be sequencing on, you'll need to follow these first four protocols to generate two pools of ZIKV amplicons

1. [Make primer pools](primer-pool-recipe.md) _Note: the first time you do this, make up multiple aliquots of the pools_
2. [Viral RNA extraction](viral-rna-extraction.md)
3. [Two-step RT-PCR for amplicon generation](amplicon-generation.md)
4. [Post-PCR clean-up and amplicon quantification](cleanup-and-amplicon-quantification.md)

#### Protocols for MinION

5. [MinION library preparation](minion-library-preparation.md)
6. [MinION flowcell priming and library loading](minion-flowcell-loading.md)

#### Protocols for Illumina MiSeq

7. [MiSeq library preparation](miseq-library-preparation.md)

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

##### MiSeq sequencing

* KAPA Hyper Library Prep (Roche, cat. no. 07962363001)
* SureSelectxt2 indexes, MSQ, 16 (Agilent, cat. no. G9622A)
* MiSeq Reagent Kit v2 (500 cycle) (Illumina, cat. no MS-102-2003) _Note: this might be provided by your Core facility_
* D1000 ScreenTape (Agilent, cat. no. 5067-5582)
* D1000 Reagents (Agilent, cat. no. 5067-5583)
* KAPA Library Quantification Kit for Illumina platforms (Roche, cat. no 07960140001) _may not be necessary if you will be submitting the library to a core facility that does dilutions, pooling, and loading for you. Check with your core facility regarding what services they provide._
