# Protocol for multiplex PCR and sequencing of Zika virus from clinical isolates

Version 1.0 (Jan 5 2016)

Allison Black (black.alli@gmail.com)

_please note that I did not develop this protocol; all credit goes towards the individuals listed below._

* Multiplex PCR: Josh Quick and Nick Loman, University of Birmingham
* MinION sequencing protocol: Josh Quick and Nick Loman, University of Birmingham
* MiSeq sequencing protocol: Nathan Grubaugh, The Scripps Research Institute

_additional thanks to Josh and Nate for lots of useful advice, which I've tried to work in to the protocol as helpful hints_

## Equipment and Reagent Checklist

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

## Protocol

## Viral RNA extraction

_We use the QiaAMP Viral RNA mini kit, and the protocol below is based off of the manufacturer's instructions. The protocol below is written out assuming a clinical sample volume of 200 uL. Note that you need 4x the amount of Buffer AVL as serum, so set Buffer AVL and carrier RNA amounts based off of your available sample volume. The volume of carrier RNA that should be used is 1/100th the volume of buffer AVL added._

_We've also run this protocol on the QIAcube. To make the protocol compatible we use 140 uL of serum/urine and elution volume of 50 uL_

_Other equivalent extraction methods can also be used (e.g. Trizol, Omega etc.)_

_In addition to the samples, remember to have an extraction negative control (add nuclease-free water instead of serum/urine)._

1.	Remove serum samples from freezer and allow to come up to room temperature.
2.	Prepare stock solution of Buffer AVL + carrier RNA. _Per sample_ add:

  *	800 uL Buffer AVL
  * 8 uL Carrier RNA (previously resuspended in Buffer AVE)

3.	Pipette 800 uL of Buffer AVL containing carrier RNA into 2.0 ml tube.
4.	Add 200 uL of serum or urine to the tube with Buffer AVL.
5.	Mix sample and buffer AVL by pulse-vortexing for 15 seconds.
6.	Incubate mixture at room temperature for 10 minutes. Spin down.
8.	Add 800 uL of 100% Ethanol to the sample + Buffer AVL mixture.
9.	Mix by pulse vortexing for 15 seconds.
10.	Spin down fluid from tube walls.
11.	Take 2ml collection tubes and put the columns in the tubes.
12. Perform the following 3 times in order to have run all 1800 uL of ethanol+sample+Buffer AVL mixture through the filter column.

  * Pipette 600 uL ethanol+sample+Buffer AVL mixture into the column directly onto the filter. _do not touch tip to filter._
  * Centrifuge at room temperature at 8000 rpm for 1 minute. Take column out and place in a new 2 mL collection tube.

13.	Add 500 uL of Buffer AW1 to the column. _this volume does not need to be increased even if the sample volume was greater than 200 uL._
14. Spin sample at 8000 rpm for 1 minute. Transfer column to a new 2 mL collection tube.
15.	Add 500 uL of Buffer AW2.
16.	Spin at full speed (14,000 RPM) for 3 minutes. Place column in a new collection tube and spin at full speed for 1 minute.
17.	Place column in clean, labelled 1.5 mL tube. Add 50 uL of Buffer AVE at room temperature to column.
18.	Incubate at room temp for 1 minute.
19.	Centrifuge at room temp at 8000 rpm for 1 minute.
20.	Place extracted RNA in eluate on ice if proceeding directly to PCR or store in freezer.

## Two-step RT-PCR for amplicon generation

_Perform the following in a hood or in a pre-PCR designated area_

#### Reverse Transcription

1. For each sample, mix the following in a PCR tube (preferably individually-hinged strip tubes):

  * 7uL RNA
  * 1 uL random hexamers (50 ng/uL)

2. Mix by inversion then spin down.
3. Heat to 70 degrees Celsius for 7 minutes, then place on ice to prevent secondary structure from re-forming.
4. Add the following to each tube:

  * 10 uL ProtoScript II Reaction Mix
  * 2 uL ProtoScript II Enzyme Mix

5. Place on the thermocycler and run the following program:

  * 25 degrees Celsius for 5 minutes _**note: this step can be swapped for a 5 minute hold at room temperature_
  * 48 degrees Celsius for 15 minutes  
  * 80 degrees Celsius for 5 minutes
  * Hold at 10 degrees Celsius

#### Amplification PCR

6. Make the following mastermix in a 1.5 mL tube. You will need both a Pool 1 mastermix and a Pool 2 mastermix.

  * 8 uL Q5 Reaction Buffer
  * 0.8 uL 10 mM dNTPs
  * 0.4 uL Q5 DNA polymerase
  * 24.4 uL Nuclease-free water
  * 2.7 uL of primer pool (either pool 1 **OR** pool 2)

7. Label 2 0.2mL PCR tubes for each sample.
8. Pipette 36 uL of mastermix into the corresponding tubes.
9. Add 4 uL of cDNA into the corresponding tubes.
10. Mix by inversion and spin down.
11. Place on the thermocycler and run the following program:

  * Step 1: 98 degrees Celsius for 30 seconds (initial denaturation)
  * Step 2: _perform 45 cycles of_
    * 98 degrees Celsius for 15 seconds
    * 65 degrees Celsius for 5 minutes
  * Step 3: Hold at 4 degrees Celsius.

## Post-PCR clean-up and amplicon quantification

#### Clean-up

1. Label two sets of DNA Lo-Bind tubes (one set for bead clean-up and the other for eluate with cleaned up amplicons)
2. Allow AMPure XP beads to come up to room temperature, and homogenize by vortexing.
3. Add 72 uL of AMPure XP beads to each tube in one set of the Lo-Bind tubes. _this is a 1.8x bead clean-up_
4. Pipette 40 uL of PCR product into correspond tube with beads.
5. Incubate beads and PCR products at room temperature on a hula mixer for 5 minutes. _if no access to a hula mixer, you can flick mix gently throughout the incubation_
6. Place tubes on magnetic rack and incubate until the solution is clear.
7. Discard the supernatant being careful to not disturb pellet.
8. Add 200 uL of 80% EtOH to each tube to wash pellet, incubate 30 seconds, then discard EtOH wash.
9. Repeat previous wash again, pipetting off as much EtOH as possible.
10. Spin tubes down, replace on magnetic rack, and pipette off any additional EtOH. Leave tubes open to dry.
11. Allow pellets to dry (roughly 5 minutes). The pellet should appear matte but not dry to the point where cracks form.
12. When pellet is dry add 31 uL of nuclease-free water to each tube, remove from rack and flick gently to resuspend beads. _if storing amplicons for longer periods of time, qiagen Elution Buffer can be used instead of nuclease-free water._
13. Once pellets have been resuspended incubate at room temperature on a hula mixer for 5 minutes.
14. Replace tubes on magnetic rack and incubate until solution fully clears.
15. Carefully pipette off 31 uL of supernatant without disturbing beads and place into new Lo-Bind tubes.

#### Quantification

_You should use the Qubit High Sensitivity dsDNA kit_

1. Make up Qubit Mastermix in a falcon tube :

  * 199 uL Buffer _per sample_
  * 1 uL dye _per sample_

2. Add 190 uL mastermix + 10 uL standard 1 into one Qubit tube.
3. Add 190 uL mastermix + 10 uL standard 2 into another Qubit tube.
4. Add 199 uL mastermix + 1 uL cleaned-up PCR products into Qubit tube for each pool.
5. Quantify the concentration of dsDNA in the sample.

  * Negative extraction controls usually have concentrations around 3 or 4 ng/uL.
  * Negative PCR controls should have concentrations < 1 ng/uL.
  * Properly amplified samples should have between 5 and 100 ng/uL.

## MinION Library Preparation

#### Normalize amplicon concentrations

_turn on 65 degree Celsius dry bath on now_

1. Divide 1500 ng by the number of _non-negative control_ samples you will run. This gives you the desired amount of DNA for each amplicon pool. Generally this will be:

  * 150 ng of DNA per amplicon pool if running 5 genomes (10 pools) and negative controls.
  * 125 ng of DNA per amplicon pool if running 6 genomes.

  _note that the PCR control should always be sequenced_

2. Determine how many uL of cleaned-up amplicons you need to add to reach the target amount of DNA. Bring this volume up to 30 uL with nuclease-free water if the full amount of cleaned-up amplicons is not used.

_example where target amount is 150 ng per pool_

| sample        | pool concentration | volume amplicons | volume water |
| ------------- |--------------------| -----------------|--------------|
| samp1_pool1   | 41.8 ng/uL         | 3.6 uL           | 26.4 uL      |
| samp1_pool2   | 55.4 ng/uL         | 2.7 uL           | 23.7 uL      |

#### End repair and clean-up

1. Add the following to each sample:

  * 4.2 uL Ultra II End-Prep buffer
  * 1.8 uL Ultra II End-Prep enzyme mix

2. Mix by inversion and spin down. Total reaction volume is 36 uL.
3. Incubate at 20 degrees Celsius for 5 minutes. _this incubation can also be at room temperature_
4. Incubate at 65 degrees Celsius for 5 minutes.
5. Perform 1:1 bead clean-up with AMPure beads (36 uL of beads to 36 uL of end-prepped DNA) as described above (5 minute incubations mixing throughout, two 200 uL 80% EtOH washes).

**this time resuspend pellet in 15 uL nuclease-free water**

#### Barcoding and pooled clean-up

1. Thaw a native barcode for each sample. _for 6 genomes you'll have 12 barcodes since each pool has it's own barcode_
2. Add 5 uL of one of the barcodes NB01-NB12 to the appropriate sample. Mix by inversion.
3. Add 20 uL of Blunt/TA Ligase Master Mix to each tube. Mix by inversion. Spin down.
4. Incubate at room temperature for 10 minutes.
5. Incubate at 65 degrees Celsius for 5 minutes to kill ligase. _Ensure that the temperature is correct. Samples will be pooled after this step so you need to prevent any further ligation that may cause sample-barcode mixtures._
6. Pool all of the samples together into a single Lo-Bind tube. If you have 12 samples the total volume will be 480 uL.
7. Add 480 uL AMPure XP beads to the pooled samples (1:1 clean-up). _if fewer than 12 samples, add the same volume of beads as pooled samples_
8. Incubate at room temperature on a hula mixer for 5 minutes.
9. Place tube in magnetic rack and wait until solution fully clears.
10. Pipette off and discard the supernatant. Be careful to not disturb the pellet.
11. Perform 2 washes with 80% EtOH as described above. Use sufficient EtOH to cover the pellet (this might be more than 200 uL).
12. Allow pellet to air dry to the point that the pellet looks matte but not cracked. _because of the size of the pellet this can take a while. You can speed this up a bit by incubating for short periods at 65 degrees Celsius. Ensure you keep checking the pellet periodically though so you don't over dry._
13. Elute the pooled reaction in 39 uL of nuclease-free water, resuspend pellet and incubate at room temperature on a hula mixer for 5 minutes. Remove the supernatant to a clean Lo-Bind tube.

#### Pooled library quantification

1. Prepare Qubit standards and sample reactions as described above.
2. Target is >1000 ng DNA in the remaining 38 uL of pooled library.

#### MinION adaptor ligation

1. Add 10 uL BAM to cleaned-up library. Mix by inversion.
2. Then add 2 uL BHP. Mix by inversion.
3. Then add 50 uL Blunt/TA Ligase Master Mix. Mix by inversion and spin down.
4. Incubate at room temperature for 10 minutes. _start making the cleaned MyOne C1 beads during this incubation, see directions below_
5. Add 1 uL HPT. Mix by inversion and spin down.
6. Incubate at room temperature for 10 minutes.

> MyONE C1 bead preparation

> 1. Vortex MyOne C1 beads until homogenous.

>2. Pipette 50 uL of beads and transfer to new Lo-Bind tube.

>3. Pellet beads on magnetic rack and discard supernatant.

>4. Wash beads twice with 100 uL BBB pipetting to resuspend. Ensure to scrape beads off of tube sides when resuspending pellets. Discard supernatant. Do not let the beads dry after the second wash.

>5. Resuspend cleaned beads in 100 uL BBB. These cleaned beads will be used for adaptor ligation clean-up.

#### Purify the adapted,tethered library

1. Add 100 uL washed MyOne C1 beads to the adapted,tethered library.
2. Incubate at room temperature on a hula mixer for 5 minutes.
3. Place on magnet. Once solution clears discard the supernatant.
4. Wash pellet twice with 150 uL of BBB pipetting to resuspend the pellet **after each wash**.
5. Spin down, replace tube on magnet and pipette off any residual BBB. Close lid, **do not let pellet dry**.
6. Resuspend the beads in 25 uL of ELB. _you'll need to scrape the beads off of the tube walls to ensure you fully resuspend the beads. This can take a while, but it's important to do it well to get sufficient library. The tube sides should not show any smear of beads._
7. Pellet beads on magnetic rack and transfer eluate (tethered,adapted library in ELB) to new Lo-Bind tube.

#### Adapted, tethered library quantification and normalization

1. Make up Qubit Master Mix and standards as before.
2. Quantify adapted, tethered library using 1 uL of library to 199 uL Qubit Master Mix.
3. Determine the total ng of library, target is >100 ng. _if you have >200 ng of library you have enough adapted, tethered library for a second run. You can save this leftover library by keeping it at -20 degrees Celsius._
4. Normalize the library to have 100 ng of library loaded on the flowcell. Bring library up to 37.5 uL volume with nuclease-free water. _multiply the library concentration by 24 uL (amount of library left post quantification) to determine total ng of library._

| Library       |   concentration    | total_amt_library | volume_library_needed | volume_nfw |
| ------------- |--------------------| ------------------|-----------------------|------------|
| library1      | 9.9 ng/uL          |      237.6 ng     |       10.1 uL         | 27.4 uL    |

5. Add 37.5 uL of RBF1 at room temperature to the library/nuclease-free water solution. Total library+RBF1 volume will be 75 uL. _this volume is specific to R9.4 spot on flow cells. If you're using an older non-spot on flowcell follow the protocol for your flowcell version._
6. Mix by inversion and spin down.

## Prime MinION flowcell and load the library (with lots of notes drawn from experience making mistakes)

**protocol assumes a R9.4 spot on flowcell**

_Ensure you've aleady QC'd the flowcell to determine the number of active pores._

_Because the priming process takes a bit of time, you can do the priming steps during incubation periods earlier on in the protocol._

1. Prepare the priming solution by mixing 500 uL of RBF1 with 500 uL of nuclease-free water. Mix by inversion and spin down. _use the same pipette to ensure that the ratio is correct even if there is error in the pipette._
2. Open the sample port on the flowcell (leave the SpotON port closed). Draw back a few uL of buffer to ensure there's fluid up to the lip of the sample port. _this is to prevent you from introducing bubbles when priming the flowcell. This is easiest to do by inserting the tip of the pipette into the port and rather than drawing up volume using your normal pipetting action, simply adjust the pipette to increase the volume it holds, which very slowly allows you to draw up the buffer, ensuring you don't take up too much. This will also prevent shakes where you might accidentally push more buffer than you drew up into the flowcell, introducing a bubble._
3. Load 500 uL of priming solution via the sample port. Wait 10 minutes. _again, it's more accident proof to load this by adjusting the pipette volume down rather than using a normal pipetting push down. Leave a few uL of priming solution in pipette tip to prevent introduction of a bubble._
4. Load 300 uL of priming solution via the sample port. Wait 10 minutes.
5. Lift the SpotON tab to reveal the SpotON port. Load 200 uL of priming solution via the SpotON port. _you want to avoid hitting the port with the pipette tip. To load slowly push a droplet of fluid out of the tip and touch the edge of the droplet to the port which should pull the droplet off the tip and into the port._
6. Load the library+RBF1 mixture via the SpotON port using the same technique. _IMPORTANT: sometimes the capillary action doesn't work, and your first droplet of library just sits on the SpotON port without getting sucked in. If this happens don't panic, just make up another 200 uL of priming solution (100 uL of RBF1 + 100 uL of nuclease-free water) and repeat your 200 uL prime via the SpotON port. Then load the library after this second spot on prime._
7. Close the sample port and SpotON port, making sure the bung goes into the SpotON port hole. Close the MinION lid.
8. In MinKNOW choose the MAP_48Hr_Sequencing_Run_SQK_LSK208_FLOMIN106.py program. _Frustratingly there are pretty regular updates for MinKNOW and the autoupdate function doesn't seem to work on Mac. So if you're using a Mac as your sequencing laptop you'll probably need to do a complete uninstall and reinstall of MinKNOW each time there is an update. Also, if the MinKNOW run isn't working once you've started the sequencing run (e.g. performs calibration but never moves along to sequencing) try uninstalling and reinstalling MinKNOW and restarting the run._
