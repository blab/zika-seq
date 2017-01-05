# Protocol for multiplex PCR and sequencing of Zika virus from clinical isolates

Version 1.0 (Jan 5 2016)
Allison Black (black.alli@gmail.com)

_please note that I did not develop this protocol; all credit goes towards the individuals listed below._
* Multiplex PCR: Josh Quick and Nick Loman, University of Birmingham
* MinION sequencing protocol: Josh Quick and Nick Loman, University of Birmingham
* MiSeq sequencing protocol: Nathan Grubaugh, The Scripps Research Institute

_additional thanks to Josh and Nate for lots of useful advice, which I've tried to work in to the protocol where appropriate_

## Equipment and Reagent Checklist

## Protocol

### Viral RNA extraction

1.

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

> MyONE CI bead preparation
1. Vortex MyOne C1 beads until homogenous.
2. Pipette 50 uL of beads and transfer to new Lo-Bind tube.
3. Pellet beads on magnetic rack and discard supernatant.
4. Wash beads twice with 100 uL BBB pipetting to resuspend. Ensure to scrape beads off of tube sides when resuspending pellets. Discard supernatant. Do not let the beads dry after the second wash.
5. Resuspend cleaned beads in 100 uL BBB. These cleaned beads will be used for adaptor ligation clean-up.

#### Purify the adapted,tethered library

1. Add 100 uL washed MyOne C1 beads to the adapted,tethered library.
2. Incubate at room temperature on a hula mixer for 5 minutes.
3. Place on magnet. Once solution clears discard the supernatant.
4. Wash pellet twice with 150 uL of BBB pipetting to resuspend the pellet **after each wash**.
5. Spin down, replace tube on magnet and pipette off any residual BBB. Close lid, **do not let pellet dry**.
6. Resuspend the beads in 25 uL of ELB. _you'll need to scrape the beads off of the tube walls to ensure you fully resuspend the beads. This can take a while, but it's important to do it well to get sufficient library. The tube sides should not show any smear of beads._
7. Pellet beads on magnetic rack and transfer eluate to new Lo-Bind tube.

#### Adapted, tethered library quantification
