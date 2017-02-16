# MiSeq library preparation

This bench protocol is our version of this [protocol developed by Nathan Grubaugh.](https://docs.google.com/document/d/1PilT4w5jHO-ROsE8TL5WBGa0wSCdTHAsNl1LIOYiTgk/edit)

#### Normalize amplicon concentrations

1. Based on the amplicon concentration you calculated post PCR-amplification and clean-up, determine the number of uL of amplicons you need to have 50 ng of DNA. Bring this volume up to 50 uL with nuclease-free water. **have 50 uL samples in 0.2mL PCR tubes** _Note: if your amplicon concentrations are sufficiently high that you would need less than 1 uL of your pool, then just add 1 uL to prevent issues with pipetting smaller volumes._

  _Once you've seen that things are working well and that you have sufficient coverage, you can run more samples at a time by combining pool 1 and pool 2 at this point. To do this combine 25 ng of amplicons from pool 1 and 25 ng from pool 2, and bring this up to 50 uL. Note that if you do this you won't be able to assess pool specific contamination._

#### End-repair using the Kapa Hyper prep kit

1. Add the following to each 50 uL sample:

  * 7 uL End Repair and A-tailing Buffer
  * 3 uL End Repair and A-tailing Enzyme mix

  _Note: it's easier just to make a master mix of these reagents and then only pipetting once. When making the mastermix budget ~10% more for pipetting error, and round up for fractions. E.g. if you have 24 samples, then you would budget 2.4 reactions for pipetting error. Round that up to 3 reactions, and then make the mastermix based off of calculations for 27 samples._

2. Run the following program on the thermocycler:

  * 20 degrees Celsius for 30 minutes
  * 65 degrees Celsius for 30 minutes
  * 4 degrees Celsius hold indefinitely

_While the samples are incubating on the thermocycler prepare the adaptor ligation mastermix. You can also stop here if you need to. Store samples at 4 degrees Celsius if you're not proceeding further with the protocol at this point._

#### Adaptor ligation and clean-up

1. Make sure that you have sufficient XT2 adaptors to have a unique adaptor for each sample.
2. Make up the following mastermix for the adaptor ligation. Volumes are per sample (budget extra mastermix for pipetting error as described above).
_Note: buffer can be brought up to temp quickly by hand warming. Keep ligase on ice. Once totally made up keep mastermix on ice._

  * 30 uL Ligation buffer (found in Kapa Hyper Prep Kit)
  * 10 uL DNA ligase (in Kapa kit)
  * 5 uL nuclease-free water

3. Once the thermocycler program has ended, take samples off the machine. Put on ice or in a cooled rack.
4. Add 45 uL of mastermix to each sample. Keep on ice.
5. Spin down adaptors. Open lids one at a time to prevent cross contamination.
6. Add 5 uL of unique adaptor to your sample. Be sure to keep track of which adaptor goes with which sample. Flick mix and spin down.
7. Incubate in a thermocycler at 20 degrees Celsius for 15 minutes. _Note you can hold at 4 degrees after this step. To save time aliquot out AMPure beads during this incubation._

> Directions for beads preparation. This is a 0.8:1 bead clean-up:

> 1. Vortex beads until homogenized.

> 2. Pipette 88 uL of beads into a Lo-Bind tube for each sample.

> 3. Allow the beads to come to room temperature. This happens pretty quickly when they are aliquotted out.

> 4. Label 0.2mL PCR tubes for the cleaned up products.

8. Add all 110 uL of ligation product to the tube with 88 uL of beads, pipetting up and down a couple of times to mix.
9. Incubate at room temperature for 10 minutes.
10. Place tubes on magnetic rack. Once the solution clears pipette off the supernatant without disturbing the beads.
11. Keeping tubes on the rack wash beads with 200 uL 80% EtOH. Resuspend the beads briefly by pipetting.
12. Discard the EtOH wash. Then repeat the wash a second time.
13. When removing the second wash remove as much EtOH as you can with the 200uL pipette, then go back with a 10uL pipette and remove any remaining EtOH.
14. Allow the beads to air dry for 5 minutes. _When you go to resuspend the beads make sure that you don't get any residual ethanol that is on the sidewalls of the tube mixed in with the water._
15. Resuspend the beads by pipetting in 25 uL of nuclease-free water. Return the tubes to the rack and allow the solution to clear.
16. Pipette off 20 uL of eluate placing in 0.2mL PCR tubes. _Ensure that you don't suck up any beads into the eluate that this point._

#### Library amplification

1. Make up the following mastermix for library amplification. Volumes are per sample. Budget a couple of extra reactions to account for pipetting error.

  * 25 uL 2X KAPA HiFi HotStart ReadyMix
  * 5 uL Illumina primer mix

2. Add 30 uL of mastermix to each adaptor-ligated library. Flick mix and spin down.
3. Run the following program on the thermocycler:

  i. 98 degrees Celsius for 45 seconds
  ii. 98 degrees Celsius for 15 seconds
  iii. 60 degrees Celsius for 30 seconds
  iv. 72 degrees Celsius for 30 seconds
   Repeat steps 2 through 4 for a **total** of 6 cycles
  v. 72 degrees Celsius for 1 minute
  vi. 4 degrees Celsius hold indefinitely

_you can stop here and store samples at 4 degrees Celsius or continue on with the clean-up_

_Note: amplification will only occur if adaptors were properly ligated on. So as an additional QC step you can quantify the library pre-pcr and post-pcr. If the post-pcr concentrations are higher then adaptor ligation of your library was successful._

4. Allow AMPure XP beads to come to room temperature and vortex until homogenized.
5. Add 40 uL of beads to each 50 uL samples (this is a 0.8:1 bead clean-up). Pipette to mix and incubate at room temperature for 10 minutes.
6. Place tubes on magnetic rack and wait until solution clears. Then pipette off and discard supernatant.
7. Wash pellet with 200 uL of 80% EtOH pipetting gently to resuspend beads, incubate for 30 seconds and discard wash.
8. Repeat 80% EtOH wash a second time. Remove as much EtOH as possible without disturbing beads by first removing wash with 200 uL pipette and then going through again with a 10 uL pipette.
9. Leave tubes on magnetic rack and allow to dry at room temperature for 5 minutes.
10. Remove tube from rack and add 25 uL of nuclease-free water to the beads. Resuspend well by pipetting and replace tube on the rack.
11. When solution has cleared pipette off 20 uL of supernatant without disturbing beads and place into new tubes.

#### Library quantification and pooling

_Before you sequence your library you should QC it on a TapeStation. You can QC the individual libraries or the pooled libraries. Economically it's better to QC the pooled library, and if it looks great then continue, if it doesn't QC well then it's a good idea to QC each of the libraries and determine whether the issue is with all of the samples or only with a few._

_At Fred Hutch we use the DNA 5000 assay, however check the size ranges of the assays available to you to choose the right one. Peak fragment size in your library should be ~ 580 nt (400 bp amplicons plus the ligated adaptors). If you have ~140 bp bands (adaptor dimers) perform the post amplification cleanup again. Fragments that are too small will interfere with cluster generation._

1. Quantify the DNA concentration using the Qubit High Sensitivity DNA kit (or equivalent) from 1 µL of each library.

_At least 0.76 ng/µL is required to achieve 2 nM for library pooling. Libraries will need to be concentrated or re-amplified if less than this amount. To ensure that your quantification is accurate include a sample with known concentration as a control._

3.	Convert DNA libraries from weight to moles:

	* Molecular weight [nM] = Library concentration [ng/µL] / ((avg. library size x 650)/1,000,000)

  * Example: if avg. size of library is 580 bp and concentration is 2.5 ng/µL…

    (580 x 650) / 1,000,000 = 0.377

  	2.5 / 0.377 = 6.6 nM

_The concentration that you are aiming will depend on how your core facility wants to receive samples or whether you will load the sample yourself. At Fred Hutch, the genomics core wants at least 30 uL of pooled library at a 6nM concentration. If you've done serial dilutions to get down to this concentration keep the more concentrated library available in case the core finds your library to be too dilute when they re-QC it prior to loading._

_Leftover pooled library can be stored at 4 degrees Celsius for short periods of time or at -20 degrees Celsius for longer periods._

_Once the MiSeq is running have a look at the cluster generation statistics and the quality score. You want the cluster generation to be around 1000-1200, but it's better to be a bit below that range than above that range. You want the quality score to be >80%._
