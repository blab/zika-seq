# MinION library preparation (1D chemistry)

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
5. Incubate at 65 degrees Celsius for 10 minutes to kill ligase. _Ensure that the temperature is correct. Samples will be pooled after this step so you need to prevent any further ligation that may cause sample-barcode mixtures._
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

1. Add 20 uL BAM to cleaned-up library. Mix by inversion.
3. Then add 50 uL Blunt/TA Ligase Master Mix. Mix by inversion and spin down.
4. Incubate at room temperature for 10 minutes.

#### Purify the adapted,tethered library

1. Add 40 uL of resuspended, room-temperature AMPure XP beads to the library.
2. Incubate at room temperature on a hula mixer for 5 minutes.
3. Place on magnet. Once solution clears discard the supernatant.
4. Add 140 uL of ABB to the pellet. Resuspend the pellet by flicking. Then spin down, replace on magnet, and pipette off the supernatant.
5. Repeat step 4 (perform another resuspension wash with 140 uL of ABB).
6. Spin down, replace tube on magnet and pipette off any residual ABB. Close lid, **do not let pellet dry**.
6. Resuspend the beads in 12 uL of ELB.
7. Incubate at room temperature for 10 minutes.
8. Pellet beads on magnetic rack and transfer eluate (tethered,adapted library in ELB) to new Lo-Bind tube.
9. Quantify adapted, tethered library using 1 uL of library to 199 uL Qubit Master Mix.

_Library can now be stored at -20 degrees Celsius if necessary, but you will have better results if you load the library and sequence at this point._

#### Adapted, tethered library quantification and normalization

Combine the following to make 75uL total of library and buffer mixture that you will load onto the flowcell.

* 35uL of RBF
* 3.5uL of nuclease free water
* 25.5 uL of LLB
* 11 uL of library (full amount from above after qubiting)
