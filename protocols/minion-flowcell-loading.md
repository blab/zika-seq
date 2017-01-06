# Prime MinION flowcell and load the library

### _with lots of notes drawn from experience making mistakes_

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
