Python powerloss model of 2or3 level buck or 'buckboost' topology

The model is designed to work for these configurations:
-2lvl and 3lvl buck
-Q4 FET for conduction loss only (bb and boost mode not yet implemented).
-single, parallel or series cyntec inductor
-single or parallel Q1 and Q2 FETs


The hs and ls mosfet models are based on the paper found in the 'references' folder.

This model utilizes estimated  package inductances along with datasheet derived mosfet non-linear capacitances to model
the turn-on and turn-off switching vgate, ids and vds waveforms.  Switching losses are calculated from these waveforms.


