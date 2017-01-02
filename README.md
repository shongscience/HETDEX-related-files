# HETDEX-related-files
Codes related to the HETDEX survey

1. the "falltest" folder is for tracking the progress of the current fall field observations. A trivial python script. 

2. the "optimizingplan" folder is for investigating how to minimize the effect of time-dependent sampling patterns on power spectrum measurement. A relatively long serious cpp program. 
  
  2.1 Problem Description: 
  
  The Hetdex survey consists of 4000+ shots to cover a large amount of sky for detecting BAOs from spectroscopic LAEs. If we have N shots, theoretically we have N!(factorial) kinds of combinations of "temportal sequence" to fill out the scheduled shots; "temporal observing patterns". 
  
  Some terrestrial weather patterns leave their traits as "short-term unexpected correlated errors" on powerspectrum measurements and lunar phases also as "periodic correlated errors". 
  
  2.2 Ruler Model:
  
  This program estimates these two errors for each temporal filling pattern using my "ruler model". You can contact me for the details about this rule model. 
