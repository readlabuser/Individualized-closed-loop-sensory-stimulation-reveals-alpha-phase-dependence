# Individualized-closed-loop-sensory-stimulation-reveals-alpha-phase-dependence

This repository accompanies the manuscript "Individualized closed-loop sensory stimulation reveals alpha phase-dependence of sound-driven brain activity in human electroencephalogram (EEG) recordings" (in preparation).  In this study, individualized alpha peak frequencies and P50 latencies were estimated during a resting-state recording session, and random-pahse presentation of pink noise acoustic stimuli with an ISI of ~900ms. From this individualized measures, parameters for the end-point corrected Hilbert Transform (ecHT) were selected for instantaneous estimation of indivudalized alpha phases and targetted presentation of pink noise acoustic stimuli at randomly itnerleaved peak and trough phases of ongoing alpha oscillations in real time during the third and fourth recording sessions.  The included code evaluates the physiological effects of peak and trough stimulation during these final recording sessions.

Code for this project can be broken into two segments: Utilities & Main. Utilities are functions that are necessary for processing the ecHT device output, statistical procedures, or aesthic processes like figure generation. Main are scripts that call upon data and the functions within Utilities in order to generate the analyses and figures from the manuscript. Additionally, all analyses assume user has the following toolboxes and packages:

  1. Circular Statistics (MATLAB)
  2. Chronux (MATLAB)
  3. Signal Processing (MATLAB)
  4. EEGLAB (MATLAB)
  5. fooof (python)
  6. fooof_mat (MATLAB)

In addition, the following two MATLAB functions from MathWorks are used when analyzing (fdr_bh) and plotting (al_goodplot):

  1. Credit to David Groppe for fdr_bh: https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh

  2. Credit to Antoine Legouhy for al_goodplot: https://www.mathworks.com/matlabcentral/fileexchange/91790-al_goodplot-boxblot-violin-plot
