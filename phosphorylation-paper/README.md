## Manifest

* `AMBER-wt-analysis` - contains all of the analysis scripts, processed data and figure generating ipython notebooks for the AMBER WT data presented in the paper 
* `spin-probe-analysis` - contains the scripts, processed data and figure generation ipython notebooks for the CHARMM MTSL labeled AurA simulation data 

## Folding@Home data manifest 
Projects analyzed for this manuscript:
* 11414: AurA noPhos Tpx2, pdb 1OL5 (AMBER) 
  * RUN0: WT
* 11419: AurA noPhos TPX2, pdb 1OL5 (AMBER) 
  * RUN0: WT
  * RUN1: WT
  * RUN2: WT
  * RUN3: WT
* 11418: AurA noPhos noTPX2, pdb 1OL5 (AMBER) 
  * RUN0: WT
  * RUN1: WT
  * RUN2: WT
  * RUN3: WT
  * RUN4: WT
* 11428: AurA Phos TPX2, pdb 1OL5 (AMBER) 
  * RUN0: WT
* 11429: AurA Phos noTPX2, pdb 1OL5 (AMBER) 
  * RUN0: WT
* 11431: Spin probe labeled AURKA (CHARMM forcfield) -  pdbs: 1OL7-tpx2-new_03082017.pdb, 5L8K-tpx2-new_03072017.pdb, 1ol5-prepped.pdb
(Note that each source PDB starts with a different residue number)
  * RUN0 1OL5-notpx2-nophos
  * RUN1 1OL5-notpx2-phos
  * RUN2 1OL5-tpx2-nophos
  * RUN3 1OL5-tpx2-phos
  * RUN4 1OL7-notpx2-nophos
  * RUN5 1OL7-notpx2-phos
  * RUN6 1OL7-tpx2-nophos
  * RUN7 1OL7-tpx2-phos
  * RUN8 5L8K-notpx2-nophos
  * RUN9  5L8K-notpx2-phos
  * RUN10 5L8K-tpx2-nophos
  * RUN11 5L8K-tpx2-phos
* 11432: WT CHARMM AURKA, all from pdb 1OL5 
  * RUN0 TPX2 Phos 
  * RUN1 TPX2 noPhos
  * RUN2 noTPX2 Phos
  * RUN3 noTPX2 noPhos
