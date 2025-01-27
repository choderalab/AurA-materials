# Analysis of EPR spin probes (PROJ 11431)

This section contains analysis scripts and processed datafiles for simulations with EPR spin probes (PROJ11431)

## Manifest

* `compute-distances.py` - compute various distances in parallel using multiprocessing
* `data/run*.npy` - numpy-format lists of distance arrays (in nm) generated by `compute-distances.py`; `distances[clone,index]` is distance array for `clone` in nanometers
  * index 0 : spin probe (CYR) NO group nitrogen-nitrogen (atom name NN) distances
  * index 1 : spin probe (CYR) alpha carbon distances
  * index 2 : R255 CZ - T288 CA distance
  * index 3 : F275 CZ - I193 CG2 distance (low in DFG-in)
  * index 4 : W277 CE2 - I193 CG2 distance (low in DFG-out)
  * index 5 : P282 O - R285 H distance (activation loop helical contact)
  * index 6 : S283 O - R286 H distance (activation loop helical contact)
  * index 7 : L225 CA - S284 CA distance (CYR-CYR here)
  * index 8 : Pseudotorsion for residues 282-285
  * index 9 : Pseudotorsion for residues 283-286

* `Analysis of spin probe distances.ipynb` - Jupyter notebook for analysis of distances

## Run index for project 11431
```
RUN INDEX

0: 1OL5-notpx2-nophos
1: 1OL5-notpx2-phos
2: 1OL5-tpx2-nophos
3: 1OL5-tpx2-phos
4: 1OL7-notpx2-nophos
5: 1OL7-notpx2-phos
6: 1OL7-tpx2-nophos
7: 1OL7-tpx2-phos
8: 5L8K-notpx2-nophos
9:  5L8K-notpx2-phos
10: 5L8K-tpx2-nophos
11: 5L8K-tpx2-phos
```
