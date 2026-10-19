Predicted structures for the 20 submitted Track 3 designs
=========================================================
Named by the identifier in submission.csv (tnf3_01 ... tnf3_20).

complex/<name>_complex.pdb   chain X = designed binder; chains A,B = TNF-alpha
                             protomers (157 aa, assay construct TNA-H4211)
binder_only/<name>_binder.pdb   chain X extracted
structures_index.tsv         ddG, binder-target pair ipTM, complex pLDDT,
                             binder pTM, His substitutions, cluster, length

These are the stage-06 re-folds: each design AFTER its histidine switch
substitutions were installed, which is what the reported pair ipTM refers to.

Two caveats that matter for reading them
----------------------------------------
1. The target chains were PINNED TO THE CRYSTAL via template (1TNF A+B).
   Without it, single-sequence Boltz folds each protomer well (chain pTM ~0.85)
   but does not dock them (target-internal pair ipTM 0.291), leaving the
   receptor groove -- which sits at the protomer-protomer interface -- absent
   from the prediction. With the template it is ~0.71. Target coordinates are
   therefore close to crystal geometry by construction, not predicted.
2. Confidence was scored on BINDER-TARGET chain pairs only, taking the minimum
   of the two. Global ipTM on a two-protomer crop is dominated by the
   target-internal pair and is not meaningful here.

Predictions are single-sequence (no MSA): the binders are de novo and have no
meaningful MSA, and the target's would be identical across every design.
Per-design pae/pde/plddt .npz arrays exist in the source tree and are available
on request.
