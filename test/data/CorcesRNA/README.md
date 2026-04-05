This fixture contains a minimal, biologically interpretable subset of the MrBiomics `CorcesRNA` RNA-seq recipe outputs.

Included query sets:
- `Bcell`: strong B-cell identity genes and ranked scores
- `Ery`: strong erythroid identity genes and ranked scores
- `Mono`: strong monocyte identity genes and ranked scores

Why these three:
- They represent clearly distinct hematopoietic lineages.
- They produce inspectable enrichment outputs instead of ambiguous edge cases.
- Three queries are enough to exercise group-level aggregation and summary plots without making the fixture large.

Included files:
- `*_up_features_annot.txt`: gene-set inputs for `ORA_GSEApy` and `RcisTarget`
- `*_featureScores_annot.csv`: ranked gene-score inputs for `preranked_GSEApy`
- `ALL_features_annot.txt`: shared background universe for the gene-set tests

Source:
- copied from `results/CorcesRNA/dea_limma/normCQN_OvA_cell_type/feature_lists/`
