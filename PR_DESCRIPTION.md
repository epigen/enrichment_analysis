## Summary
This PR improves group summary plotting robustness and edge-case handling, while keeping the change scope focused to plotting logic.

## What Changed

### 1. `workflow/scripts/overview_plot.R`
- Added informative fallback plot generation instead of creating empty files.
- Added `make_message_plot()` and `make_no_sig_plot()` helpers.
- Empty aggregated input now produces message PNGs:
  - `No results found\naggregated result file is empty`
- Single-query groups now produce message PNGs:
  - `Only one query, so no group plot for <query>`
- `specificTerms` now produces a clear "No significant matches found" message plot when no terms pass threshold.

### 2. Summary plot readability improvements
- Improved term label handling for long names (truncate + wrap).
- Refined visual style (titles/subtitles, point/legend styling, grid colors).
- Increased minimum dynamic plot dimensions to avoid cramped layouts.

## Robustness Fixes (Edge Cases)

### Empty BED handling in GREAT paths
- Added early-exit guards for empty inputs in both GREAT scripts:
  - `workflow/scripts/region_enrichment_analysis_GREAT.R` (see lines 30-34)
  - `workflow/scripts/region_gene_association_GREAT.R` (see lines 32-38)
- Behavior:
  - If query/background BED is empty in enrichment, the rule writes an empty result CSV and exits cleanly.
  - If query BED is empty in region-gene association, the rule creates empty `genes.txt`, `region_gene_associations.csv`, and `region_gene_associations.pdf`, then exits cleanly.
- This avoids downstream failures and keeps output contracts intact.

### GREAT region-gene association coordinate fix (BED vs GREAT indexing)
- Adjusted exported `start` coordinate in `workflow/scripts/region_gene_association_GREAT.R` by `-1` before writing `region_gene_associations.csv`.
- Code reference: `workflow/scripts/region_gene_association_GREAT.R` line 74.
- Reason:
  - BED coordinates are 0-based.
  - GREAT association coordinates are 1-based.
- This makes exported association coordinates directly comparable to BED inputs and resolves the off-by-one indexing issue in reported GREAT region-gene associations.

## Behavioral Impact
- Summary outputs are now always interpretable:
  - either a proper `topTerms/specificTerms` plot,
  - or a clear explanatory message PNG.
- No change to the intended ranking logic:
  - `topTerms`: top-ranked per query, unioned across queries.
  - `specificTerms`: significant-first, then prioritized by weaker signal in other queries.

## Validation
- Re-ran aggregated `visualize` jobs and regenerated summary outputs.
- Tested top-`n` aggregated summary plots (`summary_topTerms`) across all configured tools/software:
  - GREAT
  - LOLA
  - ORA_GSEApy
  - preranked_GSEApy
  - pycisTarget
  - RcisTarget
- Verified single-query case now renders informative message PNGs.
- Verified `RcisTarget topTerms` selection still matches expected logic (`NES` descending per query, union across queries).

## Scope
- Kept config/test file modifications out of the final patch.
- Effective code change is focused on:
  - `workflow/scripts/overview_plot.R`
