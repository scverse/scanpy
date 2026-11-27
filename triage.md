## Important: needs your attention

- **#4499**: multi-graph Leiden (from maltekuehl, discussed with ilan-gold). niklasmueboe points out squidpy already has this (scverse/squidpy#1019). You need to decide whether it belongs in scanpy. It would go in 1.14 at the earliest.

## Everything else

- **#4526**: a real bug, worth fixing. `sc.get.aggregate` fails when multi-column group labels give the same underscore-joined name (`("a_b","c")` and `("a","b_c")`). Label it a bug and add it to 1.14.0, or 1.13.1 if the fix is cheap.
- **#4541**: feature request for a Seurat `DimHeatmap` equivalent, asked for twice before (#505, #904). The author is ready to open a PR and only needs to know which plotting backend you want. Answer "v2 only" (or whichever you prefer) to unblock it.
- **#4525, #4503, #4502**: your own trackers. They have no milestones or labels (#4525 has an area label), but there's nothing to triage.

I only read the diffs for the fix PRs. For #4512, #4499 and #4493 I went by the descriptions and the size of the diffs.
