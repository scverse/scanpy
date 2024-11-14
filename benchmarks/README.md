# Scanpy Benchmarks

This directory contains code for benchmarking Scanpy using [asv][].

The functionality is checked using the [`benchmark.yml`][] workflow.
Benchmarks are run using the [benchmark bot][].

[asv]: https://asv.readthedocs.io/
[`benchmark.yml`]: ../.github/workflows/benchmark.yml
[benchmark bot]: https://github.com/apps/scverse-benchmark

## Data processing in benchmarks

Each dataset is processed so it has

- `.layers['counts']` (containing data in C/row-major format) and `.layers['counts-off-axis']` (containing data in FORTRAN/column-major format)
- `.X` and `.layers['off-axis']` with log-transformed data (formats like above)
- a `.var['mt']` boolean column indicating mitochondrial genes

The benchmarks are set up so the `layer` parameter indicates the layer that will be moved into `.X` before the benchmark.
That way, we don’t need to add `layer=layer` everywhere.
