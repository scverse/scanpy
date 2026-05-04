```{include} ../README.md
:end-before: '## Citation'
```

::::{grid} 1 2 3 3
:gutter: 2

:::{grid-item-card} Installation {octicon}`plug;1em;`
:link: installation
:link-type: doc

New to *scanpy*? Check out the installation guide.
:::

:::{grid-item-card} Tutorials {octicon}`play;1em;`
:link: tutorials/index
:link-type: doc

The tutorials walk you through real-world applications of scanpy.
:::

:::{grid-item-card} API reference {octicon}`book;1em;`
:link: api/index
:link-type: doc

The API reference contains a detailed description of
the scanpy API.
:::

:::{grid-item-card} Discussion {octicon}`megaphone;1em;`
:link: https://discourse.scverse.org

Need help? Reach out on our forum to get your questions answered!
:::

:::{grid-item-card} GitHub {octicon}`mark-github;1em;`
:link: https://github.com/scverse/scanpy

Find a bug? Interested in improving scanpy? Checkout our GitHub for the latest developments.
:::

:::{grid-item-card} GPU acceleration {octicon}`rocket;1em;`
:link: https://rapids-singlecell.readthedocs.io

Working with millions of cells and have a GPU available?
[rapids-singlecell][rsc] mirrors the scanpy API with order-of-magnitude speedups on core preprocessing, neighbors, embedding, and clustering steps.
:::
::::

[rsc]: https://rapids-singlecell.readthedocs.io

**Other resources**

* Follow changes in the {ref}`release notes <release-notes>`.
* Discover tools that build on or complement scanpy in the [scverse ecosystem](https://scverse.org/packages/#ecosystem), and follow scverse-wide news at [scverse.org/blog](https://scverse.org/blog/).
* Check out our {ref}`contribution guide <contribution-guide>` for development practices.
* Consider citing [Genome Biology (2018)] along with original {doc}`references <references>`.

% put references first so all references are resolved

% NO! there is a particular meaning to this sequence

```{toctree}
:hidden: true
:maxdepth: 1

installation
tutorials/index
usage-principles
how-to/index
api/index
external/index
release-notes/index
community
dev/index
contributors
references
```

[contribution guide]: dev/index.md
[genome biology (2018)]: https://doi.org/10.1186/s13059-017-1382-0
[github]: https://github.com/scverse/scanpy
