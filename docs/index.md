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
::::

**Other resources**

* Follow changes in the {ref}`release notes <release-notes>`.
* Find tools that harmonize well with anndata & Scanpy at [scverse.org/packages/](https://scverse.org/packages/)
* Check out our {ref}`contribution guide <contribution-guide>` for development practices.
* Consider citing [Genome Biology (2018)] along with original {doc}`references <references>`.

## News

```{include} news.md
:start-after: '<!-- marker: after prelude -->'
:end-before: '<!-- marker: before old news -->'
```

{ref}`(past news) <News>`

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
ecosystem
release-notes/index
community
news
dev/index
contributors
references
```

[contribution guide]: dev/index.md
[genome biology (2018)]: https://doi.org/10.1186/s13059-017-1382-0
[github]: https://github.com/scverse/scanpy
