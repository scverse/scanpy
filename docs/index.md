```{include} ../README.md
:end-line: 41
```

```{eval-rst}
.. role:: small
```

```{eval-rst}
.. role:: smaller
```

::::{grid} 1 2 3 3
:gutter: 2

:::{grid-item-card} Installation {octicon}`plug;1em;`
:link: installation
:link-type: doc

New to *scanpy*? Check out the installation guide.
:::

:::{grid-item-card} Tutorials {octicon}`play;1em;`
:link: _tutorials
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
* Find tools that harmonize well with anndata & Scanpy via the {doc}`external API <external/index>` and the {doc}`ecosystem page <ecosystem>`.
* Check out our {ref}`contributing guide <contribution-guide>` for development practices.
* Consider citing [Genome Biology (2018)] along with original {doc}`references <references>`.

## News

```{include} news.md
:start-line: 9
:end-line: 32
```

{ref}`(past news) <News>`

% put references first so all references are resolved

% NO! there is a particular meaning to this sequence

```{toctree}
:hidden: true
:maxdepth: 1

_tutorials
usage-principles
installation
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

[discourse]: https://discourse.scverse.org/
[genome biology (2018)]: https://doi.org/10.1186/s13059-017-1382-0
[github]: https://github.com/scverse/scanpy
