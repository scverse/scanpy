```{include} ../README.md
:end-line: 17
```

```{eval-rst}
.. include:: _key_contributors.rst
```

```{eval-rst}
.. role:: small
```

```{eval-rst}
.. role:: smaller
```

* Discuss usage on [Discourse] and development on [GitHub].
* Get started by browsing {doc}`tutorials <tutorials>`, {doc}`usage principles <usage-principles>` or the main {doc}`API <api>`.
* Follow changes in the {ref}`release notes <release-notes>`.
* Find tools that harmonize well with anndata & Scanpy via the {doc}`external API <external>` and the {doc}`ecosystem page <ecosystem>`.
* Check out our {ref}`contributing guide <contribution-guide>` for development practices.
* Consider citing [Genome Biology (2018)] along with original {doc}`references <references>`.

# News

```{include} news.md
:start-line: 9
:end-line: 32
```

{ref}`(past news) <News>`
# Latest additions

```{include} release-notes/release-latest.md
```

% put references first so all references are resolved

% NO! there is a particular meaning to this sequence

```{toctree}
:hidden: true
:maxdepth: 1

tutorials
usage-principles
installation
api
external
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
