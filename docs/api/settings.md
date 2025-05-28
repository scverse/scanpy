(settings)=

## Settings


```{eval-rst}
.. currentmodule:: scanpy
```

A convenience function for setting some default {obj}`matplotlib.rcParams` and a
high-resolution jupyter display backend useful for use in notebooks.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   set_figure_params
```

An object that allows configuring Scanpy.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   settings
```

Some selected settings are discussed in the following.

Influence the global behavior of plotting functions. In non-interactive scripts,
you'd usually want to set `settings.autoshow` to `False`.

% no :toctree: here because they are linked under the class

```{eval-rst}
.. autosummary::
   :nosignatures:

   ~settings.autoshow
   ~settings.autosave
```

IO related settings for saving figures, caching files and storing datasets.


```{eval-rst}
.. autosummary::
   :nosignatures:

   ~settings.figdir
   ~settings.cachedir
   ~settings.datasetdir
   ~settings.file_format_figs
   ~settings.file_format_data
```

The verbosity of logging output, where verbosity levels have the following
meaning: 0='error', 1='warning', 2='info', 3='hint', 4=more details, 5=even more
details, etc.

```{eval-rst}
.. autosummary::
   :nosignatures:

   ~settings.verbosity
```

Print versions of packages that might influence numerical results.

```{eval-rst}
.. autosummary::
   :nosignatures:
   :toctree: ../generated/

   logging.print_header
   logging.print_versions

```
