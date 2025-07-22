(settings)=

## Settings


```{eval-rst}
.. currentmodule:: scanpy
```

A convenience function for setting some default {obj}`matplotlib.rcParams` and a
high-resolution jupyter display backend useful for use in notebooks.

```{eval-rst}
.. autosummary::
   :signatures: none
   :toctree: ../generated/

   set_figure_params
```

An object that allows configuring Scanpy.

```{eval-rst}
.. autosummary::
   :signatures: none
   :toctree: ../generated/

   settings
```

Some selected settings are discussed in the following.

Presets allow to set the behavior of many scanpy functions at once:

```{eval-rst}
.. autosummary::
   :signatures: none
   :toctree: ../generated/

   Preset
```

Verbosity controls the amount of logging output:

```{eval-rst}
.. autosummary::
   :signatures: none
   :toctree: ../generated/

   Verbosity
```

Influence the global behavior of plotting functions. In non-interactive scripts,
you'd usually want to set {attr}`settings.autoshow` to `False`.

```{eval-rst}
.. autosummary::

   settings.autoshow
   settings.autosave
```

IO related settings for saving figures, caching files and storing datasets.


```{eval-rst}
.. autosummary::

   settings.figdir
   settings.cachedir
   settings.datasetdir
   settings.file_format_figs
   settings.file_format_data
```

Print versions of packages that might influence numerical results.

```{eval-rst}
.. autosummary::
   :signatures: none
   :toctree: ../generated/

   logging.print_header
```
