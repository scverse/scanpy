{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. note::
   Both backends are accessible as ``{{ fullname }}``.
   The active backend is chosen by :attr:`scanpy.settings.preset`:
   the default is the legacy matplotlib backend;
   set it to :attr:`scanpy.Preset.ScanpyV2Preview` to use the HoloViews backend.

.. tab-set::

   .. tab-item:: New (HoloViews)
      :sync: scanpy-preset-v2

      .. autofunction:: scanpy.plotting._v2.{{ objname }}
         :no-index:

   .. tab-item:: Legacy (matplotlib)
      :sync: scanpy-preset-legacy

      .. autofunction:: {{ objname }}
