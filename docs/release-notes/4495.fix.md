Fix {func}`scanpy.external.exporting.cellbrowser` builds with categorical colors and a relative `data_dir` when generating a new configuration.
Existing `cellbrowser.conf` files are preserved: correct affected color-file references or regenerate the configuration after saving any custom settings.
Generated color-file references are absolute and must be updated if the raw export directory is moved.
