#!/usr/bin/env bash

set -xe  # print commands and exit on error

python scripts/scanpy.py moignard15 pca --figdir examples/figs/ --savefigs png
python scripts/scanpy.py moignard15 tsne --figdir examples/figs/ --savefigs png
python scripts/scanpy.py moignard15 diffmap --figdir examples/figs/ --savefigs png --plotparams legendloc "upper left"
python scripts/scanpy.py moignard15 dpt --figdir examples/figs/ --savefigs png
python scripts/scanpy.py paul15 dpt --figdir examples/figs/ --savefigs png
python scripts/scanpy.py paul15 difftest --figdir examples/figs/ --savefigs png
python scripts/scanpy.py krumsiek11 sim --figdir examples/figs/ --savefigs png
python scripts/scanpy.py krumsiek11 dpt --figdir examples/figs/ --savefigs png --plotparams layout 3d
