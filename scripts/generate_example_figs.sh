#!/usr/bin/env bash

set -x  # print commands

./scripts/scanpy.py moignard15 pca --figdir examples/figs/ --savefigs
./scripts/scanpy.py moignard15 tsne --figdir examples/figs/ --savefigs
./scripts/scanpy.py moignard15 diffmap --figdir examples/figs/ --savefigs --plotparams legendloc "upper left"
./scripts/scanpy.py moignard15 dpt --figdir examples/figs/ --savefigs
./scripts/scanpy.py moignard15 dpt tsne --figdir examples/figs/ --savefigs
./scripts/scanpy.py paul15 dpt --figdir examples/figs/ --savefigs
./scripts/scanpy.py paul15 difftest --prev dpt --figdir examples/figs/ --savefigs
./scripts/scanpy.py paul15pca dpt --savefigs --figdir examples/figs/
./scripts/scanpy.py paul15pca difftest --prev dpt --savefigs --figdir examples/figs/
./scripts/scanpy.py krumsiek11 sim --figdir examples/figs/ --savefigs
./scripts/scanpy.py krumsiek11 dpt --figdir examples/figs/ --savefigs --plotparams layout 3d legendloc "upper left"
