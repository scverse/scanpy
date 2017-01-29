#!/usr/bin/env bash

set -x  # print commands

python scripts/scanpy.py exdata -v 0
python scripts/scanpy.py examples -v 0
python scripts/scanpy.py pca moignard15 --subsample 5 --savefigs --recompute -v 0
python scripts/scanpy.py tsne moignard15 --subsample 5 --savefigs --recompute -v 0
python scripts/scanpy.py diffmap moignard15 --subsample 5 --savefigs --recompute -v 0
python scripts/scanpy.py dpt moignard15 --subsample 5 --savefigs --recompute -v 0
python scripts/scanpy.py dpt krumsiek11 --savefigs --recompute -v 0
python scripts/scanpy.py dpt paul15 --savefigs --recompute -v 0
python scripts/scanpy.py difftest paul15 --savefigs --recompute -v 0
python scripts/scanpy.py sim toggleswitch --savefigs --recompute -v 0
