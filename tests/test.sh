#!/usr/bin/env bash

set -x  # print commands

python scripts/scanpy.py exdata 
python scripts/scanpy.py examples 
python scripts/scanpy.py pca moignard15 --subsample 5 --savefigs --recompute -v 0
python scripts/scanpy.py tsne moignard15 --subsample 5 --savefigs --recompute -v 0
python scripts/scanpy.py diffmap moignard15 --subsample 5 --savefigs --recompute -v 0
python scripts/scanpy.py dpt moignard15 --subsample 5 --savefigs --recompute -v 0
python scripts/scanpy.py dpt paul15 --savefigs --recompute -v 0
python scripts/scanpy.py difftest paul15 --savefigs --recompute -v 0
python scripts/scanpy.py toggleswitch sim --savefigs --recompute -v 0
python scripts/scanpy.py toggleswitch dpt --savefigs --recompute -v 0