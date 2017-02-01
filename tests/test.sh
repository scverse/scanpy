#!/usr/bin/env bash

set -x  # print commands

python scripts/scanpy.py exdata 
python scripts/scanpy.py examples 
python scripts/scanpy.py moignard15 pca --subsample 5 --savefigs --recompute -v 0
python scripts/scanpy.py moignard15 tsne --subsample 5 --savefigs --recompute -v 0
python scripts/scanpy.py moignard15 diffmap --subsample 5 --savefigs --recompute -v 0
python scripts/scanpy.py moignard15 dpt --subsample 5 --savefigs --recompute -v 0
python scripts/scanpy.py moignard15 dpt tsne --subsample 5 --savefigs --recompute -v 0
python scripts/scanpy.py paul15 dpt --savefigs --recompute -v 0
python scripts/scanpy.py paul15 difftest --prev dpt --savefigs --recompute -v 0
python scripts/scanpy.py sim toggleswitch --savefigs --recompute -v 0
python scripts/scanpy.py dpt toggleswitch --savefigs --recompute -v 0