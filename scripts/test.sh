#!/usr/bin/env bash

set -x  # print commands

./scripts/scanpy.py exdata 
./scripts/scanpy.py examples 
./scripts/scanpy.py moignard15 pca --subsample 5 --savefigs --recompute -v 0
./scripts/scanpy.py moignard15 tsne --subsample 5 --savefigs --recompute -v 0
./scripts/scanpy.py moignard15 diffmap --subsample 5 --savefigs --recompute -v 0
./scripts/scanpy.py moignard15 dpt --subsample 5 --savefigs --recompute -v 0
./scripts/scanpy.py moignard15 dpt tsne --subsample 5 --savefigs --recompute -v 0
./scripts/scanpy.py paul15 dpt --savefigs --recompute -v 0
./scripts/scanpy.py paul15 difftest --prev dpt --savefigs --recompute -v 0
./scripts/scanpy.py toggleswitch sim --savefigs --recompute -v 0
./scripts/scanpy.py toggleswitch dpt --savefigs --recompute -v 0