#!/usr/bin/env bash

set -x  # print commands

./scanpy.py exdata 
./scanpy.py examples 
./scanpy.py moignard15 pca --subsample 5 --savefigs -r all -v 0
./scanpy.py moignard15 tsne --subsample 5 --savefigs -r all -v 0
./scanpy.py moignard15 diffmap --subsample 5 --savefigs -r all -v 0
./scanpy.py moignard15 dpt --subsample 5 --savefigs -r all -v 0
./scanpy.py moignard15 dpt tsne --subsample 5 --savefigs -r all -v 0
./scanpy.py paul15 dpt --savefigs -r all -v 0
./scanpy.py paul15 difftest --prev dpt --savefigs -r all -v 0
./scanpy.py toggleswitch sim --savefigs -r all -v 0
./scanpy.py toggleswitch dpt --savefigs -r all -v 0