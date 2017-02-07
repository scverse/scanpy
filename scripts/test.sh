#!/usr/bin/env bash

set -x  # print commands

./scanpy.py exdata 
./scanpy.py examples 
./scanpy.py moignard15 pca --s 5 --savefigs -r all -v 0
./scanpy.py moignard15 tsne --s 5 --savefigs -r all -v 0
./scanpy.py moignard15 diffmap --s 5 --savefigs -r all -v 0
./scanpy.py moignard15 dpt --s 5 --savefigs -r all -v 0
./scanpy.py moignard15 dpt tsne --s 5 --savefigs -r all -v 0
./scanpy.py paul15pca dpt --savefigs -r all -v 0
./scanpy.py paul15pca difftest --prev dpt --savefigs -r all -v 0
./scanpy.py toggleswitch sim --savefigs -r all -v 0
./scanpy.py toggleswitch dpt --savefigs -r all -v 0