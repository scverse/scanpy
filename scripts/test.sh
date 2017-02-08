#!/usr/bin/env bash

set -x  # print commands

./scanpy.py exdata 
./scanpy.py examples 
./scanpy.py moignard15 pca -s 5 --savefigs -r read -v 0
./scanpy.py moignard15 tsne -s 5 --savefigs -r read -v 0
./scanpy.py moignard15 diffmap -s 5 --savefigs -r read -v 0
./scanpy.py moignard15 dpt -s 5 --savefigs -r read -v 0
./scanpy.py moignard15 dpt tsne -s 5 --savefigs -r read -v 0
./scanpy.py burczynski06 difftest --savefigs -r read -v 0
./scanpy.py paul15pca dpt -s 5 --savefigs -r read -v 0
./scanpy.py paul15pca difftest -s 5 --prev dpt -o names '2,3' --savefigs -r read -v 0
./scanpy.py toggleswitch sim --savefigs -r read -v 0
./scanpy.py toggleswitch dpt --savefigs -r read -v 0