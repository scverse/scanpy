#!/usr/bin/env bash

set -x  # print commands

./scanpy.py exdata 
./scanpy.py examples 
# moignard15, start with rereading
./scanpy.py moignard15 pca -s 5 --savefigs -r read -v 0
./scanpy.py moignard15 tsne -s 5 --savefigs -r -v 0
./scanpy.py moignard15 diffmap -s 5 --savefigs -r  -v 0
./scanpy.py moignard15 dpt -s 5 --savefigs -r -v 0
./scanpy.py moignard15 dpt -p basis tsne -s 5 --savefigs -r -v 0
# burczynski06
./scanpy.py burczynski06 difftest --savefigs -r read -v 0
# paul15pca
./scanpy.py paul15pca dpt -s 5 --savefigs -r read -v 0
./scanpy.py paul15pca difftest -o smp dpt_groups names '2,3' -s 5 --savefigs -r -v 0
# toggleswitch
./scanpy.py toggleswitch sim --savefigs -r read -v 0
./scanpy.py toggleswitch dpt --savefigs -r -v 0