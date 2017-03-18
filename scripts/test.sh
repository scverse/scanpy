#!/usr/bin/env bash

set -x  # print commands

./scripts/scanpy.py exdata
./scripts/scanpy.py exparams
# moignard15, start with rereading
./scripts/scanpy.py moignard15 pca -ss 5 -s -r read -v 0
./scripts/scanpy.py moignard15 tsne -ss 5 -s -r -v 0
./scripts/scanpy.py moignard15 diffmap -ss 5 -s -r  -v 0
./scripts/scanpy.py moignard15 dpt -ss 5 -s -r -v 0
./scripts/scanpy.py moignard15 dpt -p basis=tsne -ss 5 -s -r -v 0
# burczynski06
./scripts/scanpy.py burczynski06 diffrank -s -r read -v 0
# paul15pca
./scripts/scanpy.py paul15pca dpt -ss 5 -s -r read -v 0
./scripts/scanpy.py paul15pca diffrank -o smp=dpt_groups names='2,3' -ss 5 -s -r -v 0
# toggleswitch
./scripts/scanpy.py toggleswitch sim -s -r read -v 0
./scripts/scanpy.py toggleswitch dpt -s -r -v 0