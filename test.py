import pooch
import tarfile
import os


input_url = "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"


filepath = pooch.retrieve(url=input_url, known_hash=None)

print("File path: ", filepath)
tar = tarfile.open(filepath, "r:gz")
tar.extractall(path=os.path.join( os.path.dirname( filepath ), '../..' ))
tar.close()
