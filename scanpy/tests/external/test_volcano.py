import pandas as pd
import numpy as np
import scanpy.external as sce
from matplotlib.testing.compare import compare_images

"""
Test that volcano works
"""
volcano_df = pd.read_csv('.\_data\volcano_sample.csv')
volcano_cleaned = volcano_df[-volcano_df['coef'].isnull()]

log_fold_change = volcano_cleaned['coef'].values
log_p_vals = np.log10(volcano_cleaned['fdr']).values
log_p_vals[log_p_vals == -np.inf] = np.min(log_p_vals[log_p_vals != -np.inf]) 
ids = volcano_cleaned['primerid'].values

sce.tl.volcano_plot(log_fold_change, log_p_vals, ids, out_dir = '/_images/test_volcano.png')
# mast_clean = mast_file[-mast_file['coef'].isnull()]




assert compare_images('.\_images\test_volcano.png', '.\_images\master_volcano.png') is None
master_volcano_df = pd.read_csv('.\_data\master_volcano_df.csv')
sce.tl.volcano_plot(master_volcano_df, out_dir = '_images/test_volcano.png')
assert compare_images('_images/test_volcano.png', '_images/master_volcano.png', tol=5) is None

def test_volcano():
    """
    Test that volcano module option for volcano_plot() works
    """
    sce.tl.volcano('.\_data\volcano_data.csv')
    master_volcano_sample_df = pd.read_csv('_data/master_voclano.csv')
    test_volcano_sample_df = pd.read_csv('_data/test_volcano.csv')
    assert test_volcano_sample_df.equals(master_volcano_sample_df) 
