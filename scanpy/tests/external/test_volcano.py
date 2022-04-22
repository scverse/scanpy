import pandas as pd
import numpy as np
import scanpy.external as sce
from matplotlib.testing.compare import compare_images
import os

"""
Test that volcano works
"""

test_file_path = './_images/test_volcano.png'
if os.path.exists(test_file_path):
    os.remove(test_file_path)

def test_volcano():
    """
    Test that volcano module option for volcano_plot() works
    """
    volcano_df = pd.read_csv('./_data/volcano_sample.csv')
    volcano_cleaned = volcano_df[-volcano_df['coef'].isnull()]
    log_fold_change = volcano_cleaned['coef'].values
    log_p_vals = np.log10(volcano_cleaned['fdr']).values
    log_p_vals[log_p_vals == -np.inf] = np.min(log_p_vals[log_p_vals != -np.inf]) 
    ids = volcano_cleaned['primerid'].values

    sce.tl.volcano_plot(log_fold_change, log_p_vals, ids, out_dir = './_images/test_volcano.png')

    # test reproducibility
    master_volcano_sample_df = pd.read_csv('_data/volcano_sample.csv')
    test_volcano_sample_df = pd.read_csv('_data/volcano_samplecsv')
    assert compare_images('./_images/test_volcano.png', './_images/master_volcano.png', tol=5) is None
    assert test_volcano_sample_df.equals(master_volcano_sample_df) 
