import pandas as pd
"""
Test that volcano works
"""
volcano_df = pd.read_csv('.\_data\volcano_sample_df.csv')
sce.tl.volcano_plot(volcano_df, out_dir = '/_images/test_volcano.png')
assert compare_images('.\_images\test_volcano.png', '.\_images\master_volcano.png') is None
master_volcano_df = pd.read_csv('.\_data\master_volcano_df.csv')
sce.tl.volcano_plot(master_volcano_df, out_dir = '_images/test_volcano.png')
assert compare_images('_images/test_volcano.png', '_images/master_volcano.png', tol=5) is None

def test_volcano():
    """
    Test that volcano module option for volcano_plot() works
    """
    sce.tl.volcano('.\_data\volcano_data.volcano_data.rnk', hallmark_gene_sets_list='KEGG_2016', type='gseapy', out_dir='_data/test_gseapy_df.csv')
    master_volcano_sample_df = pd.read_csv('_data/master_voclano_df.csv')
    test_volcano_sample_df = pd.read_csv('_data/test_volcano_df.csv')
    assert test_volcano_sample_df.equals(master_volcano_sample_df) 
