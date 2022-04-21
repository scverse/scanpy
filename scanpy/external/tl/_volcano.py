from typing import Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def volcano_plot(
  logfoldchange : DataFrame, 
  p_vals : DataFrame, 
  gene_ids: DataFrame,
  logfold_minthres: Optional[float] = -0.6,
  logfold_maxthres: Optional[float] = 0.6,
  pval_thres: Optional[float] = 0.05,
  figsize: Optional[Tuple[float, float]] = (6.5,6),
  title: Optional['str'] = None,
  out_dir: Optional['str'] = "volcanoplot.png"
  ):
  """
  Plot ranking of genes using volano plot

  Visualizes and identifies statistifically significant gene expression changes from two different 
  experimental conditions in terms of log fold change and p value

  This requires having run :func:‘~scanpy.get.rank_gene_groups_df’

  Parameters
  ----------
  logfoldchanges
    df containing logfoldchanges to plot
  p_vals
    df containing p-values --> should we adjust before or after passing into function?
  logfold_minthres
      Value to display as minimum log2fold threshold.
  logfold_maxthres
      Value to display as maximum log2fold threshold.
  pval_thres
      Value to display as p-value threshold.
  figsize
      Output size of figure (width, height).
  title
      Title to display as caption of figure.
  output_file
      Name of output plot image (include ".png" extension)

  Returns
  -------
  Returns volcano plot highlighting which genes are more expressed (or enriched) in one group

  Notes
  ------

  Example
  -------
  >> volcano_plot(logfoldchanges, p_vals)
  """
  # default current dir, add user input option

  p_vals[p_vals == -np.inf] = np.min(p_vals[p_vals != -np.inf])
  
  fig = plt.figure(figsize = (8, 6))
  ax = fig.add_subplot(1, 1, 1)

  # find genes with significant expression
  id_sig = np.where((p_vals < -3.5) & (np.abs(logfoldchange) > 0.95))[0]
  # id_sig2 = np.where((log_p_vals < -3.5) & (np.abs(log_fold_change) > 0.95))[0]


  # 
  ax.scatter(logfoldchange, -p_vals, s = 5, c = 'lightgray', label = 'Not Significant')
  ax.scatter(logfoldchange[id_sig], -p_vals[id_sig], s = 5, c = 'r', label = 'Significant')

  # label significant points
  for i in id_sig:
    ax.text(logfoldchange[i], -p_vals[i], gene_ids[i])

  ax.legend(fontsize = 12, markerscale = 5)
  ax.set_xlabel('Log fold change coefficient', fontsize = 16)
  ax.set_ylabel('Adjusted p-values', fontsize = 16)

  ax.axhline(y = pval_thres, color = 'g') 
  ax.axvline(x = logfold_minthres, color = 'b')
  ax.axvline(x = logfold_maxthres, color = 'b')
  plt.savefig('volcanoplot.png', out_dir)
