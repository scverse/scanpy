## Wish List

### Urgent

* fix installation of scanpy / import anywhere on the system

* seaborn layout for difftest?

* speedup by *sparse* importing? seems mainly important for displaying help?

* Error 'QxcbConnection' on susi / warning raised when Lukas runs scanpy on the server

* storing figure files on http://falexwolf.de is good solution?

* make sure drawg converges for a given set of parameters!

### Less Urgent

* recheck command line help and parameters: are they intuitive?

* check sc.help() versus help()

* command line output: instead of setting default verbosity to 1, set it to 2 ->
  this provides lots of ouput and suggestions for new users, but allows to
  decrease it for more experienced users / timing information important for users? / 
  tipps marked with "-->" are very important for new users but not for experienced ones / 
  should all be managed via verbosity, maybe config file that sets the settings in 
  settings.py?

* support for reading fcs cytof files (Maren)

* support for sparse data types

* integrate exisiting unit tests

* Lukas: frequent task: map dataset into other datasets / do an analysis of one
  of the large datasets

* is there a way to get a list of available genes to choose from on the command
  line? for plotting?

### Solved Feb 3, 2017 at 5:26

* ~~data and result class that is flexible but allows for convenient way of
  defining subgroups etc. (-> the user should be able to define different
  subgroups, different colorings, and coloring by genes) / among others, mind
  that sparse data types have to be supported [one idea would be to allow
  for a class `rowgroups`, and a class `colgroups` as possible attributes of data and 
  result classes. class `rowgroups` allows storing names, ids, masks, colors. make
  sure we do not reinvent the wheel: can pandas help us?~~ -> hack with dict

* ~~merge plotting functions / add option for placing the legend outside~~ -> at least
  for visualization tools done

### Solved Feb 2, 2017 at 09:58

* ~~add more examples / structure examples according to data type and technical
  aspects, such as runtime comparisons~~ -> ok

* ~~speed-up through writing of data dictionary / option recompute more powerful~~

### Solved Feb 1, 2017 at 19:06

* ~~speed-up computation of M matrix in DPT by using symmetries~~ -> the latter
  didn't work, but a huge speed-up was possible by performing PCA

* ~~consistent way of performing preprocessing: preprocessing functions take X or
  data class as input? / way to list different preprocessing functions? -> 
  integrate discussion in readme~~ -> two different pp modules `simple and advanced`

* ~~preprocess PCA option for nearest-neighbor/graph-type tools? / do it by hand? 
  compare for paul15 / include in examples~~ -> write `Xpca` during pp

* ~~naming conventions of parameters: stick to scikit-learn conventions,
  e.g. 'n_components'?~~ -> no! we use 'nr_comps' instead, all parameter names
  involving 'number' are called 'nr_'

