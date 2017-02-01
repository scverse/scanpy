## Tests

Call `./tests/test.sh` to test and `./tests/generate_example_figs.sh` to
reproduce the figures [here](https://github.com/theislab/scanpy) and there
[there](https://github.com/theislab/scanpy/examples).

## Wishlist

### Urgent

* data and result class that is flexible but allows for convenient way of
  defining subgroups etc. (-> the user should be able to define different
  subgroups, different colorings, and coloring by genes) / among others, mind
  that sparse data types have to be supported

* seaborn layout for difftest?

* merge plotting functions / add option for placing the legend outside

* check command line help and parameters: are they intuitive?

* bugs: Error 'QxcbConnection' on susi / warning raised when Lukas runs it on the server

* installation decription / script usage optimal?

* storing figures on http://falexwolf.de, optimal solution?

* add more examples / structure examples according to data type and technical
  aspects, such as runtime comparisons

### Less Urgent

* command line output: instead of setting default verbosity to 1, set it to 2 ->
  this provides lots of ouput and suggestions for new users, but allows to
  decrease it for more experienced users

* speed increase by intelligent importing

* preprocess PCA option for nearest-neighbor/graph-type tools? / do it by hand?
  compare for paul15 / include in examples

* speed-up computation of M matrix by using symmetries

* support for reading fcs cytof files (Maren)

* support for sparse data types

* integrate exisiting unit tests

* Lukas: frequent task: map dataset into other datasets / do an analysis of one
  of the large datasets

