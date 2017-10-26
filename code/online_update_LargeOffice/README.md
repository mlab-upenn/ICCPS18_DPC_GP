## Code for optimal subset of data selection

Different `evolving_ig_*` files differ in three variables
* `n_days_init` - Number of days from June for initial data. This data is used to learn an initial GP model.
* `n_days_next` - Number of days from June for initial data. This is the new data generated online based on the month of operation. See below.
* `shift` - shift of 30 means next data is from July, shift of 30+31 means next data is from August and shift of 30+31+31 means next data is from September.

In the current setup, both `n_days_init` and `n_days_next` are generated randomly.