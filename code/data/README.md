## Initialization data ##

`init` **data from May**

* data generated using `generateInitData.m`
* `init-` contains data with all inputs following the nominal rule-based strategy.

## Training data ##

`train` **data from June**

* data generated using `generateTrainData.m`
* `train-unconstrained` contains data with `22<=ClgSP<=32`, `24<=KitchenClgSP<=32`, `22<=GuestClgSP<=26`, `12<=SupplyAirSP<=14` and `3.7<=ChwSP<=9.7`. Other inputs follow the nominal rule-based strategy.

## Testing data ##

`test` **data from July**
* data generated using `generateTestDataRandom.m`
* `test-constrained` contains data with `23<=GuestClgSP<=25`, `12.5<=SupplyAirSP<=13.5` and `5.2<=ChwSP<=8.2`. Other inputs follow the nominal rule-based strategy.
* `test-unconstrained` contains data with `22<=GuestClgSP<=26`, `12<=SupplyAirSP<=14` and `3.7<=ChwSP<=9.7`. Other inputs follow the nominal rule-based strategy.
* `test-` contains data with all inputs following the nominal rule-based strategy.