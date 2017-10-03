# LargeOffice

## Training data ##

`train` **data from June**

* data generated using `generateTrainData.m`.
* `train-nominal-` contains data with all inputs following the nominal rule-based strategy.
* `train-unconstrained-` contains data with `23<=ClgSP<=28`, `11<=SupplyAirSP<=15` and `3.7<=ChwSP<=9.7`. Other inputs follow the nominal rule-based strategy.
* `train-ramped2-` contains data with `23<=ClgSP<=28`, `11<=SupplyAirSP<=15` and `3.7<=ChwSP<=9.7` and `|ChwSP(k)-ChwSP(k-1)|<=2`. Other inputs follow the nominal rule-based strategy.

## Testing data ##

`test` **data from July**

* data generated using `generateTestData.m`.
* `test-nominal-` contains data with all inputs following the nominal rule-based strategy.
* `test-unconstrained-` contains data with `23<=ClgSP<=28`, `11<=SupplyAirSP<=15` and `3.7<=ChwSP<=9.7`. Other inputs follow the nominal rule-based strategy.
* `test-ramped2-` contains data with `23<=ClgSP<=28`, `11<=SupplyAirSP<=15` and `3.7<=ChwSP<=9.7` and `|ChwSP(k)-ChwSP(k-1)|<=2`. Other inputs follow the nominal rule-based strategy.

## Simulation data ##

`sim` **data from June-August**

* data generated using `generateSimData.m`.
* `sim-nominal-` contains data with all inputs following the nominal rule-based strategy.


# LargeHotel

## Initialization data ##

`init` **data from May**

* data generated using `generateInitData.m`.
* `init-` contains data with all inputs following the nominal rule-based strategy.

## Training data ##

`train` **data from June**

* data generated using `generateTrainData.m`.
* `unconstrained-` contains data with `22<=ClgSP<=32`, `24<=KitchenClgSP<=32`, `22<=GuestClgSP<=26`, `12<=SupplyAirSP<=14` and `3.7<=ChwSP<=9.7`. Other inputs follow the nominal rule-based strategy.
* `train-ramped2-` contains data with `22<=GuestClgSP<=26`, `12<=SupplyAirSP<=14`, `3.7<=ChwSP<=9.7` and `|ChwSP(k)-ChwSP(k-1)|<=2`. Other inputs follow the nominal rule-based strategy.

## Testing data ##

`test` **data from July**

* data generated using `generateTestData.m`.
* `test-constrained-` contains data with `23<=GuestClgSP<=25`, `12.5<=SupplyAirSP<=13.5` and `5.2<=ChwSP<=8.2`. Other inputs follow the nominal rule-based strategy.
* `test-unconstrained-` contains data with `22<=GuestClgSP<=26`, `12<=SupplyAirSP<=14` and `3.7<=ChwSP<=9.7`. Other inputs follow the nominal rule-based strategy.
* `test-ramped1-` contains data with `22<=GuestClgSP<=26`, `12<=SupplyAirSP<=14`, `3.7<=ChwSP<=9.7` and `|ChwSP(k)-ChwSP(k-1)|<=1`. Other inputs follow the nominal rule-based strategy.
* `test-ramped2-` contains data with `22<=GuestClgSP<=26`, `12<=SupplyAirSP<=14`, `3.7<=ChwSP<=9.7` and `|ChwSP(k)-ChwSP(k-1)|<=2`. Other inputs follow the nominal rule-based strategy.
* `test-` contains data with all inputs following the nominal rule-based strategy.

