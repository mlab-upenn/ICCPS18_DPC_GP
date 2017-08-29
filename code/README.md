# Important #

We use the following convention of time steps in the model:

```
y(k) = GP(u(k-n), u(k-n+1), \dots, u(k-1))
```

Note that the inputs are up to `k-1`, not `k`.

# Experiment Design #

## Input ##

## Output ##

# Training #

Two folders involved:
* `training` contains code to train the models from the data of the previous step.
* `models` contains the resulting models.

## Input ##

Data files from previous step.

## Output ##

Model files from training.

The GP models are stored in the directory `models`.
Each model is a MAT file of the following variables:
* `normalization` : a structure of normalization parameters for different variables, each field is a structure of two fields `max` and `min`. Example: `normalization.y` contains `max` and `min` for the output normalization (power).
* `stepsahead` is the number of steps of looking-ahead: 0 if the predicted output is for the current step, 1 if the predicted output is for the next step, and so on.
* `lags` is a structure that contains the lag values for different input variables, **counting backward from the stepshead horizon**. Example: `lags.u = [3 1]` with `stepsahead = 0` results in `[u(k-3), u(k-1)]` as the input; with `stepsahead = 1`, it becomes `[u(k-2), u(k)]`.
* `hyp` is the hyperparameter structure of GPML.
* `cov` is the covariance structure of GPML.
* `mean` is the mean function structure of GPML.
* `lik` is the likelihood structure of GPML.
* `X` is the matrix of row vectors of training inputs.
* `Y` is the matrix of row vectors of training outputs.
* `hypnames` is a cell array of the names of the covariance hyperparameters.
* `flogtheta` is the training results as returned by GPML.
* `validation` is the validation results, a structure of:
    * `X` are the test inputs
    * `targets` are the test outputs (unnormalized)
    * `Y` are the predicted outputs (unnormalized)
    * `S2` are the predicted variances (unnormalized)
    * `loss` is a structure of loss values (AE, SE, LPD, RMSE, SMSE, MSLL).

## How to use ##

# Control #

Folders involved:
* `control` contains the control code to run the simulation in closed-loop with OpenBuildNet.
* `controlresults` contains the results of the control simulation.

