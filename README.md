# COVID19_ImpactOfPolicies
Data and code for article "Impact of government policies on the COVID-19 pandemic unraveled by mathematical modelling"


## Data:

All data is in MATLAB format

populations_europe - populations of European countries used for modelling

pandemic_europe - confirmed cumulative and daily cases, deaths and recoveries

restrictions - normalized levels of government policies

beta_opt - values of $\beta_{opt}$ estimated using adjoint sensitivity analysis


## Code:

SEIR - implementation of the SEIR model

nlinfun, nlinbeta - helper functions used to generate matrices for estimation

main - parameter estimation for three modelling approaches, prediction of cases and error calculation
