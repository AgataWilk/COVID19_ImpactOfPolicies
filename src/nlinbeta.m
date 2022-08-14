function beta = nlinbeta(x, restrictions, bCountries, modelling)
% x - a vector of parameters of the model
% restrictions - matrix corresponding to odpowiadająca równaniom na obostrzenia - [1, -o1, -o2, ... -on]
% bCountries - matrix of zeros and ones to generate a vector of b parameters
% For common and independent apptroaches the function is basically the same
% except in the common model the computation for all countries is done at
% the same time and for the independent models each country is computed
% separately in a loop
nRestrictions = 13;
nCountries = 42;
if modelling == "individualized"
    beta = restrictions*[1;reshape(x(1:nRestrictions),[],1)].*(bCountries*...
        [reshape(x((nRestrictions+1):nRestrictions+nCountries),[],1)]);
else 
    beta = restrictions*[1;reshape(x(1:nRestrictions),[],1)]*x(nRestrictions+1);
end