function beta = nlinbeta(x, restrictions, bCountries, modelling)
% x - a vector of parameters of the model
% restrictions - matrix corresponding to odpowiadająca równaniom na obostrzenia - [1, -o1, -o2, ... -on]
% bCountries - matrix of zeros and ones to generate a vector of b parameters
% For common and independent apptroaches the function is basically the same
% except in the common model the computation for all countries is done at
% the same time and for the independent models each country is computed
% separately in a loop

if modelling == "individualized"
    beta = restrictions*[1;x(1);x(2);x(3);x(4);x(5);x(6);x(7);x(8);x(9);x(10);...
        x(11);x(12);x(13)].*(bCountries*[x(14);x(15);x(16);x(17);x(18);x(19);x(20);x(21);x(22);x(23);...
       x(24);x(25);x(26);x(27);x(28);x(29);x(30);x(31);x(32);x(33);...
       x(34);x(35);x(36);x(37);x(38);x(39);x(40);x(41);x(42);x(43);...
       x(44);x(45);x(46);x(47);x(48);x(49);x(50);x(51);x(52);x(53);x(54);x(55)]);
else 
    beta = restrictions*[1;x(1);x(2);x(3);x(4);x(5);x(6);x(7);x(8);x(9);x(10);...
        x(11);x(12);x(13)]*x(14);
end