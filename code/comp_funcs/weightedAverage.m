function outputArg = weightedAverage(I,weights)
%Calculates the weighted average
%   I=input Irradiance matrix from the CMOS Camera (must be double image)
%   outputArg= average intensity in the subregion of interest
%% Calculate the weighted average value
new_I = im2double(I).*255;
D_ponderada = weights.*new_I;
outputArg = sum(D_ponderada(new_I>0))./sum(weights(new_I>0));
end