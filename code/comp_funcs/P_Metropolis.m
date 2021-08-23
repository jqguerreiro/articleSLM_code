function outputArg = P_Metropolis(cost_old,cost_new,Tk)
%Applies the Metropolis principle to calculate the SAA probability of
%acceptance.
outputArg = exp(-(cost_new-cost_old)/Tk);