%****************************************************************************
%
% TAU_N.M
%
% programmer:	Shayne McGregor
% date:		10 March 2009
%
% subject: This function is used to calculate the "integral time scale
% determining the time period required to gain a new degree of freedom"
% (cf p252 of Davis (1976))
%
% Davis, R.E. (1976), "Predictability of sea surface temperatures and
% sea level pressure anomalies over the North Pacific Ocean"
% Journal of Physical Oceanography, 6, 249-266.
%
%****************************************************************************
% INPUT
% autocorrA = autocorrelation function for time series A
% autocorrB = autocorrelation function for time series B
% deltat = time interval between each of the N observations
%
% OUTPUT
% taun = integral time scale determining the time period required to gain
% a new degree of freedom (cf p252 of Davis (1976))
%
%****************************************************************************
%
function [taun]=tau_n(autocorrA,autocorrB,deltat)
%
% calculate taun
%
taun=0;
for i=1:length(autocorrA),
    taun=taun+autocorrA(i)*autocorrB(i)*deltat;
end
%
%****************************************************************************
