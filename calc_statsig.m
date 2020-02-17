% This script was written to calculate the statistical
% significance of the correlation coeffients using the
% reduced degrees of freedom of davis 1976.
%
% Inputs: (ts1,ts2,delt)
%  ts1 is the first time series
%  ts2 is the second time series
%  delt is the time thing (does not appear to have an effect on results)
%       delt will default to 1 if not specified
% Beware of the means of the time series influencing outputs
% Outputs: [r,sig,df]
%   r is correlation
%   sig is significance (0, 90, 95 or 99)
%   df is the degrees of freedom
% If the output is not specified, r will be the output.
% 
% written by shayne mcgregor 
% on the 27th Jan 2010
%
%-----------------------------------------------------
function[r,sig,df]=calc_statsig(ts1,ts2,delt);

rtmp=corrcoef(ts1,ts2);
r=rtmp(2,1);

if nargin < 3
    delt = 1;
end

%--------------------------------------------------------
% calculate the "effective number of degrees of freedom"
%--------------------------------------------------------

% calculate the normalised autocorrelation functions for each of the two
% time series

autocorrA=xcorr(ts1,'coeff');
autocorrB=xcorr(ts2,'coeff');

% calculate the integral time scale determining the time period required to
% gain a new "degree of freedom" (cf p252 Davis (1976))

[taun]=tau_n(autocorrA,autocorrB,delt);

N=length(ts1);
df=N*delt/taun;

%--------------------------------------------------------------------------
% calculate the t-statistic using the "effective number of degrees of
% freedom"
% (i.e. substitute "N-2" in the t-statistic equation with "df")
%--------------------------------------------------------------------------

if r == '1'
  'CORRELATION COEFFICIENT = 1.  WILL GET DIVIDE BY ZERO IN T-STATISTIC!!!'
end

Nminus2=df;
tstatistic=(r*sqrt(Nminus2))/(sqrt(1-r^2));

%--------------------------------------------------------------------------
% now calculate the statistical significance
% values come from t distribution in Schaum's stats 
%--------------------------------------------------------------------------

dfc=[(1:30),40,60,120];
sig90c=[3.08,1.89,1.64,1.53,1.48,1.44,1.42,1.40,1.38,1.37,1.36,1.36,1.35,1.34,1.34,1.34,1.33,1.33,1.33,1.32,1.32,1.32,1.32,1.32,1.32,1.32,1.31,1.31,1.31,1.31,1.30,1.29,1.28];
sig95c=[6.31,2.92,2.35,2.13,2.02,1.94,1.90,1.86,1.83,1.81,1.80,1.78,1.77,1.76,1.75,1.75,1.74,1.73,1.73,1.72,1.72,1.72,1.71,1.71,1.71,1.71,1.70,1.70,1.70,1.70,1.68,1.67,1.66];
sig99c=[31.82,6.96,4.54,3.75,3.36,3.14,3.00,2.90,2.82,2.76,2.72,2.68,2.65,2.62,2.60,2.58,2.57,2.55,2.54,2.53,2.52,2.51,2.50,2.49,2.48,2.48,2.47,2.47,2.46,2.46,2.42,2.39,2.36];

ii=find(dfc<=df);
%if length(ii)>=1
sig90=sig90c(ii(end));
sig95=sig95c(ii(end));
sig99=sig99c(ii(end));
%else
%keyboard
%end

if abs(tstatistic)>=sig99
sig=99;
elseif abs(tstatistic)>=sig95
sig=95;
elseif abs(tstatistic)>=sig90
sig=90;
else
sig=0;
end



