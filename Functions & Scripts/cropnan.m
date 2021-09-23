function xcrop = cropnan(x)
% removes NaN and inf periods from a panel-form time series

nant = any(isnan(x),2);
inft = any(isinf(x),2);

xcrop = x(~nant & ~inft,:);