function [BinData, xR]  = data_binning1(x,y,xRange,varargin)
% Function for binning data into specified bins
% x -> x range
% y -> variable value
% xRange -> range of x bins
% specify fourth argument to find mean ignoring NaN
[a,b]   = histc(x,xRange);
nbin    = length(xRange)-1;
BinData = zeros(3,nbin);

argsin  = nargin;

if(argsin == 3)
    for i = 1:nbin
        WhichBin    = (b == i);
        BinData(1,i)= mean(y(WhichBin));   
        
        BinData(2,i)= std(y(WhichBin));
        BinData(3,i)= sum(WhichBin);
        BinData(4,i)= std(x(WhichBin));
    end
else
    for i = 1:nbin
        WhichBin    = (b == i);
        BinData(1,i)= mean(y(WhichBin),'omitnan');
        BinData(2,i)= std(y(WhichBin),'omitnan');
        tt          = y(WhichBin);
        ss          = isnan(tt);
        BinData(3,i)= sum(tt) - sum(ss);
        BinData(4,i)= std(x(WhichBin),'omitnan');
    end
end
dx      = (xRange(2) - xRange(1))/2;
xR      = xRange(1:end-1) + dx;

end
    
