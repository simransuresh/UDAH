
function [noisevariance]=noise(a,xlat,xlong,xdepth,config)

% function [noisevariance]=noise(a,xlat,xlong,xdepth,config)
%
% for float calibration!!
%
% INPUT:
% xlat = a column vector (nx1)
% xlong = a column vector (nx1)
% a = (nx1) vector consisting of n casts
% - this can contain NaN's.
%
% OUTPUT:
% noisevariance = variance of noise between closest casts
% - an estimate of the noise due to variations between casts,
%   one variance value per level.
%
% NOTE:
% noise = measured value - true value.
% 2*noisevariance = variance between casts. Here, the variance is
% based on closest casts.
%
% modified by L.Boehme, IfM Kiel
% lboehme@ifm.uni-kiel.de
%
% last modification: 10/2003


n=length(a);
nn = length(unique(xlat));
diff=NaN*ones(n,1);
noisevariance=zeros(1);

if nn<2, noisevariance(1:n) = nan; return, end

% for each measurement
for i=1:n
    % find closest measurement
    [wei]=find_form(xlong,xlat,xdepth,xlong(i),xlat(i),xdepth(i),str2num(config.MAPSCALE_RADIUS_LARGE),str2num(config.MAPSCALE_PHI_LARGE)) ; 
    index=find(wei<1);
    [tmp,j]=max(wei(index));
%		if isempty(index), keyboard; end
	diff(i)=a(i)-a(index(j));
end

ii=find(isnan(diff)==0);
noisevariance=sum(diff(ii).^2)/(2*length(ii));
