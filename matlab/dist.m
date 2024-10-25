
function [rng] = dist(lat1,lon1,lat2,lon2)

%  [rng] = dist(lat1,lon1,lat2,lon2)
%  Calculates great circle distance between points 
%  on a spherical geoid
%  The output is in km!
%
%  lat1,lon1 = nx1 in degrees
%  lat2,lon2 = nx1 or 1x1 in degrees
%  
% by L.Boehme, IfM Kiel
% lboehme@ifm.uni-kiel.de
%
% last modification: 03/2004
%*****************************************************

%  Dimension tests

if ~isequal(size(lat1),size(lon1),size(lat2),size(lon2))
    
    if (size(lon2)==[1 1] & size(lon2)==[1 1])
        
        lat2 = ones(size(lat1)).*lat2;
        lon2 = ones(size(lon1)).*lon2;

    else
    
      msg = 'Inconsistent dimensions for latitude and longitude';
	  error(msg); 
	  return
      
    end
    
end

% unit conversion from degrees to radians

lat1 = lat1.*pi./180;
lon1 = lon1.*pi./180;
lat2 = lat2.*pi./180;
lon2 = lon2.*pi./180;


% calculate the distance

dd1 = sin(lat1).*sin(lat2);
dd2 = cos(lat1).*cos(lat2).*cos(lon2-lon1);
dd3 = dd1+dd2;
rad   = acos(dd3).*180./pi;  % in degree


rng = rad.*60.*1.853;
