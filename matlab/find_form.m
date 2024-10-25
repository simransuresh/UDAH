    function [wei]=find_form(lon_his,lat_his,z_his,lon_float,lat_float,z_float,Lx,Lz) 
%   function [wei]=find_form(lon_his,lat_his,z_his,lon_float,lat_float,z_float,Lx,Lz) 
%
%   lon_his, lat_his       == coordinates of historical datapoints (vectors)
%   z_his                  == depths at historical profiles
%   lon_float, lat_float   == coordinates of float profile
%   z_float                == depth at float location
%
%   Lx                     == characteristical length scale in km (150)
%                             at which the weights becom 0.5 (due to distance)
%   Lz                     == scale factor for f/H difference (.3)
%
%   latter 2 arguments are optional.   
%   on return weights (wei) are given
%
%   weights datapoints to the float position using the topographic-following
%   mapping scheme from Davies (1998)
%   scales set for the north atlantic from Rhein & Fischer (2002)
%
% 
%  last modified: 10/2003 by L. Boehme
%
%	Modified 2009.03.23 by Benjamin Rabe

  if nargin<6;
    Lx=680; Lz=.24;
  end
  
  %Lx=(-1)*log(.5)/(Lxx^2);

   
    
%  derive contributions to generalized distance
% ********************************************************************************   

%  units in km for distance and depth
%  now calculate distance between float position and all other data

  %xd=(lon_float-lon_his).*cos((lat_float+lat_his)./2*pi/180)*60*1.852; 
  %yd=(lat_float-lat_his)*60*1.852;  
  %dxy_sq=(xd.^2+yd.^2);   % distance from float to data  squared
% Use accurate distance function -- BR
  %% distance from float to data  squared
  dxy_sq=(m_idist(lon_his,lat_his,lon_float,lat_float)/1000).^2;
	%% m_idist returns nan for equal positions !
  dxy_sq(lon_his+lat_his*i==lon_float+lat_float*i) = 0;

%	if ~isempty(dxy_sq), keyboard; end

  R_dis=1/Lx.*sqrt(dxy_sq);   % isotropic distance weighting 

% derive planetary vorticity at each data point
   Plv=2*7.292e-5.*sin(lat_his.*pi/180);
   Plv_float=2*7.292e-5.*sin(lat_float.*pi/180);

% now calculate f/H scaling accordingly
   ii=find(z_his==0);z_his(ii)=0.001;              % avoid division by zero
   zd_sq=((Plv_float/z_float-Plv./z_his)).^2;  
   znom=(Plv_float/z_float).^2+(Plv./z_his).^2;
   R_depth=(1/Lz).*sqrt(zd_sq./znom);  % weight for depth dependence

%  generalized distance, that is a function of distance and waterdepth; f/H on larger scale   
   RR=R_dis.^2+R_depth.^2;

%  WEIGHTS
   wei=exp(-RR);                % exists for each data point
%
%	[exp(-R_dis.^2), exp(-R_depth.^2)]'
%	[xd,yd,sqrt(dxy_sq)]'
%	tt = sqrt(dxy_sq);
%		if any(tt<10000)
%			tt, lbdist(lat_his,lon_his,lat_float,lon_float)'
%			input('');
%		end
%	wei'
%	input('')
%
