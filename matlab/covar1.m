% function [a] =covar1(x1,x2,lam,phi,q)
%
% Output:
% a 	m*n  unscaled Gaussian covariance function
%
% Input:
% x1	m*3  model grid, col 1+2+3 are lat, long, depth respectively
% x2	n*3  data grid, col 1+2+3 are lat, long, depth respectively
% lam      scalar spatial scale (in km)
% phi      scalar cross_isobaric scale 
% q        scalar time scale (in days)
%
%	Modified 2009.03.23 by Benjamin Rabe

function  [a] =covarxy(x1,x2,lam,phi,q)

m=size(x1);
n=size(x2);
a=zeros(m(1),n(1));
one=ones(n(1),1);

%lam3=3*lam;

Zx1=x1(:,4); % depth at x1 datapoint
Zx2=x2(:,4); % depth at x2 datapoint

% derive planetary vorticity at each data point
PV_x1=(2*7.292e-5.*sin(x1(:,1).*pi/180))./Zx1;
PV_x2=(2*7.292e-5.*sin(x2(:,1).*pi/180))./Zx2;

   
   
for i=1:m(1)
    
    %  units in km for distance and depth
    %  now calculate distance between float position and all other data
    %xd=(x1(i,2)-x2(:,2)).*cos((x1(i,1)+x2(:,1))./2*pi/180)*60*1.852; 
    %yd=(x1(i,1)-x2(:,1))*60*1.852;  
    %dxy_sq=(xd.^2+yd.^2);   % distance from float to data  squared
% Use accurate distance function -- BR
  %% distance from float to data  squared
  dxy_sq=(m_idist(x2(:,2),x2(:,1),x1(i,2),x1(i,1))/1000).^2;
	%% m_idist returns nan for equal positions !
  dxy_sq(x1(i,1)+x1(i,2)*i==x2(:,1)+x2(:,2)*i) = 0;
    %Lx=(log(.5)*(-1))./(lam3^2);
    R_dis=sqrt(dxy_sq')./lam;   % isotropic distance weighting 

    % now calculate f/H scaling accordingly
    zd_sq=((PV_x1(i)-PV_x2)).^2;  
    znom=(PV_x1(i)).^2+(PV_x2).^2;
    R_depth=((1/phi).*sqrt(zd_sq./znom))';  % weight for depth dependence
    
    %  generalized distance, that is a function of distance and waterdepth; f/H on larger scale   
    RR=R_dis.^2+R_depth.^2;


%  WEIGHTS
   a(i,:)=exp(-RR); 
    
end

%%%if any(abs(a)>100), disp('covar1.m 62'), keyboard; end
%
