function [vgrid,vgriderror,vdata,vdataerror]=mapping22(v,posfloat,...
    posdata,lam,phi,q,s,e)


% Inputs:
%
% v=data   m1*1 vector of historical data. This does not take NaN's.
% posfloat 1*3  vector of float location [latitude,longitude,dates]
% posdata  m1*3 matrix of data locations [latitude,longitude,dates]
% lam      scalar spatial scale (in km)
% phi      scalar cross_isobaric scale 
% q        scalar time scale (in days)
% s        scalar signal level (in variance terms)
% e        scalar noise level (in variance terms)
%
% Outputs:
%
% vgrid    m2*rank matrix of mapped fields (n is the number of repeated data)
% vgriderror  m2*rank matrix of error estimates of the mapped data
%             (assumed statistics)
%
% Modified on 2009.07.31 by Benjamin Rabe
%	Mean of data (v) is now calculated using method by Bretherton, etal
%
[m,n]=size(v);
%
% Set aside the memory that I need for this mapping
%
[nnn,i]=size(posfloat);
[mmm,i]=size(posdata);
vgrid=zeros(nnn,1);
vgriderror=vgrid;
vdata=zeros(mmm,1);
vdataerror=vdata;
%
% Create the data-data covariance function.  
%    
[Cdd] = s*covar2(posdata,posdata,lam,phi,q);
Cdd(1:mmm+1:mmm^2+mmm)=diag(Cdd)+e*ones(mmm,1);    % data-data cov. matrix
Cdd=inv(Cdd);
%
% calculate the objectively mapped fields on the data and grid fields.
% note the mean field is removed to ensure that estimate is unbiased
%
%
% try to estimate the mean field using the OI weights
% (eqn 21 of Bretherton, etal) -- as in Owens and Wong (2009; DSRI)
sum_Cdd = sum(sum(Cdd));
%%meanv=mean(v);
meanv = sum(Cdd*v)/sum_Cdd;
%%wght=Cdd*(v-meanv*ones(length(v),1));
wght = Cdd*( v - meanv );

% Create the grid-data covariance function.  
%

for j=1:mmm,
    [Cmd] = s*covar2(posdata(j,:),posdata,lam,phi,q);
    [vdata(j,1)] = Cmd*wght+ meanv;
    %%[vdataerror(j,1)]=sqrt(s-Cmd*Cdd*Cmd');
    % include error in mean (eqn 24 of Bretherton, etal)
    [vdataerror(j,1)] = sqrt(s - Cmd*Cdd*Cmd' + ...
		(1- sum(Cmd*Cdd,2)).^2/sum_Cdd);
end %for

[Cmd] = s*covar2(posfloat,posdata,lam,phi,q);
[vgrid] = Cmd*wght+ meanv;
%%[vgriderror]=sqrt(s-Cmd*Cdd*Cmd');
[vgriderror] = sqrt(s - Cmd*Cdd*Cmd' + ...
	(1- sum(Cmd*Cdd,2)).^2/sum_Cdd);

%disp('mapping fct 2nd stage...')

%	keyboard


return

