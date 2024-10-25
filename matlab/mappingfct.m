function [vgrid,vgriderror,varargout] = ...
	mappingfct(lons,lats,julds,invar,newlons,newlats,newjulds, ...
	xyscale,zscale,tscale,Hgrid,hlon,hlat,Hdat,varargin)
%function [vgrid,vgriderror,vgrid2,vgriderror,vdatanos] = ...
%	mappingfct(lats,lons,julds,invar,newlons,newlat,newjulds, ...
%	xyscale,zscale,tscale,H,hlon,hlat,Hdat,mincasts)
%
% Modified by Benjamin Rabe on 2009.04.03
%	Now use 2-stage mapping:
%		1. Large scale (smooth)
%		2. Small scale (also individual stations)
%	All points that are not mapped in stage 2. just use value from 1.
%	If more than config2.CONFIG_MIN_CASTS # stations available for mapping,
%	choose 1/3 random, 1/3 large scale, 1/3 small scale
% Modified by Benjamin Rabe on 2009.07.24
%	Use time scale in mapping
% Modified by Benjamin Rabe on 2009.09.29
%	Add no of data points for each grid point as optional output
%
config.MAPSCALE_RADIUS_LARGE = num2str(xyscale(1));
config.MAPSCALE_PHI_LARGE = num2str(zscale(1));
config.MAPSCALE_RADIUS_SMALL = num2str(xyscale(2));
config.MAPSCALE_PHI_SMALL = num2str(zscale(2));
config.MAPSCALE_TIME = num2str(tscale);
%config.CONFIG_MIN_CASTS = varargin{1};
config.CONFIG_MIN_CASTS = '1';
config2 = config; config2.CONFIG_MIN_CASTS = '20'%'20'
%config2 = config; config2.CONFIG_MIN_CASTS = '100'%'40'
config.CONFIG_MIN_CASTS_ALT = '1'
iindex = find(~isnan(invar));
%
Hdat = abs(Hdat); % depth always positive !!!!!!!!!
Hgrid = abs(Hgrid);
%	xyscale = [300 100], zscale = [0.2 0.1]
%% Add random number to lon/lat/depth to avoid singular covariance matrices
%% -- 10^-4 deg (m) is only 1.1 km in latitude (1 mm in depth) -- negligible!
% lats = lats+rand(size(lats,1),1)*10^-2;
% lons = lons+rand(size(lons,1),1)*10^-2;
% Hdat = Hdat+rand(size(Hdat,1),1)*10^-2;
% check if topography passed on -- call mapping functions accordingly
  H_DATA.LON = lons(iindex);
  H_DATA.LAT = lats(iindex);
  H_DATA.Z = Hdat(iindex);
  H_DATA.JULD = julds(iindex);
  invar = invar(iindex); % eliminate nan in invar as well !!!
  posxyzt = [lats(iindex) lons(iindex) julds(iindex) Hdat(iindex)];
% posxyzt(:,[1 2 4]) = posxyzt(:,[1 2 4])+rand(size(posxyzt,1),3)*10^-4;
%% Add random number to mapping variable ...
% invar = invar+rand(size(invar,1),1)*std(invar)*10^-5;
%%%		whos
  % run mapping loop
  iilength = length(newlons); icount = -1;
%	ilind = find(newlons<-140&newlons>-145&newlats>75)
%  for ii = ilind(:)'  %1:iilength
%  for ii = floor(iilength/1.7):iilength
  for ii = 1:iilength
	%% select only data within certain radius of mapping position
	f.lon = newlons(ii); f.lat = newlats(ii); f.z = Hgrid(ii);
	f.day = newjulds(ii);
	[iindex,DO] = getpos2(H_DATA,f,config,1);
	if DO~=0
	 %% if many reference profiles, reduce -- otherwise, use all
	 if length(iindex) > str2num(config2.CONFIG_MIN_CASTS)*3
		%[iindex2,dummy] = getbest2(H_DATA,f,config2,1,iindex);
		[iindex2,dummy] = getbest2_randompts(H_DATA,f,config2,1,iindex);
		iindex = iindex(iindex2);

				if ii > 250 & ii < 280, iindex', end

	 end % if length
	 siglev = signal(invar(iindex)');
	 errlev = noise(invar(iindex),lats(iindex),lons(iindex),...
		Hdat(iindex),config);
	 %[vgrid,vgriderror,vdata,vdataerror]=...
	 %	mapping1(v,posfloat,posdata,lam,phi,q,s,e)
	 %Hgrid = interp2(hlon,hlat,H,newlon(ii),newlat(jj),'linear');
	 %% Measurement depth has to be regular!!! not true, but almost...(;
	 newposxyzt = [newlats(ii) newlons(ii) newjulds(ii) Hgrid(ii)];
	 if ( length(errlev)==1 & all(isnan(errlev)) ) | ...
		sum(sum(diff(posxyzt(iindex,:),1,1))) == 0
	  %% If only one data point, just assign value !
%		ii,iindex,size(invar),size(vgrid)
		%keyboard
	  vgrid(ii) = meanmiss(invar(iindex)); vgriderror(ii) = nan;
	 else
	  %% Stage 1: Large scale
		%mapping1(invar(iindex),newposxyzt,posxyzt(iindex,:),...
	  [tvgrid,tvgriderror,vdata,vdataerror] = ...
		mapping11(invar(iindex),newposxyzt,posxyzt(iindex,:),...
		xyscale(1),zscale(1),[],siglev,errlev);

%		if any(abs(tvgrid)>100), disp('mapping11.m 81'), keyboard; end
	  %%% Add random number to lon/lat/depth to avoid singular covariance
	  %%% matrices -- 10^-4 deg (m) is only 1.1 km in latitude (1 mm in
	  %%% depth) -- negligible!
	  %%% Iterate this procedure until grid value within range of data
	  %%% If this doesn't work, need condition for sing. matrix in covar1.m!
	  iitermax = 10; iiter = 0; tposxyzt = posxyzt(iindex,:);
	  tinvar = invar(iindex);
	  while ( tvgrid > max(tinvar) | tvgrid < min(tinvar) ...
			| isnan(tvgrid)) & ( iiter < iitermax )
		iiter = iiter+1;
		tposxyzt = posxyzt(iindex,:);
	  	tposxyzt(:,[1 2 4]) = posxyzt(iindex,[1 2 4]) - ...
			rand(length(iindex),3)*10^-2;
%%%		if any(abs(tposxyzt(:,1))>90), keyboard; end
		errlev = noise(tinvar,tposxyzt(:,1),tposxyzt(:,2),...
			tposxyzt(:,4),config);
		[tvgrid,tvgriderror,vdata,vdataerror] = ...
			mapping11(tinvar,newposxyzt,tposxyzt,...
			xyscale(1),zscale(1),[],siglev,errlev);
	  end % if tvgrid
	  if iiter > 0
		iiter,tvgrid,max(tinvar),min(tinvar)
		disp('Matrix nearly singular -- add random noise to lat etc.')
	  end % if iiter
	  if tvgrid<min(tinvar) | tvgrid>max(tinvar) | isnan(tvgrid)
	   %%% Adding random noise didn't work -- average nearby data... TODO
	   disp('Matrix still singular -- average nearby values')
	   [dummy,pind] = sort(posxyzt(iindex,1)+posxyzt(iindex,2)*i);
	   iindex = iindex(pind);
	   [unipos,pind2,pind] = unique(round(( posxyzt(iindex,1) + ...
		posxyzt(iindex,2)*i)*10^3));
           tposxyzt = posxyzt(iindex(pind2),:);
	   tinvar = tinvar(pind2);
	   for ji = 1:length(unipos)
		idp = find(round(( posxyzt(iindex,1) + ...
			posxyzt(iindex,2)*i)*10^3) == unipos(ji));
		if length(idp) > 1
			tposxyzt(ji,:)=meanmiss(posxyzt(iindex(idp),:));
			tinvar(ji,:)=meanmiss(invar(iindex(idp)));
		end % if length(idp
	   end % for ji
	   errlev = noise(tinvar,tposxyzt(:,1),tposxyzt(:,2),...
		tposxyzt(:,4),config);
	   siglev = signal(tinvar');
	   [tvgrid,tvgriderror,vdata,vdataerror] = ...
		mapping11(tinvar,newposxyzt,tposxyzt,...
		xyscale(1),zscale(1),[],siglev,errlev);
	  end % if tvgrid<min
%	  if ( tvgrid<min(tinvar) | tvgrid>max(tinvar) | isnan(tvgrid) ) &
%		range(tinvar) < std(invar)*10^-5
%	  end % if tvgrid<min
	  if tvgrid>min(tinvar) & tvgrid<max(tinvar)
	   %%% residual at data points...
	   %resinvar = invar(iindex) - vdata ;
	   resinvar = tinvar - vdata ;
	   %  nindex = abs(resinvar) > 2*std(invar(iindex));
	   %if tvgrid > 40, keyboard; end
	   %  if sum(nindex) > 1
	   %%%% If outliers in data, redo stage 1 mapping
	   %iindex = iindex(~nindex);
	   %siglev = signal(invar(iindex)');
	   %errlev = noise(invar(iindex),lats(iindex),lons(iindex),...
	   %Hdat(iindex),config);
	   %[tvgrid,tvgriderror,vdata,vdataerror] = ...
	   %mapping1(invar(iindex),newposxyzt,posxyzt(iindex,:),...
	   %xyscale(1),zscale(1),tscale,siglev,errlev);
	   %  resinvar = invar(iindex) - vdata ;
	   %  end % if sum
	   %keyboard
	   %sigresinvar = signal(vdata'); % WRONG!!! signal of residuals, not
	   % from mapped data from 1st stage!!!
	   sigresinvar = signal(resinvar');
	   %% Stage 2: Small scale
	   %whos
	   %mapping2(resinvar,newposxyzt,posxyzt(iindex,:),...
	   [tvgrid2,tvgriderror2,tvdata,tvdataerror] = ...
		mapping22(resinvar,newposxyzt,tposxyzt,...
		xyscale(2),zscale(2),tscale,sigresinvar,errlev);
	   %mapping1nopv(invar(iindex),newposxyzt,posxyzt(iindex,:),...
	   %xyscale(2),siglev,errlev);
	   %% If stage 2 mapping didn't work, use just stage 1 !
	   if tvgrid+tvgrid2<min(tinvar) | tvgrid+tvgrid2>max(tinvar) | ...
		isnan(tvgrid2)
		tvgrid2 = 0; % newposxyzt has one row, so one mapped value!
	   end % if tvgrid+tvgrid2
%	   tvgrid2(isnan(tvgrid2)) = 0;
	   tvgriderror2(~isnan(tvgrid2)) = tvgriderror2(~isnan(tvgrid2))+i*0;
	   tvgriderror2(isnan(tvgrid2)) = 0+tvgriderror*i;
	   vgrid(ii) = tvgrid+tvgrid2;
	   vgriderror(ii) = tvgriderror2; %Error from stage 2, if existing !!!
	   vgrid2(ii) = tvgrid2; % 2nd stage residual map
%		if vgrid(ii)<min(tinvar)|vgrid(ii)>max(tinvar)
%			disp('mappingfct.m 180'), keyboard;
%		end % if vgrid2
	   %%%if vgrid(ii) < 0, keyboard; end
	  else
		disp('Matrix still singular -- set grid value to nan')
			tvgrid, tvgriderror
			tvgrid,max(tinvar),min(tinvar)
%				keyboard;
	  	vgrid(ii) = nan; vgriderror(ii) = nan;
		vgrid2(ii) = nan;
	  end % if tvgrid
	 end % if length
%	whos, iindex, vgrid(ii), input('...3rd mappingfct')
	 vdatanos(ii) = length(iindex); % number of data points used
	else
	 vgrid(ii) = nan; vgriderror(ii) = nan;
	 vgrid2(ii) = nan; vdatanos(ii) = nan;
	 %vdata(ii) = nan; vdataerror(ii) = nan;
	end % if DO~=0
	icount = progcounter(iilength,ii,1,icount);
	%sprintf('%d(%d)',ii,iilength)
  end %for ii
%
if nargout > 3
	varargout{1} = vgrid2; varargout{2} = vgriderror;
	if nargout > 4
		varargout{3} = vdatanos;
	end % if nargout > 4
end % if nargout > 3
%whos
if ~exist('vgrid','var'), vgrid = []; vgriderror = []; end

%	keyboard
%
