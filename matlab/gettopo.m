function [wdeps,wlons,wlats] = gettopo(inlons,inlats)
% function [wdeps,wlons,wlats] = gettopo(inlons,inlats)
%
% get waters depths at each data location (lon/lat)
disp('get topography for each data location...')
%% convert to Lambert Azimuthal coords
%% associate closest point in topographic database
whos
odstr = now
if isempty(dir('/opt/matlab/bin/matlab'))
	topolim = 10^5;
	llsecint = 1;
else
	topolim = 5*10^5;
	llsecint = 2.5;
	disp('you are on linsrv4!')%, input('...')
end % if isempty
% lon limits for loop -- assure that exclusive limits include all data!
lonsecs = [-180:llsecint:180-llsecint; -180+llsecint:llsecint:180]; 
lonsecs(1,:) = lonsecs(1,:)-10^-10; lonsecs(2,end) = lonsecs(2,end)+10^-10;
illength = size(lonsecs,2)
%wdeps = []; wlons = []; wlats = [];
wdeps = ones(size(inlons))*nan; wlons = wdeps; wlats = wdeps;
	%oinlons = [];
twsint = 2;
%% do this loop piecewise to speed up process in matlab
for jj = 1:illength
	lonsecs(:,jj)
	[jdind,tinlons,tinlats] = varsel(lonsecs(:,jj),inlons,inlats,1==1);
	if length(tinlons>1) & all(diff(tinlons)==0)
		tinlons(end) = tinlons(end)+10^-10;
	end
	if sum(jdind(:))>0
	 [twdep,twlons,twlats] = ...
		m_ibcao([min(tinlons) max(tinlons) min(tinlats) max(tinlats)]);
	 twdep = twdep(:); twlons = twlons(:); twlats = twlats(:);
	 twdep = twdep(1:twsint:end); twlons = twlons(1:twsint:end);
	 twlats = twlats(1:twsint:end);
	 disp(['only use non-nan values in topo -- m_ibcao returns nan for' ...
		' lon=-180 and lat=90']);
	 [twdep,twlons,twlats] = varsel(twdep,twlons,twlats);
	 %[jdind,twlons,twlats,twdep]=varsel(lonsecs(:,jj),twlons,twlats,twdep);
	 %twdep = twdep(:); twlons = twlons(:); twlats = twlats(:);
	 % convert to Lambert Azimuthal coords
	 %[wxkm,wykm] = llkm(twlats,twlons);
	 %[lxkm,lykm] = llkm(tinlats(:)',tinlons(:)');
		whos
	%if length(tinlons>1) & all(diff(tinlons)==0)
	if length(twlons)>topolim %10^5%1.5*10^5
			disp('check topo selection -- very large!!!')
			keyboard
	end
	 [plons,plats,idind,lld,jjdind]= ...
		distall2(tinlons(:),tinlats(:),twlons,twlats);
		whos
	 clear lld
	 %idind = abs((wxkm+wykm*i)-(lxkm(ii)+lykm(ii)*i));
	 %idind= abs( (wxkm(:,ones(size(lxkm)))+wykm(:,ones(size(lxkm)))*i) ...
	 %	- (lxkm(ones(size(wxkm)),:)+lykm(ones(size(wxkm)),:)*i) );
	 %jdind = min(idind,[],1);
	 %[idind,jdind] = find( idind == jdind(ones(size(wxkm)),:) );
	 if length(idind) ~= length(tinlons)
		disp('error in size of topo and station pos! correcting...')
			keyboard
		[dummy,jjdind,dummy] = unique(jjdind);
		whos
		idind = idind(jjdind);
		whos
%%%		input('...')
	 end % if length
	 %oinlons(end+(1:length(idind))) = tinlons;
	 %wdeps(end+(1:length(idind))) = twdep(idind);
	 %wlons(end+(1:length(idind))) = twlons(idind);
	 %wlats(end+(1:length(idind))) = twlats(idind);
	 if sum(jdind) ~= length(idind) | length(jdind)~=length(wdeps) | ...
		max(idind) > length(twdep)
			keyboard;
	 end % if sum
	 wdeps(jdind) = twdep(idind);
	 wlons(jdind) = twlons(idind);
	 wlats(jdind) = twlats(idind);
	 progcounter(illength,jj)
		if any(isnan(wdeps(jdind)))
			disp('nan in wdeps -- this shouldn''t be!!!')
			keyboard
		end % if any
	end % if sum
end % for ii
datestr(odstr)
datestr(now)
%
