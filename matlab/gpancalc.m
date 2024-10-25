function [sdynh,varargout] = gpancalc(salg,ptmpg,pressg,pref,lats)
% function [sdynh,sgpan,sgvel] = gpancalc(salg,ptmpg,pressg,pref,lats)
%
% calculate specific volume anomaly and geopotential anomaly
db2pas = 10^4;
%pressg = pressg(ones(size(salg,2),1),:,ones(size(salg,3),1));
ptmpg = permute(ptmpg,[2 1 3]); salg = permute(salg,[2 1 3]);
if prod(size(pressg))~=max(size(pressg))
	pressg = permute(pressg,[2 1 3]); pq = 1==1;
else
	pq = 0==1;
end % if prod
tempg = sw_temp(salg,ptmpg,pressg,0);
svan = sw_svan(salg,tempg,pressg);
idind = abs(pressg-pref);
if pq
		disp('not yet implemented!!!'), return
	idind2 = min(idind,[],2);
	idind = find(idind2(:,ones(size(pressg,2),1))==idind);
	dpress = diff([0 pressg]);
else
	idind = find(min(idind)==idind);
	dpress = diff([0 pressg]);
	dpress = [dpress(1)+dpress(2)/2 (dpress(2:end-1)+dpress(3:end))/2 ...
		dpress(end)];
	sgpan = (summiss((svan(:,1:idind).* ...
		dpress(ones(size(svan,1),1),1:idind))')' ...
		- svan(:,idind)) * db2pas;
	sdynh = sgpan./getg(lats,0);
	if nargin == 6
	 	lons = varargin{1};
		f = sw_f(lats); f = f(:);
		gvdist = m_idist(lons(1:end-1),lats(1:end-1), ...
			lons(2:end),lats(2:end));
		sgvel = -( sgpan(2:end)-sgpan(1:end-1) ) ./ ...
			( f .* gvdist );
	end % if varargin == 6
end % if pq
% output...
if nargout > 1
	varargout{1} = sgpan;
	if nargout == 2
		varargout{2} = svel;
	end % if nargout == 2
end % if nargout > 1
%
