function  [pos,DO] = getpos2(H_DATA,f,CONFIG,DO);

% [pos,DO] = getpos2(H_DATA,f,CONFIG,DO);
%
% get historical data  within large scale radius
%
% do=0 means that an error occured
%
% by L.Boehme, IfM Kiel
% lboehme@ifm.uni-kiel.de
%
% last modification: 11/2003
% Modified on 2009.02.20 by Benjamin Rabe
%	adjust to own mapping routines...
%*****************************************************
if DO==1
    %rng = lbdist(H_DATA.LAT,H_DATA.LON,f.lat,f.lon);
    rng = m_idist(H_DATA.LON,H_DATA.LAT,f.lon,f.lat)/1000;
    II  = find(rng<=str2num(CONFIG.MAPSCALE_RADIUS_LARGE));
%	rng(II), H_DATA, f
    % use f/H coordinates to get data within large scale radius and PHI_large
    [wei]=find_form(H_DATA.LON(II),H_DATA.LAT(II),H_DATA.Z(II),f.lon,f.lat,f.z,str2num(CONFIG.MAPSCALE_RADIUS_LARGE),str2num(CONFIG.MAPSCALE_PHI_LARGE)) ;
%    dummy2=find(wei>=exp(-1)*str2num(CONFIG.MAPSCALE_PHI_LARGE));
    dummy2=find(wei>=exp(-1));% ???no obvious reason why this factor is used???

    II=II(dummy2);

    % historical profiles relevant for this float profile
%    pos.LAT = H_DATA.LAT(II);
%    pos.LON = H_DATA.LON(II);
%    pos.JULD= H_DATA.JULD(II);
%    pos.PRES= H_DATA.PRES;
%    pos.PSAL= H_DATA.PSAL(II,:);
%    pos.PTMP= H_DATA.PTMP(II,:);
%    pos.Z   = H_DATA.Z(II);
    pos = II;
    %pos.SOURCE = H_DATA.SOURCE(II);

    if isfield(CONFIG,'CONFIG_MIN_CASTS_ALT')
	if length(pos)<abs(str2num(CONFIG.CONFIG_MIN_CASTS_ALT))
        	%disp('Not enough historical profiles')
		DO=0;
	end
    else
	if length(pos)<str2num(CONFIG.CONFIG_MIN_CASTS)*3
        	%disp('Not enough historical profiles')
		DO=0;
	end
    end % if str2num
else
        disp('Bad Input!!')
        DO=0;
        pos=nan;
    
end
