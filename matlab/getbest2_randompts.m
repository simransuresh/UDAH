function     [POS,DO] = getbest2(POS,f,CONFIG,DO,varargin);

%    [POS,DO] = getbest2(POS,f,CONFIG,DO,II);
%
% get best historical data 
%
% do=0 means that an error occured
%
% by L.Boehme, IfM Kiel
% lboehme@ifm.uni-kiel.de
%
% last modification: 11/2003
% Modified on 2009.02.20 by Benjamin Rabe
%	adjust to own mapping routines...
%	Now works even if all 
%	
%*****************************************************
if DO==1
    % select points, if index passed on  -- Benjamin Rabe
    if nargin == 5
	II = varargin{1};
	POS.LON = POS.LON(II); POS.LAT = POS.LAT(II); POS.Z = POS.Z(II);
	if ~isempty(POS.JULD)
		POS.JULD = POS.JULD(II);
	end % if ~isempty
    end % if isfield
    %  Now get random points -- BR
	%	for ijk = []
    tindex_random = [];
    while(length(tindex_random)<str2num(CONFIG.CONFIG_MIN_CASTS))
	ttr = round(rand(1, ...
		ceil(str2num(CONFIG.CONFIG_MIN_CASTS))).*length(POS.LAT)) ;
%	LL=find(index_random==0);
%	index_random(LL)=1;
	ttr = ttr(ttr~=0);
	tti=1:length(ttr);
	tindex_random(end+tti) = ttr(tti);
	[dummy,ttj,dummy] = unique(tindex_random);
	tindex_random = tindex_random(ismember(1:length(tindex_random),ttj));
	tindex_random = tindex_random(1:min(length(tindex_random), ...
		str2num(CONFIG.CONFIG_MIN_CASTS)));
    end % while
    index_random=unique(tindex_random);
    index_junk = [1:length(POS.LAT)];
    index_junk = index_junk(~ismember(index_junk,index_random));
%    ii=[];
%    for h=1:length(POS.LAT)
%        a=find(index_random==index_junk(h));
%        if(isempty(a)==0)
%            ii=[ii,h];
%        end
%    end
%    index_junk(ii)=[];
    index_remain=index_junk;
	%	end
%    index_random = [];
%    index_remain = 1:length(POS.LON);

    % sort remaining by large scale
    [wei]=find_form(POS.LON(index_remain),POS.LAT(index_remain),POS.Z(index_remain),f.lon,f.lat,f.z,str2num(CONFIG.MAPSCALE_RADIUS_LARGE),str2num(CONFIG.MAPSCALE_PHI_LARGE)) ;
    [wei,ii] = sort(wei);
    index_large=index_remain(end-str2num(CONFIG.CONFIG_MIN_CASTS)+1:end);
    
    index_junk = index_remain;
    index_junk = index_junk(~ismember(index_junk,index_large));
%    ii=[];
%    for h=1:length(index_remain)
%        a=find(index_large==index_junk(h));
%        if(isempty(a)==0)
%            ii=[ii,h];
%        end
%    end
%    index_junk(ii)=[];
    index_remain=index_junk;
    % sort remaining by small scales
    [wei]=find_form(POS.LON(index_remain),POS.LAT(index_remain),POS.Z(index_remain),f.lon,f.lat,f.z,str2num(CONFIG.MAPSCALE_RADIUS_SMALL),str2num(CONFIG.MAPSCALE_PHI_SMALL)) ;
% no time weighting, for now! -- Benjamin Rabe
%	added time weighting again, as now included in mapping!!!
   if ~isempty(str2num(CONFIG.MAPSCALE_TIME))
	wei=wei.*exp(-(POS.JULD(index_remain)-f.day).^2 / ...
		str2num(CONFIG.MAPSCALE_TIME)^2);
   end % if ~isempty
    [wei,ii] = sort(wei);
    index_remain=index_remain(ii);
    % In case small scale points all/mostly in randrom or none exist -- BR
    %if length(index_remain) < str2num(CONFIG.CONFIG_MIN_CASTS)*2
    if length(index_remain) < str2num(CONFIG.CONFIG_MIN_CASTS)
     if isempty(index_remain)
	index_small = index_junk;
     else
	index_small = index_remain;
     end % if isempty
    else
     index_small=index_remain(end-str2num(CONFIG.CONFIG_MIN_CASTS)+1:end);
     %index_small=index_remain(end-str2num(CONFIG.CONFIG_MIN_CASTS)*2+1:end);
    end % length(index_remain)...
    
    % Now we have all historical profiles needed to calibrate
    % this float profile!
    % index values must be bigger than 0 !!!
    K=find(index_random<1);
    if ~isempty(K)
        disp([num2str(lenght(K)) ' random points lost....'])
    end
    clear K
    K=find(index_large<1);
    if ~isempty(K)
        disp([num2str(lenght(K)) ' large scale points lost....' ])
    end
    clear K
    K=find(index_small<1);
    if ~isempty(K)
        disp([num2str(lenght(K)) ' small scale points lost....' ])
    end
    clear K
   
    index=([index_random ,index_large,index_small]');
    K=find(index>0);       
    index=sort(index(K));
    clear K
    
    %POS.LAT     = POS.LAT(index);
    %POS.LON     = POS.LON(index);
    %POS.JULD    = POS.JULD(index);
    %POS.PRES    = POS.PRES;
    %POS.PSAL    = POS.PSAL(index,:);
    %POS.PTMP    = POS.PTMP(index,:);
    %POS.Z       = POS.Z(index);
    %%POS.SOURCE  = POS.SOURCE(index);

%		keyboard

    POS = index;
    
end %DO==1
