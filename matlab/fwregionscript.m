%function fwregionscript()
%

%clear
global tpath

%%%hq = 1 % top 500m
%%%hq = 2 % top 150m
%%%hhq = 3 % sal < 34
hhq = 1 % top 500 m
hq = 1; % FW content salinity < 34 only for anomalies to own 91to08 mean !!!
	% and EWG summer
cdatachoice = '6'
if ~exist('ffseason','var')
	ffseason = ''%'2006to2008JAS'%1992to1999JAS'%'yearlyavgJAS'
end % if ~exist

if ~exist('vname','var'), vname = 'hfw_', end
if ~exist('testq','var'), testq = '', end
if ~exist('sref','var'), sref = 35; fsref = sref; end
if ~exist('fsref','var'), sref = 35; fsref = sref; end
if ~exist('plim','var')
	fplim = '';
elseif isreal(plim(1))
	fplim = sprintf('ignore%06.1fto%06.1fm.',plim(1),plim(2))
else
	fplim = sprintf('ml%06.1fto%06.1fm.',imag(plim(1)),imag(plim(2)))
end % if ~exist('plim
flayer = '';

if exist('layersal','var')
	if layersal ~= 34 | syear < 1992
		flayer = sprintf('layertoisohanline%04.1f.',layersal);
	end % if layersal
end % if ~exist(flayer
if ~exist('hfwlims','var')
	hfwlims = [-5 30; 0 100]; % simple outlier elimination FW inventories
end

% random number generator...
if isempty(strfind(testq,'randommaptest'))
	randreset();
else
	randrand();
end % if isempty(strfind(testq

%ggq = '2dbbins.qc' % new -- before had 10 db bins for ITP, ARK 22 & 23, POPS
tggq = '2dbbins.qc.' % new -- before had 10 db bins for ITP, ARK 22 & 23, POPS
ggq = sprintf(['%s%ssref%04.1f.%s'],tggq,fplim,fsref,flayer)

	for ijk = []

%%%ii = [3 6:20 22:30]
% put WOD09 first, so 'ctdelimdouble.m' keeps non-WOD source of duplicates!!!
% order is important: SCICEX data source 32 is preferred over 33 !
%ii = [26 3 6:20 22:23 25 27 30 32 33 31] % JAS
if ~exist('ii','var')
	ii = [26 3 6:20 21:23 25 27 30 32 33 31] % JAS -- with XCTD (2008) !!!
end
%	ii = [32 33 31] % testing duplicate station eliminator...
%ii = [3 6:20 22:23 25:27 30] % JAS
%%%	ii = ii([2:14 16 18:end]) % JAS -- only ship data (no ITP/POPS)!!!
%%%ii = [3 19 22 24 26 28:29] % MAM
%
%%%pseason = [{'monthly'} {'phc'}]
%%%cdatachoice = '3'
%%%pseason = [{'all'} {'mn91to08'}]
%%%cdatachoice = '7'
if ~exist('pseason','var')
	pseason = [{'summerJAS'} {'ewg'}]
	cdatachoice = '6'
end % if ~exist('pseason
%%%pseason = [{'winterMAM'} {'ewg'}]
%%%cdatachoice = '7'
%%%pseason = [{'JAS'} {'mn91to08'}]
%%%cdatachoice = '4'
%%%pseason = [{'MAM'} {'mn91to08'}]
%%%cdatachoice = '5'
juldlims = []
%hfwlims = [0 450; 0 450]; % simple outlier elimination halocline layer depth

%%%		end % for ijk

% get filenames for chosen season and data set indices
fnames =  climatanofilematch(ii,pseason,ggq)
if ~isempty(ggq), ggq = ['.' ggq]; end

%%%			for ijk = []

%load previously calculated freshwater content and anomalies to climatology
if ~exist('plim','var')
	% (duplicates eliminated in the process!!!)
	%%%[alldat,dexc22] = fwregionload(hq,fnames,vname,1)
	[alldat,dexc22] = fwregionload(hq,fnames,vname,1)%3);
else
	%load FW inv. without plim to select profiles
	% (use gap test from normal FW inv.)
	% duplicates not eliminated!
	[alldat,dexc22] = fwregionload(hq,fnames,vname,0)

		testalldat = alldat;

	ggqpl = sprintf(['%s%ssref%04.1f.'],tggq,'',fsref)
	fnamespl =  climatanofilematch(ii,pseason,ggqpl)
	[alldatpl,dexc22pl] = fwregionload(hq,fnamespl,vname,0)
	lljpl = alldatpl.lons+alldatpl.lats+alldatpl.julds*i;
	llj = alldat.lons+alldat.lats+alldat.julds*i;
	%lljpl = round(lljpl.*10^5*[1+i]).*10^-5*[1+i];
	%llj = round(llj.*10^5*[1+i]).*10^-5*[1+i];
	%% all locations not eliminated in FW inventories without gap and
	%% also contained in those with gap
	lljind = find(ismember(llj,lljpl));
	%% all locations not contained in FW inventories with gap but in those
	%% without gap
	lliind = find(~ismember(lljpl,llj));
	%dtno = eps(meanmiss(alldatpl.lons(:)))*10^5
	%dtno = 10^-10
	%[dummy,dummy,lljind,lld,lliind] = ...
	%	distall2(alldatpl.lons+dtno,alldatpl.lats+dtno, ...
	%	alldat.lons,alldat.lats);
	alldat.lons = [alldat.lons(lljind);alldatpl.lons(lliind)];
	alldat.lats = [alldat.lats(lljind);alldatpl.lats(lliind)];
	alldat.julds = [alldat.julds(lljind);alldatpl.julds(lliind)];
	%% set missing FW inventories in calculation with gap to zero
	alldat.hfw = [alldat.hfw(lljind);zeros(size(alldatpl.hfw(lliind)))];
	alldat.ahfw = [alldat.ahfw(lljind);nan(size(alldatpl.ahfw(lliind)))];
	alldat.wdeps = [alldat.wdeps(lljind);alldatpl.wdeps(lliind)];
	alldat.maxpres = [alldat.maxpres(lljind);alldatpl.maxpres(lliind)];
	% eliminate duplicates...
	[idind,dind,alldat.lons,alldat.lats,alldat.julds] = ...
	  ctdelimdouble(alldat.lons,alldat.lats,alldat.julds,alldat.maxpres);
	alldat.elimdoubleno = [length(alldat.maxpres)-length(idind)];
	alldat.maxpres = alldat.maxpres(idind);
	alldat.wdeps = alldat.wdeps(idind);
	eval(['alldat.' vname(vname~='_') ' = alldat.' ...
		vname(vname~='_') '(idind);']);
	eval(['alldat.a' vname(vname~='_') ...
		' = alldat.a' vname(vname~='_') '(idind);']);;
	clear idind dind llj lljpl lliind lljind alldatpl dexc22pl
end % if exist('plim


%			save ttemp.noplim.fwregionscript.matchlocs.mat
%			save ttemp.mlplim.fwregionscript.matchlocs.mat
% eliminate unrealistic FW content values by range and region in FW content
[alldat.lons,alldat.lats,alldat.julds,alldat.hfw,alldat.ahfw, ...
	alldat.wdeps] = ...
	fwelimoutliers(alldat.lons,alldat.lats,alldat.julds,alldat.hfw, ...
	alldat.ahfw,alldat.wdeps,hfwlims);

% select data by region
[fregions,fids,frlims,rlons,rlats,rjulds,rwdeps,rhfw,rahfw] = ...
	regionselect(alldat.lons,alldat.lats,alldat.julds,alldat.wdeps, ...
	alldat.hfw,alldat.ahfw,'Z');
%%% eliminate duplicate stations (using lat/lon/juld/maximum profile pressure)
%%%[idind,dind,rlons,rlats,rjulds] = ...
%%%	ctdelimdouble(rlons,rlats,rjulds,rmaxpres)
%%%% already done in 'fwregionload.m' !!!
%%%rwdeps = rwdeps(idind); rhfw = rhfw(idind); rahfw = rahfw(idind);

%	keyboard;

%%%if ~isempty(ggq), ggq = ['.' ggq]; end
%%%eval(['save fwregionts.ttemp' ggq testq 'mat -v7.3'])
eval(['save fwregionts.ttemp.' ffseason '.' int2str(hhq) '.' ...
	cdatachoice(1) '.' vname(1:end-1) ggq testq 'mat -v7.3'])

%		input('xxx')

% load and prepare climatology
%% mapped data
[alldat,juldlims,ewgfname,jsize] = ...
	climatload(juldlims,[67 90],sref,cdatachoice,dexc22);
alldat = alldat{1};
%% topography
%if strcmp(cdatachoice(1),'4')
%	load owntopoJAS.mat wdeps
%elseif strcmp(cdatachoice(1),'5')
%	load owntopoMAM.mat wdeps
%elseif strcmp(cdatachoice(1),'3')
%	if exist('phctopo.mat','file')
%		load phctopo.mat wdeps
%	end % if exist
%elseif strcmp(cdatachoice(1),'6') | strcmp(cdatachoice(1),'7')
%	load ewgtopo.mat wdeps
%else
%	disp('no other topo prepared...'); return
%end % if strcmp
%%
%if ~exist('wdeps','var') % no topo file -- get topo from database...
%	[wdeps,dummy1,dummy2] = gettopo(alldat.lons,alldat.lats);
%%	save phctopo.mat wdeps
%	clear dummy1 dummy2
%end
%alldat.wdeps = wdeps;

%eval(['save fwregionts.ttemp2.' int2str(hhq) '.' cdatachoice(1) '.' ...
%	vname(1:end-1) ggq testq '.mat'])
eval(['save fwregionts.ttemp2.' ffseason '.' int2str(hhq) '.' ...
	cdatachoice(1) '.' vname(1:end-1) ggq testq 'mat -v7.3'])

%%%		end % for ijk

disp(['loading fwregionts.ttemp2.' ffseason '.' int2str(hhq) '.' ...
	cdatachoice(1) '.' vname(1:end-1) ggq testq 'mat'])

eval(['load fwregionts.ttemp2.' ffseason '.' int2str(hhq) '.' ...
	cdatachoice(1) '.' vname(1:end-1) ggq testq 'mat'])

format compact
% mapping scales
% [large scale distance, PV; small scale distance, PV]
if ~exist('tscale','var')
	tscale = 60 % seasonal -- within +-2 month most weight!!!
end % if ~exist('tscale
%tscale = []
%%%txyscales = [600 0.2; 100 0.1] % annual averages
%%%txyscales2 = [1000 0.2; 100 0.1]
%%txyscales = [400 0.4 tscale; 50 0.05 tscale] % multi-year...
%%txyscales2 = [600 0.4 tscale; 50 0.05 tscale]
if ~exist('txyscales','var')
% used as final scales for mapping in Arctic FW paper
	txyscales = [600 1 tscale; 300 0.4 tscale] % multi-year...
	txyscales2 = [600 1 tscale; 300 0.4 tscale]
end % if ~exist('txysclaes
% used for testing effect of extrapolation-- TESTscales
%%%txyscales = [100 1 tscale; 50 0.4 tscale] % multi-year...
%%%txyscales2 = [100 1 tscale; 50 0.4 tscale]
% used for testing effect of extrapolation-- TEST2scales
%%%txyscales = [200 1 tscale; 100 0.4 tscale] % multi-year...
%%%txyscales2 = [200 1 tscale; 100 0.4 tscale]
%
%txyscales2 = [400 0.4; 50 0.05]
%txyscales2 = [300 0.4; 50 0.05]
% assign different scales for different regions...
xyscales(1) = {txyscales2};
xyscales(2:length(fids)) = {txyscales};
clear txyscales
% climatology grid subsampling (less dense close to north pole)
%	now use all points as equal area projection useful for integral...
cind = 1:length(alldat.lons);
%% current minimum grid distance is about 6km
%ttcind = cind(alldat.lats<87);
%%tcind = cind(alldat.lats>86&alldat.lats<87); ttcind = [ttcind tcind(1:2:end)];
%tcind = cind(alldat.lats>87&alldat.lats<88); ttcind = [ttcind tcind(1:2:end)];
%tcind = cind(alldat.lats>88&alldat.lats<89); ttcind = [ttcind tcind(1:4:end)];
%tcind = cind(alldat.lats>89); ttcind = [ttcind tcind(1:8:end)];
%cind = ttcind; clear ttcind tcind
% year / seasonal limits
%juldlims2{1} = [datenum(1991,1,1) ...
%	datenum(2000,1,1) 6 10]; fseason = '1991to1999JAS'; % JAS
%juldlims2{1} = [datenum(1992,1,1) ...
%	datenum(1994,1,1) 6 10]; fseason = 'testing1992to1993'; % JAS
%juldlims2{1} = [datenum(1992,1,1) ...
%	datenum(2000,1,1) 6 10]; fseason = '1992to1999JAS'; % JAS
%	datenum(2000,1,1) 2 6]; fseason = '1992to1999MAM'; % MAM
%	datenum(2000,1,1)]; fseason = '1992to1999'; % all year
%juldlims2{1} = [datenum(2003,1,1) ...%compare to Carmack et al.(2008;ASOF book)
%	datenum(2007,1,1)]; fseason = '2003to2006'; % all year
%juldlims2{1} = [datenum(2006,1,1) ...
%	datenum(2008,1,1) 6 10]; fseason = '2006to2007JAS'; % JAS
%juldlims2{1} = [datenum(2008,1,1) ...
%	datenum(2010,1,1) 6 10]; fseason = '2008to2009JAS'; % JAS
%juldlims2{1} = [datenum(2006,1,1) ...
%	datenum(2010,1,1) 6 10]; fseason = '2006to2009JAS'; % JAS
%juldlims2{1} = [datenum(2006,1,1) ...
%	datenum(2009,1,1) 6 10]; fseason = '2006to2008JAS'; % JAS
%	datenum(2010,1,1) 2 6]; fseason = '2005to2009MAM'; % MAM
%juldlims2{1} = [datenum(2001,1,1) ...
%	datenum(2006,1,1) 6 10]; fseason = '2001to2005JAS'; % JAS
%juldlims2{1} = [datenum(2001,1,1) ...
%	datenum(2009,1,1) 6 10]; fseason = '2001to2008JAS'; % JAS
%	datenum(2009,1,1) 2 6]; fseason = '2001to2008MAM'; % MAM
%	datenum(2008,1,1) 6 10]; fseason = '2001to2007JAS'; % JAS
%	datenum(2009,1,1)]; fseason = '2001to2008'; % all year
%juldlims2{1} = [datenum(1991,1,1) ...
%	datenum(2009,1,1)]; fseason = '1991to2008'; % JAS
%juldlims2{1} = [datenum(1991,1,1) ...
%	datenum(2009,1,1)]; fseason = '1991to2008all'; % JAS
for ii = 1:18
 for jj = 1:12
	%newjulds{(ii-1)*12+jj} = datenum(ii+1990,jj,15);
	%newjulds{(ii-1)*3+jj} = datenum(ii+1990,jj+6,15);
	%newjulds{jj} = datenum(2007,8,1);
 end % for jj
end % for ii
%newjulds = newjulds((ii-1)*3+[3]);
%juldlims2{1} = [datenum(1991,1,1) ...
%	datenum(2009,1,1) 6 10]; fseason = '1991to2008JAS'; % JAS
%	datenum(2010,1,1) 2 6]; fseason = '1991to2009MAM'; % MAM
%	datenum(2009,1,1)]; fseason = '1991to2008'; % all year
for ii = 1:12
	%juldlims2{ii} = [datenum(1991,1,1) ...
	%	datenum(2009,1,1) ii-1 ii+1];
	%fseason{ii} = ['1991to2009M' int2str(ii)];% M# denotes month of year
end % for ii = 1:12
for ii = 1:19 % 2001 to 2008 (2009 for MAM)
	%juldlims2{ii} = [datenum(1990+ii,1,1) ...
	%	datenum(1990+ii+1,1,1) 6 10]; fseason = 'yearlyavgJAS';
	%	datenum(1990+ii+1,1,1) 2 6]; fseason = 'yearlyavgMAM';
	%	datenum(1990+ii+1,1,1)]; fseason = 'yearlyavg'; % all year
	%juldlims2{ii} = [datenum(1990+ii,1,1) ...
	%	datenum(1990+ii+2,1,1)] ; fseason='twoyearlyavg' % all year
end % for ii = 1:18
juldlims2, fseason

%		fseason = ['TEST' fseason]


if iscell(fseason) & isempty(ffseason)
	ffseason = fseason{1}(1:11);
elseif ~iscell(fseason)
	ffseason = fseason;
end

%keyboard

%jjj = 2
%if exist('newjulds','var'), jjend = length(newjulds)+1, else, jjend = 2, end

rjulds=rjulds(:); rlons=rlons(:); rlats=rlats(:); rwdeps=rwdeps(:);
rhfw=rhfw(:); rahfw=rahfw(:);

%%%		for ijk = []

% This part is for parallel processing -- multiple time periods take long...
%clear;
%eval(['load ttemp.fwregiontsoutliers.' int2str(iiload) '.253.mat']); 
%%%ttempoutlierfname = fwregionscriptoutliertfilegen(rjulds,newjulds,juldlims2)
%%eval(['load ' ttempoutlierfname ' trahfw trhfw tridind tridind2 trjulds' ...
%%%	' trlats trlons trwdeps iiloads']);
iiloads = {length(juldlims2)*[0 1]+1}
%jjint = 48, jjit = 96 % for iiload <=3
%jjint = 6, jjit = 36 % for iiload>=25 & iiload<=34
%jjint = 2, jjit = 36 % for iiload>=97 & iiload <=102
jjit = 0
%jjint = 1
%%%if ~exist('ttemp.fwregiontsoutliers.status.mat','file')
	iiload = 1;
eval(['save ttemp.fwregiontsoutliers.status.' vname(1:end-1) '.mat iiload']);
%%%end % if ~exist
eval(['load ttemp.fwregiontsoutliers.status.' vname(1:end-1) '.mat iiload']);
iiload = iiload+1;
eval(['save ttemp.fwregiontsoutliers.status.' vname(1:end-1) '.mat iiload']);
%%iiload = 193
%%jjj = jj+1+jjint*(iiload-2), clear ii jj
%jjj = 2+jjint*(iiload-1), clear ii jj
%%jjj = 2;
%jjend = jjj+jjint-1

% select data by time intervals
disp('sort data into time intervals and eliminate outliers%%')
%for jj = 2:length(juldlims2)+1
for ii = 1:size(rjulds,1)
%%% if (iiload-1) > length(iiloads{ii}), disp('no more to process'), return, end
 jjj = iiloads{ii}(iiload-1)+1, jjend = iiloads{ii}(iiload)
%%% if jjend > length(newjulds)+1-jjit, disp('no more to process'), return, end
%%%	iiload, jjj, jjend, return
 for jj = jjj:jjend % mapping time periods
	disp('loop counters...'), ii, jj
	disp('of total '),length(juldlims2)+1, size(rjulds,1)
	[dummy,idind] = julselect(rjulds{ii,1},juldlims2{jj-1});
	tridind{ii,jj} = find(idind);
%	rlons{ii,jj}=rlons{ii,1}(idind); rlats{ii,jj}=rlats{ii,1}(idind);
%	rwdeps{ii,jj}=rwdeps{ii,1}(idind); rjulds{ii,jj}=rjulds{ii,1}(idind);
%	rhfw{ii,jj} = rhfw{ii,1}(idind); rahfw{ii,jj} = rahfw{ii,1}(idind);
	%% eliminate unrealistic FW content values by std deviation in anomalies
	[trlons{ii,jj},trlats{ii,jj},trjulds{ii,jj},trahfw{ii,jj}, ...
		trhfw{ii,jj},trwdeps{ii,jj},tridind2{ii,jj}] = ...
		fwelimoutliers2(rlons{ii,1}(idind), ...
		rlats{ii,1}(idind),rjulds{ii,1}(idind),rahfw{ii,1}(idind), ...
		rhfw{ii,1}(idind),rwdeps{ii,1}(idind));
%%%    if rem(jj,5) == 0
	eval(['save ttemp.fwregionts.procstatus22.' int2str(jjj) '.' ...
		int2str(jjend) '.' vname(1:end-1) '.mat']);
%%%    end % if rem(jj,5) == 0
 end % for jj
end % for ii
% %rlons=rlons(:,1:end-1); rlats=rlats(:,1:end-1); rjulds=rjulds(:,1:end-1);
% %rwdeps=rwdeps(:,1:end-1); rhfw=rhfw(:,1:end-1);
%end % for jj
clear idind

%%eval(['save ttemp.fwregiontsoutliers.' int2str(iiload) '.' ...
%eval(['save ttemp.fwregiontsoutliers.' int2str(jjj) '.' ...
%	int2str(jjend) '.mat']); return
%
%%%			end % for ijk
%
%load ttemp.fwregiontsoutliers.252.merged.mat

%keyboard
eval(['save fwregionts.ttemp22.' ffseason '.' int2str(hhq) '.' ...
	cdatachoice(1) '.' vname(1:end-1) ggq testq 'mat -v7.3'])

			end % for ijk

%	return % restart ollie job, as 96h not enough for all steps...

if ~isempty(ggq),
 if ggq(1) ~= '.', ggq = ['.' ggq]; end
end % if ~isempty

disp(['load fwregionts.ttemp22.' ffseason '.' int2str(hhq) '.' ...
	cdatachoice(1) '.' vname(1:end-1) testq ggq 'mat'])

eval(['load fwregionts.ttemp22.' ffseason '.' int2str(hhq) '.' ...
	cdatachoice(1) '.' vname(1:end-1) ggq testq 'mat'])

%	tscale = 120;
%	fseason = '1991to2008all2'; ffseason = fseason;
%	txyscales = [1000 1 tscale; 300 0.4 tscale] % multi-year...
%	txyscales2 = [1000 1 tscale; 300 0.4 tscale]
%	xyscales(1) = {txyscales2};
%	xyscales(2:length(fids)) = {txyscales};
%	for ii = 1:18
%	 for jj = 1:12
%		newjulds{(ii-1)*12+jj} = datenum(ii+1990,jj,15);
%		%%%newjulds{(ii-1)*3+jj} = datenum(ii+1990,jj+6,15);
%		%newjulds{jj} = datenum(2007,8,1);
%	 end % for jj
%	end % for ii
%	%%%newjulds = newjulds((ii-1)*3+[1 2 3]);

if size(alldat.hfw,1)>20
	alldat.hfw = alldat.hfw';
end % if size

% Load EWG grid for mapping with equal area around each grid point
%ewgdat = load('/csys/ocean2/brabe/climatologies/ewgclimatology.mean.summer.mat')
%ewglons = ewgdat.alldat.lons; ewglats = ewgdat.alldat.lats; clear ewgdat
load ewggrid.mat
cind = 1:length(ewglats);
cind = cind(ewglats>80 | (ewglats>68&(ewglons<=-90|ewglons>=90)) );
%% EWG topography
if exist('ewgtopo.mat','file')
	load ewgtopo.mat wdeps
else
	wdeps = [];
end % if exist
if length(wdeps) ~= length(ewglons)
	[wdeps,dummy1,dummy2] = gettopo(ewglons,ewglats);
	save ewgtopo.mat wdeps
	clear dummy1 dummy2
end % if length
ewgwdeps = wdeps;
cind = cind(ewgwdeps(cind)<0); % only ocean, not land ;)
% map data onto climatological data points
% for now ignore time, so get an unweighted time mean of all data!
% Use equal area projection grid for easier FW volume content / mean calculation

%%%			end % for ijk

%%%		return

% This part is for parallel processing -- multiple time periods take long...
%clear;
%eval(['load ttemp.fwregiontsoutliers.' int2str(iiload) '.253.mat']); 
sfac = 5+i;
%%%ttempmappingfname=fwregionscriptmappingtfilegen(trjulds,newjulds,juldlims2,sfac)
%%%eval(['load ' ttempmappingfname ' oahfw oahfwerr ohfw ohfwdatanos ohfwerr' ...
%%%	' ohfws2 iiloads alons alats']);
iiloads = {length(juldlims2)*[0 1]+1}
%jjint = 48, jjit = 96 % for iiload <=3
%jjint = 6, jjit = 36 % for iiload>=25 & iiload<=34
%jjint = 2, jjit = 36 % for iiload>=97 & iiload <=102
jjit = 0
%jjint = 1
%%%if ~exist('ttemp.fwregionts.status.mat','file')

	iiload = 1;
eval(['save ttemp.fwregionts.status.' vname(1:end-1) '.mat iiload']);
%%%end % if ~exist
eval(['load ttemp.fwregionts.status.' vname(1:end-1) '.mat iiload']);
iiload = iiload+1;
eval(['save ttemp.fwregionts.status.' vname(1:end-1) '.mat iiload']);

%%iiload = 193
%%jjj = jj+1+jjint*(iiload-2), clear ii jj
%jjj = 2+jjint*(iiload-1), clear ii jj
%%jjj = 2;
%jjend = jjj+jjint-1

% map data ...
if exist('newjulds','var')
 if length(newjulds)==size(trjulds,2)-1 % use several time-selected datasets
					% for mapping of each time (newjulds)!!!
  % jjend = length(newjulds)+1
  disp('mapping each time-selected data at one associated time (newjulds)')
  for ii = 1:size(trjulds,1) % selected regions
%%%   if (iiload-1) > length(iiloads{ii}), disp('no more to process'), return, end
   jjj = iiloads{ii}(iiload-1)+1, jjend = iiloads{ii}(iiload)
%%%   if jjend > length(newjulds)+1-jjit, disp('no more to process'), return, end
%%%	load ttemp.fwregionts.procstatus.2.253.hfw.mat
%%%		jjj = jj; clear jj
   for jj = jjj:jjend % mapping times
    disp('loop counters...'), ii, jj
    %% 2-stage mapping process all in mapping function!
    [dummy,dummy,alons{ii,jj},alats{ii,jj},ohfw(ii,jj),oahfw(ii,jj), ...
 	ohfwerr(ii,jj),oahfwerr(ii,jj), ...
 	ohfws2(ii,jj),dummy,ohfwdatanos(ii,jj)] = ...
 	climatano3(trhfw(ii,jj),trahfw(ii,jj),trlons(ii,jj),trlats(ii,jj), ...
 	zeros(size(ewglons))',zeros(size(ewglons))', ...
 	ewglons,ewglats,xyscales{ii}, ...
 	trwdeps(ii,jj),ewgwdeps,trjulds(ii,jj), ...
 	newjulds{jj-1}*ones(size(ewglons)),cind); % freshwater content
 %%%	eval(['save ttemp.fwregionts.0.' int2str(jjend) '.mat']); return
%    if rem(jj,4)== 0, save('proctemp.progress.mat','jj','vname','jjend'); end
		disp('rem(jj,5) == 0'), rem(jj,5) == 0, rem(jj,5)
    if rem(jj,5) == 0
	eval(['save ttemp.fwregionts.procstatus.' int2str(iiload) '.' ...
		int2str(jjend) '.' vname(1:end-1) '.mat']);
    end % if rem(jj,12) == 0
   end % for jj
  end % for ii
%  eval(['save ttemp.fwregionts.' int2str(iiload) '.' int2str(jjend) '.mat']);
  %% save temporary file from parallel processing for later merging
  eval(['save ttemp.fwregionts.' int2str(jjj) '.' ...
	int2str(jjend) '.' vname(1:end-1) '.mat']); %%%return
 else % use one or more time-selected datasets for mapping of
					% all times (newjulds)!!!
  % jjend = length(newjulds)+1
  disp('mapping each (regionally selected) dataset at all times (newjulds)')
  for ii = 1:size(trjulds,1)
%%%   if (iiload-1) > length(iiloads{ii}), disp('no more to process'), return, end
   jjj = iiloads{ii}(iiload-1)+1, jjend = iiloads{ii}(iiload)
%%%   if jjend > length(newjulds)+1-jjit, disp('no more to process'), return, end
   for jj = jjj:jjend
    disp('loop counters...'), ii, jj
    %% 2-stage mapping process all in mapping function!
    [dummy,dummy,alons{ii,jj},alats{ii,jj},ohfw(ii,jj),oahfw(ii,jj), ...
 	ohfwerr(ii,jj),oahfwerr(ii,jj), ...
 	ohfws2(ii,jj),dummy,ohfwdatanos(ii,jj)] = ...
 	climatano3(trhfw(ii,2),trahfw(ii,2),trlons(ii,2),trlats(ii,2), ...
 	zeros(size(ewglons))',zeros(size(ewglons))', ...
 	ewglons,ewglats,xyscales{ii}, ...
 	trwdeps(ii,2),ewgwdeps,trjulds(ii,2), ...
 	newjulds{jj-1}*ones(size(ewglons)),cind); % freshwater content
 %%%	eval(['save ttemp.fwregionts.0.' int2str(jjend) '.mat']); return
%    if rem(jj,4)== 0, save('proctemp.progress.mat','jj','vname','jjend'); end
    if rem(jj,5) == 0
	eval(['save ttemp.fwregionts.procstatus.' int2str(iiload) '.' ...
		int2str(jjend) '.' vname(1:end-1) '.mat']);
    end % if rem(jj,12) == 0
   end % for jj
  end % for ii
%  eval(['save ttemp.fwregionts.' int2str(iiload) '.' int2str(jjend) '.mat']);
  %% save temporary file from parallel processing for later merging
  eval(['save ttemp.fwregionts.' int2str(jjj) '.' ...
	int2str(jjend) '.' vname(1:end-1)  '.mat']); %%%return
 %%% return
 end % if length(newjulds)==size(rjulds,2)-1
else % use one or more time-selected datasets for mapping without time scale
 disp('mapping one or more selected datasets without timescale')
 for jj = jjj:length(juldlims2)+1
  for ii = 1:size(trjulds,1)
   %% 2-stage mapping process all in mapping function!
   [dummy,dummy,alons{ii,jj},alats{ii,jj},ohfw(ii,jj),oahfw(ii,jj), ...
	ohfwerr(ii,jj),oahfwerr(ii,jj), ...
	ohfws2(ii,jj),dummy,ohfwdatanos(ii,jj)] = ...
	climatano3(trhfw(ii,jj),trahfw(ii,jj),trlons(ii,jj),trlats(ii,jj), ...
	zeros(size(ewglons))',zeros(size(ewglons))', ...
	ewglons,ewglats,xyscales{ii}, ...
	trwdeps(ii,jj),ewgwdeps,cind); % freshwater content
  end % for ii
 end % for jj
end % if exist('newjulds


%			return
%
%%%			end % for ijk
%
%
%% load merged file from parallel processing
%%%%save tproc.tmpvar.mat ffseason, clear, load tproc.tmpvar.mat ffseason
%%%%eval(['load ttemp.fwregionts.' ffseason '.merged.mat']);
%%load ttemp.fwregionts.2.253.hfw.mat
%load ttemp.fwregionts.252.merged.mat

eval(['save fwregionts.ttemp3.' ffseason '.' int2str(hhq) '.' ...
	cdatachoice(1) '.' vname(1:end-1) ggq testq 'mat -v7.3'])
%save fwregionts.ttemp3.yearlyavg.mat
%save fwregionts.ttemp3.yearlyavgMAM.mat
%save fwregionts.ttemp3.yearlyavgJAS.mat

%%%			return

%%%		end % for ijk

disp(['load fwregionts.ttemp3.' ffseason '.' int2str(hhq) '.' ...
	cdatachoice(1) '.' vname(1:end-1) ggq testq 'mat'])

eval(['load fwregionts.ttemp3.' ffseason '.' int2str(hhq) '.' ...
	cdatachoice(1) '.' vname(1:end-1) ggq testq 'mat'])
%load fwregionts.ttemp3.yearlyavg.mat
%load fwregionts.ttemp3.yearlyavgMAM.mat
%load fwregionts.ttemp3.yearlyavgJAS.mat
	clear alldlons alldlats alldjulds alldwdeps alldohfw alldohfwerr
	clear alldoahfw alldoahfwerr alldohfws2 alldohfwdatanos

clear tcind ttcind

%%%		fids = 'ZAKE' % use different sub-regions for one mapping
%for jj = 2:length(juldlims2)+1
%for ii = 1:size(rjulds,1)
for ii = 1:length(fids)
 for jj = 2:size(ohfw,2)
  if length(fids) ~= size(ohfw,1)
	ij = 1; % use different sub-regions for one mapping
  else
	ij = ii; % use regions according to mapping (by region)
  end % if length
  % select data by region
  [tfregions{ii},tfids(ii),tfrlims{ii},alldlons{ii,jj},alldlats{ii,jj}, ...
	dummy,alldwdeps{ii,jj},alldohfw{ii,jj},alldohfwerr{ii,jj}, ...
	alldoahfw{ii,jj},alldoahfwerr{ii,jj}, ...
	alldohfws2{ii,jj},alldohfwdatanos{ii,jj}] = ...
	regionselect(ewglons(cind),ewglats(cind),zeros(size(cind)), ...
	ewgwdeps(cind),ohfw{ij,jj}(:),ohfwerr{ij,jj}(:),oahfw{ij,jj}(:), ...
	oahfwerr{ij,jj}(:),ohfws2{ij,jj}(:),ohfwdatanos{ij,jj}(:),fids(ii));
 end % for jj
end % for ii

% reassign variables if isohaline depth / layer thickness is wanted
if strcmp(vname,'ssdh_')
	alldossdh = alldohfw; alldoassdh = alldoahfw;
	alldossdhs2 = alldohfws2; alldossdhdatanos = alldohfwdatanos;
	alldossdherr = alldohfwerr; alldoassdherr = alldoahfwerr;
	ossdh = ohfw; ossdherr = ohfwerr; oassdh = oahfw; oassdherr = oahfwerr;
	rssdh = rhfw; rassdh = rahfw;
	clear alldohfw alldoahfw alldohfwerr alldoahfwerr ohfw ohfwerr ...
		oahfw oahfwerr rhfw rahfw alldohfws2 alldohfwdatanos
%			save ttemp.fwregionscript.1.mat, return
elseif strcmp(vname,'h_')
	alldoh = alldohfw; alldoah = alldoahfw;
	alldohs2 = alldohfws2; alldohdatanos = alldohfwdatanos;
	alldoherr = alldohfwerr; alldoaherr = alldoahfwerr;
	oh = ohfw; oherr = ohfwerr; oah = oahfw; oaherr = oahfwerr;
	rh = rhfw; rah = rahfw;
	clear alldohfw alldoahfw alldohfwerr alldoahfwerr ohfw ohfwerr ...
		oahfw oahfwerr rhfw rahfw alldohfws2 alldohfwdatanos
elseif strcmp(vname,'mld_')
	alldomld = alldohfw; alldoamld = alldoahfw;
	alldomlds2 = alldohfws2; alldomlddatanos = alldohfwdatanos;
	alldomlderr = alldohfwerr; alldoamlderr = alldoahfwerr;
	omld = ohfw; omlderr = ohfwerr; oamld = oahfw; oamlderr = oahfwerr;
	rmld = rhfw; ramld = rahfw;
	clear alldohfw alldoahfw alldohfwerr alldoahfwerr ohfw ohfwerr ...
		oahfw oahfwerr rhfw rahfw alldohfws2 alldohfwdatanos
%			save ttemp.fwregionscript.1.mat, return
end % if strcmp


%%%save fwregionts.ttemp4.tmpcheck2.mat
%
%	input(['fwregionts.ttemp4.' ffseason '.' int2str(hhq) '.' ...
%		cdatachoice(1) ggq testq 'mat'])

exist(['fwregionts.ttemp4.' ffseason '.' int2str(hhq) '.' cdatachoice(1) ...
        ggq testq 'mat'],'file')

if exist(['fwregionts.ttemp4.' ffseason '.' int2str(hhq) '.' cdatachoice(1) ...
	ggq testq 'mat'],'file') == 2
	eval(['save fwregionts.ttemp4.' ffseason '.' int2str(hhq) '.' ...
		cdatachoice(1) ggq testq 'mat -APPEND -v7.3'])
else
	eval(['save fwregionts.ttemp4.' ffseason '.' int2str(hhq) '.' ...
		cdatachoice(1) ggq testq 'mat -v7.3'])
end % if exist
%save fwregionts.ttemp4.yearlyavg.mat
%save fwregionts.ttemp4.MAMyearlyavg.mat
%save fwregionts.ttemp4.JASyearlyavg.mat
%save fwregionts.ttemp4.JAStwoyearlyavg.mat

%%%	end % for ijk

%disp(['load fwregionts.ttemp4.' ffseason '.' int2str(hhq) '.' cdatachoice(1) 'mat'])

% reset random number generator to default
randreset();

	return

%eval(['load fwregionts.ttemp4.' ffseason '.' int2str(hhq) '.' cdatachoice(1) 'mat'])
%load fwregionts.ttemp4.yearlyavg.mat
%load fwregionts.ttemp4.MAMyearlyavg.mat
%load fwregionts.ttemp4.JASyearlyavg.mat
%load fwregionts.ttemp4.twoyearlyavg.mat
%load fwregionts.ttemp4.1992to2011MM.3.6.2dbbins.qc.sref35.0..xy600and300km.mat
load fwregionts.ttemp4.1992to2011MM.seasonalcycle.3.6.2dbbins.qc.sref35.0..xy600and300km.mat

% calculate regional timeseries of FW content
% currently this is the median of the FW inventory (layer thickness of pure FW)
if exist('alldohfw','var')
	[hfwts,hfwerrts,ahfwts,ahfwerrts,mhfwts,mahfwts] = ...
		fwregionts(alldohfw,alldohfwerr,alldoahfw,alldoahfwerr);
end % if exist
if strcmp(vname,'h_')
	[hts,herrts,ahts,aherrts,mhts,mahts] = ...
		fwregionts(alldoh,alldoherr,alldoah,alldoaherr);
end % if strcmp
whos
%
