% fwregionscriptscript()
%

global tpath
%tpath = fullpathget('~/working/lfw/') % for use on ollie, with /tmp
tpath = fullpathget('~/work/lfw/') % for use on ollie, with /tmp
%%%		for ijk = []
clear
%%%juldlims2{1} = [datenum(1992,1,1) ...
%%%	datenum(2000,1,1) 6 10]; fseason = '1992to1999JAS'; % JAS
%%%ffseason = '1992to1999JAS'
%%%juldlims2{1} = [datenum(2006,1,1) ...
%%%	datenum(2009,1,1) 6 10]; fseason = '2006to2008JAS'; % JAS
%%%ffseason = '2006to2008JAS'
%%%juldlims2{1} = [datenum(2006,1,1) ...
%%%	datenum(2011,1,1) 6 10]; fseason = '2006to2011JAS'; % JAS
%%%ffseason = '2006to2011JAS'
%%%juldlims2{1} = [datenum(2009,1,1) ...
%%%	datenum(2012,1,1) 6 10]; fseason = '2009to2012JAS'; % JAS
%%%ffseason = '2009to2012JAS'
%%%ffseasons(end-2:end) = 'MAM', fseasons(end-2:end) = 'MAM'
%%%juldlims2{1}(end-1:end) = [2 7] % MAM
%%%for ii = 1:12
%%%	juldlims2{ii} = [datenum(2006,1,1) ...
%%%		datenum(2015,1,1) ii-1 ii+1];
%%%	fseason{ii} = ['2006to2014M' int2str(ii)];% M# denotes month of year
%%%end % for ii = 1:12
%%%ffseason = fseason{1}(1:11);
%%for ii = 1:12
%%	juldlims2{ii} = [datenum(2006,1,1) ...
%%		datenum(2015,1,1) ii-2 ii+2];
%%	fseason{ii} = ['2006to2014M' int2str(ii) ...
%%		'PM1'];% M# denotes month of year -- +- 1 month !!!
%%end % for ii = 1:12
%%ffseason = fseason{1}(1:14);
%%%jyears = 6
%%%for ii = 1:12 % months
%%%  for jj = 1:jyears % years
%%%	fseason{(ii-1)*jyears+jj} = [int2str(2006+(jj-1)) 'M' int2str(ii)];
%%%	%% M# denotes month of year
%%%	newjulds{(ii-1)*jyears+jj} = datenum(2006+(jj-1),ii,15);
%%%  end % for jj = 1:jyears
%%%end % for ii = 1:12
%%%juldlims2{1} = [datenum(2006,1,1) datenum(2012,1,1)];
%%%ffseason = '2006to2011MM'
%%%jyears = 20, ww = 3
%%%for jj = 1:jyears % years
%%%  for ii = 1:12 % months
%%%	fseason{(jj-1)*12+ii} = [int2str(1992+(jj-1)) 'M' int2str(ii)];
%%%	%% M# denotes month of year
%%%	newjulds{(jj-1)*12+ii} = datenum(1992+(jj-1),ii,15);
%%%	startyr = 1992+jj-ww-1; endyr = 1992+jj+ww-1;
%%%	startmnth = ii; endmnth = ii;
%%%	if startyr < 1992, startyr = 1992; endyr = startyr+ww*2; end
%%%	if endyr > 2011, endyr = 2011; startyr = endyr-ww*2; end
%%%	juldlims2{(jj-1)*12+ii} = ...
%%%		[datenum(startyr,endmnth,1) datenum(endyr,endmnth,1)];
%%%  end % for ii = 1:12
%%%end % for jj = 1:jyears
%%%ffseason = '1992to2011MMwindow2'
%%%%jyears = 2, ww = 3, syear = 2013 % only process 2013/2014 !!!
%%%%%%jyears = 24, ww = 1, syear = 1992
%jyears = 24, ww = 3, syear = 1992
%jyears = 39, ww = 3, syear = 1980
jyears = 27, ww = 3, syear = 1992
%%%jyears = 9, ww = 1, syear = 2006 % for mld
for jj = 1:jyears % years
  for ii = 1:12 % months
	fseason{(jj-1)*12+ii} = [int2str(syear+(jj-1)) 'M' int2str(ii)];
	%% M# denotes month of year
	newjulds{(jj-1)*12+ii} = datenum(syear+(jj-1),ii,15);
	startyr = syear+jj-ww-1; endyr = syear+jj+ww-1;
	startmnth = ii; endmnth = ii;
	if startyr < syear, startyr = syear; endyr = startyr+ww*2; end
	if endyr > 2018, endyr = 2018; startyr = endyr-ww*2; end
	juldlims2{(jj-1)*12+ii} = ...
		[datenum(startyr,endmnth,1) datenum(endyr,endmnth,1)];
  end % for ii = 1:12
end % for jj = 1:jyears
%%%ffseason = '2006to2014MMwindow5'
%%%ffseason = '1992to2015MMwindow5'
%ffseason = '1980to2018MMwindow4new'
ffseason = '1992to2018MMwindow4new'
%%%%ffseason = '2013to2014MMwindow4'
%%%fseason{1} = ['2006to20010MM'];% M# denotes month of year
%juldlims2{1} = [datenum(2006,1,1) ...
%	datenum(2011,1,1) 2 4]; fseason = '2006to2010Mar'; % monthly...
%%%ffseason = '2006to2010Mar'
%%%vname = 'hfw_', sref = 34.8, fsref = sref, layersal = sref
vname = 'hfw_', sref = 35, fsref = sref%, layersal = 34
%%%vname = 'h_', sref = 0, fsref = 35, layersal = 34
%%%vname = 'ssdh_', sref = -1, fsref = 35
%%%vname = 'mld_', sref = [{{[]}} {[35]}] % MLD testing !!!
%%%testq = '.negativehfwtest.'
%testq = '.randommaptest4.'
%%%testq = '.xy100and50km.'
%%%testq = '.xy400and200km.'
testq = '.xy600and300km.'
%%%testq = '.xy200and100km.'
%
tscale = 30 % 30 days -- within +-1/2 month most weight!!! seasonal cycle!!!
%%%txyscales = [100 1 tscale; 50 0.4 tscale] % for seasonal cycle
%%%txyscales2 = [100 1 tscale; 50 0.4 tscale]
%txyscales = [100 1; 50 0.4] % for MLD (avg. of time periods)
%txyscales2 = [100 1; 50 0.4]
% used for testing effect of extrapolation
%%%tscale = 120 % 30 days -- within +-1/2 month most weight!!!
%%%txyscales = [100 1 tscale; 50 0.4 tscale] % multi-year...
%%%txyscales2 = [100 1 tscale; 50 0.4 tscale]
%%%txyscales = [200 1 tscale; 100 0.4 tscale] % multi-year...
%%%txyscales2 = [200 1 tscale; 100 0.4 tscale]
%%%txyscales = [400 1 tscale; 200 0.4 tscale]
%%%txyscales2 = [400 1 tscale; 200 0.4 tscale]
txyscales = [600 1 tscale; 300 0.4 tscale]
txyscales2 = [600 1 tscale; 300 0.4 tscale]
%
%hfwlims = [-inf inf; -inf inf]; % simple outlier elimination FW inventories
			%% leave out for now, as second elimination by std. dev.
%% put back in, as extremes in EB and BG not elminated!
hfwlims = [0 40; 0 20]; % simple outlier elimination FW inventories
%%%hfwlims = [0 400; 0 200]; % simple outlier elimination interface depth
%%%hfwlims = [0 2; 0 1]; % simple outlier elimination for ssdh
%%%hfwlims = [0 100; 0 100]; % simple outlier elimination for mld
%
%%%plim = [-0.1 10.2]
%
%%%ii = [6:18 20:21 23:38] % all year
%%%ii = [26 3 6:20 21:23 25 27 30 32 33 31] % JAS -- with XCTD (2008) !!!
%%%ii = [45 6:18 20:21 23:24 25 27:30 33 31:32 34 35:44 46:47]% all year -- only skip WOD09 sel,
ii = [52]% all year -- AO_phys dataset (Ext. UDASH / Myriel Vredenborg)
%%%ii = [45 6:18 20:21 23:24 25 27:30 33 31:32 34 35:44 46:47 49:51]% all year NEW!
%%%ii = [45 6:18 20:21 23:24 25 27:30 33 31 34 35:44 46:47 49:51]% all year -- for MLD!!!
%	skip WOD09 sel as only JAS! currently use WOD09 with auto-qc
%	skip SCICEX 97/98 as minimum pressure > 20 dbar (i.e. not for MLD !!!)
%%%ii = [6:18 20:21 23:36] % all year -- for seasonal cycle
%%%	ii = ii(ii~=34)
%
pseason = 'mn91to08'
pseason = [{'all'} {pseason}]
cdatachoice = '9'
fwregionscript
		return

%clear
%%%		end


%%%		for ijk = []
juldlims2{1} = [datenum(1992,1,1) ...
	datenum(2000,1,1) 6 10]; fseason = '1992to1999JAS'; % JAS
ffseason = '1992to1999JAS'
vname = 'h_', sref = 0, fsref = 35
%%testq = '.xy200and100km.'
hfwlims = [0 450; 0 450]; % simple outlier elimination halocline layer depth
fwregionscript
clear
%%%		end
%%%		for ijk = []
juldlims2{1} = [datenum(2006,1,1) ...
	datenum(2009,1,1) 6 10]; fseason = '2006to2008JAS'; % JAS
ffseason = '2006to2008JAS'
vname = 'hfw_', sref = 35, fsref = sref
%testq = '.randommaptest4.'
%%%testq = '.xy100and50km.'
% used for testing effect of extrapolation
%%%tscale = 120 % 30 days -- within +-1/2 month most weight!!!
%%%txyscales = [100 1 tscale; 50 0.4 tscale] % multi-year...
%%%txyscales2 = [100 1 tscale; 50 0.4 tscale]
%%%txyscales = [200 1 tscale; 100 0.4 tscale] % multi-year...
%%%txyscales2 = [200 1 tscale; 100 0.4 tscale]
hfwlims = [-5 30; 0 100]; % simple outlier elimination FW inventories
%%%plim = [-0.1 10.2]
fwregionscript
clear
%%%		end
		for ijk = []
juldlims2{1} = [datenum(2006,1,1) ...
	datenum(2009,1,1) 6 10]; fseason = '2006to2008JAS'; % JAS
ffseason = '2006to2008JAS'
vname = 'h_', sref = 0, fsref = 35
%%testq = '.xy200and100km.'
hfwlims = [0 450; 0 450]; % simple outlier elimination halocline layer depth
fwregionscript
		end
		return
clear
vname = 'ssdh_', sref = -1, fsref = 35
hfwlims = [0 5; 0 5]; % simple outlier elimination sea surface dynamic height
fwregionscript
%
