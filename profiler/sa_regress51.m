function [slp,range,ukn1,stats,movernset,ukn,uknhigh,uknlow,min_area,max_area,h2,h3,h4,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18,h19,h20,chandata] = sa_regress51(chandata,ida,islope,lbarea,lbslope,idfd,idfm,ielev,ichi,movernset,name,arc_workdir,mat_workdir,rmspike,wind,cont_intv,ks_window)

% sa_regress51.m performs the primary regression(s) of slope-area data 
% based on the channel parameters.  It is written to work with data from 
% profile51.m and associated functions, movavg51.m, sa_analysis51.m.  
% It takes a minimum and maximum drainage area to regress slope-area data.

% USAGE: 
%     [slp,range,ukn1,stats,movernset,ukn,uknhigh,uknlow,min_area,max_area,
%     h2,h3,h4,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18,h19,h20,chandata] = sa_regress51(chandata,ida,islope,lbarea,lbslope,idfd,idfm,ielev,ichi,movernset,name,arc_workdir,mat_workdir,rmspike,wind,cont_intv,ks_window)
%
% OUTPUT:
%    SLP, output concavity from slope-area data
%    RANGE, output error on concavity
%    UKN1, output intercept, ks
%    STATS, error on intercept
%    UKNHIGH, UKNLOW, error estimates on steepness index
%    MIN_AREA, MAX_AREA, minimum and maximum drainage area for regression
%    h#, handles for ArcGIS rendering
%    CHANDATA, primary array of channel data
% 
% INPUT 
%    DFD, vector of distance from divide
%    PELEV, elevation array
%    SMOOTH_PELEV, smoothed elevation array
%    DRAINAREA, array of drainage areas
%    IDA, sampled drainage areas (output of SA_ANALYSIS)
%    ISLOPE, sampled channel gradients (output of SA_ANALYSIS)
%    LBAREA, log-binned drainage area array
%    LBSLOPE, log-binned slopes
%    IDFD, sampled distance from divide array (output of SA_ANALYSIS)
%    IDFM, sampled distance from mouth array (output of SA_ANALYSIS)
%    IELEV, sampled elevation array (output of SA_ANALYSIS)
%    ICHI, sampled chi array (see profile for chi definition)
%    MOVERNSET, value of channel concavity, usually -0.4 to -0.6
%    NAME, stream name, defined in profile.m
%    ARC_WORKDIR, ArcGIS workspace
%    MAT_WORKDIR, MATLAB workspace
%    RMSPIKE, 1 if spikes were removed, 0 if not. set in profile.m
%    WIND, smoothing window size, zero if no smoothing, else value comes from movavg.m
%    KS_WINDOW, size of window over which to calculate steepness index
%    CONT_INTV, contour interval

% set global varible (tribsection).  This is used to number the
% files associated with the trib sections, which are read in ArcMap

% global TRIBSECTION

dfd = chandata(:,1);
pelev = chandata(:,2);
drainarea = chandata(:,3);
smooth_pelev = chandata(:,4);

disp('From which plot would you like to pick regression limits? (or enter "d" for manual input)') 
regress_opt = input('a) logS-logA (fig2 plot3); b) long profile (fig2 plot1); c) dist-log(gradient) (fig3 plot2):  ','s');

while ~strcmp(regress_opt,'a') & ~strcmp(regress_opt,'b') & ~strcmp(regress_opt,'c') & ~strcmp(regress_opt,'d'),
    %case where you didn't enter a, b, c, or d:
    disp('Follow the directions, you fool!  Enter either a, b, c or d:')
    regress_opt = input('a) logS-logA; b) long profile; c) dist-log(gradient); d) manual: ','s');
end

% OPTION 1: pick regression limits from logS-logA
if regress_opt == 'a'
    disp('Click on minimum THEN maximum bounds for drainage area from LOG(S)-LOG(A) PLOT (fig2, plot 3)')
    disp('Include at least 3 data points -- crosses on LOG(S)-LOG(A) PLOT')
    figure (2)
    subplot(3,1,3)
    [area,slope] = ginput(2);
    while area(1,1) >= area(2,1),       %loop to make sure min, max values entered in correct order
        disp('Follow the directions!')
        disp('Click on MINIMUM (left) then MAXIMUM (right) bounds for drainage area from LOG(S)-LOG(A) PLOT (fig2, plot 3)')
        disp('Include at least 3 data points -- crosses on LOG(S)-LOG(A) PLOT')
        [area,slope] = ginput(2);
    end
    min_area = area(1,1);
    max_area = area(2,1);
end

% OPTION 2: pick regression limits off Long Profile.  NOTE: also saves targi, targj 
if regress_opt == 'b'
    disp('Click on LEFT (max dfm) then RIGHT (min dfm) bounds for regression from STREAM PROFILE (fig2, plot 1)')
    disp('Regress bounds must include at least 3 data points -- crosses on LOG(S)-LOG(A) PLOT')
    figure (2)
    subplot(3,1,1)
    [distfm,elevat] = ginput(2);
    while distfm(2,1) >= distfm(1,1),       %loop to make sure max, min dfm values entered in correct order
        disp('Follow the directions!')
        disp('Click on LEFT THEN RIGHT bounds for regression from STREAM PROFILE (fig2, plot 1)')
        disp('Regress bounds must include at least 3 data points -- crosses on LOG(S)-LOG(A) PLOT')
        [distfm,elevat] = ginput(2);
    end
%keyboard
    distfm=distfm*1000;             % convert to meters
    dfm=chandata(:,7);
    mindistfm = abs(dfm-distfm(1,1));      % find closest point in chandata matrix
    maxdistfm = abs(dfm-distfm(2,1));
    
    [xx1,minflag] = min(mindistfm);
    [xx2,maxflag] = min(maxdistfm);
    
    min_area = chandata(minflag,3);
    max_area = chandata(maxflag,3);
    
    % Save coordinates of regression limits in x,y space.  Might be useful to
    % be able to export these to the ArcGIS project to see where things are.
    min_mapx = (chandata(minflag,6));
    min_mapy = (chandata(minflag,5));
    max_mapx = (chandata(maxflag,6));
    max_mapy = (chandata(maxflag,5));
end

% OPTION 3: pick regression limits from logS-dfm plot (NOTE: also
% saves targi, targj)
if regress_opt == 'c'
    disp('Click on LEFT (max dfm) THEN RIGHT (min dfm) bounds for regression from DIST-LOG(GRADIENT) PLOT (fig3, plot 2)')
    disp('Include at least 3 data points -- crosses on DIST-LOG(GRADIENT) PLOT')
    figure (3)
    subplot(3,1,2)
    [distfm,gradientpick] = ginput(2);
    while distfm(2,1) >= distfm(1,1),       %loop to make sure min, max values entered in correct order
        disp('Follow the directions!')
        disp('Click on LEFT THEN RIGHT bounds for regression from DIST-LOG(GRADIENT) PLOT (fig3, plot 2)');
        [distfm,gradientpick] = ginput(2);
    end
    distfm=distfm*1000;                     % convert to meters
    dfm=chandata(:,7);
    mindistfm = abs(dfm-distfm(1,1));      % find closest point in chandata matrix
    maxdistfm = abs(dfm-distfm(2,1));
    
    [xx1,minflag] = min(mindistfm);
    [xx2,maxflag] = min(maxdistfm);
    
    min_area = chandata(minflag,3);
    max_area = chandata(maxflag,3);
    
    % Save coordinates of regression limits in x,y space.  Might be useful to
    % be able to export these to the ArcGIS project to see where things are.
    min_mapx = (chandata(minflag,6));
    min_mapy = (chandata(minflag,5));
    max_mapx = (chandata(maxflag,6));
    max_mapy = (chandata(maxflag,5));   
end

% OPTION 4: Type in regression limits manually (the old fashioned way)
if regress_opt == 'd'
    min_area = input('Enter the minimum drainage area to consider...');
    max_area = input('Enter the maximum drainage area to consider...');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%start commented replotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % replot all original data
% % convert dfd to km
dfm_km = chandata(:,7)./1000;
idfm_km = idfm./1000;
% %
% % Plot on new figure log slope vs. log distance from mouth.
 figure(3)
% % Plot on new figure log slope vs. distance from mouth.
 subplot(3,1,2)

ax = axis;
% % Plot the longitudinal profile
 figure(2)
% hold off
 subplot(3,1,1)
v = axis;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Plot slope vs. drainage area in log space.
figure(2)
subplot(3,1,3)
w = axis;

% Select range of drainage area to fit.
ind = find(ida>min_area & ida<max_area);
sslope = islope(ind);
sarea = ida(ind);
sdfd = idfd(ind);
selev = ielev(ind);
schi = ichi(ind);
if max(size(selev)) <= 2
    disp('ERROR: Fewer than 2 contour crossings in regression: uncertainty undefined')
    disp('ERROR: Will not import into arcgis -- select new regression bounds')
end 
% Take the log of slope and area, set up a two-column area matrix.
logsarea=log(sarea');
logsslope=log(sslope');
areamatrix=[ones(size(logsslope)) logsarea];

% Calculate the least-squares linear fit and 95% confidence intervals.
[b,bint,r,rint,stats] = regress(logsslope,areamatrix,0.05);

% Calculate y2, the slopes on the best-fit line.
y2 = b(2).*logsarea + b(1);
exy2 = exp(y2);
ukn1 = exp(b(1));

% Calculate the error on the m/n.
range= (bint(2)-bint(4))/(-2);

% Calculate the m/n (slp).
slp= b(2)*(-1);
movernvect = [slp range -1*range];

%
% Now regress on chi-elev data, for accurate ksn (for theta_ref) and for
% uncertainty bounds
%
% set up a two-column chi matrix.
		chimatrix=[ones(size(selev')) schi];

	% Calculate the least-squares linear fit and 95% confidence intervals.
		[bc,bintc,rc,rintc,statsc] = regress(selev',chimatrix,0.05);
        
	% Calculate ks, the slope of the best-fit line.
        zo_c = bc(1);
        ks_chi = bc(2);
        ukn = ks_chi;
        statsc(1);
        logukn = log(ukn);
        logyhat = logukn+(movernset*logsarea);
        yhat = exp(logyhat);
        uknvect = [ukn1 ukn];

	% Calculate the error on ks.
		range_chi = (bintc(2)-bintc(4))/(-2);
        uknlow = bintc(2);
        uknhigh = bintc(4);
        
    % Calculate the fit segment for later plotting
        calc_elev = ks_chi.*schi + zo_c;
       
	% calculate errors on residuals and plot
        figure(4);
		r_err = rc - rintc(:,1);
		%subplot(4,1,4)
		errorbar(rc,r_err,'-om')
		w1 = axis;
		%axis ([0 w1(2) -3 3])

% Plot the best fits to the slope-area data, add text to the plot.
figure(2)
subplot(3,1,3)
wrange = abs(w(3)-w(4));
%put in variable_color to have random colors on plots for fits; I think easier to read plots with colors always the same.
variable_color=0;   %change to 1 to have random colors generated each time!
if variable_color,      %try random colors for the fit lines, to be able to differentiate one from another:
    color1=rand(1,3).*0.9;
    color2=rand(1,3).*0.9;
    h8=loglog(sarea,exy2,'color',color1);
    h9=loglog(sarea,yhat,'--','color',color2);
else
    h8=loglog([sarea(1) sarea(length(sarea))],[exy2(1) exy2(length(exy2))],'b.-');
    h9=loglog([sarea(1) sarea(length(sarea))],[yhat(1) yhat(length(yhat))],'c.-');
end

%place theta, ksn on this plot, above the fit section:
minareapt=min(sarea);
maxareapt=max(sarea);
xloc=exp((log(minareapt)+log(maxareapt))/2);
%placement of auto text, above line first, then below if it will be > 10^0
%keyboard
yloc1=exp(b(2).*log(xloc) + b(1) + 1.2);
yloc2=exp(b(2).*log(xloc) + b(1) + 2);
if yloc2 > 0.4,
    yloc1=exp(b(2).*log(xloc) + b(1) - 1.2);
    yloc2=exp(b(2).*log(xloc) + b(1) - 2);
end
h13=text(xloc,yloc2,['\theta=',num2str(slp,2),' \pm ',num2str(range,2)],'fontsize',10,'horizontalalignment','center');
h14=text(xloc,yloc1,['k_s_n=',num2str(ukn,3)],'fontsize',10,'horizontalalignment','center');

disp('Click the upper left corner to locate temporary parameter text.')
h3=gtext({['\theta=',num2str(slp,2),' \pm ',num2str(range,2)],['k_s= ',num2str(ukn1,3)],['R^2: ',num2str(stats(1),2)],['\theta_r_e_f=',num2str(-1*movernset),' k_s_n=',num2str(ukn,3),' (',num2str(uknlow,3),'-',num2str(uknhigh,3),')'],['Fit betw ',num2str(min_area,2)],['   and ',num2str(max_area,2)]},'fontsize',10);

% Compute the steady-state profile, based on derived ukn, and m/n (b(2)).
% Slope (Se) and elevation (zss):
Se = ukn1*sarea.^b(2);
[i,j] = find(ida==sarea(1));
zss(1) = ielev(i,j);
for a = 2:length(sarea)
    zss(a) = Se(a-1).*(sdfd(a-1)-sdfd(a))+zss(a-1);
end

% Compute the steady-state profile, based on forced ukn, and m/n.
% Slope (Se) and elevation (zss):
fSe = ukn*sarea.^movernset;
[i,j] = find(ida==sarea(1));
zfss(1) = ielev(i,j);
for a = 2:length(sarea)
    zfss(a) = fSe(a-1).*(sdfd(a-1)-sdfd(a))+zfss(a-1);
end        

% Plot the steady-state and forced profiles in blue and green.
% convert dfd to km
sdfd_km = sdfd./1000;
sdfm = max(dfm_km)-sdfd_km;
figure(3)
hold on
    subplot(3,1,1)
    if variable_color,
        h10=plot(schi,calc_elev,'color',color1);
%             h11=loglog(sdfm,fSe,'--','color',color2);
        h11=1;
        h15=1;
        h16=1;
	else
        h10=plot([schi(1) schi(length(schi))],[calc_elev(1) calc_elev(length(calc_elev))],'b.');
%             h11=loglog([sdfm(1) sdfm(length(sdfm))],[fSe(1) Se(length(fSe))],'c.');
        h11=1;
        h15=plot(schi,calc_elev,'b');
%             h16=loglog(sdfm,fSe,'c-');
        h16=1;
    end
%
subplot(3,1,2)
if variable_color,
    h2=semilogy(sdfm,Se,'color',color1);
    h4=semilogy(sdfm,fSe,'--','color',color2);
    h17=0;
    h18=0;
else 
    h2=semilogy([sdfm(1) sdfm(length(sdfm))],[Se(1) Se(length(Se))],'b.');
    h4=semilogy([sdfm(1) sdfm(length(sdfm))],[fSe(1) Se(length(fSe))],'c.');
    h17=semilogy(sdfm,Se,'b-');
    h18=semilogy(sdfm,fSe,'c-');    
end
figure (2)
hold on
subplot(3,1,1)
if variable_color,
    h12=plot(sdfm,zss,'color',color1);
    h7=plot(sdfm,zfss,'--','color',color2);
    h19=0;  %just need these to export to not get errors
    h20=0;
else
    h12=plot([sdfm(1) sdfm(length(sdfm))],[zss(1) zss(length(zss))],'b.');
    h7=plot([sdfm(1) sdfm(length(sdfm))],[zfss(1) zfss(length(zfss))],'c.');
    h20=plot(sdfm,zss,'b-');
    h19=plot(sdfm,zfss,'c-');    
end

%**************************************************************************
% Routine to plot ks vs dfm along the profile.  User provides a moving window 
% length; steepness indices are calculated for midpoints of these windows along
% the profile.  Figure 3 subplot (3,1,3) plots ks vs distance from mouth.
% NOTE TO USER: If moving window is too small, it is possible that there
% will be no data in the log(S) and log(A) arrays.  In this case, you will
% get a "Divide by zero" warning.  Experiment with the best moving window
% length to suit your needs.
%**************************************************************************
% ks_window = input('Enter moving window length (km) for plotting ks vs dfm: ');
% while ischar(ks_window) | isempty(ks_window),   %attempt at error checking
%     disp('Enter number! in Km!  Don''t make it to small or it''ll crash.')
%     ks_window = input('Enter moving window length (km) for plotting ks vs dfm: ');
% end

interval = ks_window*1000;     % convert to meters
halfint = 0.5*interval;         % divide in half
streamlength = max(chandata(:,1));
help = round(streamlength/halfint);

% find midpoints of regression segments; save locations in new matrix
% called 'midpts', whose dimensions depend on the length of the window
% and stream.  These will eventually give dfd locations to plot
% against calculated ks values for that segment.

%diff = zeros(help,3);
newdfm = zeros(help,1);
newks = zeros(help,1);
midptsx = zeros(help,1);
midptsy = zeros(help,1);

for x = 0:help
    diff1 = abs(chandata(:,1)-(x*halfint));            % find index of lower regression limit
    diff2 = abs(chandata(:,1)-((x+2)*halfint));        % find index of upper regression limit
    diffmid = abs(chandata(:,1)-((x+1)*halfint));      % find index of regression midpoint
    [uu,vv] = min(diff1);              % Note: first # will give values; second will give indices.
    [ww,yy] = min(diff2);
    [aa,bb] = min(diffmid);
    newdfm(x+1,1) = chandata(bb,7);
    midptsx(x+1,1) = chandata(bb,6);
    midptsy(x+1,1) = chandata(bb,5);
    min_ar = chandata(vv,3);
    max_ar = chandata(yy,3);
    
    % Select range of drainage area to fit, calc ks_n only for non-zero elements
    ind = find(ida>min_ar & ida<max_ar);
    k = size(nonzeros(ind));
    if k(1) > 0
    sslope = islope(ind);
    sarea = ida(ind);
    sdfd = idfd(ind);
    selev = ielev(ind);
    % Take the log of slope and area, set up a two-column area matrix.
    logsarea=log(sarea');
    logsslope=log(sslope');
    areamatrix=[ones(size(logsslope)) logsarea];
    
    % Calculate steepness for each interval [(U/k)^(1/n)], based on
    % movernset.  Save value in matrix along with distance from divide.
    logsareamean = mean(logsarea);
    logsslopemean = mean(logsslope);
    logukn = (logsslopemean-movernset*logsareamean);
    uknstream = exp(logukn);

    chandata(yy:vv,8) = uknstream;
    newks(x+1,1) = uknstream;
    end

end  
%
% Extract non-zero elements for plotting
%
newdfm = newdfm(find(newks));
newks = newks(find(newks));

figure(3)
subplot(3,1,3)
hold off
plot(newdfm/1000,newks,'o')
hold on
set(gca,'Xdir','reverse')
xlabel('Distance from mouth (km)')
ylabel('Averaged ks')
% Use axis limits from subplot (3,1,2)
try
    axis ([ax(1) ax(2) 0 max(newks)])
catch
    warning('problem with plotting fig3, plot 3.')
end
%********************************************************************
% END of ks vs DFD plotting subroutine
%********************************************************************




