function [ida,islope,lbarea,lbslope,chandata, ans2] = sa_analysis51(chandata,movernset,name,arc_workdir,mat_workdir,rmspike,wind,no_step,ks_window,cont_intv)

% sa_analysis51.m contains most of the analysis of slope/area data
% extracted from the profile.  It calculates channel gradient on fixed 
% vertical interval.  It is written to work with dem data extracted from a channel
% profile.  It is associated with profile.m and profileplot.m.  The function
% requires sa_regress.m to preform the slope-area regression

% USAGE: 
%   [ida,islope,lbarea,lbslope,chandata,ans2] = sa_analysis51(chandata,movernset,name,arc_workdir,mat_workdir,rmspike,wind,no_step,ks_window,cont_intv)
% Output:
%     IDA, output vector of drainage areas interpolated at vertical interval
%     ISLOPE, output vector of gradients interpolated at vertical interval
%     LBAREA, log-binned drainage area array
%     LBSLOPE, log-binned slopes
%     CHANDATA, 10-column matrix of channel data (see profile.m)
%     ANS2, flag to save slope-area data and regression
% Input:
%     CHANDATA, input matrix of distance, elevation, drainage area, smoothed elevation and more (from profile.m)
%     MOVERNSET, value of channel concavity, usually -0.4 to -0.6
%     NAME, stream name, defined in profile.m
%     ARC_WORKDIR, ArcGIS workspace
%     MAT_WORKDIR, MATLAB workspace
%     RMSPIKE, 1 if spikes were removed, 0 if not. set in profile.m
%     WIND, smoothing window size, zero if no smoothing, else value comes from movavg.m
%     NO_STEP, flag for whether step-remover was used in profile.m
%     KS_WINDOW, size of window over which to calculate steepness index
%     CONT_INTV, contour interval
%  
% The step option refers to DEMs where you need to extract and use the raw 
% elevation and position data from the contour crossings to run analysis 
% independent of introduced artifacts, which we have found to be specific 
% to some DEMs (e.g. the LandcareResearch NZ 25m dem).

% NOTE:
% This code includes both step removing options as well as no-step
% removal.  In the step option, the code uses both the complete Chandata 
% and the filtered Chandata_ns where 'ns' is for 'no_step'.  The 'enter 
% contour interval' prompt is removed, the smoothed profile is 
% not plotted, and the contour crossing points are plotted with a circle.

global TRIBSECTION

dfd = chandata(:,1);
pelev = chandata(:,2);
drainarea = chandata(:,3);
smooth_pelev = chandata(:,4);
dfm = chandata(:,7);

if no_step,  %case where step remover is used.
    %%%%%%%%%%%INTERACTIVE PARAMETER--Comment out interactive, uncomment set value, if desired
    yesno_usgs = answer_yn('STEP REMOVER:  Enter "y" for USGS 10m DEM processing, "n" for user-defined step height (e.g. 20m for NZ): ');
    %yesno_usgs=1;  %1 to use usgscontour.m, 0 to use waipaoa vertical step
    %remover.
    
    if yesno_usgs,   
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [usgsind,cont_intv]=usgscontour(chandata);
        chandata_ns=chandata(usgsind,:);    %only temporarily is this correct; it then get changed for consistency so that chandata_ns = chandata.  No, now chandata_ns is removed, not used anywhere else.
        ida=chandata_ns(:,3)';
        idfd=chandata_ns(:,1)';
        dfd_ns = chandata_ns(:,1);
        pelev_ns = chandata_ns(:,2);
        ielev = chandata_ns(:,2)';    
        imax = length(ielev);
        idfm = max(dfd)-idfd;
%        chandata_ns = chandata;  %see above comment wrt chandata_ns.  Or remove all chandata_ns because it is now obsolete.
    else   
        no_step=2;
%        return
        [ind,cont_intv]=stepremover(chandata);
        chandata_ns=chandata(ind,:);    %only temporarily is this correct; it then get changed for consistency so that chandata_ns = chandata.  No, now chandata_ns is removed, not used anywhere else.
        ida=chandata_ns(:,3)';
        idfd=chandata_ns(:,1)';
        dfd_ns = chandata_ns(:,1);
        pelev_ns = chandata_ns(:,2);
        ielev = chandata_ns(:,2)';    
        imax = length(ielev);
        idfm = max(dfd)-idfd;
    end
 
    %calculate slope from STEP REMOVED elev data: 
    for a = 1:length(ielev)-1
        islope(a) = (ielev(a+1)-ielev(a))/(idfd(a)-idfd(a+1));
    end
    for a = imax
        islope(a) = (ielev(a)-ielev(a-1))/(idfd(a-1)-idfd(a));
    end
    islope = abs(islope);   %not sure this is good to have in here
    
%
% set variables for computing logbin averages, ida2, ielev2, etc
%
    ielev2 = ielev;
    ida2 = ida;
    idfd2 = idfd;               
    idfm2 = idfm;                       
    
else    %case where step remover is NOT used; just subsample/interpolate at entered increment:
    
    % Make contoured slope array from smoothed data.
    pmax = length(dfd);
    ielevmin = 10*ceil(.1*smooth_pelev(1));
    ielevmax = 10*floor(.1*smooth_pelev(pmax));
    telev = smooth_pelev;
    for a = 2:length(telev)
        if telev(a) <= telev(a-1)
            telev(a)=telev(a-1)+0.000001;
        end
    end
    ielev = ielevmin:cont_intv:ielevmax;
    imax=length(ielev);
    ida = interp1(telev,drainarea,ielev);
    idfd = interp1(telev,dfd,ielev);
    idfm = max(dfd)-idfd;  
    
    %calculate slope from SMOOTHED elev data: 
    for a = 1:length(ielev)-1
        islope(a) = (ielev(a+1)-ielev(a))/(idfd(a)-idfd(a+1));
    end
    for a = imax
        islope(a) = (ielev(a)-ielev(a-1))/(idfd(a-1)-idfd(a));
    end
    islope = abs(islope);   %not sure this is good to have in here
    
    dfm = chandata(:,7);
%
% Before exit loop make another contoured slope array from unsmoothed data for use
% below to determine logbin average slopes
%       
    pmax = length(dfd);
    ielevmin2 = 10*ceil(.1*pelev(1));
    ielevmax2 = 10*floor(.1*pelev(pmax));
    telev = pelev;
    
    for a = 2:length(telev)
        if telev(a) <= telev(a-1)
            telev(a)=telev(a-1)+0.000001;
        end
    end
    ielev2 = ielevmin2:cont_intv:ielevmax2;
    imax2=length(ielev2);
    ida2 = interp1(telev,drainarea,ielev2);
    idfd2 = interp1(telev,dfd,ielev2);
    idfm2 = max(idfd2)-idfd2;

%
% END IF LOOP for step remover vs. smoothing option loop
%
end

figure(3), clf
figure(2), clf

%
% Determine logbin average slopes, but from raw, unsmoothed data in both
% cases -- step removed or smoothed data, use ielev2, ida2, idfd2
%

%calculate slope for unsmoothed data: 

for a = 1:length(ielev2)-1
    islope2(a) = (ielev2(a+1)-ielev2(a))/(idfd2(a)-idfd2(a+1));
end
for a = imax2
    islope2(a) = (ielev2(a)-ielev2(a-1))/(idfd2(a-1)-idfd2(a));
end
islope2 = abs(islope2);   
%
temparea = zeros(1,40);
tempslope = zeros(1,40);
for a = 1:80
    logamin = 3+(a-1)*.1;
    logamax = logamin+.1;
    temparea(a) = 10^(logamin+.05);
    j = find(ida2 >= 10^logamin & ida2 < 10^logamax);
    k = size(nonzeros(j));
    if k(1) > 0
        tempslope(a) = 10^(mean(log10(islope2(j))));
    end
    %
end
%
% Extract non-zero elements
%
lbarea = temparea(find(tempslope));
lbslope = tempslope(find(tempslope));

%
% End Log-Bin Average Computation
%

%
% Begin ichi Calculation - for the Integral Method
%
istop = length(idfm);

ichi(1)=0;
for i = 2:istop
ichi(i)=ichi(i-1)+((ida(i).^(movernset)+ida(i-1).^(movernset))./2).*(idfm(i)-idfm(i-1));
end

ichi = ichi';			%transpose matrix

%
% End ichi Calculation - for the Integral Method
%

% convert dfd to km
dfd_km = dfd./1000;
dfm_km = dfm./1000;
dfm_km = chandata(:,7)./1000;
idfm_km = idfm./1000;

% Plot on new figure elevation vs. chi (Wiki Plot).
figure(3)
orient tall
hold off
subplot(3,1,1)
plot(ichi,ielev,'m')
hold on
set(gca,'Xdir','reverse')
%ax = axis;
title(name)
xlabel ('chi')
ylabel ('elevation (m)')
%axis ([ax(1) ax(2) 1e-4 1e0])

% Plot on new figure log slope vs. distance from mouth.
subplot(3,1,2)
hold off
semilogy(idfm_km,islope,'m+')
hold on
set(gca,'Xdir','reverse')
ax = axis;
xlabel ('distance from mouth (km)')
ylabel ('gradient')
axis ([ax(1) ax(2) 1e-4 1e0])

% Plot the longitudinal profile (elevation vs. distance from mouth).
figure(2)
orient tall
hold off
subplot(3,1,1)
plot(dfm_km,pelev,'g')
hold on
%   plot(dfm_ns_km,pelev_ns,'mo','MarkerSize',2)
set(gca,'Xdir','reverse')
%loop to plot either smoothed data, or step removed data
if no_step,      %case for plotting only chosen step points:
    plot(idfm./1000,ielev,'m');
else    %case for plotting smoothed data channel 
    plot(dfm_km,smooth_pelev,'m');
end
v = axis;

%title(name)
%Title of figure 2 has lots of info on processing routine:
if rmspike,
    str0 = ' spikes removed;';
else
    str0 = '';
end

if no_step == 1,
    str1 = ' USGS step removal;';
elseif no_step == 2,
    str1 = ' step removal;';
elseif no_step == 0,
    str1 = '';
end
if wind,  %case where smoothing done; no step removal
    str1 = ([name,';',str0,' smoothing window=',num2str(wind),'m; contour=',num2str(cont_intv),'m']);
else    %case with step removal most likely, no smoothing allowed
    str1 = ([name,';',str0,str1,' No smoothing; contour=',num2str(cont_intv),'m']);
end
title(str1)

xlabel ('distance from mouth (km)')
ylabel ('elevation (m)')

% Plot the drainage area vs. distance from mouth.
subplot(3,1,2)
hold off
semilogy(dfm_km,chandata(:,3),'m')
set(gca,'Xdir','reverse')
hold on
xlabel ('distance from mouth (km)')
ylabel ('drainage area (m^2)')
axis ([v(1) v(2) 1e5 1e11])

% Plot slope vs. drainage area in log space.
subplot(3,1,3)
hold off
loglog (ida,islope,'m+')
hold on
loglog (lbarea,lbslope,'rs','MarkerFaceColor','r','MarkerSize',3)
xlabel ('drainage area (m^2)')
ylabel ('gradient')
w = axis;
axis ([1e3 1e11 1e-4 1e0])
hold on    


SA_fits = zeros(1,14);    
answer3 = 1;
while answer3,
    
    %   first write out data to file for reading in Arcview.
    %   This is added to the base stream table
    %   call function sa_regress to run slope-area regressions
    
    [slp,range,ukn1,stats,movernset,ukn,uknhigh,uknlow,min_area,max_area,h2,h3,h4,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18,h19,h20,chandata] = sa_regress51(chandata,ida,islope,lbarea,lbslope,idfd,idfm,ielev,ichi,movernset,name,arc_workdir,mat_workdir,rmspike,wind,cont_intv,ks_window);
   % keyboard
    answer5 = answer_yn('Do you want to remember this fit (readable into ArcMap)?');
    if answer5,
        newdata = [slp,range,ukn1,stats(1),movernset,ukn,uknhigh,uknlow,min_area,max_area,rmspike,wind,cont_intv,no_step];
        SA_fits = [SA_fits; newdata];
        
        theta_ref = -1*movernset;
        toGisData3 = [theta_ref ukn uknhigh uknlow slp range ukn1 min_area max_area rmspike wind cont_intv];
%       toGisData3 = [slp range ukn1 ukn uknhigh uknlow min_area max_area theta_ref rmspike wind cont_intv]; 
        [v,d]=version;
            if str2num(v(1))==6
                dlmwrite([[arc_workdir,name],num2str(TRIBSECTION,0),'.txt'], toGisData3,' ');
            else
                dlmwrite([[arc_workdir,name],num2str(TRIBSECTION,0),'.txt'], toGisData3,'delimiter',' ','precision','%1.14g');
            end
%    
        TRIBSECTION = TRIBSECTION + 1;
        delete(h3)    %h3 is handle to the large text on the plots
    else
%        keyboard
        delete([h2,h3,h4,h7,h8,h9,h10,h12,h13,h14,h15,h17,h18,h19,h20])   %delete all of the fit lines on the plot--h# are handles (h11, h16 removed - not needed unless variable_color=1 in sa_regress
    end
    
    % ask if more fits are necessary.
    answer3 = answer_yn('Select another range to fit?');
end
%add ref concavity label to lower left of fig2, subplot(3,1,3)--SA
%assume axis goes from 1e-4 to 1e0
theta_ref = -1*movernset;
figure(2), subplot(3,1,3)
text(2e3,2e-4,sprintf('\\theta_r_e_f = %s',num2str(theta_ref)));

% *************************************************************************
% Optional loop to mark knickpoints from long profile on DEM.  Moved
% from profile3.m so that variables ida, islope, idfm can be used...
% *************************************************************************

knick_data_arc = zeros(1,min(size(chandata))+1);
knick_data = zeros(1,min(size(chandata))+1);
%%%%%%%%%%%INTERACTIVE PARAMETER--Comment out interactive, uncomment set value, if desired
knickpt = answer_yn('Mark points on long profile?');
%knickpt=0;   %1 to mark points on profile, 0 to not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while knickpt,
    disp('From which plot would you like to mark points? ')
    knick_opt = input('a) logS-logA (fig2 plot 3); b) long profile (fig2 plot1); c) logS-dfm (fig3 plot 2): ','s');
    while ~strcmp(knick_opt,'a') & ~strcmp(knick_opt,'b') & ~strcmp(knick_opt,'c'),
        %case where you didn't enter a, b, or c:
        disp('Nope, try again.  Enter either a, b, or c:')
        knick_opt = input('a) logS-logA (fig2 plot 3); b) long profile (fig2 plot1); c) logS-dfm (fig3 plot 2): ','s');
    end
    
    % A few needed variables regardless of choice
    dfm = chandata(:,7);
    dfm_km = dfm./1000;
    pelev = chandata(:,2);
    
    if knick_opt == 'a'
        disp('Select point on Slope-Area plot (fig2, plot 3)')
        figure (2)
        subplot(3,1,3)
        [area,slope] = ginput(1);
        % Find closest point to user input in ida vector
        resarea = abs(ida-area);
        [xx,newflag] = min(resarea);
        % Find closest point in the chandata array
        drainarea = chandata(:,3);
        resarea = abs(drainarea-area);
        [x,flag] = min(resarea);
    end
    
    if knick_opt == 'b'
        disp('Select point on longitudinal profile (fig2, plot 1)' )
        figure(2)
        subplot(3,1,1)        
        [dist,elev] = ginput(1);
        % Find closest point to user input in dfm vector
        dist = dist*1000;
        resdfm = abs(dfm-dist);
        [x,flag] = min(resdfm);
        % find closest point in the idfm array
        residfm = abs(idfm-dist);
        [xx,newflag] = min(residfm);     
    end
    
    if knick_opt == 'c'
        disp('Select point on logS-DFM plot (fig3, plot 2)')
        figure (3)
        subplot(3,1,2)
        [dist,slope] = ginput(1);
        % Find closest point to user input in ida vector
        dist = dist*1000;
        resdfm = abs(dfm-dist);
        [x,flag] = min(resdfm);
        % find closest point in the idfm array
        residfm = abs(idfm-dist);
        [xx,newflag] = min(residfm);
    end
    
    figure (2)
    % Mark point on long profile
    subplot(3,1,1)
    plot(dfm_km(flag,1),pelev(flag,1),'+')
    % Mark point on logS-logA plot
    subplot(3,1,3)
    loglog(ida(1,newflag),islope(1,newflag),'o') 
    figure (3)
    % Mark point on chi-elev plot
    subplot(3,1,1)
    %semilogy(idfm_km(1,newflag),islope(1,newflag),'o')
    plot(ichi(newflag,1),ielev(1,newflag),'o')
    % Mark point on logS-dist plot
    subplot(3,1,2)
    plot(idfm_km(1,newflag),islope(1,newflag),'o')

    % Classify your point, if desired.
    
    knick_name = '0';
    %%%%%%%%%%%INTERACTIVE PARAMETER--Comment out interactive, uncomment set value, if desired
    knick_classify = answer_yn('Classify this point?');
%    knick_classify = 0;  %1 to classify, 0 to not classify
    %%%%%%%%%%%%%%%%%%%%%%
    if knick_classify,
        disp('How do you want to classify it? ')
        knick_name = input('1) Major Knick; 2) Minor Knick; 3) Start of Steep Sect.; 4) End of Steep Sect.; 5) Other ? ','s');
        while ~strcmp(knick_name,'1') & ~strcmp(knick_name,'2') & ~strcmp(knick_name,'3') & ~strcmp(knick_name,'4') & ~strcmp(knick_name,'5'),
            disp('Enter a number from 1 through 5 only:')
            knick_name = input('1) Major Knick; 2) Minor Knick; 3) Start of Steep Sect.; 4) End of Steep Sect.; 5) Other ? ','s');
        end
    else
    end     

%   warning('problem?')  %put here because I'm having an odd warning message associated with the above lines
    
    % Find point in chandata matrix corresponding to row # found above,
    % append to knick_data file.
    knick_data = [knick_data; chandata(flag,:), str2num(knick_name)];
    knick_data_arc = [knick_data_arc; chandata(flag,:), str2num(knick_name)];

    % knick_data will have 11 columns of data for each point marked, the
    % first 10 are from chandata (all of it), and the eleventh is the classification
    % number; if you don't classify, the default value is zero

    % Repeat for another knickpoint?
    knickpt = answer_yn('Mark another point?');
    %
end 
% remove the first row of zeros from both data outputs
size_knick_data = size(knick_data);
knick_data = knick_data(2:size_knick_data(1),:);
size_knick_data_arc = size(knick_data_arc);
knick_data_arc = knick_data_arc(2:size_knick_data_arc(1),:);

ans2 = answer_yn(sprintf('Save slope-area fits, chandata, knickdata, and print figures?'));

check = size(knick_data);
if check(1,1) > 0 
    if ans2,
        eval ([' save ',mat_workdir,name,'_knick.mat knick_data -MAT'])
%        
        [v,d]=version;
        if str2num(v(1))==6
            dlmwrite([arc_workdir,name,'_knick.txt'], knick_data_arc, ' ');
        else
            dlmwrite([arc_workdir,name,'_knick.txt'], knick_data_arc,'delimiter',' ','precision','%1.14g');
        end
    end
end     
% End of the knickpoint loop

%   ask if the fit's data should be saved.
check = size(SA_fits);
if check(1,1) > 1  
    %answer4 = answer_yn(sprintf('Should these fits be saved (to %s)?',[name,'_SA_fits.mat']));
    %answer4=ans2;
    if ans2,
        eval ([' save ',mat_workdir,name,'_SA_fits.mat SA_fits -MAT'])
    end
end