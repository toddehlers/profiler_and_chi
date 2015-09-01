function [chandata] = auto_ks_calc(chandata,movernset,cont_intv,ks_window)
% 
% AUTO_KS_CALC.m provides the automated addition of ks values to the
% chandata matrix as an eight column.  This way channel steepness can be
% plotted as an overlay on long profiles or on network maps
% 
% USAGE
%     [chandata] = auto_ks_calc(chandata,movernset,cont_intv,ks_window)
% INPUT
%     chandata
%     movernset = theta or the reference concavity (needs negative sign)
%     cont_intv = vertical interval at which slopes are calculated (m)
% 	  ks_window = along-channel distance over which ks is calculated (km)
% OUTPUT
%     chandata
% 

% all control variables now imported, and t_interval renamed "ks_window"
%
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% Code extracted from SA_ANALYSIS51.m     (revised to use smoothed elev data)
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

%data are extracted from chandata to new variables
dfd = chandata(:,1);
pelev = chandata(:,2); %actual values pulled from DEM
drainarea = chandata(:,3);
smooth_pelev = chandata(:,4); %smoothed elevation values
dfm = chandata(:,7);
pmax = length(dfd);

%The raw channel data are not used for calculating ks values.
%Rather, the channel is interpolated between contour intervals.
%The interpolated variables begin with an "i"
ielevmin = 10*ceil(.1*smooth_pelev(1));
ielevmax = 10*floor(.1*smooth_pelev(pmax));
%telev is an elevation array in which elevation is always decreasing
%downstream.  
telev = smooth_pelev;
for a = 2:length(telev)
    if telev(a) <= telev(a-1)
        telev(a)=telev(a-1)+0.000001;
    end
end
%ielev is array of elevation points at each contour interval
%note that min(ielev) can be greater than min(pelev)
%and max(ielev) can be less than max(pelev)
ielev = ielevmin:cont_intv:ielevmax;
imax=length(ielev);
%below interprets drainage area values at each ielev value
%NMG I believe telev is used because otherwise there might not be a unique
%relationship between elevation and drainage area.
ida = interp1(telev,drainarea,ielev);
%below calculates a distance from divide for each ielev value
idfd = interp1(telev,dfd,ielev);
idfm = max(dfd)-idfd;

%calculate slope from contoured elev data:
%If the channel does not have much relief, it's possible that the
%contouring algorithm will make a channel with one point.  In this case, we ignore
%this channel.  Note that this can happen even if the total relief in the
%channel is greater than the contour interval, because of the floor and
%ceil calculations used to find the contour interval values.
%NMG made a change here to avoid calculating ks if the smoothed channel
%only has one point.  However, at the end of this code, the eighth column 
%of chandata (ks vals) will be assigned a value of -9999 to indicate that 
%no calculations were made.
if length(ielev)>1
    
    for a = 1:length(ielev)-1
        islope(a) = (ielev(a+1)-ielev(a))/(idfd(a)-idfd(a+1));
    end
    for a = imax
        islope(a) = (ielev(a)-ielev(a-1))/(idfd(a-1)-idfd(a));
    end

    
    islope = abs(islope);   %not sure this is good to have in here

    %interpolate slope values back to every point along channel, for
    %auto-ks
    ks_slope = interp1(ielev,islope,pelev,'linear','extrap');
    
    %below added by nmg to check that all ks_slope values are positive
    [minksslope,minslploc]=min(ks_slope);
    while minksslope<0 %enter loop if ks_slope has a negative value
        helper_elev=pelev(minslploc); %elevation at point with a negative slope

        %this code assumes that negative values will always be outisde the
        %elevation range of ielev, so find which extreme you are at in the
        %if statement below        
        if helper_elev<ielevmin
            %for consistency, set all points with elevations less than the
            %minimum contour to have the same ksn value as those at the
            %minimum contour point
            ksnhelper=islope(1)*(ida(1)^(-1*movernset)); 
            [helper_indeces]=find(pelev<ielevmin);
            for jj=1:length(helper_indeces)
                ks_slope(helper_indeces(jj))=ksnhelper*(drainarea(helper_indeces(jj))^movernset);
            end
            disp(sprintf('FYI - fixed negative interpolated slope values at lower end of channel'));
        else
            ksnhelper=islope(imax)*(ida(imax)^(-1*movernset));
            [helper_indeces]=find(pelev>ielevmax);
            for jj=1:length(helper_indeces)
                ks_slope(helper_indeces(jj))=ksnhelper*drainarea(helper_indeces(jj))^movernset;;
            end
            disp(sprintf('FYI - fixed negative interpolated slope values at upper end of channel'));
        end
        [minksslope,minslploc]=min(ks_slope); %check to make sure other 
        %end of data doesn't have negative values too.
    end     

    % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    % Code extracted from from SA_REGRESS51.m
    % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    %**************************************************************************
    % Routine to plot ks vs dfm along the profile.  User provides a moving window
    % length; steepness indices are calculated for midpoints of these windows along
    % the profile.
    % NOTE TO USER: If moving window is too small, it is possible that there
    % will be no data in the log(S) and log(A) arrays.  In this case, you will
    % get a "Divide by zero" warning.  Experiment with the best moving window
    % length to suit your needs.
    %**************************************************************************
    % t_interval = input('Enter moving window length (km) for plotting ks vs dfm: ');
    % while ischar(t_interval) | isempty(t_interval),   %attempt at error checking
    %     disp('Enter number! in Km!  Don''t make it to small or it''ll crash.')
    %     t_interval = input('Enter moving window length (km) for plotting ks vs dfm: ');
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
        [uu,vv] = min(diff1);                              % Note: first # will give values; second will give indices.
        [ww,yy] = min(diff2);
        [aa,bb] = min(diffmid);
        newdfm(x+1,1) = chandata(bb,7);
        midptsx(x+1,1) = chandata(bb,6);
        midptsy(x+1,1) = chandata(bb,5);
        min_ar = chandata(vv,3);
        max_ar = chandata(yy,3);

        % Select range of drainage area to fit, calc ks_n only for non-zero elements
        % Edited here to compute ks for ALL local slope values, to eliminate
        % huge number of zeroes, only auto-ks values affected
        %

        ind = find(drainarea>min_ar & drainarea<max_ar);
        k = size(nonzeros(ind));
        if k(1) > 0
            sslope = ks_slope(ind);
            sarea = drainarea(ind);
            sdfd = dfd(ind);
            selev = smooth_pelev(ind);
            % Take the log of slope and area, set up a two-column area matrix.
            logsarea=log(sarea');
            logsslope=log(sslope');
            areamatrix=[ones(size(logsslope)) logsarea];

            % Calculate steepness for each interval [(U/k)^(1/n)], based on
            % movernset.  Save value in matrix along with distance from mouth.
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
    
else %Added by NMG 
    %this is for channels in which there was only one point after smoothing
    numpts=length(chandata(:,1)); %number of points in this channel
    chandata(1:numpts,8)=-9999; %assign ks a no-data value
end