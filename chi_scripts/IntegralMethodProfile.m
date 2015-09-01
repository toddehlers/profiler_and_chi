% IntegralMethodProfile.m

% Integral method for stream profile analysis, following Perron and Royden
% (2012, ESPL, doi:10.1002/esp.3302). Plots elevation against chi, which is
% a drainage-area weighted upstream distance.  Steady-state profiles plot
% as a straight line.

% S. Miller, 5/2013
% Updated to remove hillslopes by B. Adams, 5/2015


clear all

%% Set parameters
movern=0.45;    % Observed or expected concavity
A0 = 10^6;      % Reference drainage area
Amin = 1e6;     % Minimum drainage area to plot


%% Load chandata file (i.e., output file from stream profiler with variables
% for elevation, area, and distance from mouth for an individual channel)
% and extract variables.

load F:\integral_method\moshannon_chandata.mat

%z = chandata(:,2);      % elevation (unsmoothed)
z = chandata(:,4);      % elevation (smoothed)
A = chandata(:,3);      % drainage area
dfm = chandata(:,7);    % upstream distance from stream mouth/outlet.

% Some chandata files show errors in first value of A.  Delete first values here.
z = z(2:length(z));
A = A(2:length(A));
dfm = dfm(2:length(dfm));

% Remove drainage areas smaller than Amin
z = z(A>=Amin);
A = A(A>=Amin);
dfm = dfm(A>=Amin);


%% Calculate average distance between adjacent points along profile
for i = 1:(length(dfm)-1)
    dx(i) = dfm(length(dfm)+1-i)-dfm(length(dfm)-i);
end
avgdx = mean(dx);


%% Calculate chi for various values of m/n
sumAterm = 0;
    for j = 1:length(dfm)
        %chi(k) = ((A0/A(k))^movern(i)) * avgdx;
        %chi(j) = ((A0/A(j))^movern) * avgdx + sumAterm;
        sumAterm = (A0/A(j))^movern + sumAterm;
        chi(j) = sumAterm * avgdx;
    end


%% Plot
figure(1)
plot(chi,z)
xlabel('chi (m)')
ylabel('z (m)')

