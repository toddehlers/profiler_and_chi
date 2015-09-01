% IntegralMethodProfile_findmovern.m

% Integral method for stream profile analysis, following Perron and Royden
% (2012, ESPL, doi:10.1002/esp.3302). Plots elevation against chi, which is
% a drainage-area weighted upstream distance.  Steady-state profiles plot
% as a straight line.  Iteratively finds best-fitting m/n.

% S. Miller, 5/2013
% Updated to remove hillslopes by B. Adams, 5/2015

clear all

%% Set parameters
A0 = 10^6;          % Reference drainage area (m^2)
Amin = 1e6;         % Minimum drainage area to plot
Aquery = 0;         % Exclude drainage areas above a certain value? 1 = yes, 0 = no.
Amax = 3e10;        % Maximum drainage area (m^2) used, if query == 1.


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

% Remove drainage areas greater than Amax, if this option (Aquery) was
% selected.
if Aquery == 1
    z = z(A<=Amax);
    A = A(A<=Amax);
    dfm = dfm(A<=Amax);
end


%% Calculate average distance between adjacent points along profile
for i = 1:(length(dfm)-1)
    dx(i) = dfm(length(dfm)+1-i)-dfm(length(dfm)-i);
end
avgdx = mean(dx);


%% Calculate chi for various values of m/n

movern = linspace(0,1,101);
for i = 1:length(movern)
    sumAterm = 0;
    for j = 1:length(dfm)
        sumAterm = (A0/A(j)).^movern(i) + sumAterm;
        chi(j) = sumAterm * avgdx;
    end
    chi=chi';
    X=[ones(length(chi),1) chi];
    [b,bint,r,rint,stats] = regress(z,X);
    R2(i)=stats(1);
    chimatrix(i,:)=chi;
    clear chi
end

bestindex = find(R2 == max(R2));

%% Plot

figure(1)
plot(chimatrix(bestindex,:),z)
xlabel('chi (m)')
ylabel('z (m)')

figure(2)
plot(movern,R2)
hold on
plot([movern(bestindex),movern(bestindex)],[0,1],'r')
xlabel('m/n')
ylabel('R^2')
text(movern(bestindex)+0.05,0.2,strcat('R^2 = ', num2str(R2(bestindex))));
text(movern(bestindex)+0.05,0.1,strcat('m/n = ', num2str(movern(bestindex))));

