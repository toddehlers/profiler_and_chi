function [ind,cont_intv]=stepremover(chandata)

% stepremover.m  is a generic version of usgscontour that allows extraction 
% of data from any stepped DEM.  It takes DEMs that have contour stair 
% steps and take indices of chandata only at points where contours cross.  
% The output is index numbers for rows in chandata; also the contour 
% interval used
% 
% Originally Made by Joel Johnson for 10m DEMs of Utah based on USGS quad maps.
% 
% usgscontour.m calls the following functions: closeto.m
%
% Usage:
% 	[ind,cont_intv]=stepremover(chandata)
% Output:
% 	ind - index rows
% 	cont_intv - contour interval
% Input:
% 	chandata

%optional interactive mode, to PLOT, check and change indices:  
interactive=1;
%interactive=0;

warning off MATLAB:divideByZero

if interactive,
    %loop to enter contour interval:
    cont_intv=input('Enter contour interval of steps to be removed (m):  ');
    while isempty(cont_intv),
        cont_intv=('Enter numeric contour interval! (m):  ')
    end
else    %if not interactive, just set cont_intv below, in meters!
    cont_intv=20;
end

minelev=min(chandata(:,2));
maxelev=max(chandata(:,2));
%assume min contour isn't less than zero (below sea level): cint or contour intervals:
%these are where the contour steps should be; ie these are the elevations I want data at:
mincont=0;
disp(sprintf('stepremover.m:  Contour interval %sm, starting contour %sm.',num2str(cont_intv),num2str(mincont)));
cint=mincont:cont_intv:maxelev;
%only take pts within the range of minelev to maxelev:
i=find(cint<minelev);
i=max(i);
cint=cint(i+1:length(cint));

%find the chandata indices within a range of each cint value, and pick the highest one:
% this assumes that the correct elevation is at the bottom of each stairstep
%pm (for plus or minus) is the range in meters that pts will be considered.
%pm=0.1;
pm=0.5;

ind=[];
%keyboard
for i=1:length(cint),

%first look for points with exactly the correct value:
    indrange=find(chandata(:,2)==cint(i));
    
    if isempty(indrange),
        indrange=find(chandata(:,2)>cint(i)-pm & chandata(:,2)<cint(i)+pm);
    end
    if isempty(indrange),
%       ind(i)=min(find(chandata(:,2)>cint(i)-pm))
        ind(i)=closeto(cint(i),chandata(:,2));
    else
        ind(i)=max(indrange);
    end
end

%Remove points in ind that are the same:  this happens if contours are cut out at a cliff or knickpt.
df=diff(ind);
df=[df df(length(df))];
notsameind=find(df);
ind=ind(notsameind);

%ind=indcorrected;
indcorrected=ind;

if interactive,
%    figure(110), clf
    figure(4), clf
    plot(chandata(:,7),chandata(:,2));
    hold on
    plot(chandata(ind,7),chandata(ind,2),'g.-');
    plot(chandata(indcorrected,7),chandata(indcorrected,2),'r.-');
    xlabel('distance from mouth, m');
    ylabel('elevation, m');
    title('Profile with picked usgscontour.m steps.  Red points will be used')
    disp(' ')
    disp('Check to see that stepremover.m did a good job (figure 4).')
    disp('green line:  uncorrected, original usgscontour picks')
    disp('red line:  "corrected" contour crossings, all should coincide with a local step base')
    disp(' ')
end

%keyboard

warning on MATLAB:divideByZero