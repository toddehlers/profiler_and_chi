function[atargi,atargj,atargx,atargy,label] = auto_chn_finder(fname,arc_workdir,mat_workdir,accum,network,easting,northing,crit_pix)
% 
% AUTO_CHN_FINDER.m is a simple routine to efficiently find all channel 
% heads at a specified minimum critical drainage area (crit_area).  
%
% This subroutine now called directly within profile_batch42.m
% 
% USAGE
%     [atargi,atargj,atargx,atargy,label] = auto_chn_finder(fname,arc_workdir,mat_workdir,accum,network,easting,northing,crit_pix)
% INPUT
%     fname - SEE profile.m
%     arc_workdir - SEE profile.m
%     mat_workdir - SEE profile.m
%     accum - accumulation array
%     network - zeros matrix the size of accum
%     easting - x location of DEM corner from .meta file
%     northing - y location of DEM corner from .meta file
%     crit_pix - critical pixel for critical drainage area calc 
% OUTPUT
%    -- An array of "locations" [atargi atargj atargx atargy label] where 
%    label is a number from 1 to n (number of channel heads found)
% 

% find all A > A_crit, create matrix with 1's along network, zero elsewhere
    [m,n] = size(accum);
%    network = zeros(m,n); (now passed in from outside)
    [r,s] = find(accum > crit_pix);
    network(find(accum > crit_pix)) = 1;
    network(find(accum == -9999)) = 2;
%
% index only channel elements with accum < 4*crit_pix (for efficiency)
    [i,j] = find(accum > crit_pix & accum < 4*crit_pix);
    q = length(i);
%
% trim channel elements on boundaries to avoid looking outside matrix below
    t = zeros(q,1);
    u = zeros(q,1);
    t(find(i>1&i<m))=1;
    u(find(j>1&j<n))=1;
% zero out all elements where either i or j is on a boundary
    v = t.*i;
    w = t.*j;
    vv = u.*v;
    ww = u.*w;
    ii = vv(find(vv));
    jj = ww(find(ww));
    qq = length(ii);
%
% calc 3x3 mean for all points on network (except boundaries)
%    netave = zeros(m,n);
% edited to eliminate another huge array in netave!
    chan_i = [];
    chan_j = [];
    for z = 1:qq
        k = ii(z);
        l = jj(z);
%        netave(k,l) = mean(mean(network(k-1:k+1,l-1:l+1)));
        if mean(mean(network(k-1:k+1,l-1:l+1)))==2/9
            chan_i = [chan_i;k];
            chan_j = [chan_j;l];
        end
    end
%
% create geographic coordinates for these i,j locations
    atargi = chan_i;
    atargj = chan_j;
    atargx = northing(chan_i);
    atargy = easting(chan_j);
%
% write out to location_ij.txt, for stream name simply number from 1 to p
    p = length(chan_i);
    label = [1:p];
    label = label';
%    
%   output in format of location_ij.txt [row column easting northing name]
%
    locations  = [atargi atargj atargy atargx label];
%
    [v,d]=version;
    if str2num(v(1))==6
        dlmwrite([arc_workdir,fname,'location_ij.txt'], locations,' ');
    else
        dlmwrite([arc_workdir,fname,'location_ij.txt'], locations,'delimiter',' ','precision','%1.14g');
    end
%

end
    
    
