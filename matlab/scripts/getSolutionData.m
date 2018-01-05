function [ data ] = getSolutionData( xfile, ufile, qfile )
%GETSOLUTIONDATA Summary of this function goes here
%   Detailed explanation goes here

data=[];
data.x = importdata(xfile,' ');
tmp    = importdata(ufile,' ');
data.t = tmp(:,1);
data.u = tmp(:,2:end);
tmp    = importdata(qfile,' ');
data.q = tmp(:,2:end);

end

