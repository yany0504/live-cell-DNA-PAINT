function [d, h] = densityScatter(X,param)
arguments
    X double
    param.filterFactor double = 0.3
    param.plotting double = 1
    param.colMap double = hot(200)
    param.markerSize double = 20
end
%
% function [d, h] = densityScatter(X,param)
%
% Authors:
%   (c) Matthias Chung (e-mail: mcchung@vt.edu)               in April 2016
%
%
% MATLAB Version: 8.6.0.267246 (R2015b)
%
% Description:
%   Create scatter plot of data and generates vector of density value for
%   each column of X.
%
% Input arguments:
%   X              - initial guess (column vector)
%   #param         - further options of algorithm
%     filterFactor - filters of larger distances, if f = 0 no filtering [default f = 0.3]-        to line search algorithm
%     plotting     - plot scattered data on or off [default plotting = 1]
%     colMap       - send in colormap (default colMap = hot(200)
%     markerSize   - marker size of scattered data [default markerSize = 20]
%
% Output arguments:
%   d              - density value
%   h              - scatter plot object
%
%
% Examples:
%     
%     densityScatter(rand(3,1000))
%
%     densityScatter([(randn(3,4000)-3), randn(3,4000)],{'filterFactor',20});
% 
%     d = densityScatter([(randn(2,8000)-3), randn(2,8000)],{'filterFactor',20;'markerSize',50});
%
%     mu = [2,3,-1]; sigma = [1,-1.5, 0.5; -1.5,3, -1; 0.5, -1, 2];
%     [d, h] = densityScatter(mvnrnd(mu,sigma,10000)',{'colMap',hsv(300)});
%
% default values
plotting = param.plotting;
filterFactor = param.filterFactor;
colMap = param.colMap;
markerSize = param.markerSize;

% initialize
[m,n] = size(X); d = zeros(m,1);
h = waitbar(0,'Computing densities...');
for j = 1:n % loop over every point
  waitbar(j/n);
  D = sqrt(sum(bsxfun(@minus,X(:,j),X)'.^2,2));
  filter = 1-D./(D+filterFactor); % normalized filter
  d(j) = sum(filter.*D);          % generates filtered distance measure
end
close(h)
% normalize density measure
maxd = max(d); mind = min(d); d = (d-mind)/(maxd-mind);
if plotting % plotting option
  T = linspace(0,1,size(colMap,1));
  col = interp1(T, colMap, d);
  if     m == 3
    h = scatter3(X(1,:),X(2,:),X(3,:),markerSize,col,'filled');
  elseif m == 2
    h = scatter(X(1,:),X(2,:),markerSize,col,'filled');   
  else
    fprintf('Unable to plot in %d dimensions.\n',m)
  end
end
end