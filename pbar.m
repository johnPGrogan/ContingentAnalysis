function h = pbar(pVals, varargin)
% h = pbar(pVals, varargin)
% plot a bar where p values are less than alpha
% Inputs: pVals is a vector of pvalues e.g. from tests at each time point
%         name-value pairs of optional inputs for graphical controls:
%           yVal = 0. y value to plot line along
%           alpha = .05. alpha threshold, only values below are plotted
%           plotargs = cell array of name-par to pass to plot (e.g. Color
%           or LineWidth)
% 
% Returns handle to plot
% 

names = {'yVal', 'alpha','plotargs'};
defaults = {0, .05, {'Color','k','LineWidth',5}};
[yVal, alpha, plotargs] = parsepvpairs(names, defaults, varargin{:});

pVal2 = double(pVals < alpha);
pVal2(pVal2==0) = NaN;

h = plot(pVal2 .* yVal, '-', plotargs{:});
end