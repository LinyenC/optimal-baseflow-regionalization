% Yiwen Mei (yiwen.mei@uconn.edu)
% CIRCA, University of Connecticut
% First version on 5/1/2022
% Last updated on 10/18/2022

%% Functionality:
% The script inplements the calculations of time lag/lead between the tracer
%  and streamflow time series which maximizes the correlation.

%% Input
% rc : time series of streamflow and tracer as Matlab timetable (name of the
%       date time column should be "Date");
% lag: time range from which the optimal lag/lead is identified (positive value
%       means that tracer lags streamflow and vice versa);

% flg: the type of correlation for use.
%       - default is the area of "hysteresis loop" encapsulated by streamflow
%         and tracer data (Cano-Paoli et al. 2019, Convenient use of electrical
%         conductivity measurements to investigate hydrological processes in
%         Alpine headwaters);
%       - other possible options are the three correlation coefficient options
%         in the Matlab's "corr" function, i.e., "Pearson", "Spearman", and "Kendall").

%% Output
% lt_op: optimal time lag/lead in number of time step;
% rc_o : shifted time series of streamflow and environmental tracer;
% cc_op: optimal lagged correlation coefficient;
% CC_o : correlation measures for all the time range investigated.

function [lt_op,rc_o,cc_op,CC_o]=shift_rec(rc,lag,varargin)
%% Check the inputs
narginchk(2,3);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'rc',@(x) validateattributes(x,{'timetable'},{'ncols',2},mfilename,'rc'));
addRequired(ips,'lag',@(x) validateattributes(x,{'double'},{'vector'},mfilename,'lag'));

addOptional(ips,'flg','Hysteresis',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'flg'));

parse(ips,rc,lag,varargin{:});
flg=ips.Results.flg;
clear ips varargin

%% Shift the time series
CC=nan(size(lag));
CC_o=nan(size(lag));
for l=1:length(lag)
  lt=lag(l);
  if lt>=0 % C lag Q
    Q=rc{1:end-lt,1};
    C=rc{1+lt:end,2};
  else % C lead Q
    Q=rc{1-lt:end,1};
    C=rc{1:end+lt,2};
  end

  switch flg
    case {'Spearman','Pearson','Kendall'}
      k=any(isnan([Q C]),2);
      Q=Q(~k);
      C=C(~k);
      cc=corr(Q,C,'Type',flg);
      CC_o(l)=cc;
      CC(l)=-abs(CC_o(l));

    case 'Hysteresis' % Cano-Paoli et al. (2019) not tested
      Q_i=Q(1:end-1);
      Q_i1=Q(2:end);
      C_i=C(1:end-1);
      C_i1=C(2:end);
      CC_o(l)=sum((Q_i1-Q_i).*(C_i1+C_i),'omitnan')/2;
      CC(l)=abs(CC_o(l));
%       k=all(~isnan([Q C]),2);
%       CC(l)=polyarea(Q(k),C(k));
  end
end

%% Output the shifted time series
[~,id]=min(CC);
lt_op=lag(id);
cc_op=CC_o(id);
if lt_op>=0 % C lag Q
  Q=rc{1:end-lt_op,1};
  C=rc{1+lt_op:end,2};
  t=rc.Date(1:end-lt_op);
else % C lead Q
  Q=rc{1-lt_op:end,1};
  C=rc{1:end+lt_op,2};
  t=rc.Date(1-lt_op:end);
end
rc_o=array2timetable([Q C],'RowTime',t,'VariableNames',rc.Properties.VariableNames);
rc_o.Properties.DimensionNames{1}=rc.Properties.DimensionNames{1};
end
