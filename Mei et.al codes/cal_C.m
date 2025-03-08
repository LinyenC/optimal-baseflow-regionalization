% Yiwen Mei (yiwen.mei@uconn.edu)
% CIRCA, University of Connecticut
% First version on 5/1/2022
% Last updated on 8/8/2023

%% Functionality:
% The code estimates the concentration of tracer (tested for specific electrical
%  conductivity and turbidity only) using the hybrid baseflow separation method
%  (Mei et al. 2024, Water Resour. Res. Optimal Baseflow Separation through Chemical
%   Mass Balance: Comparing the Usages of Tracers, Signature Concentration Estimation
%   Strategies, and Baseflow Filters).

%% Input
%  rc  : time series of streamflow and environmental tracer as Matlab timetable
%         (name of the date time column should be "Date"; units for the other
%         two column don't matter);
% C_md : table stores the model parameters for the Chemical Mass Balance filter
%         (possible C_md options are "CFit", "FSpl", and "Cnst");
%         - C_md of CFit contains five fields, namely "Name", "Trc", "c_str",
%           "c_f", and "n"; C_md of FSpl and Cnst contains two fields, namely
%           "Name" and "typ"
% BF_md: table stores the model parameters for the empirical baseflow filters
%         (possible BF_md options are "RDF", "SMM", "FSMM", and "BRM").
%         - BF_md of RDF contains three fields, namely "Name", "BFI_x", and "resC";
%           BF_md of SMM contains three fields, namely "Name", "N_d", and "rt";
%           BF_md of FSMM contains five fields, namely "Name", "N_d", "rt",
%           "BFI_x", and "resC"
%           BF_md of BRM contains three fields, namely "Name", "f_bmp", and "k_rise"

%% Output
%  rc_o : output time series in timetable format (from left to right are time
%          series of observed streamflow - Qtot, observed tracer concentration - C,
%          estimated baseflow - Q_bf, estimated tracer concentration in baseflow - Q_bf,
%          estimated tracer concentration in quick flow - Q_sf, flag for validation
%          time step - Vid, and estimated tracer concentration - C_tot);
%  BFI  : long-term baseflow index;
%  N_bp : number of time step with negative baseflow and baseflow larger than
%          streamflow;
% C_md_o: table stores the model parameters for the mass-balance filters with
%          the estimated parameters filled.

%% Additional note
% Require RDF.m, SMM.m, BRM.m, and Run_Length.m.

function [rc_o,BFI,N_bp,C_md_o]=cal_C(rc,C_md,BF_md)
%% Check the inputs
narginchk(3,3);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'rc',@(x) validateattributes(x,{'timetable'},{'ncols',2},mfilename,'rc'));
addRequired(ips,'C_md',@(x) validateattributes(x,{'table'},{'nonempty'},mfilename,'C_md'));
addRequired(ips,'BF_md',@(x) validateattributes(x,{'table'},{'nonempty'},...
    mfilename,'BF_md'));

parse(ips,rc,C_md,BF_md);
clear ips

%% Calculate baseflow
switch BF_md.Name{1}
  case 'RDF'
    BFI_x=BF_md.BFI_x;
    resC=BF_md.resC;
    [bf,BFI,N_bp]=RDF(rc.Qtot,BFI_x,resC);

  case 'SMM'
    sc=BF_md.sc;
    N_d=BF_md.N_d;
    rt=BF_md.rt;
    [bf,BFI,N_bp]=SMM(rc.Qtot,sc,N_d,rt);

  case 'FSMM'
    sc=BF_md.sc;
    N_d=BF_md.N_d;
    rt=BF_md.rt;
    BFI_x=BF_md.BFI_x;
    resC=BF_md.resC;
    bf0=SMM(rc.Qtot,sc,N_d,rt);
    [bf,BFI,N_bp]=RDF(rc.Qtot,BFI_x,resC,bf0);

  case 'BRM'
    f_bmp=BF_md.f_bmp;
    k_rise=BF_md.k_rise;
    [bf,BFI,N_bp]=BRM(rc.Qtot,f_bmp,k_rise);
end

rc_o=addvars(rc,bf,'NewVariableNames',{'Q_bf'});

%% Estimate the concentration
C_md_o=C_md;
switch C_md.Name{1}
  case 'FSpl'
    switch C_md.Trc{1}
      case 'EC'
        [C_bf,id]=LE_interp(rc,'max');
        C_bf(C_bf<rc.C)=rc.C(C_bf<rc.C);
        C_sf=LE_interp(rc,'min');
        C_sf(C_sf>rc.C)=rc.C(C_sf>rc.C);
      case 'Tur'
        C_sf=LE_interp(rc,'max');
        C_sf(C_sf<rc.C)=rc.C(C_sf<rc.C);
        [C_bf,id]=LE_interp(rc,'min');
        C_bf(C_bf>rc.C)=rc.C(C_bf>rc.C);
    end

    id=unique([id-1;id;id+1]); % Exclude the baseflow period
    id(id<1 | id>size(rc,1))=[];

  case 'CFit'
    sf=rc.Qtot-bf;
    q=rc.Qtot;
    q(q==0 | isnan(q))=[];
    c=rc.C;
    c(c==0 | isnan(c))=[];
    mx=max(log(c),[],'omitnan');
    mn=mean(log(c),'omitnan');
    sd=std(log(c),[],'omitnan');

    k=isnan(sf) | isnan(rc.C) | sf./rc.Qtot<.001 | rc.Qtot<=0; % Exclude the baseflow
    id=find(k);                                                % period
    switch C_md.Trc{1}
      case 'EC'
        c_str_ub=max(exp(mx+sd)./quantile(q,.99).^(-5:0)); % Unit c_b=C/Q^n, % upper bound of c_str
        c_str_0=exp(mn)/mean(q)^-1; % initial value of c_str
        n_ub=0;
        n_lb=-5;
        n_st=-.5;

        try
% ATTENTION! the order of the parameters follows the alphabetical order
          mdl=fittype('c_f*c_str*y^n*x+c_str*y^(n+1)','indep',{'x','y'},'dep','z');
          [Fu,~]=fit([sf bf],rc.C.*rc.Qtot,mdl,'Lower',[.05 0 n_lb],...
              'Upper',[.95 c_str_ub n_ub],'StartPoint',[.5 c_str_0 n_st],'Exclude',k);
        catch
          mdl=fittype('c_f*c_str*y^n*x+c_str*y^(n+1)','indep',{'x','y'},'dep','z',...
              'problem','n');
          [Fu,~]=fit([sf bf],rc.C.*rc.Qtot,mdl,'Lower',[.05 0],'Upper',[.95 c_str_ub],...
              'StartPoint',[.5 c_str_0],'Exclude',k,'problem',0);
        end
        C_sf=Fu.c_f*Fu.c_str*bf.^Fu.n;
        C_sf(rc.Qtot==0)=NaN;
        C_bf=Fu.c_str*bf.^Fu.n;
        C_bf(rc.Qtot==0)=NaN;

      case 'Tur'
        c_str_ub=max(exp(mx+sd)./quantile(q,.01).^(0:5));
        c_str_0=exp(mn)/mean(q);
        n_ub=5;
        n_lb=0;
        n_st=.5;

        try
          mdl=fittype('c_str*x^(n+1)+c_f*c_str*x^n*y','indep',{'x','y'},'dep','z');
          [Fu,~]=fit([sf bf],rc.C.*rc.Qtot,mdl,'Lower',[.05 0 n_lb],...
              'Upper',[.95 c_str_ub n_ub],'StartPoint',[.5 c_str_0 n_st],'Exclude',k);
        catch
          mdl=fittype('c_str*x^(n+1)+c_f*c_str*x^n*y','indep',{'x','y'},'dep','z',...
              'problem','n');
          [Fu,~]=fit([sf bf],rc.C.*rc.Qtot,mdl,'Lower',[.05 0],'Upper',[.95 c_str_ub],...
              'StartPoint',[.5 c_str_0],'Exclude',k,'problem',0);
        end
        C_sf=Fu.c_str*sf.^Fu.n;
        C_sf(rc.Qtot==0)=NaN;
        C_bf=Fu.c_f*Fu.c_str*sf.^Fu.n;
        C_bf(rc.Qtot==0)=NaN;
    end
    C_md_o.c_f=Fu.c_f;
    C_md_o.c_str=Fu.c_str;
    C_md_o.n=Fu.n;

  case 'Cnst'
    switch C_md.typ{1}
      case 'qt'
        C_sf=quantile(rc.C,C_md.qt_sf/100)*ones(size(rc,1),1);
        C_bf=quantile(rc.C,C_md.qt_bf/100)*ones(size(rc,1),1);
      case 'val'
        C_sf=C_md.c_sf*ones(size(rc,1),1);  
        C_bf=C_md.c_bf*ones(size(rc,1),1);
    end

    id=find(rc.C>=C_bf | rc.C<=C_sf);
end

rc_o=addvars(rc_o,C_bf,'NewVariableNames',{'C_bf'});
rc_o=addvars(rc_o,C_sf,'NewVariableNames',{'C_sf'});

%% Validation time step
k=true(size(rc,1),1);
k(id)=false;
k(isnan(rc_o.C))=false;
rc_o=addvars(rc_o,k,'NewVariableNames',{'Vid'});

%% Estimated C
C_tot=sum([C_bf.*bf C_sf.*(rc_o.Qtot-bf)],2)./rc_o.Qtot;
rc_o=addvars(rc_o,C_tot,'NewVariableNames',{'C_tot'});
end

function [C_cmp,ID]=LE_interp(Rc,mth)
%% Find the continuous periods
TSi=fillmissing(Rc.C,'constant',-999,'MaxGap',32);
TSi=~isnan(TSi);
k=Run_Length(TSi,true,Rc.C);
k=reshape(k',2,length(k)/2)';
k(k(:,1)==1,1)=1:sum(k(:,1));
TSi=Run_Length(reshape(k(:,1:2)',size(k,1)*2,1),false,[])';

%% Interpolation for each period
C_cmp=nan(size(Rc,1),1);
ID=[];
for i=1:max(TSi)
  rc=Rc(TSi==i,:);
  [Y,M,~]=datevec(rc.Date);
  C_in=accumarray([Y-min(Y)+1 M],rc.C,[],@(x) ists(x,mth,true),NaN)';
  C_in=reshape(C_in,numel(C_in),1);
  C_in(isnan(C_in))=[];
  D_bs=accumarray([Y-min(Y)+1 M],datenum(rc.Date),[],@(x) ists(x,'min',true),NaN);
  D_dur=accumarray([Y-min(Y)+1 M],rc.C,[],@(x) ists(x,mth,false),NaN);
  D_in=D_bs'+D_dur'-1;
  D_in=reshape(D_in,numel(D_in),1);
  D_in(isnan(D_in))=[];
  T_in=datetime(D_in,'ConvertFrom','datenum');

% Perform the interpolation
  if length(T_in)>1
    vq=interp1(T_in,C_in,rc.Date,'pchip',NaN);
    vq1=interp1(T_in,C_in,rc.Date,'linear',NaN);
    k=vq<0 | isnan(vq);
    vq(k)=vq1(k);
  else
    vq=repelem(C_in,size(rc,1));
  end
  C_cmp(TSi==i)=vq;

% Convert datetime to index
  k=ismember(Rc.Date,T_in);
  ID=[ID;find(k)];
end
end

function I=ists(x,mth,vflg)
if all(isnan(x))
  I=NaN;
else
  switch mth
    case 'max'
      if vflg
        I=max(x,[],'omitnan');
      else
        [~,I]=max(x,[],'omitnan');
      end
    case 'min'
      if vflg
        I=min(x,[],'omitnan');
      else
        [~,I]=min(x,[],'omitnan');
      end
  end
end
end
