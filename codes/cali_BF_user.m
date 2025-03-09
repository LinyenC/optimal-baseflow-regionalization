function cali_BF_user()
%% Model options
Name_C='CFit'; % CFit (0), FSpl (1);
Name_BF='FSMM'; % RDF (0), SMM (1), FSMM (2), BRM (3)
Trc='EC'; % EC (0), Tur (1)
em_id=2; % Objective function: MSE (1), r^2 (2), KGE (3), NSE (4)
L_thr=365; % Minimum number of tracer records
cc_thr=.3; % Correlation coefficient for filtering the sites

pth='C:\IsoHydro\RS\Baseflow Separation Optimization';
dpth=fullfile(pth,'USGS');
rpth=fullfile(dpth,sprintf('Paired_%s',Trc));
opth=fullfile(pth,sprintf('SO_%s-%s-%s',Name_C,Name_BF,Trc));
if exist(opth,'dir')~=7
  mkdir(opth);
end
ppth=fullfile(opth,'ckps');
if exist(ppth,'dir')~=7
  mkdir(ppth);
end

pflg=false;
% num_proc=str2double(getenv('SLURM_NTASKS')); % For parallel computing
% if isempty(gcp('nocreate'))
%   parpool('local',num_proc);
% end

%% Parameters
BF_thr=.5:.1:.9;
itm=cellfun(@(X) sprintf('f_BF%i',X*100),num2cell(BF_thr),'UniformOutput',false);
itm=[{'BFI_lt'} itm {'sz','mn_m','mn_o','sd_m','sd_o','bias_mn','bias_sd','rmse','crmse','CC',...
    'NSE','KGE','r2'}];

switch Name_BF
  case 'RDF'
    itm=[{'BFI_x'} itm];
    pint=[];
  case 'SMM'
    itm=[{'N_d','rt'} itm];
    pint=1;
  case 'FSMM'
    itm=[{'N_d','rt','BFI_x'} itm];
    pint=1;
  case 'BRM'
    itm=[{'f_bmp','k_rise'} itm];
    pint=[];
end
switch Name_C
  case 'CFit'
    itm=[{'c_str','c_f','n'} itm]; % c*, cf, n in Mei et al. (2024)
end
switch Trc
  case 'EC'
    itm=[{'Lcc_EC'} itm];
  case 'Tur'
    itm=[{'Lcc_Tur'} itm];
end
lag=-7:7;

%% Station records
fn=fullfile(rpth,'stn_TS.mat');
stn=matfile(fn);
stn=stn.stn;
k=stn.resC>=1 | isnan(stn.resC);
stn(k,:)=[];

Prms=nan(size(stn,1),length(itm));
% parfor i=1:size(stn,1) % For parallel computing
for i=1:size(stn,1)    
  sn=stn.Properties.RowNames{i};
  BF_md=array2table({Name_BF},'VariableNames',{'Name'});
  C_md=array2table({Name_C Trc},'VariableNames',{'Name','Trc'});

%% Read the records
  fn=fullfile(rpth,sprintf('RC-%s.mat',sn));
  rc=matfile(fn);
  rc=rc.rc;

  if sum(all(~isnan(rc{:,:}),2))>L_thr && ~isnan(stn.Area_rp(i))
%% Shift the time series
    [~,rc_s,~,~]=shift_rec(rc,lag,'Spearman'); % Spearman
    rc_s.Properties.VariableNames={'Qtot','C'};
    rc_s(isnan(rc_s.Qtot),:)=[]; % Cut the begin and end NaNs
    rc_s=retime(rc_s,'daily','fillwithmissing');

%% Flag to continue or not
    Lcc_op=corr(rc_s.Qtot,rc_s.C,'Type','Spearman','rows','complete');
    cflg=false;
    switch Trc
      case 'EC'
        if Lcc_op<-cc_thr
          cflg=true;
        end
      case 'Tur'
        if Lcc_op>cc_thr
          cflg=true;
        end
    end

    if cflg
%% Signature concentration estimation parameters
      switch C_md.Name{1}
        case 'CFit'
          C_md=addvars(C_md,NaN,'NewVariableNames',{'c_str'});
          C_md=addvars(C_md,NaN,'NewVariableNames',{'c_f'});
          C_md=addvars(C_md,NaN,'NewVariableNames',{'n'});
        case 'FSpl'
      end

%% Baseflow filter parameters
      switch BF_md.Name{1}
        case 'RDF'
          BF_md=addvars(BF_md,NaN,'NewVariableNames',{'BFI_x'});
          BF_md=addvars(BF_md,stn.resC(i),'NewVariableNames',{'resC'});
          lb=.01; % BFI_x
          ub=.99;
        case 'SMM'
          BF_md=addvars(BF_md,1,'NewVariableNames',{'sc'});
          BF_md=addvars(BF_md,NaN,'NewVariableNames',{'N_d'});
          BF_md=addvars(BF_md,NaN,'NewVariableNames',{'rt'});
          lb=[2 .5]; % N_d, rt
          ub=[max([round(3*stn.LSP(i)) 10]) 1];
        case 'FSMM'
          BF_md=addvars(BF_md,1,'NewVariableNames',{'sc'});
          BF_md=addvars(BF_md,NaN,'NewVariableNames',{'N_d'});
          BF_md=addvars(BF_md,NaN,'NewVariableNames',{'rt'});
          BF_md=addvars(BF_md,NaN,'NewVariableNames',{'BFI_x'});
          BF_md=addvars(BF_md,stn.resC(i),'NewVariableNames',{'resC'});
          lb=[2 .5 .01]; % N_d, rt, BFI_x
          ub=[max([round(3*stn.LSP(i)) 10]) 1 .99];
        case 'BRM'
          BF_md=addvars(BF_md,NaN,'NewVariableNames',{'f_bmp'});
          BF_md=addvars(BF_md,NaN,'NewVariableNames',{'k_rise'});
          dq=diff(rc.Qtot);
          dq(dq<=0)=[];
          dq(dq<=0 | isnan(dq))=[]; % k_rise in m3/s/d
          lb=[.01 quantile(dq,.01)]; % f_bmp k_rise
          ub=[.5 quantile(dq,.6)];
      end

      try
%% Execute calibrations
        rc.Properties.VariableNames={'Qtot','C'};
        rc(isnan(rc.Qtot),:)=[]; % Cut the begin and end NaNs
        rc=retime(rc,'daily','fillwithmissing');

        cpfn=fullfile(ppth,sprintf('cp-%s.mat',sn));
        prm_op=cali_OHS(rc,C_md,BF_md,lb,ub,pint,em_id,cpfn,pflg);
        [BF_md_o,~,val_B]=update_field(BF_md,prm_op);

%% Calculate Qb with the optimal parameters
        [rc_o,BFI_lt,~,C_md_o]=cal_C(rc,C_md,BF_md_o);
        [~,~,val_C]=update_field(C_md,C_md_o);
        Y=[rc_o.C_tot rc_o.C];
%         Y(any(isnan(Y),2) | rc_o.Vid==0 | rc_o.C_tot<0,:)=[];
        Y(any(isnan(Y),2) | Y(:,2)<=0,:)=[];
        [sts,ems]=errM(Y,.05,[]);

%% Calculate the baseflow day
        BFD=nan(1,length(BF_thr));
        for h=1:length(BF_thr)
          k=rc_o.Q_bf./rc_o.Qtot>BF_thr(h);
          BFD(h)=sum(k)/size(rc,1);
        end

        Prms(i,:)=[Lcc_op val_C val_B BFI_lt BFD sts(1) sts(2) sts(3) sts(4) sts(5) sts(2)/sts(3)...
            sts(4)/sts(5) ems(1:6)];

%% Save the outputs
        ofn=fullfile(opth,sprintf('RC_o-%s.mat',sn));
        parsave(ofn,rc_o,'rc');
        parsave(ofn,C_md_o,'C_md','-append');
        parsave(ofn,BF_md_o,'BF_md','-append');

      catch ME
        efn=fullfile(opth,sprintf('ER-%s.mat',sn));
        parsave(efn,ME,'ER');
        fprintf('%i - %s\n',i,sn);
      end
    end
  end
end

%% Save the optimal hyper-parameters
k=all(isnan(Prms),2);
Prms(k,:)=[];
rn=stn.Properties.RowNames;
rn(k)=[];
Prms=array2table(Prms,'VariableNames',itm);
Prms.Properties.RowNames=rn;
ofn=fullfile(opth,'Prms.Sopt.mat');
save(ofn,'Prms');

%% Zip the files
% zfn=fullfile(opth,'RC_o.zip');
% mfl=fullfile(opth,'*.mat');
% zip(zfn,mfl);
end

function prm_op=cali_OHS(rc,C_md,BF_md,lb,ub,pint,em_id,cpfn,pflg)
Opts=optimoptions('surrogateopt','Display','off','MaxFunctionEvaluations',150,...
    'CheckpointFile',cpfn,'UseParallel',pflg);
prm_op=surrogateopt(@(prm)OHS_obj(prm,rc,C_md,BF_md,em_id),lb,ub,pint,...
    [],[],[],[],Opts);
end

function EM=OHS_obj(prm,rc,C_md,BF_md,em_id)
%% Estimate C
BF_md=update_field(BF_md,prm);
rc_o=cal_C(rc,C_md,BF_md);

%% Error metric
Y=[rc_o.C_tot rc_o.C]; % C_tot is the modeled one
Y(any(isnan(Y),2) | rc_o.Vid==0 | rc_o.C_tot<0,:)=[];
[~,ems]=errM(Y,.05,[]);
if em_id==1
  EM=ems(em_id)^2; % MSE
elseif em_id==2 % r^2
  EM=1-ems(3);
elseif em_id==3 % KGE
  EM=1-ems(5);
elseif em_id==4 % NSE
  EM=1-ems(em_id);
end
end
