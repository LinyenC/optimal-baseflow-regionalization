function ResTab()
pth='<data archive>\Data';
dpth=fullfile(pth,'USGS');

C_md={'CFit','FSpl'};
BF_md={'RDF','SMM','FSMM','BRM'};
Trc={'EC','Tur'};

for t=1:length(Trc)
  for c=1:length(C_md)
    for b=1:length(BF_md)
%% Read the tables
      opth=fullfile(pth,sprintf('SO_%s-%s-%s',C_md{c},BF_md{b},Trc{t}));
      fn=fullfile(opth,'Prms.Sopt.mat');
      tab=matfile(fn);
      tab=tab.Prms;

      vn=tab.Properties.VariableNames;
      vn=cellfun(@(X) sprintf('%s_%i%i%i',X,c-1,b-1,t-1),vn(2:end),'UniformOutput',false); % Exclude Lcc_
      tab.Properties.VariableNames(2:end)=vn; % Exclude Lcc_

%       fprintf('%d\n',size(tab,1));

%% Join the tables
      if c~=1 || b~=1
        tab(:,1)=[]; % Exclude Lcc_
        Tab=outerjoin(Tab,tab,'Keys','Row','Mergekey',true);
      else
        if t==1
          Tab=tab;
        else
          Tab=outerjoin(Tab,tab,'Keys','Row','Mergekey',true);
        end
      end
    end
  end
end

%% Join with the station records
fn=fullfile(dpth,'stn_TS');
stn=matfile(fn);
stn=stn.stn;
Tab=innerjoin(stn,Tab,'Keys','Row');

%% Remove gages not in CONUS
k=Tab.Latitude>50 | Tab.Latitude<25 | Tab.Longitude<-125 | Tab.Longitude>-65;
Tab(k,:)=[];

clear stn tab

%% Remove low performance catchments
k=contains(Tab.Properties.VariableNames,'KGE_');
k=Tab(:,k);
k=all(k{:,1:8}<.3,2) | all(k{:,9:end}<.3,2); % one EC and one Tur exps' KGE>0.3
Tab(k,:)=[];

%% Remove catchments that fitting can't be performed
k=all(~isnan(Tab{:,9:188}),2) | all(~isnan(Tab{:,190:369}),2); % Exclude Lcc_
Tab(~k,:)=[]; % There is one gage can't be fitted

%% Output the merged table
ofn=fullfile(pth,'ResTab.mat');
save(ofn,'Tab');
end
