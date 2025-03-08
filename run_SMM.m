addpath('.\Mei et.al codes')

stationInfo = load('.\Mei et.al data\ResTab.mat');

Nd_model = readtable('.\predicted_N.xlsx', 'ReadVariableNames', true,'ReadRowNames', true);
Nd_model = Nd_model(~isnan(Nd_model.pred_N),:);
allStationIndex = Nd_model.Properties.RowNames;

folderPath = '.\Mei et.al data\SO_FSpl-SMM-EC\';

for i = 1:length(allStationIndex)
    % 构建完整的文件路径
    filePath = fullfile(folderPath, ['RC_o-',allStationIndex{i},'.mat']);
    % 使用load函数加载.mat文件
    loadedData = load(filePath);
    % 从加载的数据中提取实测径流量
    Qobs = loadedData.rc.Qtot; 
    
    %==========================================================
    % N_d用随机森林结果，rt用0.9
    N_d = Nd_model.pred_N(i);
    rt = 0.9;
    [Qb_m_09,~,~]=SMM(Qobs,1,N_d,rt);
    
    % N_d用随机森林结果，rt用之前率定的值
    N_d = Nd_model.pred_N(i);
    rowIndex = strcmp(stationInfo.Tab.Properties.RowNames, allStationIndex{i});
    rt = stationInfo.Tab{rowIndex, 'rt_110'};
    [Qb_m_o,~,~]=SMM(Qobs,1,N_d,rt);
    
    % N_d用随机森林结果，rt用之前率定的值的均值
    %N_d = Nd_model.pred_N(i);
    %rt = mean(stationInfo.Tab.rt_110,"omitnan"); % 是0.9128
    %[Qb_m_omean,~,~]=SMM(Qobs,1,N_d,rt);
    
    % N_d用随机森林结果，rt用之前率定的值的中位数
    N_d = Nd_model.pred_N(i);
    rt = median(stationInfo.Tab.rt_110,"omitnan"); % 是0.9688
    [Qb_m_omedian,~,~]=SMM(Qobs,1,N_d,rt);
    
    %----------------------------
    % N_d用1.6*A^0.2，rt用0.9
    rowIndex = strcmp(stationInfo.Tab.Properties.RowNames, allStationIndex{i});
    N_d = round(1.6*(stationInfo.Tab{rowIndex, 'Area_rp'}^0.2));
    rt = 0.9;
    [Qb_Af_09,~,~]=SMM(Qobs,1,N_d,rt);
    
    % N_d用1.6*A^0.2，rt用之前率定的值
    rowIndex = strcmp(stationInfo.Tab.Properties.RowNames, allStationIndex{i});
    N_d = round(1.6*(stationInfo.Tab{rowIndex, 'Area_rp'}^0.2));
    rowIndex = strcmp(stationInfo.Tab.Properties.RowNames, allStationIndex{i});
    rt = stationInfo.Tab{rowIndex, 'rt_110'};
    [Qb_Af_o,~,~]=SMM(Qobs,1,N_d,rt);
    
    % N_d用1.6*A^0.2，rt用之前率定的值的均值
    %rowIndex = strcmp(stationInfo.Tab.Properties.RowNames, allStationIndex{i});
    %N_d = round(1.6*(stationInfo.Tab{rowIndex, 'Area_rp'}^0.2));
    %rt = mean(stationInfo.Tab.rt_110,"omitnan"); % 是0.9128
    %[Qb_Af_omean,~,~]=SMM(Qobs,1,N_d,rt);
    
    % N_d用1.6*A^0.2，rt用之前率定的值的中位数
    rowIndex = strcmp(stationInfo.Tab.Properties.RowNames, allStationIndex{i});
    N_d = round(1.6*(stationInfo.Tab{rowIndex, 'Area_rp'}^0.2));
    rt = median(stationInfo.Tab.rt_110,"omitnan"); % 是0.9688
    [Qb_Af_omedian,~,~]=SMM(Qobs,1,N_d,rt);
    
    %----------------------------
    % N_d用之前率定的值，rt用0.9
    rowIndex = strcmp(stationInfo.Tab.Properties.RowNames, allStationIndex{i});
    N_d = stationInfo.Tab{rowIndex, 'N_d_110'};
    rt = 0.9;
    [Qb_o_09,~,~]=SMM(Qobs,1,N_d,rt);
    
    % N_d用之前率定的值，rt用之前率定的值的中位数
    rowIndex = strcmp(stationInfo.Tab.Properties.RowNames, allStationIndex{i});
    N_d = stationInfo.Tab{rowIndex, 'N_d_110'};
    rt = median(stationInfo.Tab.rt_110,"omitnan"); % 是0.9688
    [Qb_o_omedian,~,~]=SMM(Qobs,1,N_d,rt);
    
    %==========================================================
    Qb_o_o = loadedData.rc.Q_bf;
    date = loadedData.rc.Date;
    
    % 创建一个表格,保存表格到 CSV 文件
    dataTable = table(date, Qobs, Qb_o_o, Qb_o_09,Qb_o_omedian,Qb_m_09,Qb_m_o,Qb_m_omedian,Qb_Af_09,Qb_Af_o,Qb_Af_omedian,...
        'VariableNames', {'date','Qobs', 'Qb_o_o', 'Qb_o_09','Qb_o_omedian','Qb_m_09','Qb_m_o','Qb_m_omedian','Qb_Af_09','Qb_Af_o','Qb_Af_omedian'});
    writetable(dataTable, ['.\BF_result\',allStationIndex{i},'.csv']);
    
end

