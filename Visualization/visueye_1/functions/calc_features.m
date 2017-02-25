function [ handles, hObject ] = calc_features(handles, hObject)
%calculates Quality index as it appears in ** Baden et al 2016

% load 1D data
contents = cellstr(get(hObject,'String'));
selection = contents{get(hObject,'Value')};
sc = table2array(handles.df.cell_table(:,selection)); % single cell
dt = 1/handles.df.properties.imaging_rate;
time = (0:dt:( handles.df.properties.signalLenght-1 ) * dt);
baseline = mean( sc.mcmc_samples.Cb );

% reshape signal according to reps
C = sc.C;
S = sc.S;
numReps = size(handles.stimPartition,2);
[C_mat, S_mat] = deal(cell(numReps,1));
CS_rsp = arrayfun(@(x,y) [C(x:y); S(x:y)],handles.stimPartition(1,:),...
    handles.stimPartition(2,:),'UniformOutput',false);
%% find QI - qualitive index [Baden et al. 2016]

min_length = min(diff(handles.stimPartition))+1;
E_rc = cell2mat(cellfun(@(x) (x(1,1:min(min_length)))',CS_rsp,...
    'UniformOutput',false));
V_tc = cellfun(@(x) var(x(1,:)),CS_rsp);
qi_C = var(mean(E_rc,2))/mean(V_tc);

% find QI
E_rs = cell2mat(cellfun(@(x) (x(2,1:min(min_length)))',CS_rsp,...
    'UniformOutput',false));
V_ts = cellfun(@(x) var(x(2,:)),CS_rsp);
qi_S = var(mean(E_rs,2))/mean(V_ts);

str_2_disp = ['ca - ',num2str(qi_C),'  spk - ',num2str(qi_S)];
set(handles.qi_TXT,'String',str_2_disp);

%% find MI

stimulus = handles.stimulus;
[ mi ] = mutualInform( S, stimulus );
str_2_disp = num2str(mi);
set(handles.mi_TBX,'String',str_2_disp);
%% find XCORR
% finds xcov funciton between he actual Calcium trace and the
% stimulus at zero

maxlags = 100;
xcov_val = 1/(length(stimulus)-1) .* xcov(S,stimulus,maxlags);
xcov_val = xcov_val(maxlags + 1);
str_2_disp = num2str(xcov_val);
set(handles.xcorr_TBX,'String',str_2_disp);
%% find CSD
% % for debugging
% % create stimulus depended spiking
% dt = 1/handles.df.properties.imaging_rate;
% temp_indx = ((stimulus > 134 & stimulus < 136) | ( stimulus > 314 & stimulus < 316));
% stim_trigger = z1(stimulus .* temp_indx); 
% stim_trigger_2 = (stimulus >0 );
% altered_stim = stim_trigger;
% altered_stim(stim_trigger == 0) = 0.05;
% altered_stim =  z1(altered_stim .* abs(0.0001 * randn(length(stimulus),1)));
% test = simpp(10 .*altered_stim,dt);
% t = (0:dt:(handles.df.properties.signalLenght)*dt)';
% test_depend = histcounts(test, t);
% test_independ = histcounts(simpp(0.1 .* ones(length(stimulus),1),dt), t);
% [ csd0 ] = condStimulusDivergence( stim_trigger_2,stimulus );
% [ csd1 ] = condStimulusDivergence( stim_trigger,stimulus );
% [ csd15 ] = condStimulusDivergence( stimulus,stimulus );
% [ csd2 ] = condStimulusDivergence( test_depend,stimulus );
% [ csd3 ] = condStimulusDivergence( test_independ,stimulus );
% [ csd4 ] = condStimulusDivergence( ones(length(stimulus),1),stimulus );
% [ csd4 ] = condStimulusDivergence( zeros(length(stimulus),1),stimulus );

[ csd ] = condStimulusDivergence( S, stimulus );
str_2_disp = num2str(csd);
set(handles.csd_TBX,'String',str_2_disp);
%% update cell struct

var_name = handles.df.cell_table(1,selection).Properties.VariableNames; 
temp_struct = table2array(handles.df.cell_table(:,var_name));
temp_struct.QI_ca = qi_C; 
temp_struct.QI_spk = qi_S;
temp_struct.MI = mi;
temp_struct.xcoval = xcov_val;
temp_struct.CSD = csd;
handles.df.cell_table(:,var_name) = [];
handles.df.cell_table(:,var_name) = table(temp_struct,'VariableNames',var_name);
%% plot CSD vs. QI_C

axes(handles.axes5)
drawnow
c = reshape(linspace(1,100,100),10,[]);
c_itr = ceil([1 + qi_C*10 ,1 + csd ]);
c_itr(c_itr > 10) = 10;
c_itr(c_itr < 1) = 1;
hold all;
scatter(qi_C,csd,25,c(c_itr(1),c_itr(2)),'filled')
% axis tight
hold off;
set(gca,'FontSize',10,'Box','on','Color',[0.8863    0.8863    0.8863]);
xlabel('QI Ca^+');
ylabel('CSD');
end

