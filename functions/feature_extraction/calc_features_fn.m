function [ qi_C, mi, csd ] = calc_features_fn(df, selection, quant_size)
%calculates Quality index as it appears in ** Baden et al 2016

% load 1D data
sc = table2array(df.cell_table(:,selection)); % single cell
dt = 1/df.properties.imaging_rate;
time = (0:dt:( df.properties.signalLenght-1 ) * dt);
baseline = mean( sc.mcmc_samples.Cb );

% reshape signal according to reps
C = sc.C;
S = sc.S;
stimulus = df.stimulus;
numReps = size(df.stimPartition,2);
CS_rsp = arrayfun(@(x,y) [C(x:y); S(x:y)],df.stimPartition(1,:),...
    df.stimPartition(2,:),'UniformOutput',false);

% load quantization size
q_size = 100;
if(nargin > 2)
    q_size  = quant_size;
end
%% find QI - qualitive index [Baden et al. 2016]

min_length = min(diff(df.stimPartition))+1;
E_rc = cell2mat(cellfun(@(x) (x(1,1:min(min_length)))',CS_rsp,...
    'UniformOutput',false));
V_tc = cellfun(@(x) var(x(1,:)),CS_rsp);
qi_C = var(mean(E_rc,2))/mean(V_tc);

% find QI
% E_rs = cell2mat(cellfun(@(x) (x(2,1:min(min_length)))',CS_rsp,...
%     'UniformOutput',false));
% V_ts = cellfun(@(x) var(x(2,:)),CS_rsp);
% qi_S = var(mean(E_rs,2))/mean(V_ts);
%% find MI

[ mi ] = mutualInform( S, stimulus, q_size );
%% find CSD

[ csd ] = condStimulusDivergence( S, stimulus, q_size );
%% update cell struct

% var_name = df.cell_table(1,selection).Properties.VariableNames; 
% temp_struct = table2array(df.cell_table(:,var_name));
% temp_struct.QI_ca = qi_C; 
% temp_struct.QI_spk = qi_S;
% temp_struct.MI = mi;
% temp_struct.xcoval = xcov_val;
% temp_struct.CSD = csd;
% df.cell_table(:,var_name) = [];
% df.cell_table(:,var_name) = table(temp_struct,'VariableNames',var_name);

end

%% for CSD debugging
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

