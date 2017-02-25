function [ thres ] = findOptimalThreshold(handles,cellinfo,mode)
%findOptimalThreshold(handles,celd.data) finds the threshold giving the
%best mutual information and csd regarding the response. The optimization
%is regularized by the number of events

% defining function vars
cell_data = cellinfo.data;
sp = cell_data.S;
stim = cellinfo.stim{1}.stim;
vNames = cellinfo.Properties.VariableNames;
if(any(ismember(vNames,'properties')))
    stim_type = cellinfo.properties.stim_type;
elseif(any(ismember(vNames,'fkey'))) 
    props = handles.prop_map(cell2mat(cellinfo(1,:).fkey)); % TODO: case where cellinfo has info on more than 1 cell
    stim_type = props.stim_type;
end
[quant_size] = findNstates(stim_type);
        
% defining the optimization criterion
switch mode
    case 'csd'
        lambda = -100;
        costfun = @(thr) -(condStimulusDivergence(pruneSpikes(thr,sp), stim, quant_size) ...
            + 1/lambda.* sum((sp >= thr))); % L1 regularization induces sparsity
        thr_hat = fminbnd(costfun,0,max(sp));
    
    case 'mi'
        lambda = -1e4;
        costfun = @(thr) -(mutualInform(pruneSpikes(thr,sp), stim, quant_size) ...
            + 1/lambda.* sum((sp >= thr))); % L1 regularization induces sparsity
        thr_hat = fminbnd(costfun,0,max(sp));
end
set(handles.param_SLDR,'Value',thr_hat);
thres = thr_hat;

% costfun = @(thr) condStimulusDivergence(pruneSpikes(thr,sp),  stim, quant_size); % without regularization
% 
% costfun = @(thr) condStimulusDivergence(pruneSpikes(thr,sp), stim, quant_size) ...
%     + 1/lambda.* sum((sp >= thr).^2); % L2 regularization on the number of spikes

% test the cost function 
% out = zeros(length((0:0.01:0.5)),1);
% range = (0:0.01:0.5);
% for k = 1:numel(range)
%     out(k) = costfun2(range(k));
% end
% figure(); plot((0:0.01:0.5),out);
% 
% % some plots for debugging
% fig = figure();set(fig ,'position',[ 2000 0 500 400]);
% C = cell_data.C;
% bsl = mean(cell_data.mcmc_samples.Cb);
% S = mapminmax_nc(pruneSpikes(thr_hat,sp),bsl, max(C));
% plot(C);hold on; stem(S,'BaseValue',bsl,'ShowBaseLine','off');hold off;
% 
% % mutualInform
end

function [out] = pruneSpikes(thr, in)
% prunes spikes according to threshold

    out = in;
    out(in < thr) = 0;
end

function [q] = findNstates(stim_type)
    switch stim_type
        case 'ORI'
            q = 8;

        case 'CRP'
            q = 20;
            error('CRP not implemented yet...');
            % TODO: modify stimulus vector to have 20 states

        case 'BDN'
            % TODO: find other method for autothres
            error('BDN not implemented yet...');

        case 'CNT'
            error('CNT not implemented yet...');

        otherwise
            error('stimulus type is not recognized');    
    end
end
