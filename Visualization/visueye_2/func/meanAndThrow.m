function [ signal ]= meanAndThrow( signal, index )
% takes the samples with index of the signal. means them and then replaces the mean
% with the values of the index

sig = signal(index);
sig = reshape(sig,[],2);
sig_val = mean(sig,2);
signal(index(length(sig_val)+1:end)) = sig_val;
signal(index(1:length(sig_val))) = [];
end
