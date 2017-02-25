function [ out ] = Composite2array( in )
%Composite2array( in ) takes a composite structure and returns the
%equivilant array

cellREP = {in{:}}; % cell Array representation
matREP = cell2mat(cellREP);
out = sum(matREP,2);
end

