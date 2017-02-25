function [ ] = plotCellArr( x, y )
% multiplot( x, y ) plots cell array in same figure. y is the data and x is the independent
% parameter

index = find(~cellfun(@isempty,y));
numElements = numel(index);
factelements = factor(numElements);
if factelements == numElements % in case of a prime number
    factelements = factor(numElements+1);
end
middle = ceil(numel(factelements)/2);
[size_x, size_y] = deal(prod(factelements(1:middle)), prod(factelements(middle+1:end)));
set(gca,'FontSize', 14);
if nargin < 2 
    for n = 1:numel(index)
        s(n) = subplot(size_x,size_y,n);
        plot(y{index(n)});
        axis tight;
        title(s(n),[num2str(n),' num figure']); xlabel('[sec]');
    end
else 
    for n = 1:numel(index)
    s(n) = subplot(size_x,size_y,n);
    plot(x,y{index(n)});
    axis tight;
    title(s(n),[num2str(n),' num figure']); xlabel('[sec]');
    end
end
end


