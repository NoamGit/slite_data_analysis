function [ csd ] = condStimulusDivergence( x, y, quant_size )
% condStimulusDivergence( x, y, quant_size )
% x - response 
% y - stimulus.
% quant_size - number of states in the stimulus

% mutual information test
% x = handles.good_list(1).S; % PSTH
% r_2 = handles.good_list(2).S; 
% r_b = handles.bad_list(1).S;
% y = handles.stimulus;
% part = handles.stimPartition;

assert(numel(x) == numel(y));
if(~any(x) || all(isnan(x)))
    csd = 0;
    return;
end
q_size = 100;
if(nargin > 2)
    q_size = quant_size;
end
n = numel(x);
x = reshape(x,1,n);
y = reshape(y,1,n);

% map values to integers for sparse construction
l = min(min(x),min(y));
x = x-l+1;
y = y-l+1;
k = max(max(x),max(y));
num_val_x = numel(unique(x));
num_val_y = numel(unique(y));
try
    if(num_val_x > 1 && num_val_x>255)
        [partition, ~] = lloyds(x,255);
        [~,x] = quantiz(x,partition,(1:255));
    elseif(num_val_x > 1)
        [partition, ~] = lloyds(x,num_val_x);
        [~,x] = quantiz(x,partition,(1:num_val_x));
    end
catch exception
    disp(exception);
end

if(num_val_y > 1 && num_val_y > q_size)
    [partition, ~] = lloyds(y,q_size);
    [~,y] = quantiz(y,partition,(1:q_size));
elseif(num_val_y > 1 && ~any(rem(y,1)))
    [partition, ~] = lloyds(y,num_val_y);
    [~,y] = quantiz(y,partition,(1:num_val_y));
else
    [partition, ~] = lloyds(y,q_size);
    [~,y] = quantiz(y,partition,(1:q_size));
end

idx = 1:n;
Mx = sparse(idx,x,1,n,max(x),n);
My = sparse(idx,y,1,n,max(y),n);
Pxy = (Mx'*My/n); %joint distribution of x and y
Py = mean(My,1);
P_x_cond_y = bsxfun(@rdivide,Pxy,Py); % according to Bayes

[X,Y] = meshgrid((1:size(P_x_cond_y,2)),(1:size(P_x_cond_y,2)));
X = nonzeros(X .* ~eye(size(X))); Y = nonzeros(Y .* ~eye(size(Y)));
Dkl_xy = arrayfun(@(ii, jj) Dkl(full(P_x_cond_y(:,ii)),full(P_x_cond_y(:,jj)),[ii,jj]),X,Y); 
% this Dkl part takes most of the time because the many posssible pairs..
% consider binning the pairs into bigger junks (states), at least for the CRP
% stimulus

% TODO: simple uniform mean instead of Py normalization yields better
% seperation between stimulus specific and non-specific response. Need to
% check
% Py_repmat = repmat(Py,length(Y)/length(Py),1);
Py_repmat = ones(length(Dkl_xy),1)./sum(ones(length(Dkl_xy),1));
csd = full(dot(Py_repmat(:),Dkl_xy));
end

function out = Dkl(x1,y1, pair)
%     disp(pair);
    ind = isfinite(log2(y1+eps)-log2(x1+eps));
    if(all(~ind))
        out = 0;
    else
        z = -dot(x1+eps,log2(y1+eps)-log2(x1+eps));
        out = max(0,z);
    end
end
