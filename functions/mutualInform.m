function [ mi ] = mutualInform( x, y)
% calculates mutual information for 2 continued value vector
% Based on Mo Chen's mutInf code

assert(numel(x) == numel(y));
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

if(num_val_x>255)
    [partition, ~] = lloyds(x,255);
    [~,x] = quantiz(x,partition,(1:255));
else
    [partition, ~] = lloyds(x,num_val_x);
    [~,x] = quantiz(x,partition,(1:num_val_x));
end

if(num_val_y > 255 && ~any(rem(y,1)))
    [partition, ~] = lloyds(y,255);
    [~,y] = quantiz(x,partition,(1:255));
elseif(~any(rem(y,1)))
    [partition, ~] = lloyds(y,num_val_y);
    [~,y] = quantiz(y,partition,(1:num_val_y));
end

idx = 1:n;
Mx = sparse(idx,x,1,n,max(x),n);
My = sparse(idx,y,1,n,max(y),n);
Pxy = nonzeros(Mx'*My/n); %joint distribution of x and y
Hxy = -dot(Pxy,log2(Pxy));

Px = nonzeros(mean(Mx,1));
Py = nonzeros(mean(My,1));

% entropy of Py and Px
Hx = -dot(Px,log2(Px));
Hy = -dot(Py,log2(Py));
% mutual information
z = Hx+Hy-Hxy;
mi = max(0,z);

end

