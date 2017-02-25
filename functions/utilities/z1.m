function x = z1(y)
% linear normalize between 0 and 1
x = (y-min(y(:)))/(max(y(:))-min(y(:)))+eps;
