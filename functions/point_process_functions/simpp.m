function times=simpp(lambda, dt)

% times=simpp(lambda, dt) - generate inhomogeneous Poisson point process 
% with underlying rate process lambda 
% times - column vector with events' times
% lambda - rate process (stochastic intensity) in [spikes/second]
% dt - time resolution of lambda

% patch for a numerical problem of adding very small values to very large
lambda=lambda+2*eps(sum(lambda));

t=(0:length(lambda))'*dt;
% calculating the integral of lambda
z=[0; cumsum(lambda)*dt];       % suppose to represent the CDF / CMF of time [sec]
N=ceil(z(end));

% drawing enough exponentially distributed intervals

if N > 1e8 % conditioning size of N
    error('sum(lambda)*dt is out of memory limits');
end

zetta = round(exprnd(1, N, 1)*(1/dt))*dt; % drawing rand exp distributed points on N and switching time->sample->time
while sum(zetta)<z(end)
    zetta = [zetta; exprnd(1,N,1)];
end
% zimes are times on the homogeneous time-axis, which is, actually,
% time-rescaled relative to our real t-axis
% (with exponentially distributed intervals with mean  = 1)
zimes=cumsum(zetta);
zimes=zimes(zimes<=z(end));

% now we are projecting the homogeneous Poisson process zimes from the
% time-rescaled time axis onto our usual t-axis 
% (going through the integral of lambda)
try
    times=interp1(z, t, zimes);
catch
    disp(['z(find(diff(z)==0)) = ', num2str(z(find(diff(z)==0)))])
    disp('simpp: problem with interpolation')
end
