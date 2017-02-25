function [ varargout ] = nlextraction( CGP, spiketimes, dt, time, span, polyorder, varargin)
%[par, domain_bin, nl_est] = nlextraction( CGP, spiketimes, dt, time) estimated the non linearity in non linear poisson model
%   this function simply takes Gaussian correlated process (CGP) and its Counting
%   process (CP) and maps the ylabel of the CGP to the matching value of
%   the mean firing rate.
%         input:
%             CGP - Correlated Gaussian proccess
%             spiketimes
%             time
%             dt - CGP's time resolution
%             span - for the smoothing
%             polyorder - for Least square polynomial fitting 
%         output:
%             par - polyfit parameters
%             domain - binned domain for function
% 
%         author: Noam Cohen,2014

% test with simple examples
if nargin == 0
    demo1();
    return;
end

% input validation
if nonzeros(spiketimes < 0)
    spiketimes = spiketimes(spiketimes >= 0);
    fprintf('some spiketimes are before 0\n');
end

if nargin > 6
    num_bin = varargin{1};
    if ~isempty(varargin{2})
        nl = varargin{2};
        mu = varargin{3}; 
        sig = varargin{4};
    end
else 
    nl = []; % nl by default is empty
    mu = 1; 
    sig = 1;
end

CP = histc(spiketimes, time);

% Modifications according to Shy's code
u=CGP/std(CGP); % normalize the linear prediction
p=prctile(u',linspace(0,100,101)); % weighted binning of the functions domain
p(end)=p(end)+eps;
[N,bin]=histc(u,p); % plug CGP values to binned domain and keep the order
N=N(1:end-1);
unwant = find(bin > length(N)); % unwanted bin index
u(unwant) = [];
CP(unwant) = [];
bin(unwant) = [];
domain_bin = accumarray(bin,u,[],@mean); % the domain (X axis) of the non linearity
image_bin = accumarray(bin,CP,[],@mean); % the image (Y axis) find how many spikes for each bin
par = polyfit(domain_bin,image_bin,polyorder); % calculate polynomial fit
% nl_est  = smooth(image_bin,span);
nl_est  = polyval(par,domain_bin);

%% estimated intensity evaluation

    if ~isempty(nl)
        subplot(3,1,1);  
        plot(domain_bin, image_bin,'.')
        h = subplot(3,1,2);
        set(h,'FontSize',14)
        plot(domain_bin, nl_est,'-')
        title('estimated nl');
        xlabel('CGP'); ylabel('firing probability');
        subplot 313;
            switch nl
                case 'square'
                    plot(domain_bin, (mu+sig*domain_bin).^2);
                case 'exp'
                    plot(domain_bin, exp(mu+sig*domain_bin));
                case 'abs'
                    plot(domain_bin, abs(mu+sig*domain_bin));
            end
        
    else
        set(gca,'FontSize',14)
        plot(domain_bin, nl_est,'-b')
        hold on; plot(domain_bin, image_bin,'-*r')
        title('estimated nl');
        xlabel('CGP'); ylabel('firing probability');
        hold off;
        varargout{1} = par;
        varargout{2} = domain_bin;
        varargout{3} = nl_est;
        varargout{4} = image_bin;
    end

end
%% test for square nl

function demo1()
    GWN = 50 + 50*randn(1e6,1);
    LinKernel = normpdf((-20:30),3,10) - normpdf((-20:30),6,10);
    CGP = conv( GWN, LinKernel, 'same');
    nl = 'abs'; % abs exp square
    mu = 1; sig = 3;
    switch nl
        case 'square'
            Int = (mu+sig*CGP).^2;
        case 'exp'
            Int = exp(mu+sig*CGP);
        case 'abs'
            Int = abs(mu+sig*CGP);
    end
    dt = 1e-3; % dt = 1/N
    t = (0:length(Int)-1)'*dt;
    spiketimes = simpp(Int, dt);
    CP = histc( spiketimes, t);
    fprintf(['E(lambda) = ', num2str(length(spiketimes)/t(end)),' Hz \n']);
    num_bin = 50;
    polyorder = 8;
    span = 0.15;
    nlextraction( CGP, spiketimes, dt, t, span, polyorder, num_bin, nl, mu, sig);
end
