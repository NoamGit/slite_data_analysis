function [ y, flag ] = correctPhase(x1, x2, max_phase)
    % find optimal phase with simple xcorr to correct the artifact (source)
    % according to the signal (target)
    % determine - kernel size, binsize
    
    target=norm_nc(x2,5);
    source=norm_nc(x1,5);
    
    % compute cross correlation and find lags displacment
    [Cxy,lags]=xcorr(target,source,max_phase);
    [~,n] = max(Cxy);
    delta_i = floor(max_phase/2) + lags(n) + 1;
        if(delta_i > max_phase+1 || delta_i < 1)
            y = x1;
            flag = false;
%             a_corrected = arti.*0;
        else
            z = zeros(max_phase,1);
            z(delta_i) = 1;
            z = [0; z; 0];
            y = conv(x1,z,'same');
            flag = true;
            figure(10);
            subplot(311);plot([target source]);
            subplot(312);plot([target conv(source,z,'same')]);
            subplot(313);plot(lags,Cxy);
        end
    end