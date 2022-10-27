function [X, events ]= inject_swr_meanTemplate(X,swr_template, samplerate, nEvents, noise )
% input:
%       X            - signal, time by sensors matrix
%       swr_template - swr template data structure
%       sample rate  - sample rate
%       nEvents      - number of swr events 
%       noise        - data structure control voltage of swr and X respectively
% output:
%       X            - signal with swr events
%       events       - swr onsets and swr signal with noise
% Author: Zhibing Xiao,2022/10/26

%% preparing
% data info 
nSamples  = size(X,1);
nSensors  = size(X,2);
t = swr_template.t;
freq_low = swr_template.swrfilter.lowedge;
freq_high= swr_template.swrfilter.highedge;
% swr template
swr   =  mean(swr_template.epochData,1)-mean(mean(swr_template.epochData,1),2);
s.len = length(swr);
s.std = std(swr);
s.min = min(swr);
s.max = max(swr);
SWRs = repmat(swr,nEvents,1);
 
% control voltage of X according to noise.backgroundnoise
X  = X.*noise.backgroudnoise;

% location of SWR events, leave the start and end of the signal free
peak_loc = sort(randi([s.len+samplerate*0.05, nSamples-s.len-samplerate*0.05], nEvents, 1),'ascend');

%% set swrs
for iSensors = 1  :nSensors

    % set event domain noise
    voltage_level = unifrnd(0.81,1,nEvents,1);

    % add noise to swr
    SWRs_withnoise = SWRs.*repmat(voltage_level,1,s.len) ;% + voltage_noise;

    % Add noise to swr onsets of different sensors
    i_peak_loc  = peak_loc + randi([-samplerate*0.05 samplerate*0.05],nEvents,1);

    % swr location
    ripplePeak_epoch_str = i_peak_loc - (s.len-1)/2;
    ripplePeak_epoch_fin = i_peak_loc + (s.len-1)/2;
    epochRange = [ripplePeak_epoch_str,ripplePeak_epoch_fin];

    %% injecting SWR
    signal = X(:,iSensors);
    for i=1:size(epochRange,1)
        voltage_noise = signal(epochRange(i,1):epochRange(i,2))'  ;
        SWRs_withnoise(i,:) = SWRs_withnoise(i,:)* noise.weight_of_voltage  + voltage_noise;
        signal(epochRange(i,1):epochRange(i,2)) = SWRs_withnoise(i,:);
    end 
    X(:,iSensors) = signal;
    signal_bp  = simul_filter(signal, samplerate,freq_low,freq_high,256);

    %% extract SWR from signal
    SWRs_withnoise_real      = zeros(size(epochRange,1), s.len);
    SWRs_withnoise_real_bp   = SWRs_withnoise_real ;
    for i=1:size(epochRange,1)
        SWRs_withnoise_real(i,:)    = signal(epochRange(i,1):epochRange(i,2));
        SWRs_withnoise_real_bp(i,:) = signal_bp(epochRange(i,1):epochRange(i,2));
    end
    events.SWRs_simulated(iSensors,:,:) = SWRs_withnoise_real;
    events.SWRs_simulated_bp(iSensors,:,:) = SWRs_withnoise_real_bp;
    events.epochRange(iSensors,:,:) = epochRange;
        
end

% randomly chose one sensor's SWRs for visualization
v_simulated_swr = SWRs_withnoise_real;

%% visulization
sub_r = 3;
sub_c = 3;
figure('Position',[20,20,1000,1000]);

% plot swr template
subplot(sub_r,sub_c,1)
plot(t,swr,'Color','blue')
title('SWR template -- raw')

% plot TF of swr template
subplot(sub_r,sub_c,2)
Nx = s.len;
nsc = floor(Nx/4.5);
nov = floor(nsc/2);
nff = max(256,2^nextpow2(nsc));
spectrogram(swr,hamming(nsc),nov,nff,samplerate,'yaxis');
xticklabels( num2cell(t(get(gca,'XTick'))))
ylim([60 100])
title('SWR template -- TF')

% plot simulated swr
subplot(sub_r,sub_c,3)
plot(t,mean(v_simulated_swr,1),'','Color','r' );
hold on
plot(t,swr,'','Color','b' );
ha = plot(t,v_simulated_swr',':','Color',[166 166 166]/255 );
legend mean template all
title('simulated SWR events -- with noise ')
uistack(ha,'bottom')

% plot one example
subplot(sub_r,sub_c,4)
plot(t,v_simulated_swr(1,:)', 'Color','blue'  )
title('simulated SWR events -- one example')

% plot TF of simulated swr
subplot(sub_r,sub_c,5)
Nx = s.len;
nsc = floor(Nx/4.5);
nov = floor(nsc/2);
nff = max(256,2^nextpow2(nsc));
spectrogram(mean(v_simulated_swr,1),hamming(nsc),nov,nff,samplerate,'yaxis');
xticklabels( num2cell(t(get(gca,'XTick'))))
ylim([60 100])
title('simulated SWRs -- TF')

% plot simulated swr
subplot(sub_r,sub_c,6)
hold on
plot(t,mean(swr_template.bp_epochData,1),'Color','b')
plot(t,mean(SWRs_withnoise_real_bp,1),'Color','r')
title('band pass SWRs -- band pass')
legend template simulated 


% signal visualization
rg = 100; % mean(range(X));
subplot(sub_r,1,sub_r)
hold on
X_plt = X + repmat([1:nSensors]*rg,nSamples,1);
plot([1:nSamples]/samplerate,X_plt)
plot([max(events.epochRange(:,:,1)); max(events.epochRange(:,:,1)) ]./samplerate, [0 rg*(iSensors+1)],'Color','b')
plot([max(events.epochRange(:,:,2)); max(events.epochRange(:,:,2)) ]./samplerate, [0 rg*(iSensors+1)],'Color','r')
xlabel time(s)
yticklabels ''
xlim([epochRange(10,1)-2500,epochRange(10,2)+2500]/samplerate)
title('simulated signal (blue: starting of event     red: ending of event)  -- drag figure to see more signals')





