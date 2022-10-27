function [X, events]=inject_confound_HFB(X,swr_template, samplerate, nEvents, hfbfilter, noise, swr_events)
% input:
%       X            - signal, time by sensors matrix
%       swr_template - swr template data structure
%       sample rate  - sample rate
%       nEvents      - number of hfb events 
%       noise        - data structure control voltage of hfb and X respectively
%       swr_events   - swr onsets and swr signal with noise
% output:
%       X            - signal with hfb events
%       events       - hfb onsets and hfb signal with noise
% Author: Zhibing Xiao,2022/10/26

%%
nSamples  = size(X,1);
nSensors  = size(X,2);
t = swr_template.t;

% measure simulated swrs
swr = mean( squeeze( mean(swr_events.SWRs_simulated,1)),1);
s.len = length(swr);
s.min = min(swr) ;
s.max = max(swr);

% control voltage of X according to noise.backgroundnoise
X = X*noise.backgroudnoise;

% make rand signal
rnd_singal = randn(s.len*nEvents,1);

% make hfb signals by filter
fs = samplerate;
fl = hfbfilter.lowedge;
fh = hfbfilter.highedge;
order=256;
hfb = simul_filter(rnd_singal,fs,fl,fh,order);

hfb = reshape(hfb,[],nEvents)';

% scale voltage of hfb 
[hfb_scale, ps]= mapminmax(hfb ,s.min,s.max);
 
% location of SWR events, leave the start and end of the signal free
peak_loc = sort(randi([s.len+samplerate*0.05, nSamples-s.len-samplerate*0.05], nEvents, 1),'ascend');

%% set hfbs
for iSensors = 1  :nSensors
    
    %  set event domain noise
    voltage_level = repmat(randn(nEvents,1),1,s.len);
       
    %
    HFBs = hfb_scale.*voltage_level;

    %
    i_peak_loc  = peak_loc + randi([-samplerate*0.05 samplerate*0.05],nEvents,1);
    Peak_epoch_str = i_peak_loc-(s.len-1)/2;
    Peak_epoch_fin = i_peak_loc+(s.len-1)/2;
    
    epochRange = [Peak_epoch_str,Peak_epoch_fin];

    signal = X(:,iSensors);

    for i=1:size(epochRange,1)
        voltage_noise = signal(epochRange(i,1):epochRange(i,2))' ;
        HFBs(i,:) = HFBs(i,:)* noise.weight_of_voltage  + voltage_noise;
        signal(epochRange(i,1):epochRange(i,2)) = HFBs(i,:);
    end

    X(:,iSensors) = signal;
    HFBs_withnoise_real      = zeros(size(epochRange,1), s.len);
    for i=1:size(epochRange,1)
        HFBs_withnoise_real(i,:) = signal(epochRange(i,1):epochRange(i,2));
    end
    events.HFBs_simulated(iSensors,:,:) =  HFBs_withnoise_real;
    events.epochRange(iSensors,:,:) = epochRange;

end


%% compute spectrogram
Nx = s.len;
nsc = floor(Nx/4.5);
nov = floor(nsc/2);
nff = max(256,2^nextpow2(nsc));
for i=1:size(epochRange,1)
    HFBs_withnoise_real(i,:) = signal(epochRange(i,1):epochRange(i,2));
    [~,tf,tt,p] = spectrogram(HFBs_withnoise_real(i,:) ,hamming(nsc),nov,nff,samplerate,'yaxis','psd');
    hps(i,:,:) = 10*log10(p);
end



%% visulization
sub_r = 3;
sub_c = 2;

figure('Position',[20,20,1000,1000]);
subplot(sub_r,sub_c,1)
plot(t,hfb_scale(1,:),'Color','blue')
title('hfb template -- one example')
xlim([-500 500])

subplot(sub_r,sub_c,2)
Nx = s.len;
nsc = floor(Nx/4.5);
nov = floor(nsc/2);
nff = max(256,2^nextpow2(nsc));
[~,tf,tt,p] = spectrogram(hfb_scale(1,:),hamming(nsc),nov,nff,samplerate,'yaxis','psd');
imagesc(t,tf, 10*log10(p) )
set(gca, 'YDir', 'normal');
c=colorbar;
c.Label.String = 'PSD';
title('hfb template -- TF')
 

subplot(sub_r,sub_c,3)
plot(t,HFBs',':','Color',[166 166 166]/255 );
hold on
plot(t,mean(HFBs_withnoise_real,1),'','Color','r' );
title('hfb events -- with noise(gray:all  red:mean  blue:template)')

subplot(sub_r,sub_c,4)
imagesc(t,tf,squeeze( mean(hps,1)) )
set(gca, 'YDir', 'normal');
c=colorbar;
c.Label.String = 'PSD';
title('hfb events -- TF')

rg = 100;% mean(range(X));
subplot(3,1,3)
X_plt = X + repmat([1:nSensors]*rg,nSamples,1);
plot([1:nSamples]/samplerate,X_plt)
xlabel time(s)
xlim([epochRange(10,1)-2500,epochRange(10,2)+2500]/samplerate)
ylim([0  (nSensors+1)*rg])
hold on
plot([max(events.epochRange(:,:,1)); max(events.epochRange(:,:,1)) ]./samplerate, [0 (nSensors+1)*rg],'Color','b')
plot([max(events.epochRange(:,:,2)); max(events.epochRange(:,:,2)) ]./samplerate, [0 (nSensors+1)*rg],'Color','r')
title('simulated signal (blue: starting of event     red: ending of event)  -- drag figure to see more signals')
 



