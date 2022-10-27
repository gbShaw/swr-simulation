%

clear
close all


%% parameters
nSensors = 5;
nSeconds = 600; % 600 seconds of data
% nSubj = 2; % number of subjects to simulate
samplerate = 1000; % 1000 time points per second
nSamples = nSeconds*samplerate;
swr_template_p = 'SWRtemplate.mat';
n_SWR_Events  = 500; % number of SWR events injecting to simulated data
n_HFB_Events  = 500; % number of high frequency band events injecting to simulated data
swrfilter.lowedge  = 70; % hz,  
swrfilter.highedge = 180; % hz
swr_noise.weight_of_voltage  = 1 ;  % control voltage of swr
swr_noise.backgroudnoise = 3 ;      % control voltage of signal
hfb_noise.weight_of_voltage = 0.5;  % control voltage of hfb
hfb_noise.backgroudnoise = 1;       % control voltage of signal
hfbfilter = swrfilter; % frequency range of hfb, same as swr

%% load template SWR
swr_template = load(swr_template_p);
swr_template.swrfilter = swrfilter;


%% generate dependence of the sensors
A = randn(nSensors);
[U,~] = eig((A+A')/2);
covMat = U*diag(abs(randn(nSensors,1)))*U';

%% make data
X = nan(nSamples, nSensors);
X(1,:) = randn([1 nSensors]);

for iT=2:nSamples
    X(iT,:) = 0.95*(X(iT-1,:) + mvnrnd(zeros(1,nSensors), covMat));% add dependence of the sensors
end


%% Injecting SWR events
[X1, swr_events ] = inject_swr_meanTemplate(X, swr_template, samplerate, n_SWR_Events, swr_noise  );
sgtitle('Step 1: Adding SWRs')

%% Injecting confounds
%% Injecting high frequency band events
[X2, hfb_events]=inject_confound_HFB(X1,swr_template, samplerate, n_HFB_Events, hfbfilter, hfb_noise, swr_events);
X_simulated = X2;
sgtitle('Step 2: Adding HFBs')

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualization for simulated SWR
Nx = length(swr_template.t);
nsc = floor(Nx/4.5);
nov = floor(nsc/2);
nff = max(256,2^nextpow2(nsc));

% extract SWRs according to events onset
epochRange = swr_events.epochRange;
SWRs=[];
for iSensors = 1:nSensors
    for i=1:size(epochRange,2)
        signal = X_simulated(:,nSensors);        
        swr = signal(epochRange(iSensors,i,1):epochRange(iSensors,i,2));        
        [~,tf,tt,p] = spectrogram(swr ,hamming(nsc),nov,nff,samplerate,'yaxis','psd');
        lp = 10*log10(p);

        % Drop the swr that overlaps with hfb
        if ~isempty(intersect(epochRange(iSensors,i,1):epochRange(iSensors,i,2),...
                hfb_events.epochRange(iSensors,i,1):hfb_events.epochRange(iSensors,i,2)))
            swr = swr*nan;
            lp = lp*nan;
        end

        SWRs(iSensors,i,:) = swr;
        SWRs_tf(iSensors, i,:,:) = lp;
    end
end

swr = squeeze (nanmean(SWRs(iSensors,:,:),[1 2]));
figure('Position',[100 100 600 200])
subplot(1,2,1)
plot(swr_template.t,mean(swr_template.epochData,1)-mean(mean(swr_template.epochData,1),2),...
    'color','b') 
hold on 
plot(swr_template.t,swr,'Color','r')
title('mean SWRs and template')
legend template simulated

subplot(1,2,2)
imagesc(swr_template.t,tf,  squeeze( nanmean(SWRs_tf,[1 2])))
set(gca, 'YDir', 'normal');
ylim([60 100])
c=colorbar;
c.Label.String = 'PSD';
title('simulated SWRs -- TF')
sgtitle('Simulated SWRs')


