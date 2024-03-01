%%
clear all
clc
addpath('../../matlab/')
% addpath('D:\Dropbox\projects\matlab\pulseq\matlab\seq_girf\slice_select')
addpath('C:\Users\mloecher\Dropbox\projects\matlab\pulseq\matlab\seq_girf\slice_select')

% addpath('../../../bloch/')
% addpath('../../../girf_stuff/waveforms/')
% addpath('../seq_girf/slice_select/')

%%

% set system limits
lims = mr.opts('MaxGrad', 32, 'GradUnit', 'mT/m', ...
               'MaxSlew', 160, 'SlewUnit', 'T/m/s', ...
               'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 40e-6);
             

TR = 1000e-3;

dt_adc = 8e-6;
N_adc = 8000;
t_adc_dur = N_adc*dt_adc;  % This is just now for reference, make sure it has no decimal I think? 
adc = mr.makeAdc(N_adc, 'Duration', t_adc_dur, 'system', lims);

% Skope trigger
mr_trig = mr.makeDigitalOutputPulse('ext1','duration', lims.gradRasterTime);
t_trigger = 400e-6;

delay_TR = TR - mr.calcDuration(adc) - t_trigger;

%%

dirs = 'xyz';
all_type = {'blip', 'biblip', 'chirp', 'rand'};

seq=mr.Sequence(lims);

N_sync = 20;
[g_blip, g_blip_re] = get_test_wave('blip', 0, lims, 150, 1);
g_blip.delay = 100e-6;
for i = 1:N_sync
    seq.addBlock(mr_trig, mr.makeDelay(t_trigger));
    seq.addBlock(g_blip, adc, mr.makeDelay(mr.calcDuration(adc)+20e-6));
    seq.addBlock(mr.makeDelay(delay_TR));
end

seq.addBlock(mr.makeDelay(3.0));


for i_type = 1:numel(all_type)
    type = all_type{i_type};
for polarity = [-1, 1]
for idx = 0:8
    [g_blip, g_blip_re] = get_test_wave(type, idx, lims, 150, polarity);
    g_blip.delay = 100e-6;
for i_dir = 1:3 
    g_blip.channel = dirs(i_dir);
    
    seq.addBlock(mr_trig, mr.makeDelay(t_trigger));
    seq.addBlock(g_blip, adc, mr.makeDelay(mr.calcDuration(adc)+20e-6));
    seq.addBlock(mr.makeDelay(delay_TR));
    
end  % i_dir
end  % idx
end  % polarity
end  % i_type

%%

[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%%
[grad_waveforms, t_rf, t_adc] = seq.gradient_waveforms_mwl();

%%
%-------  Get TR start times and the peak rf times (excitation)
ts_rf = cell2mat(t_rf).';
ts_adc = cell2mat(t_adc.');

start_idx = round(squeeze(ts_adc(:,1) - t_trigger)/lims.gradRasterTime);

fprintf('N per TR = %d \n', size(grad_waveforms, 2)/size(ts_rf,1));
fprintf('start_idx(1:10) =  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d \n', start_idx(1:10))
fprintf('Scan time = %.2f min\n',size(grad_waveforms, 2)*lims.gradRasterTime/60);


%%

figure; hold on;
i_tr = 220;
plot(grad_waveforms(1,start_idx(i_tr)+1:start_idx(i_tr+2)), 'LineWidth', 2)
plot(grad_waveforms(2,start_idx(i_tr)+1:start_idx(i_tr+2)), 'LineWidth', 2)
plot(grad_waveforms(3,start_idx(i_tr)+1:start_idx(i_tr+2)), 'LineWidth', 2)

%%

seq.setDefinition('FOV', [100e-3 100e-3 100e-3]);
seq_name = sprintf('girf_skope2');
seq.setDefinition('Name', seq_name);

seq.write([seq_name '.seq'])       % Write to pulseq file

%%

%%

load('../../../girf_stuff/waveforms/wave_save_v2.mat')