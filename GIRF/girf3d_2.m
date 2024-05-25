%%
clear all
clc

% Addpath to pulseq if needed
% addpath('D:\Dropbox\projects\matlab2\pulseq142\matlab')

%%

% We make a separate sequence for 
for fix_i_type = 1:4

% set system limits
lims = mr.opts('MaxGrad', 32, 'GradUnit', 'mT/m', ...
               'MaxSlew', 160, 'SlewUnit', 'T/m/s', ...
               'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 40e-6);

app = 0.4;
tbw = 4;
thickness = 1.0e-3;
flip = 30;

[rf, gss] = mr.makeSincPulse(flip * pi / 180,'system',lims,'Duration',4e-3,...
    'SliceThickness',thickness,'apodization',app,'timeBwProduct',tbw);

max_rf = max(rf.signal);

gss_re = mr.makeTrapezoid('z','Area',-gss.area/2,'system',lims);


t_grad_delay = 4e-3;
t_adc_start = 10e-3;
t_adc_dur = 64e-3;
TR = 200e-3;

dt_adc = 8e-6;
N_adc = round(t_adc_dur/dt_adc);
adc = mr.makeAdc(N_adc, 'Duration', t_adc_dur, 'system', lims);


% ------------------
% PE

fov=300e-3; 
Npe=5;  % 5x5 phase encoding

deltak=1/fov;
phaseAreas = ((0:Npe-1)-Npe/2+0.5)*deltak;

g_pe_1 = mr.makeTrapezoid('x', 'Area', phaseAreas(1), 'system', lims);
g_pe_2 = mr.makeTrapezoid('y', 'Area', phaseAreas(1), 'system', lims);
t_pe = mr.calcDuration(g_pe_1);
pe_delay = 0.1e-3;

% ---------------
% Seq

all_type = {'blip', 'biblip', 'chirp', 'rand', 'flow80'};


fix_type = all_type{fix_i_type}

slew = lims.maxSlew/lims.gamma;
all_dir = {'zxy', 'xyz', 'yxz'};

all_spos = 20e-3 .* [-1  -0.5  0  0.5  1];  % 5 slices with 1cm spacing
% all_spos = 20e-3 .* [1];

Nslice = numel(all_spos);

seq=mr.Sequence(lims);


% x and y encoding loops
for ipe2 = 1:Npe
    i_g_pe_2 = mr.makeTrapezoid('y','Area',phaseAreas(ipe2),...
                             'Duration',t_pe, 'system', lims, 'delay', pe_delay);
                         
%     label_pe2 = mr.makeLabel('SET','PAR', ipe2-1);
                         
for ipe1 = 1:Npe
    i_g_pe_1 = mr.makeTrapezoid('x','Area',phaseAreas(ipe1),...
                             'Duration',t_pe, 'system', lims, 'delay', pe_delay);
                         
%     label_pe1 = mr.makeLabel('SET','LIN', ipe1-1);
    

for polarity = [-1, 1]
    
for i_type = fix_i_type:fix_i_type
    type = all_type{i_type};
for idx = 0:8
if strcmpi(type, 'flow80') && (idx > 4)
    continue  % There are fewer of these to measure
end

    % label_wave_type = mr.makeLabel('SET','ECO', i_type-1);
    % label_wave_idx = mr.makeLabel('SET','SET', idx);

    [g_blip, g_blip_re] = get_test_wave(type, idx, lims, slew, polarity);
    g_blip.delay = t_grad_delay;
    
%     label_polarity = mr.makeLabel('SET','PHS', round((polarity+1)/2));
    
for i_dir = 1:3  % Gradient exis to play on
    
%     label_dir = mr.makeLabel('SET','AVG', i_dir-1);
    
    gss.channel = all_dir{i_dir}(1);
    gss_re.channel = all_dir{i_dir}(1);
    
    g_blip.channel = all_dir{i_dir}(1);
    g_blip_re.channel = all_dir{i_dir}(1);
    
    i_g_pe_1.channel = all_dir{i_dir}(2);
    i_g_pe_2.channel = all_dir{i_dir}(3);

for i_slice = 1:Nslice
    
%     label_slice = mr.makeLabel('SET','SLC', i_slice-1);

    spos = all_spos(i_slice);
    rf.freqOffset = gss.amplitude * spos;

    cur_time = 0;
    
    seq.addBlock(gss_re);
%     seq.addBlock(gss_re, mr.makeLabel('SET','LIN', ii));
%     ii = ii + 1;
    cur_time = cur_time + mr.calcDuration(gss_re);
     
    % == Slice Select:
    seq.addBlock(rf, gss);
%     seq.addBlock(rf, gss, mr.makeLabel('SET','LIN', ii));
%     ii = ii + 1;
    
    cur_time = cur_time + mr.calcDuration(rf, gss);
    
    % == Phase Encode and SS refocus:
    seq.addBlock(gss_re, i_g_pe_1, i_g_pe_2);
    cur_time = cur_time + mr.calcDuration(gss_re, i_g_pe_1, i_g_pe_2);
    
%     seq.addBlock(label_pe2);
%     seq.addBlock(label_pe2, label_pe1);
%     seq.addBlock(label_pe2, label_pe1, label_wave_type, label_wave_idx, label_polarity, ...
%                label_dir, label_slice);
    
    delay_adc_start = mr.makeDelay(t_adc_start - cur_time);
    seq.addBlock(delay_adc_start);
    cur_time = cur_time + mr.calcDuration(delay_adc_start);
    
     % == ADC and test waveform:
    seq.addBlock(adc, g_blip);
    cur_time = cur_time + mr.calcDuration(adc, g_blip);
    
    delay_refocus = mr.makeDelay(1e-3);
    seq.addBlock(delay_refocus);
    cur_time = cur_time + mr.calcDuration(delay_refocus);
    
    % == Refocus the test wave and the phase encode
    i_g_pe_1.amplitude = -i_g_pe_1.amplitude;
    i_g_pe_2.amplitude = -i_g_pe_2.amplitude;
    seq.addBlock(g_blip_re, i_g_pe_1, i_g_pe_2);
    cur_time = cur_time + mr.calcDuration(g_blip_re, i_g_pe_1, i_g_pe_2);
    
    % == Revert the PE lobes back to normal for next loop
    i_g_pe_1.amplitude = -i_g_pe_1.amplitude;
    i_g_pe_2.amplitude = -i_g_pe_2.amplitude;
    
    delay_TR = mr.makeDelay(TR - cur_time);
    seq.addBlock(delay_TR);
    cur_time = cur_time + mr.calcDuration(delay_TR);

    
end  % i_slice

end  % i_dir

end  % polarity 

end  % wave idx
end  % wave type

end  % ipe1
end  % ipe2



%---------- check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
    
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%---------- prepare sequence export

seq.setDefinition('FOV', [20e-3 20e-3 20e-3]);
seq_name = sprintf('girf3d2_%s',fix_type);
seq.setDefinition('Name', seq_name);

seq.write([seq_name '.seq'])       % Write to pulseq file

fprintf('SAVED %s !!!\n', seq_name);

end


