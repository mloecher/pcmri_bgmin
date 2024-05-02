%%
clear all
clc
addpath("..\..\matlab\")

%%
% load("./linord_skip_3200.mat");
load("./linord_skip_600.mat");
linord = linord_skip;
NN = numel(linord);

%%
% rootdir = "D:\Dropbox\dev\gropt3\python\notebooks\round2_invivo\waves_2023_0826_5b\";
rootdir = "C:\Users\mloecher\Dropbox\dev\gropt3\python\notebooks\round2_invivo\waves_2023_0826_5b\";
filelist = dir(fullfile(rootdir, 'v*'));

% ii = 1
% loop goes here
for ii = 1:numel(filelist)

venc = str2num(filelist(ii).name(2:4));
ro_asym = str2num(filelist(ii).name(10:13));
fov_in = str2num(filelist(ii).name(18:20));

seq_name = sprintf("v%03d_e%03d_fov%03d", venc, 100*ro_asym, fov_in);
fprintf('==== Seq Name: %s ====\n', seq_name);

%----------------------------------------
% set system limits
lims = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
               'MaxSlew', 80, 'SlewUnit', 'T/m/s', ...
               'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
               'adcDeadTime', 20e-6);


fov = fov_in * 1e-3;

Nx=128; Ny=128;      % Define FOV and resolution

Navg = 1;
gap_time = 0e-3;
slice_pos = 0;

ro_time = 2.5e-3;

flip = 10;
thickness=6.0e-3;  % slice  

rfSpoilingInc=117;  
deltak=1/fov;


[rf, gss] = mr.makeSincPulse(flip * pi/180,'system',lims,'Duration',1.2e-3,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);

gx_full = Nx*deltak;
Nadc = round(ro_asym*2*Nx);
gx = mr.makeTrapezoid('x','FlatArea',ro_asym*gx_full,'FlatTime',Nadc*10e-6,'system',lims);
adc = mr.makeAdc(Nadc,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',lims);

phaseAreas = ((0:Ny-1)-Ny/2)*deltak;

%-----------------------------------

seq=mr.Sequence(lims);
all_g_out = {};

for i_waves = 0:3
    
    mat = load(rootdir + filelist(ii).name + "\waves" + num2str(i_waves) + ".mat");
    all_g_out{i_waves+1} = mat.g_out;

    g_bipolar0 = mr.makeArbitraryGrad('z', lims.gamma .* mat.bp_out0, 'system',lims);
    g_bipolar1 = mr.makeArbitraryGrad('z', lims.gamma .* mat.bp_out1, 'system',lims);
    
    g_spoiler0 = mr.makeArbitraryGrad('z', lims.gamma .* mat.sp_out0, 'system',lims);
    g_spoiler1 = mr.makeArbitraryGrad('z', lims.gamma .* mat.sp_out1, 'system',lims);
    
    %-----------------------------------------
    % Readout prewinder
    pre_dur = 1.0e-3;
    
    gx_pre = mr.makeTrapezoid('x', 'Area', -gx_full * (ro_asym - 1/2), 'duration', pre_dur, 'system', lims);
    t_ro_mid = mr.calcDuration(gx_pre);
    
    pe_duration = pre_dur;
    % Phase encodes
    gy_pe = mr.makeTrapezoid('y', 'Area', phaseAreas(1), 'duration', pe_duration, 'system', lims);
    t_pe_mid = mr.calcDuration(gy_pe);
    
    % Slice select rephaser
    gss_re = mr.makeTrapezoid('z','Area',-gss.area/2,'system',lims,'duration', pre_dur);
    t_ss_mid = mr.calcDuration(gss_re) + mr.calcDuration(g_bipolar0);
    
%     fprintf('Middle Timings:  %.2f  %.2f  %.2f \n', t_ro_mid*1000, t_pe_mid*1000, t_ss_mid*1000)
    
    %-----------------------------------------
    
    % Calc end duration
    
    % Readout spoiler
    gx_end0 = mr.makeTrapezoid('x','Area',2*gx_full,'system',lims);
    
    % Phase encodes
    gy_end0 = mr.makeTrapezoid('y','Area',phaseAreas(1),'system',lims);
    
    
    % SS
    if isstruct(g_spoiler0)
        gz_end0 = g_spoiler0;
        gz_end1 = g_spoiler1;
    else
        gz_end0 = mr.makeTrapezoid('z','Area',4/thickness,'system',lims);
        gz_end1 = gz_end0;
    end
    
    
%     fprintf('End Timings:  %.2f  %.2f  %.2f \n', mr.calcDuration(gx_end0)*1000, mr.calcDuration(gy_end0)*1000, mr.calcDuration(gz_end0)*1000)
    
    end_dur = mr.calcDuration(gx_end0, gy_end0, gz_end0);
    
%     fprintf('End Duration:  %.2f \n', end_dur*1000)
    
    gx_end = mr.makeTrapezoid('x','Area',2*gx_full,'Duration',end_dur,'system',lims);
    gy_end = mr.makeTrapezoid('y','Area',phaseAreas(1),'Duration',end_dur, 'system',lims);
    if ~isstruct(g_spoiler0)
        gz_end0 = mr.makeTrapezoid('z','Area',4/thickness,'Duration',end_dur,'system',lims);
        gz_end1 = gz_end0;
    end
    
    %-----------------------------------------
    
    % Single actual sequence
    rf_phase=0;
    rf_inc=0;
    for i=1:NN
        for j = 1:2
            % RF and ADC phase cycling
            rf.phaseOffset=rf_phase/180*pi;
            adc.phaseOffset=rf_phase/180*pi;
            rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase=mod(rf_phase+rf_inc, 360.0);
            
            % Slice off iso-center
            rf.freqOffset = gss.amplitude * slice_pos;
    
            % RF and slice select gradient
            seq.addBlock(rf, gss);
            
            % Make bipolars
            if j == 1
                g_bp = g_bipolar0;
            elseif j == 2
                g_bp = g_bipolar1;
            end
            
            i_pe = linord(i);
    
            % Calc phase encodes
            i_gy_pe = mr.makeTrapezoid('y','Area',phaseAreas(i_pe),...
                                     'Duration',pe_duration, 'system',lims);
                                 
            i_gy_end = mr.makeTrapezoid('y','Area',-phaseAreas(i_pe),...
                                     'Duration',end_dur, 'system',lims);
    
            gy_pe_d = i_gy_pe;
            gx_pre_d = gx_pre;
            
            % Add all middle gradients
            seq.addBlock(gx_pre_d, gy_pe_d, gss_re);
            
            seq.addBlock(g_bp);
            
            % gap for eddy current testing
            if gap_time > 0
                seq.addBlock(mr.makeDelay(gap_time));
            end
            
            % ADC event
            seq.addBlock(gx, adc);
            
            if j == 1
                gz_end = gz_end0;
            elseif j == 2
                gz_end = gz_end1;
            end
            
            % Spoilers and refocusing
            seq.addBlock(gx_end, i_gy_end, gz_end);
            
            seq.addBlock(mr.makeDelay(40e-6));
    
        end
    end
    
    seq.addBlock(mr.makeDelay(2));

end

%----------------------------------------------

[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

[grad_waveforms, t_rf, t_adc] = seq.gradient_waveforms_mwl();
ts_rf = cell2mat(t_rf).';
ts_adc = cell2mat(t_adc.');
start_idx = round(squeeze(ts_rf(:,1) - ts_rf(1,1))/lims.gradRasterTime);


fprintf('Test grad diff RMSE =');
for jj = 0:3
    g_test = grad_waveforms(3,start_idx(1+jj*numel(start_idx)/4)+1:start_idx(3+jj*numel(start_idx)/4));
    diff_waves = g_test/lims.gamma - all_g_out{jj+1};
    test_grad_norm = norm(diff_waves);
    fprintf('  %.1e',test_grad_norm);
    % figure; hold on;
    % plot(g_test/lims.gamma)
    % plot(all_g_out{jj+1})
end
fprintf('\n',test_grad_norm);

fprintf('Acquisition time = %.1f \n', size(grad_waveforms,2) * 10e-6);

seq.setDefinition('FOV', [fov fov thickness]);
seq.setDefinition('Name', seq_name);

% seq.write('./seqs_single_2023_0827/' + seq_name + '.seq')
seq.write('./seqs_single_2023_0827q/' + seq_name + '.seq')

% end loop here
end % ii filelist loop
%%
[1000*pe_duration, 1000*mr.calcDuration(adc), 1000*mr.calcDuration(gx), 1000*end_dur] 
 [start_idx(1), start_idx(2), start_idx(3)]


%%

fprintf('Test grad diff RMSE =');
for jj = 0:3
g_test = grad_waveforms(3,start_idx(1+jj*numel(start_idx)/4)+1:start_idx(3+jj*numel(start_idx)/4));
diff_waves = g_test/lims.gamma - all_g_out{jj+1};
test_grad_norm = norm(diff_waves);
fprintf('  %.1e',test_grad_norm);
% figure; hold on;
% plot(g_test/lims.gamma)
% plot(all_g_out{jj+1})
end
fprintf('\n',test_grad_norm);

%%
size(mat.bp_out1)

%%

figure;
hold on;
plot(mat.bp_out0)
plot(mat.bp_out1)

figure;
hold on;
plot(diff(mat.bp_out0)/10e-6)
plot(diff(mat.bp_out1)/10e-6)

%%
figure;
hold on;
plot(mat.sp_out0)
plot(mat.sp_out1)

%%
% g_test = grad_waveforms(3,start_idx(end-3)+1:start_idx(end-1));
figure; hold on;
plot(g_test/lims.gamma)
plot(all_g_out{4})
%%



%%

seq_name = sprintf("v%03d_e%03d_fov%03d", venc, 100*ro_asym, fov_in)
