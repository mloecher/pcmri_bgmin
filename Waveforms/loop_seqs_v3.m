%%
clear all
clc
addpath("..\..\matlab\")

%%
% all_fnames = dir("D:\Dropbox\dev\gropt3\python\notebooks\best_waves_2023_0423\*.mat");
% all_fnames = dir("C:\Users\mloecher\Dropbox\dev\gropt3\python\notebooks\best_waves_2023_0423\*.mat");
% all_fnames = dir("C:\Users\mloecher\Dropbox\dev\gropt3\python\notebooks\round2_slidingte_waves\waves_2023_0719_4\*.mat");
all_fnames = dir("D:\Dropbox\dev\gropt3\python\notebooks\round2_slidingte_waves\waves_2023_0719_5\*.mat")


for ii = 1:numel(all_fnames)

wave_dir = all_fnames(ii).name(3);
suffix = all_fnames(ii).name(5:end-4);
seq_name = all_fnames(ii).name(1:end-4)


% set system limits
lims = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
               'MaxSlew', 80, 'SlewUnit', 'T/m/s', ...
               'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
               'adcDeadTime', 20e-6);


fov=320e-3; Nx=128; Ny=128;      % Define FOV and resolution

Navg = 1;
gap_time = 0e-3;
slice_pos = 0;

TR = 12e-3;

ro_asym = 3/4;
ro_time = 2.5e-3;

flip = 10;
thickness=6.0e-3;  % slice  

rfSpoilingInc=117;  
deltak=1/fov;


[rf, gss] = mr.makeSincPulse(flip * pi/180,'system',lims,'Duration',1.2e-3,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4);

gx_full = Nx*deltak;
gx = mr.makeTrapezoid('x','FlatArea',ro_asym*gx_full,'FlatTime',ro_asym*2*Nx*10e-6,'system',lims);
adc = mr.makeAdc(ro_asym*2*Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',lims);

phaseAreas = ((0:Ny-1)-Ny/2)*deltak;


mat = load([all_fnames(ii).folder '\' all_fnames(ii).name]);

g_bipolar0 = mr.makeArbitraryGrad('z', lims.gamma .* mat.bp_out0, 'system',lims);
g_bipolar1 = mr.makeArbitraryGrad('z', lims.gamma .* mat.bp_out1, 'system',lims);

g_spoiler0 = mr.makeArbitraryGrad('z', lims.gamma .* mat.sp_out0, 'system',lims);
g_spoiler1 = mr.makeArbitraryGrad('z', lims.gamma .* mat.sp_out1, 'system',lims);


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

fprintf('Middle Timings:  %.2f  %.2f  %.2f \n', t_ro_mid*1000, t_pe_mid*1000, t_ss_mid*1000)



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


fprintf('End Timings:  %.2f  %.2f  %.2f \n', mr.calcDuration(gx_end0)*1000, mr.calcDuration(gy_end0)*1000, mr.calcDuration(gz_end0)*1000)

end_dur = mr.calcDuration(gx_end0, gy_end0, gz_end0);
if end_dur < 1.1e-3
    end_dur = 1.1e-3;
end

fprintf('End Duration:  %.2f \n', end_dur*1000)

gx_end = mr.makeTrapezoid('x','Area',2*gx_full,'Duration',end_dur,'system',lims);
gy_end = mr.makeTrapezoid('y','Area',phaseAreas(1),'Duration',end_dur, 'system',lims);
if ~isstruct(g_spoiler0)
    gz_end0 = mr.makeTrapezoid('z','Area',4/thickness,'Duration',end_dur,'system',lims);
    gz_end1 = gz_end0;
end


seq=mr.Sequence(lims);

% for slice_pos = [0, 15e-3, 30e-3]
for gap_time = ((0:40) * 40e-6)
for ia = 1:Navg
    
    % Single actual sequence
    rf_phase=0;
    rf_inc=0;
    for i=1:Ny
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
            
            % Calc phase encodes
            i_gy_pe = mr.makeTrapezoid('y','Area',phaseAreas(i),...
                                     'Duration',pe_duration, 'system',lims);
                                 
            i_gy_end = mr.makeTrapezoid('y','Area',-phaseAreas(i),...
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
            
            TR_delay = TR - mr.calcDuration(gss) - 1e-3 - mr.calcDuration(g_bp) - gap_time - mr.calcDuration(gx) - mr.calcDuration(gx_end);
            
            end_delay = 10e-6*(double(mat.N/2) - double(mat.sp_idx1));

            if TR_delay > end_delay
                seq.addBlock(mr.makeDelay(TR_delay-end_delay));
                TR_delay = end_delay;
            end
            
            if j == 1
                gz_end = gz_end0;
            elseif j == 2
                gz_end = gz_end1;
            end
            
            % Spoilers and refocusing
            seq.addBlock(gx_end, i_gy_end, gz_end);
            
            if TR_delay > 0
                seq.addBlock(mr.makeDelay(TR_delay));
            end

            
            
            
        end
    end
end

end  % end gap time
% end


[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

seq.setDefinition('FOV', [fov fov thickness]);
seq.setDefinition('Name', seq_name);

seq.write(['./seqs_2023_0719_5/' seq_name '.seq'])


[grad_waveforms, t_rf, t_adc] = seq.gradient_waveforms_mwl();
ts_rf = cell2mat(t_rf).';
ts_adc = cell2mat(t_adc.');
start_idx = round(squeeze(ts_rf(:,1) - ts_rf(1,1))/lims.gradRasterTime);
g_test = grad_waveforms(3,start_idx(1)+1:start_idx(3));
diff_waves = g_test/lims.gamma - mat.g_out;
test_grad_norm = norm(diff_waves);

fprintf('Test grad diff RMSE = %.2e \n',test_grad_norm);
    
if (test_grad_norm > 1e-4)
   fprintf('*************************************************\n')
   fprintf('*************************************************\n')
   fprintf('       Test grad norm error ! ! ! \n')
   fprintf('*************************************************\n')
   fprintf('*************************************************\n')
end


end % End file loop