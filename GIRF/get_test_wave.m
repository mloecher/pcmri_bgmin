function [g_blip, g_blip_re] = get_test_wave(type, idx, lims, slew, polarity)

load('./wave_save_v2.mat')
load('./bipolars.mat')


if idx == 0
    
    g_blip = mr.makeTrapezoid('z','system',lims,'flatTime',0, 'riseTime', 10e-6, 'fallTime', 10e-6, 'amplitude', 0);
    g_blip_re = mr.makeTrapezoid('z','system',lims,'flatTime',0, 'riseTime', 10e-6, 'fallTime', 10e-6, 'amplitude', 0);
    
else

    if strcmpi(type, 'blip')

        t_ramp = 60e-6 + 10e-6*idx;
        g_blip = mr.makeTrapezoid('z','system',lims,'flatTime',0, 'riseTime', t_ramp, 'fallTime', t_ramp, ...
                             'amplitude', polarity*lims.gamma*slew*t_ramp);

        g_blip_re = g_blip;
        g_blip_re.amplitude = -g_blip_re.amplitude;


    elseif strcmpi(type, 'biblip')

        t_ramp = 60e-6 + 10e-6*idx;
        g_blip = mr.makeTrapezoid('z','system',lims,'flatTime',0, 'riseTime', t_ramp, 'fallTime', t_ramp, ...
                             'amplitude', polarity*lims.gamma*slew*t_ramp);  

        g_blip_neg = g_blip;
        g_blip_neg.amplitude = -g_blip_neg.amplitude;
        g_blip_neg.delay = mr.calcDuration(g_blip);

        g_blip = mr.addGradients({g_blip, g_blip_neg}, lims);

        g_blip_re = mr.makeTrapezoid('z','system',lims,'flatTime',0, 'riseTime', 10e-6, 'fallTime', 10e-6, 'amplitude', 0);

    elseif strcmpi(type, 'chirp')

        blip = all_blip_chirp(idx,:);
        g_blip = mr.makeArbitraryGrad('z', polarity*lims.gamma .* blip, 'system',lims);

        g_blip_re = mr.makeTrapezoid('z','Area',-polarity*lims.gamma.*sum(blip)*lims.gradRasterTime,'system',lims);

    elseif strcmpi(type, 'rand')

        blip = all_blip_rand(idx,:);
        g_blip = mr.makeArbitraryGrad('z',polarity*lims.gamma .* blip, 'system',lims);

        g_blip_re = mr.makeTrapezoid('z','Area',-polarity*lims.gamma.*sum(blip)*lims.gradRasterTime,'system',lims);



    elseif strcmpi(type,'flow80')
        
        if idx == 1
            vel_res = 80;  % venc of scan
            vel_res = vel_res * 2;  % For two sided encoding

            T_kv = 0.5e-3;
            d_T = 0.01e-3;
            good_grad = 0;

            while good_grad < 1
                try
                    max_area_kv = 100 / vel_res / T_kv / 2;
                    g_bp = mr.makeTrapezoid('z', lims, 'Area', -max_area_kv, ...
                    'Duration', T_kv, 'MaxGrad', lims.gamma*32, 'MaxSlew', lims.gamma*100);
                    good_grad = 1;
                catch
                    T_kv = T_kv + d_T;
                end
            end

            g_blip = g_bp;
            g_blip.amplitude = polarity * g_blip.amplitude;

            g_blip_neg = g_blip;
            g_blip_neg.amplitude = -g_blip_neg.amplitude;
            g_blip_neg.delay = mr.calcDuration(g_blip);


            g_blip = mr.addGradients({g_blip, g_blip_neg}, lims);
        elseif idx == 2
            
            blip = g_save_reg;
            g_blip = mr.makeArbitraryGrad('z',polarity*lims.gamma .* blip, 'system',lims);

        elseif idx == 3
            blip = g_save_150;
            g_blip = mr.makeArbitraryGrad('z',polarity*lims.gamma .* blip, 'system',lims);
        elseif idx == 4
            blip = g_save_75;
            g_blip = mr.makeArbitraryGrad('z',polarity*lims.gamma .* blip, 'system',lims);
        end

        g_blip_re = mr.makeTrapezoid('z','system',lims,'flatTime',0, 'riseTime', 10e-6, 'fallTime', 10e-6, 'amplitude', 0);



    end
    
end

end

