def get_skope_girf(skope_res, wave_types = [0,], shift = 358, weiner_reg = 1e-8, offset = 3.55805e-04, dt_gg = 10e-6, dt_cam = 1e-6, 
                   do_filt=0, wave_dir = 0, fit_dir = 3, new_vals = False):

    all_gg = skope_res['all_grad'][:,0,:,0,2]
    all_cam = skope_res['all_cam'][:,:,:,wave_dir,:]

    all_X = []
    all_Y0 = []
    all_Y1 = []

    all_gg_out = []
    all_y0 = []
    all_y1 = []

    if np.ndim(wave_types) == 0:
        wave_types = [wave_types,]
    
    for i_type in wave_types:
        for i_wave in range(all_gg.shape[1]):
            
            gg = 1e3*all_gg[i_type, i_wave]
            gg = np.hstack([np.zeros(shift), gg])
            gg = fft_offset(gg, dt_gg, -offset)

            wave1a = 267e6*all_cam[i_type,0,i_wave,0]
            wave1b = 267e6*all_cam[i_type,1,i_wave,0] 
            y0 = (wave1a - wave1b)/2
            y0 = np.hstack([np.zeros(10*shift), y0])

            wave1a = 267e6*all_cam[i_type,0,i_wave,fit_dir]
            wave1b = 267e6*all_cam[i_type,1,i_wave,fit_dir] 
            y1 = (wave1a - wave1b)/2
            y1 = np.hstack([np.zeros(10*shift), y1])

            # The saved gradients are too long
            gg = gg[:y0.size//10]

            y0 += get_b0_comp(gg*1e-3, y0, wave_dir=wave_dir, new_vals=new_vals)
            if fit_dir == 0:
                y1 += get_b0_comp(gg*1e-3, y1, wave_dir, new_vals=new_vals)


            # # Filtering
            if do_filt:
                Nf0 = 1000
                Nf1 = int(do_filt)
                filt0 = kaiser(2*Nf0, 12)
                filt1 = kaiser(2*Nf1, 12)

                filt0 = filt0[:Nf0]
                filt1 = filt1[Nf1:]

                filt = np.hstack([filt0, np.ones(y1.size-Nf0-Nf1), filt1])
                y0 *= filt
                y1 *= filt

            # gg = np.hstack([gg, 0])
            # y0 = np.hstack([y0, np.zeros(10)])
            # y1 = np.hstack([y1, np.zeros(10)])

            X = np.fft.ifftshift(np.fft.fft(gg))
            Y0 = np.fft.ifftshift(np.fft.fft(y0))
            Y1 = np.fft.ifftshift(np.fft.fft(y1))

            all_gg_out.append(gg)
            all_y0.append(y0)
            all_y1.append(y1)

            dw = (Y0.size - X.size)/2
            Y0 = Y0[int(np.floor(dw)):-int(np.ceil(dw))] * X.size / Y0.size
            Y1 = Y1[int(np.floor(dw)):-int(np.ceil(dw))] * X.size / Y1.size

            all_X.append(X)
            all_Y0.append(Y0)
            all_Y1.append(Y1)

    all_X = np.array(all_X)
    all_Y0 = np.array(all_Y0)
    all_Y1 = np.array(all_Y1)
    all_gg_out = np.array(all_gg_out)
    all_y0 = np.array(all_y0)
    all_y1 = np.array(all_y1)

    H0 = (np.conj(all_X)*all_Y0).sum(0) / ( (np.abs(all_X)**2.0).sum(0) + weiner_reg)
    H1 = (np.conj(all_X)*all_Y1).sum(0) / ( (np.abs(all_X)**2.0).sum(0) + weiner_reg)
    f_H = np.fft.fftshift(np.fft.fftfreq(H0.size, dt_gg))

    all_waves = {'all_gg_out':all_gg_out,
                 'all_y0':all_y0,
                 'all_y1':all_y1,
                 'all_X':all_X,
                 'all_Y0':all_Y0,
                 'all_Y1':all_Y1
                 }

    return f_H, H0, H1, all_waves