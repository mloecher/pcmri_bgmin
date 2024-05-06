import numpy as np
from scipy.signal.windows import kaiser
import scipy.signal

def fft_offset(v, dt, offset = None):
    f_meas = np.fft.fftfreq(len(v), dt)
    v2 = np.fft.fft(v)
    v2 *= np.exp(-1j * f_meas * offset * 2 * np.pi)
    v = np.real(np.fft.ifft(v2))
    return v

def get_b0_comp(gg, y, dt_gg = 10e-6, dt_cam = 1e-6, wave_dir = 0, new_vals = False):

    if new_vals:
        amps = [[0.004017448984090,0.107122369111000,0.014783694408800,],
            [0.139922320843000,0.054929979145500,0.023704057559400,],
            [0.361955434084000,0.060874473303600,-0.124086134136000,]]

        t_constants = [[1.999913334850000,0.153517320752000,0.033818759024100,],
            [0.139710530639000,0.011686773970700,0.000876913836692,],
            [0.166305392981000,0.054901216179100,0.028627114370500,]]

        # print(amps)

    else:
        amps = [[0.00250313989818,0.0984266996384,0.0248604007065,],
        [ 0.142157003284,0.0598868988454,0.0427879989147,],
        [0.0585108995438,0.339473009109,-0.0912802964449,]]

        t_constants = [[1.99994003773,0.167437002063,0.0499899014831,],
        [0.141018003225, 0.0111469998956,0.000501247006468,],
        [0.262499988079, 0.148469001055,0.0262141004205,]]

        # print(amps)

    Nconvolve = 200000
    tt = np.arange(Nconvolve) * dt_gg

    b0_resp_x = np.zeros_like(tt)
    b0_resp_y = np.zeros_like(tt)
    b0_resp_z = np.zeros_like(tt)

    for i in range(3):
        b0_resp_x += amps[0][i] * np.exp(-tt/t_constants[0][i])
        b0_resp_y += amps[1][i] * np.exp(-tt/t_constants[1][i])
        b0_resp_z += amps[2][i] * np.exp(-tt/t_constants[2][i])

    tt_gg = np.arange(gg.size) * dt_gg
    tt_y = np.arange(y.size) * dt_cam

    # if new_vals:
    #     if wave_dir == 0:
    #         b0_resp = b0_resp_x
    #     elif wave_dir == 1:
    #         b0_resp = b0_resp_y
    #     elif wave_dir == 2:
    #         b0_resp = b0_resp_z
    #     else:
    #         print("ERROR: wavedir not right")
    # else:
    if wave_dir == 0:
        b0_resp = b0_resp_z
    elif wave_dir == 1:
        b0_resp = b0_resp_x
    elif wave_dir == 2:
        b0_resp = b0_resp_y
    else:
        print("ERROR: wavedir not right")
    

    b0z = scipy.signal.convolve(np.gradient(gg)/dt_gg, b0_resp, 'full')[:gg.size]
    b0z_i = np.interp(tt_y, tt_gg, b0z)

    return 2*b0z_i

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