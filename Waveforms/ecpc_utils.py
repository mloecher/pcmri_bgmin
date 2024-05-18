import sys 
sys.path.append('D:\\Dropbox\\dev\\gropt2_proj\\python\\')
import gropt2

import cv2

from scipy.io import loadmat
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import windows

def wave_plotter(g_out, res, g_out_ref = None, res_ref = None, return_im = False, title = None):

    Nax = 3
    if 'H1' in res:
        Nax += 2
        
    fig, axs = plt.subplots(Nax, 1, figsize=(12,8))
    
    ii = 0
    axs[ii].axhline(0, lw=1, alpha=0.3, color='k')
    axs[ii].plot(g_out)
    axs[ii].spines['right'].set_visible(False)
    axs[ii].spines['top'].set_visible(False)
    axs[ii].axvspan(res['win0_idx'][0], res['win1_idx'][0], facecolor='r', alpha=0.2)


    resp = gropt2.get_girf_ec_pc_response(g_out, 1000*res['H'], res['excite_idx_in'][0])
    # plt.figure(figsize = (12,4))

    ii = 1
    axs[ii].axhline(0, color='k', alpha=0.2)
    axs[ii].axhline(0+res['tol'], color='r', alpha=0.7, ls=":")
    axs[ii].axhline(0-res['tol'], color='r', alpha=0.7, ls=":")
    axs[ii].plot(resp)
    axs[ii].axvspan(res['win0_idx'][0], res['win1_idx'][0], facecolor='r', alpha=0.2)
    axs[ii].spines['right'].set_visible(False)
    axs[ii].spines['top'].set_visible(False)

    # print('resp norm:', 100*np.linalg.norm(resp[res['win0_idx'][0]:res['win1_idx'][0]]))

    if g_out_ref is not None:
        resp_ref = gropt2.get_girf_ec_pc_response(g_out_ref, 1000*res['H'], res['excite_idx_in'][0])
        axs[ii].plot(resp_ref)
        axs[ii].axvspan(res_ref['win0_idx'][0], res_ref['win1_idx'][0], facecolor='g', alpha=0.2)
        # print('resp_ref norm:', 100*np.linalg.norm(resp_ref[res_ref['win0_idx'][0]:res_ref['win1_idx'][0]]))
    
    axs[ii].set_ylim(-.02, .02)
    # plt.xlim(win0_idx[0]-100,win1_idx[0]+100)

    ii = 2
    axs[ii].axhline(0, color='k', alpha=0.2)
    axs[ii].axhline(0+res['tol'], color='r', alpha=0.7, ls=":")
    axs[ii].axhline(0-res['tol'], color='r', alpha=0.7, ls=":")
    axs[ii].plot(resp)
    axs[ii].axvspan(res['win0_idx'][0], res['win1_idx'][0], facecolor='r', alpha=0.2)
    axs[ii].set_ylim(-5*res['tol'], 5*res['tol'])
    axs[ii].set_xlim(res['win0_idx'][0]-100//res['shrink'],res['win1_idx'][0]+100/res['shrink'])
    axs[ii].spines['right'].set_visible(False)
    axs[ii].spines['top'].set_visible(False)

    if 'H1' in res:
        resp = gropt2.get_girf_ec_pc_response(g_out, 1000*res['H1'], res['excite_idx_in'][0])
        # plt.figure(figsize = (12,4))
        ii = 3
        axs[ii].axhline(0, color='k', alpha=0.2)
        axs[ii].axhline(res['H1_target']+res['H1_tol'], color='r', alpha=0.7, ls=":")
        axs[ii].axhline(res['H1_target']-res['H1_tol'], color='r', alpha=0.7, ls=":")
        axs[ii].plot(resp)
        axs[ii].axvspan(res['win0_idx'][0], res['win1_idx'][0], facecolor='r', alpha=0.2)

        if g_out_ref is not None:
            resp_ref = gropt2.get_girf_ec_pc_response(g_out_ref, 1000*res['H1'], res['excite_idx_in'][0])
            axs[ii].plot(resp_ref)
            axs[ii].axvspan(res_ref['win0_idx'][0], res_ref['win1_idx'][0], facecolor='g', alpha=0.2)
        
        axs[ii].set_ylim(-.1, .1)
        # plt.xlim(win0_idx[0]-100,win1_idx[0]+100)

        # plt.figure(figsize = (12,4))
        ii = 4
        axs[ii].axhline(0, color='k', alpha=0.2)
        axs[ii].axhline(res['H1_target']+res['H1_tol'], color='r', alpha=0.7, ls=":")
        axs[ii].axhline(res['H1_target']-res['H1_tol'], color='r', alpha=0.7, ls=":")
        axs[ii].plot(resp)
        axs[ii].axvspan(res['win0_idx'][0], res['win1_idx'][0], facecolor='r', alpha=0.2)
        axs[ii].set_ylim(res['H1_target']-3*res['H1_tol'], res['H1_target']+3*res['H1_tol'])
        axs[ii].set_xlim(res['win0_idx'][0]-100//res['shrink'],res['win1_idx'][0]+100/res['shrink'])

    if title is not None:
        plt.suptitle(title, fontsize=14)
        axs[ii].text(0.95, 0.01, 'Norm = {:.2f}'.format(100*np.linalg.norm(resp[res['win0_idx'][0]:res['win1_idx'][0]])),
        verticalalignment='bottom', horizontalalignment='right',
        transform=axs[ii].transAxes, fontsize=15)


    if return_im:
        fig.canvas.draw()

        # convert canvas to image
        img = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
        img  = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))

        # img is rgb, convert to opencv's default bgr
        img = cv2.cvtColor(img,cv2.COLOR_RGB2BGR)

        plt.close(fig)

        return img, 100*np.linalg.norm(resp[res['win0_idx'][0]:res['win1_idx'][0]])



def run_all(fname, H0_orig, T_bp = 1.8e-3, T_spoil = 1.0e-3, 
            offset_te = None, offset = 80, dW = 10, tol = .0001, 
            shrink = 1, do_ecpc = True, l2_weight = None,
            cg_niter=80, N_feval=20000, verbose=0, d_obj_thresh = 1e-5,
            H1_orig = None, H1_tol = .001, H1_target = 0, nonconvex=False, print_optimize=False,
            smooth_weight = 0, kaiser_sigma = 300, girfopt_mode = 1, add_tik = False, minTR = False, 
            solver_type = 'minresqlp', cg_rtol = 1e-8, cg_shift = 1e-6,
            **kwargs):
    
    if minTR:
        gparams, params = prepare_bgpc_seq_minTR(fname, do_init_plot= False, shrink = shrink, T_bp=T_bp, T_spoil=T_spoil, **kwargs)
    else:
        gparams, params = prepare_bgpc_seq(fname, do_init_plot= False, do_g0_plot=False, shrink = shrink, T_bp=T_bp, T_spoil=T_spoil, **kwargs)

    H = H0_orig.copy()

    if shrink > 1:
        dH_shrink = int(round(H.size*(1-1/shrink)/2))
        if (H[dH_shrink:-dH_shrink].size%2 != 0):
            H = H[dH_shrink:-dH_shrink-1]
        else:
            H = H[dH_shrink:-dH_shrink]

    N = params['N']

    offset = offset // shrink
    dW = dW // shrink

    if offset_te is not None:
        win_start = params['new_allidx_te'][int(offset_te)] - dW//2
        win_stop = params['new_allidx_te'][int(offset_te)] + dW//2
    else:
        win_start = params['new_idx_te'] - dW//2
        win_stop = params['new_idx_te'] + dW//2
    # print('[win_start, win_stop] =', [win_start, win_stop])

    excite_idx_in = params['allidx_excite'][:2].copy()
    excite_idx_in[1] = excite_idx_in[0] + N//2   # this is a quick fix for the minTR version
    win0_idx = np.array([win_start, win_start+N//2])
    win1_idx = np.array([win_stop, win_stop+N//2])

    # print('excite_idx_in:', excite_idx_in)
    

    target = 0
    weight = 1
    if do_ecpc:
        if H1_orig is not None:
            H1 = H1_orig.copy()

            if shrink > 1:
                dH_shrink = int(round(H1.size*(1-1/shrink)/2))
                if (H1[dH_shrink:-dH_shrink].size%2 != 0):
                    H1 = H1[dH_shrink:-dH_shrink-1]
                else:
                    H1 = H1[dH_shrink:-dH_shrink]

            gparams.add_girf_ec_pc(1000*H1, excite_idx_in, win0_idx, win1_idx, H1_tol, weight)
            gparams.set_constraint_vals('GirfEC_PC', H1_target, H1_tol, 1e-2)

        gparams.add_girf_ec_pc(1000*H, excite_idx_in, win0_idx, win1_idx, tol, weight, mode=girfopt_mode)
    
    else:
        H1 = H1_orig
        
    gparams.solver_type = solver_type
    gparams.cg_rtol = cg_rtol
    gparams.cg_shift = cg_shift
    gparams.cg_niter = cg_niter
    
    if smooth_weight > 0:
        mod = 1-windows.kaiser(N, kaiser_sigma)
        mod += 0.1*np.abs(np.linspace(-1,1,N))
        mod = smooth_weight*np.fft.fftshift(mod)
        gparams.add_fft_sqdiff_mod(mod)
    
    if l2_weight is not None:
        gparams.add_obj_duty(l2_weight)

    if add_tik:
        gparams.add_sqdiff_ref(np.zeros(N))

    gparams.set_verbose(verbose)
    gparams.d_obj_thresh = d_obj_thresh
    gparams.N_feval = N_feval
    
    if nonconvex:
        gparams.nonconvex_settings()
    
   
    if print_optimize:
        from coutcatcher import capture
        with capture() as c:
            gparams.optimize()
        print("Printed to C-stdout:\n", c.cout)
        print("Printed to C-stderr:\n", c.cerr)
    else:
        gparams.optimize()
        

    g_out = gparams.get_out()
    
    bp_out0 = g_out[params['bp_idx0']:params['bp_idx1']]
    bp_out1 = g_out[params['bp_idx0']+N//2:params['bp_idx1']+N//2]
    
    sp_out0 = g_out[params['sp_idx0']:params['sp_idx1']]
    sp_out1 = g_out[params['sp_idx0']+N//2:params['sp_idx1']+N//2]
    
    res = {'g_out': g_out,
           'final_good': gparams.final_good,
           'excite_idx_in': excite_idx_in,
           'win0_idx': win0_idx,
           'win1_idx': win1_idx,
           'H': H,
           "H1_tol": H1_tol,
           "H1_target": H1_target,
           'tol':tol,
           'shrink': shrink,
           'bp_out0': bp_out0,
           'bp_out1': bp_out1,
           'sp_out0': sp_out0,
           'sp_out1': sp_out1,
           'sp_idx0': params['sp_idx0'],
           'sp_idx1': params['sp_idx1'],
           'N': N,
           'params_in':params,
              }
    
    if H1_orig is not None:
        res['H1'] = H1
    
    
    return g_out, res, gparams

def prepare_bgpc_seq(fname, dt = 10e-6, venc = 160, do_init_plot = False, do_g0_plot = False,
                     gmax = .029, smax = 79, do_bipolar = True, do_spoiler = True,
                     T_bp = 1.8e-3, T_spoil = 1.0e-3, shrink = 1):
    
    gamma = 42576000

    mat = loadmat(fname)

    gg = mat['grad_save']/gamma
    # print('gg.shape = {:}'.format(gg.shape))

    idx_ss_start = int(mat['pidx_ss_start'])
    idx_excite = int(mat['pidx_excite'])
    idx_ss_end = int(mat['pidx_ss_end'])
    idx_ssre_end = int(mat['pidx_ssre_end'])
    idx_ro_start = int(mat['pidx_ro_start'])
    idx_te = int(mat['pidx_te'])
    idx_ro_end = int(mat['pidx_ro_end'])
    idx_spoil_start = int(mat['pidx_spoil_start'])
    idx_spoil_end = int(mat['pidx_spoil_end'])

    all_idx_times = [idx_ss_start, idx_excite, idx_ss_end,
                    idx_ssre_end, idx_ro_start, idx_te, 
                    idx_ro_end, idx_spoil_start, idx_spoil_end]

    # print('all_idx_times = {:}'.format(all_idx_times))

    if do_init_plot:
        plt.figure(figsize = (12,3))
        plt.axhline(0, lw=1, alpha=0.3, color='k')
        for tt in all_idx_times:
            plt.axvline(tt, lw=1, alpha=0.5, color='r', ls=':')
        plt.plot(gg.T)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)

    allidx_excite = mat['p_allidx_excite'].astype(int).squeeze()
    allidx_start = mat['p_allidx_start'].astype(int).squeeze()
    allidx_te = mat['p_allidx_te'].astype(int).squeeze()

    spoil_M0 = mat['spoil_m0']

    venc_u = venc * 1e-3 * 1e-2 # m/ms

    gamma = 267.5153151e6  # rad/s/T
    gamma_u = gamma * 1e-6  # rad/ms/mT

    M1 = np.pi / gamma_u / venc_u  # mT * ms^2 / m

    # print('M1 = {:.3f}  M1/2 = {:.3f}'.format(M1, M1/2))

    

    N = 2*(allidx_start[1] - allidx_start[0])
    # print('N = {:}'.format(N))

    g0 = gg[2].copy()

    N_dTE = idx_te-idx_ro_start  # This is time(N) from readout start to TE
    # print('N_dTE = {:}'.format(N_dTE))

    if do_bipolar:
        g0[idx_ssre_end:idx_ro_start+10] = 0  # Zero out the previous bipolar
        g0[idx_ssre_end+N//2:idx_ro_start+10+N//2] = 0

        # This is the bipolar region
        bp_idx0 = idx_ssre_end
        bp_idx1 = int(idx_ssre_end+T_bp/dt)
        g0[bp_idx0+1:bp_idx1-1] = np.nan
        g0[bp_idx0+1+N//2:bp_idx1-1+N//2] = np.nan
        # print('bp_idx0, bp_idx1 =', [bp_idx0, bp_idx1])

        # This is the earliest possible TE given the readout
        new_idx_te = bp_idx1 + N_dTE
    else:
        new_idx_te = idx_te
        bp_idx0 = idx_ssre_end
        bp_idx1 = idx_ro_start


    # Set region for spoiler optimization
    T_spoil0 = dt*(idx_spoil_end - idx_spoil_start)
    if do_spoiler:
        if T_spoil > T_spoil0:
            didx_sp = (T_spoil-T_spoil0)/dt
            didx_sp1 = int(min(int(np.floor(didx_sp/2)), 80))
            didx_sp0 = int(didx_sp - didx_sp1)
            sp_idx0 = idx_spoil_start - didx_sp0
            sp_idx1 = idx_spoil_end + didx_sp1
        else:
            sp_idx0 = idx_spoil_start
            sp_idx1 = idx_spoil_end
        g0[sp_idx0:sp_idx1] = np.nan
        g0[sp_idx0+N//2:sp_idx1+N//2] = np.nan
    else:
        sp_idx0 = idx_spoil_start
        sp_idx1 = idx_spoil_end


    if shrink > 1:
        g0 = g0[shrink//2::shrink]

        dt = dt*shrink
        N = int(round(N/shrink))
        idx_ss_start = int(round(idx_ss_start/shrink))
        idx_excite = int(round(idx_excite/shrink))
        idx_ss_end = int(round(idx_ss_end/shrink))
        idx_ssre_end = int(round(idx_ssre_end/shrink))
        idx_ro_start = int(round(idx_ro_start/shrink))
        idx_te = int(round(idx_te/shrink))
        idx_ro_end = int(round(idx_ro_end/shrink))
        idx_spoil_start = int(round(idx_spoil_start/shrink))
        idx_spoil_end = int(round(idx_spoil_end/shrink))

        new_idx_te = int(round(new_idx_te/shrink))
        bp_idx0 = int(round(bp_idx0/shrink))
        bp_idx1 = int(round(bp_idx1/shrink))
        sp_idx0 = int(round(sp_idx0/shrink))
        sp_idx1 = int(round(sp_idx1/shrink))

        allidx_start = np.round(allidx_start/shrink).astype(int)
        allidx_excite = np.round(allidx_excite/shrink).astype(int)
        allidx_te = np.round(allidx_te/shrink).astype(int)

        # print('shrink g0.shape', g0.shape, N, dt)


    # g0 = np.hstack([g0, g0])
    N = g0.size
    # print('Final g0.shape =', g0.shape, '  N =', N)

    if do_g0_plot:
        plt.figure(figsize = (12,3))
        plt.axhline(0, lw=1, alpha=0.3, color='k')
        plt.plot(g0)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)

    gparams = gropt2.GroptParams()
    gparams.init_N(N, dt)
    gparams.init_set_vals(g0)

    gparams.add_gmax(gmax)
    gparams.add_smax(smax)

    if do_bipolar:
        gparams.add_moment(0, 0, 0, idx_ssre_end, new_idx_te, 0.0, 1e-6)
        gparams.add_moment(0, 1, 0, idx_ssre_end, new_idx_te, M1/2, 1e-6)
        gparams.add_moment(0, 0, 0, idx_ssre_end+N//2, new_idx_te+N//2, 0.0, 1e-6)
        gparams.add_moment(0, 1, 0, idx_ssre_end+N//2, new_idx_te+N//2, -M1/2, 1e-6)
    if do_spoiler:
        gparams.add_moment(0, 0, 0, sp_idx0, sp_idx1, spoil_M0, 1e-6)
        gparams.add_moment(0, 0, 0, sp_idx0+N//2, sp_idx1+N//2, spoil_M0, 1e-6)

    gparams.cg_niter = 40
    # gparams.cg_resid_tol = 1e-2
    gparams.N_feval = 20000

    new_allidx_te = allidx_te + new_idx_te - idx_te

    params = {"idx_ss_start": idx_ss_start, 
               "idx_excite": idx_excite, 
               "idx_ss_end": idx_ss_end,
               "idx_ssre_end": idx_ssre_end, 
               "idx_ro_start": idx_ro_start, 
               "idx_te": idx_te, 
               "idx_ro_end": idx_ro_end, 
               "idx_spoil_start": idx_spoil_start, 
               "idx_spoil_end": idx_spoil_end,
               "new_idx_te": new_idx_te,
               "bp_idx0": bp_idx0,
               "bp_idx1": bp_idx1,
               "sp_idx0": sp_idx0,
               "sp_idx1": sp_idx1,
               "T_bp": T_bp,
               "T_spoil": T_spoil,
               "gmax": gmax,
               "smax": smax,
               "venc": venc,
               "spoil_M0": spoil_M0,
               "M1": M1,
               "allidx_start": allidx_start,
               "allidx_excite": allidx_excite,
               "allidx_te": allidx_te,
               "new_allidx_te": new_allidx_te,
               "N": N,
               "g0": g0,
               }
    
    return gparams, params


def prepare_bgpc_seq_minTR(fname, dt = 10e-6, venc = 160, do_init_plot = False, do_g0_plot = False,
                     gmax = .029, smax = 79, do_bipolar = True, do_spoiler = True,
                     T_bp = 1.8e-3, T_spoil = 1.0e-3, shrink = 1):
    
    gamma = 42576000

    mat = loadmat(fname)

    gg = mat['grad_save']/gamma
    # print('gg.shape = {:}'.format(gg.shape))

    idx_ss_start = int(mat['pidx_ss_start'])
    idx_excite = int(mat['pidx_excite'])
    idx_ss_end = int(mat['pidx_ss_end'])
    idx_ssre_end = int(mat['pidx_ssre_end'])
    idx_ro_start = int(mat['pidx_ro_start'])
    idx_te = int(mat['pidx_te'])
    idx_ro_end = int(mat['pidx_ro_end'])
    idx_spoil_start = int(mat['pidx_spoil_start'])
    idx_spoil_end = int(mat['pidx_spoil_end'])

    all_idx_times = [idx_ss_start, idx_excite, idx_ss_end,
                    idx_ssre_end, idx_ro_start, idx_te, 
                    idx_ro_end, idx_spoil_start, idx_spoil_end]

    # print('all_idx_times = {:}'.format(all_idx_times))

    if do_init_plot:
        plt.figure(figsize = (12,3))
        plt.axhline(0, lw=1, alpha=0.3, color='k')
        for tt in all_idx_times:
            plt.axvline(tt, lw=1, alpha=0.5, color='r', ls=':')
        plt.plot(gg.T)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)

    allidx_excite = mat['p_allidx_excite'].astype(int).squeeze()
    allidx_start = mat['p_allidx_start'].astype(int).squeeze()
    allidx_te = mat['p_allidx_te'].astype(int).squeeze()

    

    spoil_M0 = mat['spoil_m0']

    venc_u = venc * 1e-3 * 1e-2 # m/ms

    gamma = 267.5153151e6  # rad/s/T
    gamma_u = gamma * 1e-6  # rad/ms/mT

    M1 = np.pi / gamma_u / venc_u  # mT * ms^2 / m

    # print('M1 = {:.3f}  M1/2 = {:.3f}'.format(M1, M1/2))

    

    N = 2*(allidx_start[1] - allidx_start[0])
    
    # print('N = {:}'.format(N))

    g0 = gg[2].copy()

    # print('INIT::  allidx_excite:', allidx_excite, '  allidx_te:', allidx_te, '  N:', N, '  g0.size:', g0.size)

    N_dTE = idx_te-idx_ro_start  # This is time(N) from readout start to TE
    # print('N_dTE = {:}'.format(N_dTE))
    
    T_bp0 = dt*(idx_ro_start - idx_ssre_end)
    T_spoil0 = dt*(idx_spoil_end - idx_spoil_start)

    # print('T_bp0 = {:}  T_spoil0 = {:}'.format(1000*T_bp0, 1000*T_spoil0))

    if do_bipolar:
        g0[idx_ssre_end:idx_ro_start+10] = 0  # Zero out the previous bipolar
        g0[idx_ssre_end+N//2:idx_ro_start+10+N//2] = 0
    if do_spoiler:
        g0[idx_spoil_start-1:N//2] = 0  # Zero out the previous bipolar
        g0[idx_spoil_start-1+N//2:] = 0

    didx_bp = 0
    if do_bipolar:
        # This is the bipolar region
        bp_idx0 = idx_ssre_end
        bp_idx1 = idx_ro_start
        g0[bp_idx0+1:bp_idx1-1] = np.nan
        g0[bp_idx0+1+N//2:bp_idx1-1+N//2] = np.nan
        # print('bp_idx0, bp_idx1 =', [bp_idx0, bp_idx1])

        if T_bp > T_bp0:
            didx_bp = int(np.round((T_bp-T_bp0)/dt))
            # print('didx_bp', didx_bp)
            g0 = np.insert(g0, bp_idx1-1, np.nan*np.ones(didx_bp))
            g0 = np.insert(g0, bp_idx1+N//2-1+didx_bp, np.nan*np.ones(didx_bp))
            N = g0.size

        bp_idx1 += didx_bp
        # This is the earliest possible TE given the readout
        new_idx_te = bp_idx1 + N_dTE
    else:
        new_idx_te = idx_te
        bp_idx0 = idx_ssre_end
        bp_idx1 = idx_ro_start

    idx_spoil_start += didx_bp
    idx_spoil_end += didx_bp

    # Set region for spoiler optimization
    didx_spoil = 0
    if do_spoiler:
        sp_idx0 = idx_spoil_start
        sp_idx1 = idx_spoil_end
        g0[sp_idx0:sp_idx1] = np.nan
        g0[sp_idx0+N//2:sp_idx1+N//2] = np.nan

        if T_spoil > T_spoil0:
            didx_spoil = int(np.round((T_spoil-T_spoil0)/dt))
            # print('didx_spoil', didx_spoil)
            g0 = np.insert(g0, sp_idx1-1, np.nan*np.ones(didx_spoil))
            g0 = np.insert(g0, sp_idx1+N//2-1+didx_spoil, np.nan*np.ones(didx_spoil))
            N = g0.size

        sp_idx1 += didx_spoil

    else:
        sp_idx0 = idx_spoil_start
        sp_idx1 = idx_spoil_end

    if shrink > 1:
        g0 = g0[shrink//2::shrink]

        dt = dt*shrink
        N = int(round(N/shrink))
        idx_ss_start = int(round(idx_ss_start/shrink))
        idx_excite = int(round(idx_excite/shrink))
        idx_ss_end = int(round(idx_ss_end/shrink))
        idx_ssre_end = int(round(idx_ssre_end/shrink))
        idx_ro_start = int(round(idx_ro_start/shrink))
        idx_te = int(round(idx_te/shrink))
        idx_ro_end = int(round(idx_ro_end/shrink))
        idx_spoil_start = int(round(idx_spoil_start/shrink))
        idx_spoil_end = int(round(idx_spoil_end/shrink))

        new_idx_te = int(round(new_idx_te/shrink))
        bp_idx0 = int(round(bp_idx0/shrink))
        bp_idx1 = int(round(bp_idx1/shrink))
        sp_idx0 = int(round(sp_idx0/shrink))
        sp_idx1 = int(round(sp_idx1/shrink))

        allidx_start = np.round(allidx_start/shrink).astype(int)
        allidx_excite = np.round(allidx_excite/shrink).astype(int)
        allidx_te = np.round(allidx_te/shrink).astype(int)

        # print('shrink g0.shape', g0.shape, N, dt)


    # g0 = np.hstack([g0, g0])
    N = g0.size
    # print('Final g0.shape =', g0.shape, '  N =', N)

    if do_g0_plot:
        plt.figure(figsize = (12,3))
        plt.axhline(0, lw=1, alpha=0.3, color='k')
        plt.plot(g0)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)

        plt.axvline(idx_ssre_end, ls=':', color='r')
        plt.axvline(new_idx_te, ls=':', color='r')

        plt.axvline(idx_ssre_end+N//2, ls=':', color='r')
        plt.axvline(new_idx_te+N//2, ls=':', color='r')

        plt.axvline(sp_idx0, ls=':', color='g')
        plt.axvline(sp_idx1, ls=':', color='g')

        plt.axvline(sp_idx0+N//2, ls=':', color='g')
        plt.axvline(sp_idx1+N//2, ls=':', color='g')

    gparams = gropt2.PyGroptParams()
    gparams.init_N(N, dt)
    gparams.init_set_vals(g0)

    gparams.add_gmax(gmax)
    gparams.add_smax(smax)

    if do_bipolar:
        gparams.add_moment(0, 0, 0, idx_ssre_end, new_idx_te, 0.0, 1e-6)
        gparams.add_moment(0, 1, 0, idx_ssre_end, new_idx_te, M1/2, 1e-6)
        gparams.add_moment(0, 0, 0, idx_ssre_end+N//2, new_idx_te+N//2, 0.0, 1e-6)
        gparams.add_moment(0, 1, 0, idx_ssre_end+N//2, new_idx_te+N//2, -M1/2, 1e-6)
    if do_spoiler:
        gparams.add_moment(0, 0, 0, sp_idx0, sp_idx1, spoil_M0, 1e-6)
        gparams.add_moment(0, 0, 0, sp_idx0+N//2, sp_idx1+N//2, spoil_M0, 1e-6)

    gparams.cg_niter = 40
    # gparams.cg_resid_tol = 1e-2
    gparams.N_feval = 20000

    new_allidx_te = allidx_te + new_idx_te - idx_te

    # print('END::  allidx_excite:', allidx_excite, '  allidx_te:', allidx_te, '  N:', N, '  g0.size:', g0.size)

    params = {"idx_ss_start": idx_ss_start, 
               "idx_excite": idx_excite, 
               "idx_ss_end": idx_ss_end,
               "idx_ssre_end": idx_ssre_end, 
               "idx_ro_start": idx_ro_start, 
               "idx_te": idx_te, 
               "idx_ro_end": idx_ro_end, 
               "idx_spoil_start": idx_spoil_start, 
               "idx_spoil_end": idx_spoil_end,
               "new_idx_te": new_idx_te,
               "bp_idx0": bp_idx0,
               "bp_idx1": bp_idx1,
               "sp_idx0": sp_idx0,
               "sp_idx1": sp_idx1,
               "T_bp": T_bp,
               "T_spoil": T_spoil,
               "T_bp0": T_bp0,
               "T_spoil0": T_spoil0,
               "gmax": gmax,
               "smax": smax,
               "venc": venc,
               "spoil_M0": spoil_M0,
               "M1": M1,
               "allidx_start": allidx_start,
               "allidx_excite": allidx_excite,
               "allidx_te": allidx_te,
               "new_allidx_te": new_allidx_te,
               "N": N,
               "g0": g0,
               }
    
    return gparams, params