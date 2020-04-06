import numpy as np
from reader import read_singlepulse, get_triggers
from tools3 import dm_range

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

# filename = 'trigger/FRB_bk.trigger'
def load_trigger_file(filename = '../data/trigger/FRB190925_CB07.trigger',
                        verbose = False,
                        read_data = True,
                        read_beam=True,
                        replace = False):
    # Set search parameters
    sig_thresh=5.0
    dm_min=0
    dm_max=np.inf
    t_window=0.5
    max_rows=None
    t_max=np.inf
    sig_max=np.inf
    dt=2*40.96
    delta_nu_MHz=300./1536
    nu_GHz=1.4
    fnout=False
    tab=None

    dm_width_filter=False

    if read_data:
        dm, sig, tt, downsample, beam = read_singlepulse(filename, beam='all')[:5]

        ntrig_orig = len(dm)

        bad_sig_ind = np.where((sig < sig_thresh) | (sig > sig_max))[0]
        sig = np.delete(sig, bad_sig_ind)
        tt = np.delete(tt, bad_sig_ind)
        dm = np.delete(dm, bad_sig_ind)
        downsample = np.delete(downsample, bad_sig_ind)
        sig_cut, dm_cut, tt_cut, ds_cut = [], [], [], []
        if read_beam:
            beam = np.delete(beam, bad_sig_ind)
            beam_cut = []

        tduration = tt.max() - tt.min()
        ntime = int(tduration / t_window)

        # Make dm windows between 90% of the lowest trigger and
        # 10% of the largest trigger
        if dm_min == 0:
            dm_min = 0.9*dm.min()
        if dm_max > 1.1*dm.max():
            dm_max = 1.1*dm.max()

        # Can either do the DM selection here, or after the loop
        dm_list = dm_range(dm_max, dm_min=dm_min)

        print("\nGrouping in window of %.2f sec" % np.round(t_window,2))
        print("DMs:", dm_list)

        tt_start = tt.min() - .5*t_window
        ind_full = []

        beams = []
        sigs = []
        times = []
        dms_ = []
        j_idx = 0
        j = 0

        # might wanna make this a search in (dm,t,width) cubes
        for dms in dm_list:
            for ii in range(ntime + 2):
                try:
                    # step through windows of t_window seconds, starting from tt.min()
                    # and find max S/N trigger in each DM/time box
                    t0, tm = t_window*ii + tt_start, t_window*(ii+1) + tt_start
                    ind = np.where((dm<dms[1]) & (dm>dms[0]) & (tt<tm) & (tt>t0))[0]

                    # Check for unique synthesized beam triggers (remove multiple instances)
                    uniques = np.unique(beam[ind])#, return_counts=True)
                    new_beam = [int(u) for u in uniques]
                    new_sig = [sig[ind][np.where(beam[ind] == u)][np.argmax(sig[ind][np.where(beam[ind] == u)])] for u in uniques]

                    # filter_sbs(new_)

                    # if 11 in scipy.signal.find_peaks(np.absolute(np.fft.fft(new_sig)))[0]:
    #                 if stderr > 0.1:
    # #                 if peaks_at_max == 1:
    #                     plt.plot(new_beam, new_sig)
    # #                     abline(slope, intercept)
    #                     print("idx", j_idx)
    #                     print("slope", slope, )
    #                     print("intercept", intercept, )
    #                     print("rvalue", rvalue, )
    #                     print("pvalue", pvalue, )
    #                     print("stderr", stderr)
    #                     print("stdev", np.std(top_sigs))
    #                     j += 1
    #                 j_idx +=1

                    if len(new_beam) > 0:
                        beams.append(np.asarray((new_beam)))
                        sigs.append(np.asarray((new_sig)))
                        times.append(np.asarray([tt[ind][np.where(beam[ind] == u)][np.argmax(sig[ind][np.where(beam[ind] == u)])] for u in uniques]))
                        dms_.append(np.asarray([dm[ind][np.where(beam[ind] == u)][np.argmax(sig[ind][np.where(beam[ind] == u)])] for u in uniques]))

                    # ind_maxsnr = ind[np.argmax(sig[ind])]
                    #
                    # sig_cut.append(sig[ind_maxsnr])
                    # dm_cut.append(dm[ind_maxsnr])
                    # tt_cut.append(tt[ind_maxsnr])
                    # ds_cut.append(downsample[ind_maxsnr])
                    # if read_beam:
                    #     beam_cut.append(beam[ind_maxsnr])
                    # ind_full.append(ind_maxsnr)
                except ValueError:
                    continue

    return beams, sigs, times, dms_

    print ("number of triggers", j)
