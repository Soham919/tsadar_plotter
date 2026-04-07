### PLOT LEARNED PARAMETERS FOR ELECTRONS ###
import numpy as np
import matplotlib.pyplot as plt
from TS_aux import shot_data

def plot_117837_EPW(data):
    lp = data['learned_parameters.csv']
    r = shot_data.loc[117837,"pointing"] - lp["Radius (\\mum)"]/10**3
    bv = (r*10**2)/shot_data.loc[117837,"timing"]
    ## Ti, Va, fract plot for all 3 species 
    plt.rcParams.update({
    #    "font.family": "serif",             # pick font family
    #    "font.serif": ["Arial"],  # pick specific font
        "font.size": 17,                     # set size
        "axes.labelsize": 17,    # font size for axis labels
        "axes.titlesize": 17,     # font size for titles
        # Tick direction
        "xtick.direction": "in",
        "ytick.direction": "in",
        # Tick width
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        # Tick length
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
    })

    # Error parameters
    window = 5  # rolling window
    Te = lp["Te_electron"].rolling(window=window).mean()
    ne = lp["ne_electron"].rolling(window=window).mean()
    std_Te = lp["Te_electron"].rolling(window=window).std()
    std_ne = lp["ne_electron"].rolling(window=window).std()

    fig1, ax1  = plt.subplots(figsize=(4,4))
    c1 = "tab:red"
    c2 = "tab:blue"
    l1, = ax1.plot(r,Te,color=c1,linewidth=2,linestyle='-',label=r'$T_{e}$')
    ax1.fill_between(r,Te-std_Te,Te+std_Te,
                    color='LightCoral',        # fill color
                    alpha=0.2,           # transparency
                    edgecolor='IndianRed',    # outline color
                    linewidth=1.5,       # outline width
                    linestyle='-',      # outline style (optional)
                    label='±1σ band')

    #ax1[2].plot(r, lp["Va_ion-2"],color='mediumblue',linewidth=2,linestyle='-')
    #ax1[2].plot(r, bv,color='k',linewidth=2,linestyle='--')
    ax1.set_title('Electron density and temperature')
    ax1.set_xlabel(r"$x (mm)$")
    ax1.set_xlim([4.8, 6.25])
    ax1.set_ylim([0.15, 1])
    ax1.invert_xaxis()
    ax1.set_ylabel(r"$T_{e} (keV)$", color = c1)
    ax1.tick_params(axis='y', colors=c1)
    ax1.spines['left'].set_color(c1)
    ax1.spines['left'].set_linewidth(1.5)
    ax1.spines['right'].set_visible(False)
    ax1.grid(True)
    #ax1[1,1].legend()

    ax2 = ax1.twinx()
    l2, = ax2.plot(r,ne,color=c2,linewidth=2,linestyle='-',label=r'$n_{e}$')
    ax2.fill_between(r,ne-std_ne,ne+std_ne,color='DodgerBlue',        # fill color
                    alpha=0.2,           # transparency
                    edgecolor='RoyalBlue',    # outline color
                    linewidth=1.5,       # outline width
                    linestyle='-',      # outline style (optional)
                    label='±1σ band')
    ax2.set_ylabel(r"$n_{e} (\times 10^{20} cm^{-3})$", color = c2)
    ax2.tick_params(axis='y', colors=c2)
    ax2.set_xlim([4.8, 6.25])
    ax2.set_ylim([0.1, 0.3])
    ax2.spines['right'].set_color(c2)
    ax2.spines['right'].set_linewidth(1.5)
    ax2.spines['left'].set_visible(False)
    #ax2.legend()

    lines = [l1, l2]
    labels = [line.get_label() for line in lines]
    ax1.legend(lines, labels, loc='upper right', frameon=False)

def plot_117837_EPW_loss(data):
    lp_1 = data['117837-5']['learned_parameters.csv']   # learned parameters 1
    lp_2 = data['117837-8']['learned_parameters.csv']   # learned parameters 2
    r1 = shot_data.loc[117837,"pointing"] - lp_1["Radius (\\mum)"]/10**3   # radius 1
    r2 = shot_data.loc[117837,"pointing"] - lp_2["Radius (\\mum)"]/10**3   # radius 1
    loss_1 = data['117837-5']['losses.csv']   # loss 1
    loss_2 = data['117837-8']['losses.csv']   # loss 2
    cl1 = np.array([(np.array(loss_1["initial_losses"][0:x+1])).sum() for x in range(len(loss_1["initial_losses"]))])  # cumulative loss 1
    cl2 = np.array([(np.array(loss_2["initial_losses"][0:x+1])).sum() for x in range(len(loss_2["initial_losses"]))])  # cumulative loss 2
    #arr = np.array([1,2,3,4,5,6,7,8,9,10])
    #cum_arr = np.array([arr[0:x+1].sum() for x in range(len(arr))])
    #print(cum_arr)
    fig2, ax2 = plt.subplots(1,2,figsize=(12,6))
    ax2[0].scatter(lp_1["Radius (\\mum)"],loss_1["initial_losses"],facecolors='none', edgecolors='blue',label='117837-5')
    ax2[0].scatter(lp_2["Radius (\\mum)"],loss_2["initial_losses"],facecolors='none', edgecolors='red',label='117837-8')
    ax2[0].set_title("loss v.s r")
    ax2[0].set_xlabel(r'$r (\mu m)$')
    ax2[0].set_ylabel('loss')
    ax2[0].legend()

    ax2[1].plot(lp_1["Radius (\\mum)"],cl1,color='blue',linestyle='--',label='117837-5')
    ax2[1].plot(lp_2["Radius (\\mum)"],cl2,color='red',linestyle='--',label='117837-8')
    ax2[1].set_title("cumulative loss v.s r")
    ax2[1].set_xlabel(r'$r (\mu m)$')
    ax2[1].set_ylabel('cum.loss')
    ax2[1].legend()