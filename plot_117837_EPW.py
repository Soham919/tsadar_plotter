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
    lp_1 = data['117837-17']['learned_parameters.csv']   # learned parameters 1
    lp_2 = data['117837-21']['learned_parameters.csv']   # learned parameters 2
    r1 = shot_data.loc[117837,"pointing"] - lp_1["Radius (\\mum)"]/10**3   # radius 1
    r2 = shot_data.loc[117837,"pointing"] - lp_2["Radius (\\mum)"]/10**3   # radius 1
    loss_1 = data['117837-17']['losses.csv']   # loss 1
    loss_2 = data['117837-21']['losses.csv']   # loss 2
    cl1 = np.array([(np.array(loss_1["initial_losses"][0:x+1])).sum() for x in range(len(loss_1["initial_losses"]))])  # cumulative loss 1
    cl2 = np.array([(np.array(loss_2["initial_losses"][0:x+1])).sum() for x in range(len(loss_2["initial_losses"]))])  # cumulative loss 2
    #arr = np.array([1,2,3,4,5,6,7,8,9,10])
    #cum_arr = np.array([arr[0:x+1].sum() for x in range(len(arr))])
    #print(cum_arr)
    fig2, ax2 = plt.subplots(1,2,figsize=(12,6))
    ax2[0].scatter(lp_1["Radius (\\mum)"],loss_1["initial_losses"],facecolors='none',s=12, edgecolors='blue',label='117837-17, Z=2')
    ax2[0].scatter(lp_2["Radius (\\mum)"],loss_2["initial_losses"],facecolors='none',s=12, edgecolors='red',label='117837-21, Z=1')
    ax2[0].set_title("loss v.s r")
    ax2[0].set_xlabel(r'$r (\mu m)$')
    ax2[0].set_ylabel('loss')
    ax2[0].legend()

    ax2[1].plot(lp_1["Radius (\\mum)"],cl1,color='blue',linestyle='--',label='117837-17, Z=2')
    ax2[1].plot(lp_2["Radius (\\mum)"],cl2,color='red',linestyle='--',label='117837-21, Z=1')
    ax2[1].set_title("cumulative loss v.s r")
    ax2[1].set_xlabel(r'$r (\mu m)$')
    ax2[1].set_ylabel('cum.loss')
    ax2[1].legend()

def plot_117837_IAW(data):
    ### PLOT LEARNED PARAMETERS for H2 AND HE SHOCKS ###
    im1 = data['117837']['117837-21']['learned_parameters.csv']
    #epw = data['117830']['92530_8_nersc.csv']

    # Convert ion fractions to number denisties




    r = shot_data.loc[117837,"pointing"] - im1["Radius (\\mum)"]/10**3
    #bv = (r*10**2)/shot_data.loc[92537,"timing"]
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

    fig1, ax1  = plt.subplots(3,1,figsize=(10,20))

    window = 2  # rolling window
    f1_1 = im1["fract_ion-1"].rolling(window=window).mean()
    f2_1 = im1["fract_ion-2"].rolling(window=window).mean()
    #f3_1 = im1["fract_ion-3"].rolling(window=window).mean()
    std_f1_1 = im1["fract_ion-1"].rolling(window=window).std()
    std_f2_1 = im1["fract_ion-2"].rolling(window=window).std()
    #std_f3_1 = im1["fract_ion-3"].rolling(window=window).std()

    t1_1 = im1["Ti_ion-1"].rolling(window=window).mean()
    t2_1 = im1["Ti_ion-2"].rolling(window=window).mean()
    #t3_1 = im1["Ti_ion-3"].rolling(window=window).mean()
    std_t1_1 = im1["Ti_ion-1"].rolling(window=window).std()
    std_t2_1 = im1["Ti_ion-2"].rolling(window=window).std()
    #std_t3_1 = im1["Ti_ion-3"].rolling(window=window).std()

    v1_1 = im1["Va_ion-1"].rolling(window=window).mean()
    v2_1 = im1["Va_ion-2"].rolling(window=window).mean()
    #v3_1 = im1["Va_ion-3"].rolling(window=window).mean()
    std_v1_1 = im1["Va_ion-1"].rolling(window=window).std()
    std_v2_1 = im1["Va_ion-2"].rolling(window=window).std()
    #std_v3_1 = im1["Va_ion-3"].rolling(window=window).std()


    ax1[0].plot(r, f1_1,color='orangered',linewidth=2, linestyle='-', label='H')
    ax1[0].fill_between(r,f1_1-std_f1_1,f1_1+std_f1_1,
                color='coral',        # fill color
                alpha=0.2,           # transparency
                edgecolor='tomato',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    #ax1[0].plot(r, lp["fract_ion-2"],color='mediumblue',linewidth=2, linestyle='-', label='Cold H')
    ax1[0].plot(r, f2_1,color='MediumSeaGreen',linewidth=2, linestyle='-', label='He')
    ax1[0].fill_between(r,f2_1-std_f2_1,f2_1+std_f2_1,
                color='MediumAquamarine',        # fill color
                alpha=0.2,           # transparency
                edgecolor='MediumSeaGreen',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )

    ax1[0].set_title('Fractions Z = 1')
    ax1[0].set_xlabel(r"$x (mm)$")
    ax1[0].set_xlim([4.75, 6.05])
    ax1[0].invert_xaxis()
    ax1[0].set_ylabel(r"$Fraction$")
    ax1[0].grid(True)
    ax1[0].legend()

    ax1[1].plot(r, t1_1,color='orangered',linewidth=2,linestyle='-',label='H')
    ax1[1].fill_between(r,t1_1-std_t1_1,t1_1+std_t1_1,
                color='coral',        # fill color
                alpha=0.2,           # transparency
                edgecolor='tomato',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    #ax1[1].plot(r, lp["Ti_ion-2"],color='mediumblue',linewidth=2,linestyle='-')
    ax1[1].plot(r, t2_1,color='MediumSeaGreen',linewidth=2,linestyle='-',label='He')
    ax1[1].fill_between(r,t2_1-std_t2_1,t2_1+std_t2_1,
                color='MediumAquamarine',        # fill color
                alpha=0.2,           # transparency
                edgecolor='MediumSeaGreen',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[1].set_title('Temperatures')
    ax1[1].set_xlabel(r"$x (mm)$")
    ax1[1].set_xlim([4.75, 6.05])
    ax1[1].set_ylim([0.0, 1.1])
    ax1[1].invert_xaxis()
    ax1[1].set_ylabel(r"$T_{i} (keV)$")
    ax1[1].grid(True)
    ax1[1].legend()

    ax1[2].plot(r, v1_1,color='orangered',linewidth=2,linestyle='-',label='H')
    ax1[2].fill_between(r,v1_1-std_v1_1,v1_1+std_v1_1,
                color='coral',        # fill color
                alpha=0.2,           # transparency
                edgecolor='tomato',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    #ax1[2].plot(r, lp["Va_ion-2"],color='mediumblue',linewidth=2,linestyle='-')
    ax1[2].plot(r, v2_1,color='MediumSeaGreen',linewidth=2,linestyle='-',label='He')
    ax1[2].fill_between(r,v2_1-std_v2_1,v2_1+std_v2_1,
                color='MediumAquamarine',        # fill color
                alpha=0.2,           # transparency
                edgecolor='MediumSeaGreen',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    #ax1[2].plot(r, bv,color='k',linewidth=2,linestyle='--')
    ax1[2].set_title('Velocities')
    ax1[2].set_xlabel(r"$x (mm)$")
    ax1[2].set_xlim([4.75, 6.05])
    ax1[2].invert_xaxis()
    ax1[2].set_ylabel(r"$V_{i} (10^{6} cm/s)$")
    ax1[2].grid(True)
    ax1[2].legend()