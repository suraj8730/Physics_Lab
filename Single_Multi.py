import os

import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib import colors
from plotly.subplots import make_subplots

import Single
import Single_h
import Multi
import Multi_height


########################################################################################################################
"""...............................................Single............................................................."""

def rough_vs_temp(N_row, N_col, Time, folder, Name):
    Time = np.arange(Time, dtype=float)
    Time = np.delete(Time, 0)
    Time = np.insert(Time, 0, 0.1)
    print(Time)
    Temp = np.arange(500, 801, 100)
    print(Temp)
    height = np.zeros([len(Temp), len(Time), N_row, N_col])
    avg = np.zeros([len(Temp), len(Time) + 1])
    std = np.zeros([len(Temp), len(Time) + 1])
    for i in range(len(Temp)):
        for j in range(len(Time)):
            avg[i][j], std[i][j], height[i][j] = Single.main_kmc_code(Time[j], N_row, N_col, Temp[i], 5)

    fig = go.Figure()
    for i in range(len(Temp)):
        fig.add_trace(go.Scatter(x=Time, y=std[i],
                                 mode='lines+markers',
                                 name=f'{Temp[i]} K'))
    fig.update_layout(
        title=f'For {N_row}*{N_col} Lattice at adsorption rate {5} ML/site-sec',
        xaxis_title='Time(sec.)',
        yaxis_title='Roughness(Monolayer)')
    fig.write_image(f"{folder}/{Name}/Roughness vs Time.svg")
    #fig.show()

    fig2 = go.Figure()
    for i in range(len(Temp)):
        fig2.add_trace(go.Scatter(x=Time, y=avg[i],
                                  mode='lines+markers',
                                  name=f'{Temp[i]} K'))
    fig2.update_layout(
        title=f'For {N_row}*{N_col} Lattice at adsorption rate {5} ML/site-sec',
        xaxis_title='Time(sec.)',
        yaxis_title='Average Height(Monolayer)')
    fig2.write_image(f"{folder}/{Name}/Avg_Height vs Time.svg")
    #fig2.show()

    """fig3 = go.Figure()
    for i in range(len(Time)):
        fig3.add_trace(go.Scatter(x=Temp, y=std.T[i],
                mode='lines+markers',
                name=f'{Time[i]} sec'))
    fig3.update_layout(
                   title = f'For {N_row}*{N_col} Lattice at adsorption rate {5} Monolayer/sec',
                  xaxis_title='Temp',
                  yaxis_title='Roughness')
    fig3.write_image(f"{folder}/{Name}/Roughness vs Temp.svg")
    fig3.show()
"""
    fig4 = make_subplots(
        rows=2, cols=2,
        specs=[[{'type': 'surface'}, {'type': 'surface'}],
               [{'type': 'surface'}, {'type': 'surface'}]],
        subplot_titles=(f"At {Temp[0]} K", f"At {Temp[1]} K", f"At {Temp[2]} K", f"At {Temp[3]} K"))
    fig4.add_trace(
        go.Surface(z=height[0][-1], colorscale='Viridis', showscale=False),
        row=1, col=1)

    fig4.add_trace(
        go.Surface(z=height[1][-1], colorscale='RdBu', showscale=False),
        row=1, col=2)

    fig4.add_trace(
        go.Surface(z=height[2][-1], colorscale='YlOrRd', showscale=False),
        row=2, col=1)

    fig4.add_trace(
        go.Surface(z=height[3][-1], colorscale='YlGnBu', showscale=False),
        row=2, col=2)

    fig4.update_layout(
        title_text=f'For {N_row}*{N_col} Lattice, at adsorption rate {5} ML/site-sec',
        height=1000,
        width=1500
    )
    fig4.write_image(f"{folder}/{Name}/Single_Surface_plot.svg")
    #fig4.show()


def ads_rate_variation(N_row, N_col, Max_level, folder, Name):
    Temp = np.arange(500, 800, 100)
    Ads_rate = np.arange(10)
    Ads_rate = np.delete(Ads_rate, 0)
    Time = np.full([len(Temp), len(Ads_rate)], 1.0)
    for i in range(len(Temp)):
        for j in range(len(Ads_rate)):
            avg, std, h_matrix, de_count, time_ = Single_h.main_kmc_code(Time[i][j], N_row, N_col, Max_level, Temp[i],
                                                                         Ads_rate[j])
            Time[i][j] = time_

    fig = go.Figure()
    for i in range(len(Temp)):
        fig.add_trace(go.Scatter(x=Ads_rate, y=Time[i],
                                 mode='lines+markers',
                                 name=f'{Temp[i]} K'))
    fig.update_layout(
        title=f'For {N_row}*{N_col} Lattice',
        xaxis_title='Adsorption Rate (ML/site-sec)',
        yaxis_title='Time(sec.)')
    fig.write_image(f"{folder}/{Name}/Time_vs_rates.svg")
    #fig.show()


def rough_vs_rates(N_row, N_col, Max_levl, folder, Name):
    Temp = np.arange(500, 900, 100)
    Ads_rate = np.arange(10)
    Ads_rate = np.delete(Ads_rate, 0)
    Rough = np.zeros([len(Temp), len(Ads_rate)])
    Timef = 1
    for i in range(len(Temp)):
        for j in range(len(Ads_rate)):
            avg, Rough[i][j], h_matrix, de_count, time_ = Single_h.main_kmc_code(Timef, N_row, N_col, Max_levl, Temp[i],
                                                                                 Ads_rate[j])

    fig2 = go.Figure()
    for i in range(len(Temp)):
        fig2.add_trace(go.Scatter(x=Ads_rate, y=Rough[i],
                                  mode='lines+markers',
                                  name=f'{Temp[i]} K'))
    fig2.update_layout(
        title=f'For {N_row}*{N_col} Lattice',
        xaxis_title='Adsorption Rate (ML/site-sec)',
        yaxis_title='Roughness(ML)')
    fig2.write_image(f"{folder}/{Name}/roughness_vs_rate.svg")
    #fig2.show()


def avg_roughness(nsitesrow, nsitescol,Time, folder, Name):
    Time = np.arange(Time, dtype=float)
    Time = np.delete(Time, 0)
    Time = np.insert(Time, 0, 0.1)
    rate_adsorption=5
    T=500
    Roughness = np.zeros([len(Time), 3])
    Avg_Rn = np.zeros([len(Time)])
    std_Rn = np.zeros([len(Time)])
    for i in range(len(Time)):
        for j in range(3):
            avg, Roughness[i][j], h = Single.main_kmc_code(Time[i], nsitesrow, nsitescol, T, rate_adsorption)
        Avg_Rn[i] = np.mean(Roughness[i])
        std_Rn[i] = np.std(Roughness[i])

    fig = go.Figure(data=go.Scatter(
        x=Time,
        y=Avg_Rn,
        error_y=dict(
            type='data',
            array=std_Rn)
    ))
    fig.update_layout(
        title=f'For {nsitesrow}*{nsitescol} Lattice at adsorption rate {rate_adsorption} ML/site-sec and {T} K temperature',
        xaxis_title='Time(sec.)',
        yaxis_title='Roughness(ML)')
    fig.write_image(f"{folder}/{Name}/{nsitescol}_Error_Roughness.svg")
    #fig.show()


def flux_vs_desorption(timef, h_max, nsitesrow, nsitescol, folder, Name):
    for j in range(700, 900, 50):
        T = j
        adsorb = []
        no_deso = []
        for i in range(1, 10, 2):
            adsorb.append(i)
            mean, roughness, height, de_count, time_ = Single_h.main_kmc_code(timef, nsitesrow, nsitescol, h_max, T, i)
            no_deso.append(de_count)
        plt.plot(adsorb, no_deso, label=f'{j} K')
    plt.xlabel('Adsorption rate(ML/site-sec)')
    plt.ylabel('No. of desorbed particle')
    plt.title('Effect of flux on number of desorbed particle at varying temperatures')
    plt.legend()
    plt.savefig(f"{folder}/{Name}/Flux_desorption.png", bbox_inches='tight')
    plt.close()


########################################################################################################################
# Multi

def rough_vs_temp_m(N_row, N_col, Time, folder, Name):
    A_comp = 0.5
    B_comp = 0.5
    Time = np.arange(Time, dtype=float)
    Time = np.delete(Time, 0)
    Time = np.insert(Time, 0, 0.1)
    print(Time)
    Temp = np.arange(500, 801, 100)
    print(Temp)
    height = np.zeros([len(Temp), len(Time), N_row, N_col])
    particle = np.zeros([len(Temp), len(Time), N_row, N_col])
    avg = np.zeros([len(Temp), len(Time) + 1])
    std = np.zeros([len(Temp), len(Time) + 1])
    for i in range(len(Temp)):
        for j in range(len(Time)):
            height[i][j], particle[i][j], p_matrix = Multi.main_kmc_code(Time[j], A_comp, B_comp, N_row, N_col, Temp[i],
                                                                         5)
            avg[i][j] = np.mean(height[i][j])
            std[i][j] = np.std(height[i][j])

    fig = go.Figure()
    for i in range(len(Temp)):
        fig.add_trace(go.Scatter(x=Time, y=std[i],
                                 mode='lines+markers',
                                 name=f'{Temp[i]} K'))
    fig.update_layout(
        title=f'For {N_row}*{N_col} Lattice at adsorption rate {5} ML/site-sec',
        xaxis_title='Time(sec.)',
        yaxis_title='Roughness(Monolayer)')
    fig.write_image(f"{folder}/{Name}/Roughness vs Time.svg")
    #fig.show()

    fig2 = go.Figure()
    for i in range(len(Temp)):
        fig2.add_trace(go.Scatter(x=Time, y=avg[i],
                                  mode='lines+markers',
                                  name=f'{Temp[i]} K'))
    fig2.update_layout(
        title=f'For {N_row}*{N_col} Lattice at adsorption rate {5} ML/site-sec',
        xaxis_title='Time(sec.)',
        yaxis_title='Average Height(Monolayer)')
    fig2.write_image(f"{folder}/{Name}/Avg_Height vs Time.svg")
    #fig2.show()

    """fig3 = go.Figure()
    for i in range(len(Time)):
        fig3.add_trace(go.Scatter(x=Temp, y=std.T[i],
                mode='lines+markers',
                name=f'{Time[i]} sec'))
    fig3.update_layout(
                   title = f'For {N_row}*{N_col} Lattice at adsorption rate {5} Monolayer/sec',
                  xaxis_title='Temp',
                  yaxis_title='Roughness')
    fig3.write_image(f"{Name}/Roughness vs Temp.svg")
    fig3.show()
    """
    fig4 = make_subplots(
        rows=2, cols=2,
        specs=[[{'type': 'surface'}, {'type': 'surface'}],
               [{'type': 'surface'}, {'type': 'surface'}]],
        subplot_titles=(f"At {Temp[0]} K", f"At {Temp[1]} K", f"At {Temp[2]} K", f"At {Temp[3]} K"))
    fig4.add_trace(
        go.Surface(z=height[0][-1], colorscale='Viridis', showscale=False),
        row=1, col=1)

    fig4.add_trace(
        go.Surface(z=height[1][-1], colorscale='RdBu', showscale=False),
        row=1, col=2)

    fig4.add_trace(
        go.Surface(z=height[2][-1], colorscale='YlOrRd', showscale=False),
        row=2, col=1)

    fig4.add_trace(
        go.Surface(z=height[3][-1], colorscale='YlGnBu', showscale=False),
        row=2, col=2)

    fig4.update_layout(
        title_text=f'For {N_row}*{N_col} Lattice, at adsorption rate {5} ML/site-sec',
        height=1080,
        width=1080
    )
    fig4.write_image(f"{folder}/{Name}/Multi_Surface_plot.svg")
    #fig4.show()


def ads_rate_variation_m(N_row, N_col, Max_level, folder, Name):
    Temp = np.arange(500, 800, 100)
    Ads_rate = np.arange(10)
    Ads_rate = np.delete(Ads_rate, 0)
    Time = np.full([len(Temp), len(Ads_rate)], 1.0)
    for i in range(len(Temp)):
        for j in range(len(Ads_rate)):
            mean, roughness, height, de_count, time_ = Multi_height.main_kmc_code(Time[i][j], N_row, N_col, Max_level,
                                                                                  Temp[i],
                                                                                  Ads_rate[j])
            Time[i][j] = time_

    fig = go.Figure()
    for i in range(len(Temp)):
        fig.add_trace(go.Scatter(x=Ads_rate, y=Time[i],
                                 mode='lines+markers',
                                 name=f'{Temp[i]} K'))
    fig.update_layout(
        title=f'For {N_row}*{N_col} Lattice',
        xaxis_title='Adsorption Rate (ML/site-sec)',
        yaxis_title='Time(sec.)')
    fig.write_image(f"{folder}/{Name}/Time_vs_rates.svg")
    #fig.show()


def avg_roughness_m(nsitesrow, nsitescol, Time, folder, Name):
    A_comp = 0.5
    B_comp = 0.5
    T = 500
    rate_adsorption = 5
    Time = np.arange(Time, dtype=float)
    Time = np.delete(Time, 0)
    Time = np.insert(Time, 0, 0.1)

    height = np.zeros([len(Time), 3])
    Avg_Rn = np.zeros([len(Time)])
    std_Rn = np.zeros([len(Time)])
    for i in range(len(Time)):
        for j in range(3):
            h, p, p_mat = Multi.main_kmc_code(Time[i], A_comp, B_comp, nsitesrow, nsitescol, T, rate_adsorption)
            height[i][j] = np.mean(h)
        Avg_Rn[i] = np.mean(height[i])
        std_Rn[i] = np.std(height[i])

    fig = go.Figure(data=go.Scatter(
        x=Time,
        y=Avg_Rn,
        error_y=dict(
            type='data',
            array=std_Rn)
    ))
    fig.update_layout(
        title=f'For {nsitesrow}*{nsitescol} Lattice at adsorption rate {rate_adsorption} ML/site-sec and {T} K temprature',
        xaxis_title='Time(sec.)',
        yaxis_title='Roughness(ML)')
    fig.write_image(f"{folder}/{Name}/Error_Roughness.svg")
    #fig.show()


def flux_vs_desorption_m(timef, h_max, nsitesrow, nsitescol, folder, Name):
    for j in range(700, 900, 50):
        T = j
        adsorb = []
        no_deso = []
        for i in range(1, 10, 2):
            adsorb.append(i)
            mean, roughness, height, de_count, time_ = Multi_height.main_kmc_code(timef, nsitesrow, nsitescol, h_max, T,
                                                                                  i)
            no_deso.append(de_count)
        plt.plot(adsorb, no_deso, label=f'{j} K')
    plt.xlabel('Adsorption rate(ML/site-sec)')
    plt.ylabel('No. of desorbed particle')
    plt.legend()
    plt.savefig(f"{folder}/{Name}/Flux_desorption.png", bbox_inches='tight')
    plt.close()


def multi_2d(nsitesrow, nsitescol, timef, folder, Name):
    A_comp = 0.5
    B_comp = 0.5
    rate_adsorption = 2
    y = 0
    fig, axs = plt.subplots(2, 3, figsize=(10, 10))
    order = []
    for i in range(300, 850, 100):
        T = i
        if y < 3:
            ax = axs[0, y]
        else:
            ax = axs[1, y - 3]
        height, particle, particl_mat = Multi.main_kmc_code(timef, A_comp, B_comp, nsitesrow, nsitescol, T,
                                                            rate_adsorption)
        # cmap = colors.ListedColormap(['white','Blue', 'red'])
        cmap = colors.ListedColormap(['red', 'Blue'])
        ax.pcolormesh(particle[::-1], cmap=cmap, edgecolors='k', linewidths=1)
        ax.set_title(f'{i} K')
        y += 1
    #plt.suptitle('Time')
    plt.savefig(f"{folder}/{Name}/Multi_2D.png", bbox_inches='tight')
    plt.close()


def multi_2d_flux(nsitesrow, nsitescol, timef, folder, Name):
    T = 700
    A_comp = 0.3
    B_comp = 0.7
    y = 0
    fig, axs = plt.subplots(1, 2, figsize=(10, 10))
    for i in range(1, 10, 8):
        rate_adsorption = i
        ax = axs[y]
        height, particle, particl_mat = Multi.main_kmc_code(timef, A_comp, B_comp, nsitesrow, nsitescol, T,
                                                            rate_adsorption)
        cmap = colors.ListedColormap(['red', 'Blue'])
        ax.pcolormesh(particle[::-1], cmap=cmap, edgecolors='k', linewidths=1)
        ax.set_title(f'{i} ML/site-sec')
        y += 1
    #plt.suptitle('Time')
    plt.savefig(f"{folder}/{Name}/2D_flux.png", bbox_inches='tight')
    plt.close()


def multi_2D_gascomposition(nsitesrow, nsitescol, timef, folder, Name):
    rate_adsorption = 7
    T = 700
    y = 0
    fig, axs = plt.subplots(1, 3, figsize=(10, 10))
    for i in range(10, 55, 20):
        A_comp = i / 100
        B_comp = 1 - A_comp
        ax = axs[y]
        height, particle, particl_mat = Multi.main_kmc_code(timef, A_comp, B_comp, nsitesrow, nsitescol, T,
                                                            rate_adsorption)
        cmap = colors.ListedColormap(['red', 'Blue'])
        ax.pcolormesh(particle[::-1], cmap=cmap, edgecolors='k', linewidths=1)
        ax.set_title(f'{i}% of A composition')
        y += 1
    #plt.suptitle('Time')
    plt.savefig(f"{folder}/{Name}/2D_A_composition.png", bbox_inches='tight')
    plt.close()


def grain_size_vs_temp(timef, nsiterow, nsitecol, folder, Name):
    rate_adsorption = 1
    Temp = np.arange(500, 950, 200)
    comp = np.arange(10, 55, 10)
    order = np.zeros((len(comp), len(Temp)))
    for i in range(len(comp)):
        for j in range(len(Temp)):
            A_comp = comp[i] / 100
            B_comp = 1 - A_comp
            height, particle, particl_mat = Multi.main_kmc_code(timef, A_comp, B_comp, nsiterow, nsitecol, Temp[j],
                                                                rate_adsorption)
            order[i][j] = (Multi.get_grain_size(nsitecol, nsiterow, 1, particle))
    fig = go.Figure()
    for i in range(len(comp)):
        fig.add_trace(go.Scatter(x=Temp, y=order[i],
                                 mode='lines+markers',
                                 name=f'{comp[i]} %'))
    fig.update_layout(
        title=f'For {nsiterow}*{nsitecol} Lattice',
        xaxis_title='Temperature (K)',
        yaxis_title='Grain size')
    fig.write_image(f"{folder}/{Name}/Grainsize_vs_temp.svg")
    #fig.show()


timef = 60
h_max = 100
nsitesrow = 25
nsitescol = 25
nsitemax = nsitescol * nsitescol
time_rough = 60
Name = ['Suraj', 'Shaswat', 'Himanshu', 'Debankit']
Sub_folder = ['Single', 'Multi']
for j in Sub_folder:
    for i in Name:
        if not os.path.exists(f'{j}/{i}/'):
            os.makedirs(f'{j}/{i}/')
        if j == 'Multi':
            rough_vs_temp_m(nsitesrow, nsitescol, timef,j, i)
            ads_rate_variation_m(nsitesrow, nsitescol, h_max,j, i)
            #rough_vs_rates_m(nsitesrow, nsitescol, h_max,j, i)
            flux_vs_desorption_m(timef, h_max, nsitesrow, nsitescol,j, i)
            avg_roughness_m(nsitesrow, nsitescol, time_rough,j, i)
            multi_2d(nsitesrow, nsitescol, timef, j, i)
            multi_2d_flux(nsitesrow, nsitescol, timef,j,i)
            multi_2D_gascomposition(nsitesrow, nsitescol, timef,j,i)
            grain_size_vs_temp(timef, nsitesrow, nsitescol, j, i)
        elif j == 'Single':
            rough_vs_temp(nsitesrow, nsitescol, timef,j, i)
            ads_rate_variation(nsitesrow, nsitescol, h_max,j, i)
            rough_vs_rates(nsitesrow, nsitescol, h_max,j, i)
            flux_vs_desorption(timef, h_max, nsitesrow, nsitescol,j, i)
            avg_roughness(nsitesrow, nsitescol,time_rough,j, i)
            avg_roughness(100,100,time_rough,j,i)

