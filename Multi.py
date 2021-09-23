import math
import random

import numpy as np


# import order_parameter

# Get Rate Hop
def get_rate_hop(site, hop_site, height, nnlist, particle_mat, T):
    site_self = 0
    site_other = 0
    neighbour_self = 0
    neighbour_other = 0
    for i in range(4):
        last_pos = len(particle_mat[site]) - 1
        particle = particle_mat[site][last_pos]
        if height[nnlist[i, site]] >= height[site]:
            if particle == particle_mat[nnlist[i, site]][last_pos]:
                site_self = site_self + 1
            else:
                site_other = site_other + 1
    for i in range(4):
        particle = particle_mat[site][len(particle_mat[site]) - 1]
        ne = nnlist[i, hop_site]
        # particle = particle_mat[hop_site][len(particle_mat[hop_site]) - 1]
        if nnlist[i, hop_site] == site:
            if height[nnlist[i, hop_site]][0] > height[hop_site][0] + 1:
                pos_nei = len(particle_mat[site]) - 1
                if particle == particle_mat[ne][pos_nei]:
                    neighbour_self = neighbour_self + 1
                else:
                    neighbour_other = neighbour_other + 1
        else:
            if height[nnlist[i, hop_site]] > height[hop_site]:
                pos_nei = len(particle_mat[nnlist[i, hop_site]]) - 1
                if particle == particle_mat[nnlist[i, hop_site]][pos_nei]:
                    neighbour_self = neighbour_self + 1
                else:
                    neighbour_other = neighbour_other + 1

    count = site_self + 0.333 * site_other
    countf = neighbour_self + 0.333 * neighbour_other
    count_ = countf - count
    k = (10 ** 13) / 4
    E_d = 1.58
    E_n = 0.27
    new_rate = k * math.exp((-E_d - 0.5 * -count_ * E_n) / (8.617343 * 10 ** (-5) * T))
    return new_rate


# Get Rate Desorption
def get_rate_desorption(site, height, nnlist, T):
    count = 0
    for i in range(4):
        if height[nnlist[i, site]] >= height[site]:
            count = count + 1
    k = 10 ** 13
    E_s = 2.32
    E_n = 0.27
    new_rate = k * math.exp(-(E_s + count * E_n) / (8.617343 * 10 ** (-5) * T))
    return new_rate


def initiate_nnlist(nnlist, nsitecol, nsiterow, nsite_max):
    for i in range(nsite_max):
        for j in range(4):
            if j == 0:
                if i < nsitecol:
                    nnlist[j, i] = nsite_max
                else:
                    nnlist[j, i] = i - nsitecol
            elif j == 1:
                if i >= nsitecol * (nsiterow - 1):
                    nnlist[j, i] = nsite_max
                else:
                    nnlist[j, i] = i + nsitecol
            elif j == 2:
                if ((i + 1) % nsitecol) == 0:
                    nnlist[j, i] = nsite_max
                else:
                    nnlist[j, i] = i + 1
            elif j == 3:
                if (i % nsitecol) == 0:
                    nnlist[j, i] = nsite_max
                else:
                    nnlist[j, i] = i - 1
    return nnlist


def update_list(rate, site, hopdir, wait_time, newrate):
    rate[hopdir][site] = newrate
    if newrate != 0:
        tau = -(math.log(random.random())) / newrate
        wait_time[hopdir][site] = tau
    else:
        wait_time[hopdir][site] = 10 ** 23


# eml function
def get_rate(site, A_comp, B_comp, height, nnlist, rate, wait_time, particle_mat, T, rate_adsorption):
    for hopdir in range(4):
        if height[site][0] != 0 and (height[nnlist[hopdir, site]][0] < height[site][0] + 1) and \
                height[nnlist[hopdir, site]][0] != -10 ** 23:
            if height[site][0] - height[nnlist[hopdir, site]][0] < 6:
                hop_site = nnlist[hopdir, site]
                new_rate = get_rate_hop(site, hop_site, height, nnlist, particle_mat, T)
            else:
                new_rate = get_rate_desorption(site, height, nnlist, T)
        else:
            new_rate = 0.0
        update_list(rate, site, hopdir, wait_time, new_rate)
    hopdir = 4
    if height[site][0] != 0:
        new_rate = get_rate_desorption(site, height, nnlist, T)
    else:
        new_rate = 0
    update_list(rate, site, hopdir, wait_time, new_rate)
    hopdir = 5
    new_rate = A_comp * rate_adsorption
    update_list(rate, site, hopdir, wait_time, new_rate)
    hopdir = 6
    new_rate = B_comp * rate_adsorption
    update_list(rate, site, hopdir, wait_time, new_rate)


def event(minimum, height, nsites_max, nnlist, particle_mat):
    event_dir = (minimum + 1) // nsites_max
    reminder = (minimum + 1) % nsites_max
    if reminder != 0:
        site = reminder - 1
    else:
        site = nsites_max - 1
        event_dir = event_dir - 1
    if event_dir == 0:
        height[site][0] = height[site][0] - 1
        height[nnlist[event_dir, site]][0] = height[nnlist[event_dir, site]][0] + 1
        particle = particle_mat[site][len(particle_mat[site]) - 1]
        del particle_mat[site][len(particle_mat[site]) - 1]
        particle_mat[nnlist[event_dir, site]].append(particle)
    elif event_dir == 1:
        height[site][0] = height[site][0] - 1
        height[nnlist[event_dir, site]][0] = height[nnlist[event_dir, site]][0] + 1
        particle = particle_mat[site][len(particle_mat[site]) - 1]
        del particle_mat[site][len(particle_mat[site]) - 1]
        particle_mat[nnlist[event_dir, site]].append(particle)
    elif event_dir == 2:
        height[site][0] = height[site][0] - 1
        height[nnlist[event_dir, site]][0] = height[nnlist[event_dir, site]][0] + 1
        particle = particle_mat[site][len(particle_mat[site]) - 1]
        del particle_mat[site][len(particle_mat[site]) - 1]
        particle_mat[nnlist[event_dir, site]].append(particle)
    elif event_dir == 3:
        height[site][0] = height[site][0] - 1
        height[nnlist[event_dir, site]][0] = height[nnlist[event_dir, site]][0] + 1
        particle = particle_mat[site][len(particle_mat[site]) - 1]
        del particle_mat[site][len(particle_mat[site]) - 1]
        particle_mat[nnlist[event_dir, site]].append(particle)
    elif event_dir == 4:
        height[site][0] = height[site][0] - 1
        del particle_mat[site][len(particle_mat[site]) - 1]
    elif event_dir == 5:
        height[site][0] = height[site][0] + 1
        particle_mat[site].append(1.0)
    elif event_dir == 6:
        height[site][0] = height[site][0] + 1
        particle_mat[site].append(2.0)
    return height, particle_mat


def get_grain_size(nsite_col, nsite_row, particle, particle_):
    order = np.zeros(nsite_row)
    for p in range(1, nsite_row + 1):
        count = 0
        for i in range(nsite_row - p + 1):
            for j in range(nsite_col - p + 1):
                flag = True
                for m in range(i, i + p):
                    for n in range(j, j + p):
                        if particle_[m][n] != particle:
                            flag = False
                if flag:
                    count += 1
        order[p - 1] = count
    #print(order)
    sum_ = np.sum(order)
    for i in range(len(order)):
        order[i] = order[i] * (i + 1)
    order_parameter = np.sum(order) / sum_
    return order_parameter


def main_kmc_code(timef, A_comp, B_comp, nsitesrow, nsitescol, T, rate_adsorption):
    nsites_max = nsitesrow * nsitescol
    events = 7
    height = np.zeros((nsitescol * nsitesrow, 1))
    nnlist = np.zeros((4, nsites_max), int)
    particle_mat = np.ones((nsites_max, 1)) * 2
    nnlist = initiate_nnlist(nnlist, nsitescol, nsitesrow, nsites_max)
    # ...MAIN LOOP (until final system time is reached)
    height = height.tolist()
    particle_mat = particle_mat.tolist()
    height.append([-10 ** 23])
    time_ = 0.0
    y = 0
    while time_ < timef:
        rate = np.full([events, nsites_max], 0.0)
        wait_time = np.full([events, nsites_max], 10 ** 23)
        for site in range(nsites_max):
            get_rate(site, A_comp, B_comp, height, nnlist, rate, wait_time, particle_mat, T, rate_adsorption)
        minimum = np.argmin(wait_time)
        min_time = np.min(wait_time)
        height, particle_mat = event(minimum, height, nsites_max, nnlist, particle_mat)
        time_ += min_time
        y += 1
    height_ = np.zeros((nsitesrow, nsitescol))
    particle_ = np.zeros((nsitesrow, nsitescol))
    for i in range(nsitesrow):
        for j in range(nsitescol):
            site = i * nsitesrow + j
            height_[i][j] = height[site][0]
            particle_ty = particle_mat[site][len(particle_mat[site]) - 1]
            particle_[i][j] = particle_ty
    return height_, particle_, particle_mat


def make_list_mat(particle_mat, nsiterow, nsitecol):
    length_1 = len(max(particle_mat, key=len))
    for i in range(nsiterow * nsitecol):
        length_2 = len(particle_mat[i])
        for j in range(length_2, length_1):
            particle_mat[i].append(0.0)
    # print(particle_mat)
