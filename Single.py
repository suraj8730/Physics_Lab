import math
import random

import numpy as np


# Get Rate Hop
def get_rate_hop(site, height, nnlist, T):
    count = 0
    for i in range(4):
        if height[nnlist[i, site]] >= height[site]:
            count = count + 1
    k = (10 ** 13) / 4
    E_d = 1.58
    E_n = 0.27
    new_rate = k * math.exp((-E_d - 0.5 * count * E_n) / (8.617343 * 10 ** (-5) * T))
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


# Make nnlist
def initiate_nnlist(nnlist, nsitescol, nsitesrow, nsite_max):
    for i in range(nsite_max):
        for j in range(4):
            if j == 0:
                if i < nsitescol:
                    nnlist[j, i] = nsite_max
                else:
                    nnlist[j, i] = i - nsitescol
            elif j == 1:
                if i >= nsitescol * (nsitesrow - 1):
                    nnlist[j, i] = nsite_max
                else:
                    nnlist[j, i] = i + nsitescol
            elif j == 2:
                if ((i + 1) % nsitescol) == 0:
                    nnlist[j, i] = nsite_max
                else:
                    nnlist[j, i] = i + 1
            elif j == 3:
                if (i % nsitescol) == 0:
                    nnlist[j, i] = nsite_max
                else:
                    nnlist[j, i] = i - 1
    return nnlist


# Update rate and wait_time matrix after rate calculation
def update_list(rate, site, hopdir, wait_time, newrate):
    rate[hopdir, site] = newrate
    if newrate != 0:
        tau = -(math.log(random.random())) / newrate
        wait_time[hopdir][site] = tau
    else:
        wait_time[hopdir][site] = 10 ** 23


# rate function
def get_rate(site, height, nnlist, rate, wait_time, T, rate_adsorption):
    for hopdir in range(4):
        neighbour = nnlist[hopdir, site]
        if height[site][0] != 0 and (height[neighbour][0] < height[site][0] + 1) and height[neighbour][0] != -10 ** 23:
            if height[site][0] - height[nnlist[hopdir, site]][0] < 6:
                new_rate = get_rate_hop(site, height, nnlist, T)
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
    new_rate = rate_adsorption
    update_list(rate, site, hopdir, wait_time, new_rate)


# Intiate the event
def event(minimum, height, nsites_max, nnlist):
    event_ = (minimum + 1) // nsites_max
    reminder = (minimum + 1) % nsites_max
    if reminder != 0:
        site = reminder - 1
    else:
        site = nsites_max - 1
        event_ = event_ - 1
    # print("Site %i and Event %i"%(site,event_))
    if event_ == 0:
        height[site][0] = height[site][0] - 1
        height[nnlist[event_, site]][0] = height[nnlist[event_, site]][0] + 1
    elif event_ == 1:
        height[site][0] = height[site][0] - 1
        height[nnlist[event_, site]][0] = height[nnlist[event_, site]][0] + 1
    elif event_ == 2:
        height[site][0] = height[site][0] - 1
        height[nnlist[event_, site]][0] = height[nnlist[event_, site]][0] + 1
    elif event_ == 3:
        height[site][0] = height[site][0] - 1
        height[nnlist[event_, site]][0] = height[nnlist[event_, site]][0] + 1
    elif event_ == 4:
        height[site][0] = height[site][0] - 1
    elif event_ == 5:
        height[site][0] = height[site][0] + 1
    return height


# main KMC code
def main_kmc_code(timef, nsitesrow, nsitescol, T, rate_adsorption):
    nsites_max = nsitesrow * nsitescol
    events = 6
    row_wait = events
    height = np.zeros((nsitescol * nsitesrow, 1))
    nnlist = np.zeros((4, nsites_max), int)
    nnlist = initiate_nnlist(nnlist, nsitescol, nsitesrow, nsites_max)
    # ...MAIN LOOP (until final system time is reached)
    height = height.tolist()
    height.append([-10 ** 23])
    time_ = 0.0
    while time_ < timef:
        rate = np.full([row_wait, nsites_max], 0.0)
        wait_time = np.full([row_wait, nsites_max], 10 ** 23)
        for site in range(nsites_max):
            get_rate(site, height, nnlist, rate, wait_time, T, rate_adsorption)
        minimum = np.argmin(wait_time)
        min_time = np.min(wait_time)
        height = event(minimum, height, nsites_max, nnlist)
        time_ += min_time
        if np.mean(height) > 357:
            time_ = timef + 1

    height_ = np.zeros((nsitesrow, nsitescol))
    for i in range(nsitesrow):
        for j in range(nsitescol):
            height_[i][j] = height[i * nsitesrow + j][0]
    return np.mean(height_), np.std(height_), height_
