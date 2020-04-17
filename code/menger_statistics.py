#!/bin/python
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from lorenz import get_lorenz
from curvature import get_menger, get_menger_centered_avg, get_menger_centered_avg_embed
from pytisean import tiseanio
from datasets import *
from plots import *
import pandas as pd
import sys
import os
plt.switch_backend('agg');
plt.rcParams.update({'figure.max_open_warning': 0})

tau_max = 100;

dataset_dict = {'lorenz': data_lorenz, 'pendulum': data_pendulum, 'pendulum_omega': data_pendulum_omega, 'rossler': data_rossler, 'musical': data_musical} 
title_dict = {'lorenz': 'Lorenz', 'pendulum': 'Pendulum', 'pendulum_omega': 'Pendulum Omega', 'rossler': 'Rossler', 'musical': 'Piano'}

def main():  

    dataset_label = sys.argv[1]; 

    try:
        dataset_func = dataset_dict[sys.argv[1]]
        plot_title = title_dict[sys.argv[1]];
    except KeyError:
              raise ValueError('Dataset not found')

    tau_avg = int(sys.argv[2]);
    n_neighbors = int(sys.argv[3]); 
    noise_level = float(sys.argv[4]);
    m = int(sys.argv[5]);

    prop_label = 'noise_level_' + str(noise_level) + '_tau_avg_' + str(tau_avg) \
                 + '_n_neighbors_' + str(n_neighbors) + '_m_' + str(m);
    
    dirname = dataset_label + '_' + prop_label
    if not os.path.exists(dirname):
        os.makedirs(dirname)    

    series = dataset_func(n=500000,timestep=0.01);
    
    np.savetxt(dirname + '/series.txt', series);
    window = 1000;
    mean = [];
    var = [];
    mad = [];
    sample_skew = [];
    pearson2 = [];

    #for start in np.arange(0,len(series) - window, window):
    for tau in range(1,tau_max):
        embed, err = tiseanio('delay', '-m', m, '-d', tau, data=series, silent=True);
        np.savetxt(dirname + '/embed_%02d.txt' % (tau), embed);
     
        menger = get_menger_centered_avg_embed(embed, tau_avg, n_neighbors)
        menger = menger[~np.isnan(menger)];
        
        plot_batman(menger, tau, dirname);
          
        mean.append(np.mean(menger));
        var.append(np.var(menger));
        sample_skew.append(np.mean((menger - np.mean(menger)) ** 3) / np.std(menger) ** 3)
        pearson2.append((np.mean(menger) - np.median(menger)) / np.std(menger))
        np.savetxt(dirname + '/menger_%02d.txt' % (tau), menger);
        
        plot_histogram(menger, tau, dirname); 
        plot_reconstruction(series, tau, dirname);

    mutual, err = tiseanio('mutual', '-D', tau_max, data=series, silent=True);
    mutual = mutual[:,1]

    var = np.asarray(var);
    mean = np.asarray(mean);
    mutual = np.asarray(mutual);
   
    plt.figure();
    #plt.title(plot_title + ' Variance'); 
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax1.plot(var, color = 'b');
    ax2.plot(mutual, color = 'r')
    #plt.axvline(x = 1, label = '$tau = 1$', color = 'b', linestyle = 'dashed')
    #plt.axvline(x = 18, label = '$tau = 18$', color = 'g', linestyle = 'dashed')
    #plt.axvline(x = 30, label = '$tau = 30$', color = 'r', linestyle = 'dashed')
    ax1.set_ylabel('Variance of curvature', fontsize=16, color = 'b');
    ax2.set_ylabel('Mutual information', fontsize=16, color = 'r')
    #plt.ylabel('Metric');
    plt.xlabel('Tau');
    #plt.legend();
    plt.savefig(dirname + '/' + 'property_variance_' + prop_label + '.png');
    np.savetxt(dirname + '/' + 'property_variance_' + prop_label + '.txt', var);

    plt.figure();
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    #plt.title(plot_title + ' Mean'); 
    ax1.plot(mean, color = 'b');
    ax2.plot(mutual, color = 'r');
    #plt.axvline(x = 1, label = '$tau = 1$', color = 'b', linestyle = 'dashed')
    #plt.axvline(x = 18, label = '$tau = 18$', color = 'g', linestyle = 'dashed')
    #plt.axvline(x = 30, label = '$tau = 30$', color = 'r', linestyle = 'dashed')
    ax1.set_ylabel('Mean curvature', fontsize=16, color = 'b');
    ax2.set_ylabel('Mutual information', fontsize=16, color = 'r')
    plt.xlabel('Tau');
    #plt.legend();
    plt.savefig(dirname + '/' + 'property_mean_' + prop_label + '.png');
    np.savetxt(dirname + '/' + 'property_mean_' + prop_label + '.txt', mean);

    plt.figure();
    plt.title(plot_title + ' Pearson 2'); 
    plt.plot(pearson2, label='menger pearson2');
    #plt.axvline(x = 1, label = '$tau = 1$', color = 'b', linestyle = 'dashed')
    #plt.axvline(x = 18, label = '$tau = 18$', color = 'g', linestyle = 'dashed')
    #plt.axvline(x = 30, label = '$tau = 30$', color = 'r', linestyle = 'dashed')
    plt.ylabel('Metric');
    plt.xlabel('Tau');
    plt.legend();
    plt.savefig(dirname + '/' + 'property_pearson2_' + prop_label + '.png');

    plt.figure();
    plt.title(plot_title + ' Sample Skewness'); 
    plt.plot(sample_skew, label='menger sample skewness');
    #plt.axvline(x = 1, label = '$tau = 1$', color = 'b', linestyle = 'dashed')
    #plt.axvline(x = 18, label = '$tau = 18$', color = 'g', linestyle = 'dashed')
    #plt.axvline(x = 30, label = '$tau = 30$', color = 'r', linestyle = 'dashed')
    plt.ylabel('Metric');
    plt.xlabel('Tau');
    plt.legend();
    plt.savefig(dirname + '/' + 'property_sample_skew_' + prop_label + '.png');

    plt.figure();
    plt.title(plot_title + ' Mutual'); 
    plt.plot(mutual, label='menger mutual');
    #plt.axvline(x = 1, label = '$tau = 1$', color = 'b', linestyle = 'dashed')
    #plt.axvline(x = 18, label = '$tau = 18$', color = 'g', linestyle = 'dashed')
    #plt.axvline(x = 30, label = '$tau = 30$', color = 'r', linestyle = 'dashed')
    plt.ylabel('Metric');
    plt.xlabel('Tau');
    plt.legend();
    plt.savefig(dirname + '/' + 'property_mutual_' + prop_label + '.png');

    plt.figure();
    #plt.title(plot_title); 
    plt.plot(var / var[0], label='curvature variance');
    plt.plot(mean / mean[0], label='curvature mean');
    #plt.plot(pearson2 / pearson2[0], label='menger pearson2');
    #plt.plot(sample_skew/ sample_skew[0], label='menger sample skew');
    plt.plot(mutual / mutual[0], label='mutual information', color = 'r');
    #plt.axvline(x = 1, label = '$tau = 1$', color = 'b', linestyle = 'dashed')
    #plt.axvline(x = 18, label = '$tau = 18$', color = 'g', linestyle = 'dashed')
    #plt.axvline(x = 30, label = '$tau = 30$', color = 'r', linestyle = 'dashed')
    plt.ylabel('Metric', fontsize=16);
    plt.xlabel('Tau', fontsize=16);
    plt.legend();
    plt.savefig(dirname + '/' + 'property_all_' + prop_label + '.png');

main();
