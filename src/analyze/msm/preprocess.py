import matplotlib.pyplot as plt
import numpy as np
import pickle
import torch
import torch.nn as nn
import mdshare  # for trajectory data
from tqdm import tqdm  # progress bar
from deeptime.util.data import TrajectoryDataset, TrajectoriesDataset
from deeptime.util.validation import implied_timescales, ck_test
from deeptime.plots import plot_implied_timescales, plot_ck_test
from deeptime.clustering import KMeans
from deeptime.markov import TransitionCountEstimator
from deeptime.markov.msm import BayesianMSM
import deeptime.markov as markov
import matplotlib.pyplot as plt

def _main(datafile, 
            numclusts,
            out_dtrajs,
            out_data,
            Timescale_image,
            discrete_image,
            ck_image,
            device="0"):

    if torch.cuda.is_available():
        device = torch.device(f"cuda:{device}")
        torch.backends.cudnn.benchmark = True
    else:
        device = torch.device("cpu")
    torch.set_num_threads(12)
    
    data = np.load(datafile)#"./run-4/selection_runs/selection-0/sd4_projection.npy")
    data = [data.astype(np.float32)]
    
    dataset = TrajectoriesDataset.from_numpy(1, data)


    cluster = KMeans(numclusts, progress=tqdm).fit_fetch(data)
    dtrajs = [cluster.transform(x) for x in data]
    msm_estimator = markov.msm.MaximumLikelihoodMSM(
                    reversible=True,
                    stationary_distribution_constraint=None
                    )
    
    lagtimes = np.arange(1, 30, dtype=np.int32)
    models = [msm_estimator.fit(dtrajs, lagtime=lag).fetch_model() for lag in tqdm(lagtimes)]
    
    ax = plot_implied_timescales(implied_timescales(models))
    ax.set_yscale('log')
    ax.set_xlabel('lagtime')
    ax.set_ylabel('timescale')
    plt.savefig(Timescale_image)
    plt.close()
    
    
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    np.savetxt(out_dtrajs, np.array(dtrajs).T, fmt='%i')
    np.savetxt(out_data, data[0])
    #counts,ybins,xbins,image = np.histogram2d(projections[0][:,0], projections[0][:,1], bins=100, norm=LogNorm())
    cb = ax2.hist2d(data[0][:,0], data[0][:,1], bins=20)
    #plt.colorbar(cb, ax=ax2)
    ax1.plot(dtrajs)
    ax1.set_xlim(0,40000)
    ax1.set_xlabel('t')
    ax1.set_ylabel('state')
    ax2.set_xlabel('SD4-1')
    ax2.set_ylabel('SD4-2')
    plt.savefig(discrete_image)
    plt.close()
    
    dtrajs = np.loadtxt(out_dtrajs)
    dtrajs = dtrajs.astype(np.int32)
    bmsms = [BayesianMSM(lagtime=lag).fit_fetch(dtrajs) for lag in tqdm([10, 20, 25, 30, 40, 50])]#, 60, 70, 80, 90, 100])]
    ck_test = bmsms[5].ck_test(bmsms, 6)
    plot_ck_test(ck_test)
    plt.savefig(ck_image)
    plt.close()
    return

def main(datafilelist, 
            numclusts,
            out_dtrajs,
            out_data,
            Timescale_image,
            discrete_image,
            ck_image,
            device="0"):

    if torch.cuda.is_available():
        device = torch.device(f"cuda:{device}")
        torch.backends.cudnn.benchmark = True
    else:
        device = torch.device("cpu")
    torch.set_num_threads(12)
    
    data = [np.load(datafilelist[0])] 
    print(data)
    for datafile in datafilelist[1:]:
        data_it = np.load(datafile)#"./run-4/selection_runs/selection-0/sd4_projection.npy")
        data.extend([data_it.astype(np.float32)])
    
    dataset = TrajectoriesDataset.from_numpy(1, data)

    cluster = KMeans(numclusts, progress=tqdm).fit_fetch(data)
    dtrajs = [cluster.transform(x) for x in data]
    msm_estimator = markov.msm.MaximumLikelihoodMSM(
                    reversible=True,
                    stationary_distribution_constraint=None
                    )
    
    lagtimes = np.arange(1, 300, dtype=np.int32)
    models = [msm_estimator.fit(dtrajs, lagtime=lag).fetch_model() for lag in tqdm(lagtimes)]
    
    ax = plot_implied_timescales(implied_timescales(models))
    ax.set_yscale('log')
    ax.set_xlabel('lagtime')
    ax.set_ylabel('timescale')
    plt.savefig(Timescale_image)
    plt.close()
    
    
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    np.savetxt(out_dtrajs, np.array(dtrajs).T, fmt='%i')
    np.savetxt(out_data, data[0])
    #counts,ybins,xbins,image = np.histogram2d(projections[0][:,0], projections[0][:,1], bins=100, norm=LogNorm())
    cb = ax2.hist2d(data[0][:,0], data[0][:,1], bins=20)
    #plt.colorbar(cb, ax=ax2)
    ax1.plot(dtrajs)
    ax1.set_xlim(0,40000)
    ax1.set_xlabel('t')
    ax1.set_ylabel('state')
    ax2.set_xlabel('SD4-1')
    ax2.set_ylabel('SD4-2')
    plt.savefig(discrete_image)
    plt.close()
    
    dtrajs = np.loadtxt(out_dtrajs)
    dtrajs = np.array(dtrajs.astype(np.int32)).T
    bmsms = [BayesianMSM(lagtime=lag).fit_fetch(dtrajs) for lag in tqdm([10, 20, 25, 30, 40, 50])]#, 60, 70, 80, 90, 100])]
    ck_test = bmsms[5].ck_test(bmsms, 6)
    plot_ck_test(ck_test)
    plt.savefig(ck_image)
    plt.close()
    return



if __name__ == "__main__":
    # Define a custom argument type for a list of strings

    import argparse
    import os
    def list_of_strings(arg):
        return arg.split(',')


    parser = argparse.ArgumentParser()

    parser.add_argument('-D',
                        '--datafilelist',
                        type=list_of_strings,
                        help='datafile list from heiranca script')

    parser.add_argument('-c',
                        '--numclusts',
                        type=int,
                        help='number of clusters')

    parser.add_argument('-o',
                        '--outdir',
                        type=str,
                        help='directory for output data')

    parser.add_argument('-I',
                        '--imdir',
                        type=str,
                        help='directory for output images')

    parser.add_argument('-d',
                        '--device',
                        type=str,
                        required=False,
                        default="0",
                        help='cuda device')
    
    #parser.add_argument('--str-list',
    #                    type=list_of_strings)

    

    args = parser.parse_args()


    try:
        os.mkdir(f'{args.outdir}')
    except:
        print(f' cant make {args.outdir}')
        pass

    try:
        os.mkdir(f'{args.imdir}')
    except:
        pass

    main(args.datafilelist, 
            args.numclusts,
            f'{args.outdir}/dtrajs.dat',
            f'{args.outdir}/data.dat',
            f'{args.imdir}/Timescale.png',
            f'{args.imdir}/discrete.png',
            f'{args.imdir}/ck.png',
            args.device)
