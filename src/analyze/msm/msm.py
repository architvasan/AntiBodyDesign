import numpy as np
from deeptime.markov import TransitionCountEstimator
import deeptime.markov as markov
import matplotlib.pyplot as plt
import networkx as nx

def msm_construct(dtrajs,
             lag,
             out_tmat,
             out_stat,
             logfile,
                ):

    count_estimator = TransitionCountEstimator(
        lagtime=lag,
        count_mode='sliding'
        )
    
    msm_estimator = markov.msm.MaximumLikelihoodMSM(
                    reversible=True,
                    stationary_distribution_constraint=None
                    )
    
    trajectory = np.loadtxt(dtrajs).astype(int).T
    msm = msm_estimator.fit(trajectory, lagtime=lag).fetch_model()
    
    np.savetxt(f'{out_tmat}', msm.transition_matrix)
    np.savetxt(f'{out_stat}', msm.stationary_distribution)

    f = open(f"{logfile}", "a")
    f.write(f"Number of states: {msm.n_states}")
    f.close()
    return msm


def draw_tmat(msm, outimage):
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    threshold = 1e-4
    title = f"Transition matrix with connectivity threshold {threshold:.0e}"
    G = nx.DiGraph()
    ax.set_title(title)
    for i in range(msm.n_states):
        G.add_node(i, title=f"{i+1}")
    for i in range(msm.n_states):
        for j in range(msm.n_states):
            if msm.transition_matrix[i, j] > threshold:
                G.add_edge(i, j, title=f"{msm.transition_matrix[i, j]:.3e}")
    
    edge_labels = nx.get_edge_attributes(G, 'title')
    pos = nx.fruchterman_reingold_layout(G)
    nx.draw_networkx_nodes(G, pos, ax=ax)
    nx.draw_networkx_labels(G, pos, ax=ax, labels=nx.get_node_attributes(G, 'title'));
    nx.draw_networkx_edges(G, pos, ax=ax, arrowstyle='-|>',
                           connectionstyle='arc3, rad=0.3')
    plt.savefig(outimage)
    plt.close()
    return


def pcca_construct(msm, 
                    n_metas,
                    probsfil,
                    assignfil,
                    tmatfil,
                    logfile):

    pcca = msm.pcca(n_metastable_sets=n_metas)

    f = open(logfile, "a")
    f.write(f"Memberships: {pcca.memberships.shape} \n")
    f.write("Stationary probability: \n")
    f.write(f"{pcca.coarse_grained_stationary_probability} \n")
    f.write("Coarse grained transition matrix \n")
    f.write(f"{pcca.coarse_grained_transition_matrix}")
    f.write(f"Metastable distributions shape: {pcca.metastable_distributions.shape}")
    f.close()

    np.savetxt(probsfil, pcca.coarse_grained_stationary_probability)
    np.savetxt(assignfil, np.array([int(p_as) for p_as in pcca.assignments]).T)
    np.savetxt(tmatfil, pcca.coarse_grained_transition_matrix)
    return pcca


def draw_pccatmat(pcca, imagefil):
    import networkx as nx
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    threshold = 1e-9
    title = f"PCCA Transition matrix with connectivity threshold {threshold:.0e}"
    G = nx.DiGraph()
    ax.set_title(title)
    for i in range(4):
        G.add_node(i, title=f"{i+1}")
    for i in range(4):
        for j in range(4):
            print(pcca.coarse_grained_transition_matrix[i, j])
            if pcca.coarse_grained_transition_matrix[i, j] > threshold:
                G.add_edge(i, j, title=f"{pcca.coarse_grained_transition_matrix[i, j]:.3e}")
    
    edge_labels = nx.get_edge_attributes(G, 'title')
    pos = nx.fruchterman_reingold_layout(G)
    nx.draw_networkx_nodes(G, pos, ax=ax)
    nx.draw_networkx_labels(G, pos, ax=ax, labels=nx.get_node_attributes(G, 'title'));
    nx.draw_networkx_edges(G, pos, ax=ax, arrowstyle='-|>',
                           connectionstyle='arc3, rad=0.3')
    
    plt.savefig(imagefil)
    plt.close()
    return

def assign_pccafrmes(dtrajs, 
                    pcca,
                    metafrmesout,
                    pcca_assignfil='none'):

    dtrajs = np.loadtxt(dtrajs)
    for t in range(len(dtrajs[0])):
        print(dtrajs[:,t])
        dtrajs_t = np.array([int(d) for d in dtrajs[:,t]])

        if pcca_assignfil != 'none':
            pcca_assign = np.loadtxt(pcca_assignfil)
            pcca_assign = np.array([int(pd) for pd in pcca_assign])

        else:
            pcca_assign = np.array([int(p_as) for p_as in pcca.assignments]).T

        pcca_clusters = {}
        pcca_frames = {}
        
        for pd in range(int(np.max(pcca_assign))+1):
            pcca_clusters[pd]=np.where(pcca_assign==pd)[0]
            pcca_frames[pd]=[]
            for c in pcca_clusters[pd]:
                pcca_frames[pd].extend(list(np.where(dtrajs_t==int(c))[0]))
            np.savetxt(f'{metafrmesout}_{pd}_{t}.dat', pcca_frames[pd], fmt='%i')
    return pcca_clusters, pcca_frames


def main(dtrajs,
        lag,
        n_metas,
        outdir,
        logfile,
        imdir='none',
        assign_pccafrms='none', 
        pcca_assignfil='none',
         ):

    f = open(logfile, "w")
    f.write(f"MSM logfile")
    f.close()

    msm = msm_construct(dtrajs,
                        lag,
                        f'{outdir}/tmat.dat',
                        f'{outdir}/stationary.dat',
                        logfile,
                        )



    pcca = pcca_construct(msm, 
                    n_metas,
                    f'{outdir}/pcca_probs.dat',
                    f'{outdir}/pcca_assignments.dat',
                    f'{outdir}/pcca_tmat.dat',
                    logfile)

    if imdir!='none':
        draw_tmat(msm, 
                f'{imdir}/tmat.png')
        
        draw_pccatmat(pcca, 
                    f'{imdir}/tmat_pcca.png')

    if assign_pccafrms!='none':
        import os
        try:
            os.mkdir(f'{outdir}/{assign_pccafrms}')
        except:
            pass

        pcca_clusters, pcca_frames = assign_pccafrmes(dtrajs, 
                                        pcca,
                                        f'{outdir}/{assign_pccafrms}/metafrmes',
                                        pcca_assignfil)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-D',
                        '--dtrajs',
                        type=str,
                        help='datafile from heiranca script')

    parser.add_argument('-L',
                        '--lag',
                        type=int,
                        help='lagtime')

    parser.add_argument('-nm',
                        '--n_metas',
                        type=int,
                        help='number of pcca metastable states')

    parser.add_argument('-O',
                        '--outdir',
                        type=str,
                        help='directory for output data')

    parser.add_argument('-I',
                        '--imdir',
                        type=str,
                        required=False,
                        default='none',
                        help='directory for output images')

    parser.add_argument('-ap',
                        '--assign_pccafrms',
                        type=str,
                        required=False,
                        default='none',
                        help='where to assign pccaframes in the output dir')

    parser.add_argument('-af',
                        '--pcca_assignfil',
                        type=str,
                        required=False,
                        default='none',
                        help='file that has the pcca assignments')

    parser.add_argument('-l',
                        '--logfile',
                        type=str,
                        required=False,
                        default='none',
                        help='logfile')

    args = parser.parse_args()

    try:
        os.mkdir(f'{args.outdir}')
    except:
        pass

    try:
        os.mkdir(f'{args.imdir}')
    except:
        pass

    main(args.dtrajs,
         args.lag,
         args.n_metas,
         args.outdir,
         #args.imdir,
         args.logfile,
         imdir=args.imdir,
         assign_pccafrms=args.assign_pccafrms, 
         pcca_assignfil=args.pcca_assignfil,
         )
