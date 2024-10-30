import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os

def show_scatter(data: np.ndarray, color: np.ndarray, output: str):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("white")
    ff = ax.scatter(data[:, 0], data[:, 1], data[:, 2], c=color)
    #ax.view_init(90, 90, 90)
    plt.colorbar(ff)
    plt.savefig(output, dpi=300)
    plt.close()


def main(dihdata, Sdata, ermsd, ZPrj4, zauto, out_image):
    dihedrals = np.load(dihdata)
    
    for i in range(dihedrals.shape[1]):
        dihedral_test = dihedrals[:,i,0]
    
    plt.plot(dihedral_test[::100])
    plt.savefig(f'{out_image}/cos_dihedral_test.png', dpi=300)
    plt.close()
    
    S = np.load(Sdata)
    plt.semilogy(S, "ro-")
    plt.savefig(f'{out_image}/S.png', dpi=300)
    plt.close()
    
    e_rmsd = np.load(ermsd)
    ZPrj4 = np.load(ZPrj4)
    
    show_scatter(ZPrj4[:, :3], e_rmsd[-1], f'{out_image}/ZPrj4.png')
    
    z = np.load(zauto)#'run-4/selection_runs/selection-0/z.npy')
    show_scatter(z, e_rmsd[-1], f'{out_image}/z_ermsds.png')



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d',
                        '--dihdata',
                        type=str,
                        help='Dihedral Data')
    
    parser.add_argument('-S',
                        '--Sdata',
                        type=Path,
                        help='S Data')
    
    parser.add_argument('-e',
                        '--ermsd',
                        type=Path,
                        help='e-rmsd data',
                        )
    
    parser.add_argument('-Z',
                        '--ZPrj4',
                        type=Path,
                        help='ZProj sd4 data',
                        )
    
    parser.add_argument('-z',
                        '--zauto',
                        type=Path,
                        help='Autoencoder projection',
                        )
    
    
    parser.add_argument('-o',
                        '--output',
                        type=Path,
                        help='Output path for Images',
                        )
    
    args = parser.parse_args()
    
    try:
        os.mkdir(args.output)
    except OSError as error:
        print(error)
    
    main(args.dihdata, args.Sdata, args.ermsd, args.ZPrj4, args.zauto, args.output)
    
    
