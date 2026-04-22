import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.io import loadmat


CASE = 1 # [1: square, 2: fork, 3: torus]
EXAMPLE = 1
x_ticks = [[0,0.5,1],[0,0.5,1],[-1,0,1]]
SENSORS = 25

folder_path = '../Case' + str(CASE) + '/Comparison/Example' + str(EXAMPLE) + '/'
sensor_path = '../Case' + str(CASE) + '/meshes_sensors/'

if __name__ == '__main__':
    means = []
    vars = []
    timers = []

    # Load meshes
    K_true = loadmat(folder_path +'K_true.mat')['K_true']
    p_fwd = np.array(loadmat(sensor_path + 'p_fwd.mat')['p_fwd'])
    e_fwd = np.array(loadmat(sensor_path + 'e_fwd.mat')['e_fwd'])
    t_fwd = np.array(loadmat(sensor_path + 't_fwd.mat')['t_fwd'])
    p_inv = np.array(loadmat(sensor_path + 'p_inv.mat')['p_inv'])
    e_inv = np.array(loadmat(sensor_path + 'e_inv.mat')['e_inv'])
    t_inv = np.array(loadmat(sensor_path + 't_inv.mat')['t_inv'])

    # Load sensors and means/vars/timers
    sensor_locs = loadmat(sensor_path + 'sensor_locs_' + str(SENSORS) + '.mat')['sensor_locs']

    eki_mean500 = loadmat(folder_path + 'eki_mean500.mat')['eki_mean500']
    eki_mean1000 = loadmat(folder_path + 'eki_mean1000.mat')['eki_mean1000']
    eki_mean5000 = loadmat(folder_path + 'eki_mean5000.mat')['eki_mean5000']
    lmap_mean = loadmat(folder_path + 'lmap_mean.mat')['lmap_mean']

    eki_var500 = loadmat(folder_path + 'eki_var500.mat')['eki_var500']
    eki_var1000 = loadmat(folder_path + 'eki_var1000.mat')['eki_var1000']
    eki_var5000 = loadmat(folder_path + 'eki_var5000.mat')['eki_var5000']
    lmap_var = loadmat(folder_path + 'lmap_var.mat')['lmap_var']


    # Create triangulations
    triang_fwd = tri.Triangulation(p_fwd[0], p_fwd[1],t_fwd.T[:,:3]-1)
    triang_inv = tri.Triangulation(p_inv[0], p_inv[1],t_inv.T[:,:3]-1)

    # Plot truth
    global_min = np.min(np.log(K_true.ravel()))
    global_max = np.max(np.log(K_true.ravel()))

    rows = 2
    columns = 5
    figsizex = 12
    figsizey = 5
    fig,axes = plt.subplots(rows,columns,figsize=(figsizex,figsizey),constrained_layout=True)

    # Plot truth
    im1 = axes[0,0].tripcolor(triang_fwd, facecolors=np.log(K_true.ravel()),vmin=global_min,vmax=global_max,cmap='turbo', shading='flat')
    axes[0,0].set_title('$u^\dagger$')
    axes[0,0].scatter(sensor_locs[:,0],sensor_locs[:,1],s=6,c='black')
    axes[0,0].scatter(sensor_locs[:,0],sensor_locs[:,1],s=3,c='white')
    axes[0,0].set_xlim([0,1])
    axes[0,0].set_ylim([0,1])
    axes[0,0].set_xticks(x_ticks[CASE-1])
    axes[0,0].set_yticks(x_ticks[CASE-1])
    axes[1,0].axis('off')

    # Plot means
    im3 = axes[0,1].tripcolor(triang_inv, facecolors=eki_mean500.ravel(),vmin=global_min,edgecolors='none',vmax=global_max,cmap='turbo', shading='flat')
    im4 = axes[0,2].tripcolor(triang_inv, facecolors=eki_mean1000.ravel(),vmin=global_min,edgecolors='none',vmax=global_max,cmap='turbo', shading='flat')
    im5 = axes[0,3].tripcolor(triang_inv, facecolors=eki_mean5000.ravel(),vmin=global_min,edgecolors='none',vmax=global_max,cmap='turbo', shading='flat')
    im6 = axes[0,4].tripcolor(triang_inv, facecolors=lmap_mean.ravel(),vmin=global_min,edgecolors='none',vmax=global_max,cmap='turbo', shading='flat')

    # Plot vars
    im7 = axes[1,1].tripcolor(triang_inv, facecolors=eki_var500.ravel(),vmin=0,edgecolors='none',vmax=0.25,cmap='turbo', shading='flat')
    im8 = axes[1,2].tripcolor(triang_inv, facecolors=eki_var1000.ravel(),vmin=0,edgecolors='none',vmax=0.25,cmap='turbo', shading='flat')
    im9 = axes[1,3].tripcolor(triang_inv, facecolors=eki_var5000.ravel(),vmin=0,edgecolors='none',vmax=0.25,cmap='turbo', shading='flat')
    im10 = axes[1,4].tripcolor(triang_inv, facecolors=lmap_var.ravel(),vmin=0,edgecolors='none',vmax=0.25,cmap='turbo', shading='flat')

    [axes[i,j].set_xticks([]) for i in range(2) for j in range(1,5)]
    [axes[i,j].set_yticks([]) for i in range(2) for j in range(1,5)]
    [axes[i,j].set_xlim([0,1]) for i in range(2) for j in range(1,5)]
    [axes[i,j].set_ylim([0,1]) for i in range(2) for j in range(1,5)]
    [axes[i,j].set_title('$u^\dagger$') for i in range(2) for j in range(1,5)]

    axes[0, 1].set_title(r'$\bar{{u}}_{{EKI}}^{{(500)}}$')
    axes[0, 2].set_title(r'$\bar{{u}}_{{EKI}}^{{(1000)}}$')
    axes[0, 3].set_title(r'$\bar{{u}}_{{EKI}}^{{(5000)}}$')
    axes[0, 4].set_title(r'$\bar{{u}}_{{MAP}}^{{(5)}}$')

    axes[1, 1].set_title(r'$\mathcal{{C}}_{{EKI}}^{{(500)}}$')
    axes[1, 2].set_title(r'$\mathcal{{C}}_{{EKI}}^{{(1000)}}$')
    axes[1, 3].set_title(r'$\mathcal{{C}}_{{EKI}}^{{(5000)}}$')
    axes[1, 4].set_title(r'$\mathcal{{C}}_{{MAP}}^{{(5)}}$')

    fig.colorbar(im1, ax=axes[0,4],fraction=1)
    fig.colorbar(im10, ax=axes[1,4],fraction=1)

    # plt.savefig(folder_path + 'LMAP_vs_EKI.eps', format='eps', bbox_inches='tight')
    # plt.savefig(folder_path + 'LMAP_vs_EKI.jpg', format='jpg', bbox_inches='tight',dpi=400)
    plt.show()
