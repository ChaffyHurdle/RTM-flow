import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.io import loadmat

CASE = 1 # [1: square, 2: fork, 3: torus]
EXAMPLE = 1
x_ticks = [[0,0.5,1],[0,0.5,1],[-1,0,1]]
SENSORS = 100

folder_path = '../Case' + str(CASE) + '/Example' + str(EXAMPLE) + '/'
sensor_path = '../Case' + str(CASE) + '/meshes_sensors/'


if __name__ == '__main__':
    true_boundary = []
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

    # Load sensors and means/vars/timers (72,100)
    sensor_locs = loadmat(sensor_path + 'sensor_locs_' + str(SENSORS) + '.mat')['sensor_locs']

    lmap_means = loadmat(folder_path + 'lmap_means.mat')['lmap_means']
    lmap_vars = loadmat(folder_path + 'lmap_vars.mat')['lmap_vars']
    # lmap_times= loadmat(folder_path + 'lmap_times.mat')['lmap_times']
    # matern_args = loadmat(folder_path + 'matern_args.mat')['matern_args']
    # exp_args= loadmat(folder_path + 'exp_args.mat')['exp_args']
    # print(lmap_times)
    # print(matern_args)
    # print(exp_args)


    # Load true front inds
    for j in range(1,6):
        true_boundary.append(loadmat(folder_path + 'front' + str(j) + '.mat')['nodes_t'])

    # Create triangulations
    triang_fwd = tri.Triangulation(p_fwd[0], p_fwd[1],t_fwd.T[:,:3]-1)
    triang_inv = tri.Triangulation(p_inv[0], p_inv[1],t_inv.T[:,:3]-1)

    # Plot truth
    global_min = np.min(np.log(K_true.ravel()))
    global_max = np.max(np.log(K_true.ravel()))

    rows = 2
    columns = 6
    figsizex = 14
    figsizey = 5
    fig,axes = plt.subplots(rows,columns,figsize=(figsizex,figsizey),constrained_layout=True)

    # Plot truth
    im1 = axes[0,0].tripcolor(triang_fwd, facecolors=np.log(K_true.ravel()),vmin=global_min,vmax=global_max,cmap='turbo', shading='flat')
    axes[0,0].set_xlim([np.min(p_fwd[0]),np.max(p_fwd[0])])
    axes[0,0].set_ylim([np.min(p_fwd[1]),np.max(p_fwd[1])])
    axes[0,0].set_title('$u^\dagger$')
    axes[0,0].scatter(sensor_locs[:,0],sensor_locs[:,1],s=6,c='black')
    axes[0,0].scatter(sensor_locs[:,0],sensor_locs[:,1],s=3,c='white')
    axes[0,0].set_xticks(x_ticks[CASE-1])
    axes[0,0].set_yticks(x_ticks[CASE-1])

    axes[1,0].axis('off')

    # Plot means
    [axes[0,i+1].tripcolor(triang_inv, facecolors=lmap_means[:,i],vmin=global_min,vmax=global_max,edgecolors='none',cmap='turbo', shading='flat') for i in range(5)]
    im2 = [axes[1,i+1].tripcolor(triang_inv, facecolors=lmap_vars[:,i],vmin=0,vmax=np.max(lmap_vars),edgecolors='none',cmap='turbo', shading='flat') for i in range(5)]
    [axes[j,i+1].set_xticks([]) for i in range(5) for j in range(2)]
    [axes[j,i+1].set_yticks([]) for i in range(5) for j in range(2)]
    [axes[j,i+1].set_xlim([np.min(p_fwd[0]),np.max(p_fwd[0])]) for i in range(5) for j in range(2)]
    [axes[j,i+1].set_ylim([np.min(p_fwd[1]),np.max(p_fwd[1])]) for i in range(5) for j in range(2)]

    [axes[0, j+1].set_title(r'$\bar{{u}}_{{MAP}}^{{({})}}$'.format(j+1)) for j in range(5)]
    [axes[0, j+1].axis('off') for j in range(5)]
    [axes[1, j+1].set_title(r'$\mathcal{{C}}_{{MAP}}^{{({})}}$'.format(j+1)) for j in range(5)]
    [axes[1, j+1].axis('off') for j in range(5)]


    # Plot fronts
    for i in range(5):
        bndry_locs_t = true_boundary[i]
        for k in range(len(bndry_locs_t)):
            axes[0,i+1].plot(
                np.array([p_fwd.T[bndry_locs_t[k,0]-1,0], p_fwd.T[bndry_locs_t[k,1]-1,0]]),
                np.array([p_fwd.T[bndry_locs_t[k,0]-1,1], p_fwd.T[bndry_locs_t[k,1]-1,1]]),
                color='white'
            )
            axes[1,i+1].plot(
                np.array([p_fwd.T[bndry_locs_t[k,0]-1,0], p_fwd.T[bndry_locs_t[k,1]-1,0]]),
                np.array([p_fwd.T[bndry_locs_t[k,0]-1,1], p_fwd.T[bndry_locs_t[k,1]-1,1]]),
                color='white'
            )

    fig.colorbar(im1, ax=axes[0,5],fraction=1)
    fig.colorbar(im2[0], ax=axes[1,5],fraction=1)

    # plt.savefig(folder_path + 'lmap_seq.eps', format='eps', bbox_inches='tight')
    # plt.savefig(folder_path + 'lmap_seq.jpg', format='jpg', bbox_inches='tight',dpi=400)
    plt.show()