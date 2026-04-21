import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.io import loadmat

folder_path = '../Case1/SensorDensity/Example2/'

CASE = 1 # [1: square, 2: fork, 3: torus]
EXAMPLE = 1
x_ticks = [[0,0.5,1],[0,0.5,1],[-1,0,1]]

folder_path = '../Case' + str(CASE) + '/SensorDensity/Example' + str(EXAMPLE) + '/'
sensor_path = '../Case' + str(CASE) + '/meshes_sensors/'

if __name__ == '__main__':
    true_boundary = []
    sensor_locs = []
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
    sensor_locs = []
    for s in [3,5,7,9,15]:
        sensors = int(s**2)
        sensor_locs.append(loadmat(sensor_path + 'sensor_locs_' + str(sensors) + '.mat')['sensor_locs'])

        means.append(loadmat(folder_path + 'means' + str(sensors) + '.mat')['means'])
        vars.append(loadmat(folder_path + 'vars' + str(sensors) + '.mat')['vars'])
        timers.append(loadmat(folder_path + 'timers' + str(sensors) + '.mat')['timers'])

    # Load true front inds
    for j in range(1,6):
        true_boundary.append(loadmat(folder_path + 'front' + str(j) + '.mat')['nodes_t'])

    rows = len(sensor_locs) + 1
    columns = 6
    figsizex = 10
    figsizey = 10
    fig,axes = plt.subplots(rows,columns,figsize=(figsizex,figsizey),constrained_layout=True)

    # Create triangulations
    triang_fwd = tri.Triangulation(p_fwd[0], p_fwd[1],t_fwd.T[:,:3]-1)
    triang_inv = tri.Triangulation(p_inv[0], p_inv[1],t_inv.T[:,:3]-1)

    # Plot truth
    global_min = np.min(np.log(K_true.ravel()))
    global_max = np.max(np.log(K_true.ravel()))
    [axes[0,i].axis('off') for i in range(5)]
    im1 = axes[0,5].tripcolor(triang_fwd, facecolors=np.log(K_true.ravel()),vmin=global_min,vmax=global_max,cmap='turbo', shading='flat')
    axes[0,5].set_xticks([])
    axes[0,5].set_yticks([])
    axes[0,5].set_xlim([0,1])
    axes[0,5].set_ylim([0,1])
    axes[0,5].set_title('$u^\dagger$')
    
    # Plot sensor locations
    [axes[i,0].scatter(sensor_locs[i-1][:,0],sensor_locs[i-1][:,1],s=6,c='black') for i in range(1,len(sensor_locs)+1)]
    [axes[i,0].scatter(sensor_locs[i-1][:,0],sensor_locs[i-1][:,1],s=3,c='white') for i in range(1,len(sensor_locs)+1)]
    [axes[i,0].set_xlim([0,1]) for i in range(1,len(sensor_locs)+1)]
    [axes[i,0].set_ylim([0,1]) for i in range(1,len(sensor_locs)+1)]
    [axes[i,0].set_xticks(x_ticks[CASE-1]) for i in range(1,len(sensor_locs)+1)]
    [axes[i,0].set_yticks(x_ticks[CASE-1]) for i in range(1,len(sensor_locs)+1)]
    axes[1,0].set_title('Sensors')

    # Plot means
    [axes[i,j+1].tripcolor(triang_inv, facecolors=means[i-1][:,j],alpha=1-vars[i-1][:,j]/0.25,
                           antialiased=True,linewidth=0.1,edgecolors='face',vmin=global_min,vmax=global_max,
                           cmap='turbo', shading='flat') for i in range(1,len(sensor_locs)+1) for j in range(5)]
    [axes[i,j+1].set_xticks([]) for i in range(1,len(sensor_locs)+1) for j in range(5)]
    [axes[i,j+1].set_yticks([]) for i in range(1,len(sensor_locs)+1) for j in range(5)]
    [axes[i,j+1].set_xlim([0,1]) for i in range(1,len(sensor_locs)+1) for j in range(5)]
    [axes[i,j+1].set_ylim([0,1]) for i in range(1,len(sensor_locs)+1) for j in range(5)]

    [axes[1, j+1].set_title(r'$\bar{{u}}_{{MAP}}^{{({})}}$'.format(j+1)) for j in range(5)]

    
    # Plot fronts
    for i in range(1,len(sensor_locs)+1):
        for j in range(5):
            bndry_locs_t = true_boundary[j]
            for k in range(len(bndry_locs_t)):
                axes[i,j+1].plot(
                    np.array([p_fwd.T[bndry_locs_t[k,0]-1,0], p_fwd.T[bndry_locs_t[k,1]-1,0]]),
                    np.array([p_fwd.T[bndry_locs_t[k,0]-1,1], p_fwd.T[bndry_locs_t[k,1]-1,1]]),
                    color='black'
                )


    fig.colorbar(im1, ax=axes[:,5],fraction=1, aspect=80)
    # plt.savefig(folder_path + 'sensor_density.jpg', format='jpg', bbox_inches='tight',dpi=400)
    plt.show()
    

    

   