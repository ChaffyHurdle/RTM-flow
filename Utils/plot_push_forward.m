function p = plot_push_forward(pressures, flow_fronts, physics_class, mesh_class, ...
    true_RTMflow)


sensor_locs_inds = [1,7,13,19,25];
figure(9)

%% Plot pressures

for i=1:5
    subplot(2,5,i)
    s = sensor_locs_inds(i);
    pressures_sensor_i = pressures(:,[5*(s-1)+1:5*s]);
    b = boxchart(pressures_sensor_i,'MarkerStyle','none');
    b.BoxFaceColor = [0 0.5 0.5];
    hold on
    scatter([1,2,3,4,5], true_RTMflow.pressure_data(s,:),'ro','filled')
    hold off
    %ylim([0.95,2])
    title(sprintf('$p(%.1f, %.1f)+\\eta$', physics_class.sensor_locs(s,1), physics_class.sensor_locs(s,2)), 'Interpreter', 'latex')
end


%% Plot front locations
for i=1:5
    t_index = find(true_RTMflow.times > physics_class.observation_times(i),1)-1;
    subplot(2,5,5+i)

    front_count = zeros(1,mesh_class.num_nodes);
    for j=1:1000
        front_nodes = unique(reshape(flow_fronts{j,i}(:,2:3),[],1));
        front_count(front_nodes) = front_count(front_nodes) + 1;
    end

    
    pdeplot(mesh_class.nodes',mesh_class.elements', ...
        XYData = front_count, XYStyle='interp', ...
        ColorMap="parula",Mesh="off",ColorBar="off")
    xlim([0,1])
    ylim([0,1])
    clim([0,1000])

    hold on
    plot_front(true_RTMflow,t_index);
    hold off
    title(sprintf('$\\Upsilon(t_{%d})$', i), 'Interpreter', 'latex')
end
ax = subplot(2,5,10);
colormap(ax,'parula');
pos = get(ax,'Position');
h = colorbar('Position', [pos(1)+(1+1/20)*pos(3)  pos(2)  pos(3)/10  pos(4)]);



end