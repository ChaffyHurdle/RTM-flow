function p = plot_push_forward(pressures, flow_fronts, physics_class, mesh_class, ...
    true_RTMflow)


figure(7)
subplot(2,5,[1,2,3,4,5])

[nObs, ~] = size(pressures);
M = physics_class.nobservations;
N = physics_class.nsensors;

% Flatten to vector for plotting
data_vector = pressures(:);

% Group IDs: 1 to M repeated N times for each column group
group_ids = repelem(1:N, M);
group_ids = repmat(group_ids, nObs, 1);
group_ids = group_ids(:);

% Within-group position (1–N)
within_group_pos = repmat(1:M, 1, N);
within_group_pos = repmat(within_group_pos, nObs, 1);
within_group_pos = within_group_pos(:);

% Plot with boxchart
boxchart(categorical(group_ids), data_vector, 'GroupByColor', within_group_pos);
hold on
scatter(0.6:0.2:25.4,reshape(transpose(true_RTMflow.pressure_data),[],1),'r')
hold off
xlabel('Sensor');
ylabel('Pressure');
legend('t_1','t_2','t_3','t_4','t_5');


for i=1:5
    t_index = find(true_RTMflow.times > physics_class.observation_times(i),1)-1;
    subplot(2,5,5+i)
    plot_front(true_RTMflow,t_index);
    hold on
    for j = 1:size(pressures,1)
        for k = 1:size(flow_fronts{j,i},1)
            plot([mesh_class.nodes(flow_fronts{j,i}(k,2),1), ...
                  mesh_class.nodes(flow_fronts{j,i}(k,3),1)],...
                 [mesh_class.nodes(flow_fronts{j,i}(k,2),2), ...
                  mesh_class.nodes(flow_fronts{j,i}(k,3),2)], ...
                 'Color', [0, 0, 0, 0.1]);
        end
    end
    hold off
    axis square;
    xlim([0,1])
    ylim([0,1])
    title('$u^{\dagger}$','interpreter','latex')
end


end