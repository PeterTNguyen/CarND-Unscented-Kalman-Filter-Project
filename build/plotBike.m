gt = dlmread("GroundTruth.csv");
est = dlmread("Estimates.csv");
lidar = dlmread("lidar_data.csv");
radar = dlmread("radar_data.csv");
nis_lidar = dlmread("NIS_lidar.csv");
nis_radar = dlmread("NIS_radar.csv");

figure(1); clf;
plot(gt(:,1),gt(:,2),'b','linewidth', 2);
hold on;grid on;
plot(est(:,1), est(:,2),'-*r');
plot(lidar(:,1), lidar(:,2),'k','marker','d','markersize',10,'linestyle','none')
plot(radar(:,1).*cos(radar(:,2)), radar(:,1).*sin(radar(:,2)),'m','marker','o','markersize',10,'linestyle','none')


figure(2);clf; hold on;grid on;
nis_threshold = 7.815;
plot(nis_radar, 'r-*');
plot(nis_lidar, 'b-o');
plot(ones(size(nis_radar))*nis_threshold);
legend('NIS Radar', 'NIS Lidar');