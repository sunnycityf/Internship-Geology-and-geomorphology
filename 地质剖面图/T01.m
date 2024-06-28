clc
clear
File = uigetfile('*.csv', '选择一个表格文件', 'MultiSelect', 'off');
Data = readmatrix(File);
Height = Data(:,5);
Distance = Data(:,6);
Lon = Data(:,3);
Lat = Data(:,4);
Title = inputdlg('请为剖面图起名：');
StartLon = Lon(1,:);
StartLat = Lat(1,:);
LengthLon = length(Lon);
LengthLat = length(Lat);
EndLon = Lon(LengthLon,:);
EndLat = Lat(LengthLat,:);
StartCoordinate = [StartLon,StartLat];
EndCoordinate = [EndLon,EndLat];
DistanceE = max(Distance);
SubTitle = sprintf('起始坐标：%f,%f， 终点坐标：%f,%f， 距离：%f 千米', StartCoordinate(1), StartCoordinate(2), EndCoordinate(1), EndCoordinate(2), DistanceE);
fig = figure('Position', [100, 100, 1000, 500]);
figure(1)
plot(Distance,Height,'-','LineWidth',1,'Color',[0,0,0])
title(Title,'FontWeight','bold')
subtitle(SubTitle, 'FontSize', 10, 'FontWeight', 'normal')
xlim([0, max(Distance)]);
ylim([(min(Height)-5),(max(Height))+5])
xlabel("长度（千米）")
ylabel("海拔（米）")
