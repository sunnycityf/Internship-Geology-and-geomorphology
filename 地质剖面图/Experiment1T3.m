clc
clear
% 程序交互
lon1 = 103.7618173;
lat1 = 31.189372;
lon2 = input('请输入第二个点的经度(度): ');
lat2 = input('请输入第二个点的纬度(度): ');
% 显示输入的经纬度并确认
disp("---------------------------------------------")
disp('您输入的经纬度为:');
disp(['第一个点经度: ', num2str(lon1)]);
disp(['第一个点纬度: ', num2str(lat1)]);
disp(['第二个点经度: ', num2str(lon2)]);
disp(['第二个点纬度: ', num2str(lat2)]);
% 核心部分
% 将经纬度转换为弧度
lon1 = deg2rad(lon1);
lat1 = deg2rad(lat1);
lon2 = deg2rad(lon2);
lat2 = deg2rad(lat2);
% 统计经纬度
geo1 = [lon1; lat1];
geo2 = [lon2; lat2];
% 调用函数计算距离(Haversine)
Distance_H = H(geo1, geo2);
% 调用函数(Vincenty)
Distance_V = V(geo1,geo2);
% 调用函数(球极坐标法)
Distance_SC = SC(geo1,geo2);
% 显示结果
disp("---------------------------------------------")
fprintf('两点之间的距离为(以Haversine公式计算) %.2f 千米\n', Distance_H);
fprintf('两点之间的距离为(以Vincenty算法计算) %.2f 千米\n', Distance_V);
fprintf('两点之间的距离为(以球极坐标法计算) %.2f 千米\n', Distance_SC);
% 创建函数(Haversine公式）
function Haversine = H(geo1,geo2)
    lon_h1 = geo1(1);
    lat_h1 = geo1(2);
    lon_h2 = geo2(1);
    lat_h2 = geo2(2);

    R = 6371.01;
    delta_lon = lon_h2 - lon_h1;
    delta_lat = lat_h2 - lat_h1;
    a = sin(delta_lat/2)^2 + cos(lat_h1) * cos(lat_h2) * sin(delta_lon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));
    Haversine = R * c;
end

% 创建函数(Vincenty算法)
function Vincenty = V(geo1,geo2)
    lon_v1 = geo1(1);
    lat_v1 = geo1(2);
    lon_v2 = geo2(1);
    lat_v2 = geo2(2);

    RL = 6378.1370000;
    RS = 6356.7253142;
    F = 1 / 298.257223563;

    L = lon_v2 - lon_v1;

    U1 = atan((1 - F) * tan(lat_v1));
    U2 = atan((1 - F) * tan(lat_v2));
    sinU1 = sin(U1);
    cosU1 = cos(U1);
    sinU2 = sin(U2);
    cosU2 = cos(U2);

    lambda = L;
    lambdaP = 2 * pi;

    iterLimit = 100;
    while abs(lambda - lambdaP) > 1e-12 && iterLimit > 0
        sinLambda = sin(lambda);
        cosLambda = cos(lambda);
        sinSigma = sqrt((cosU2 * sinLambda) ^ 2 + (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) ^ 2);
        cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
        sigma = atan2(sinSigma, cosSigma);
        sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
        cosSqAlpha = 1 - sinAlpha ^ 2;
        cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;

        C = F / 16 * cosSqAlpha * (4 + F * (4 - 3 * cosSqAlpha));
        lambdaP = lambda;
        lambda = L + (1 - C) * F * sinAlpha * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM ^ 2)));

        iterLimit = iterLimit - 1;
    end

    uSq = cosSqAlpha * (RL ^ 2 - RS ^ 2) / (RS ^ 2);
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
    deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 * (cosSigma * (-1 + 2 * cos2SigmaM ^ 2) - B / 6 * cos2SigmaM * (-3 + 4 * sinSigma ^ 2) * (-3 + 4 * cos2SigmaM ^ 2)));

    Vincenty = RS * A * (sigma - deltaSigma);
end

%创建函数(球极坐标法)
function Spherical_Coordinates = SC(geo1,geo2)

    R = 6371.01;

    lon_sc1 = geo1(1);
    lat_sc1 = geo1(2);
    lon_sc2 = geo2(1);
    lat_sc2 = geo2(2);
    
    r1 = [cos(lon_sc1)*cos(lat_sc1);sin(lon_sc1)*cos(lat_sc1);sin(lat_sc1)];
    r2 = [cos(lon_sc2)*cos(lat_sc2);sin(lon_sc2)*cos(lat_sc2);sin(lat_sc2)];
    T = acos(dot(r1,r2));
    Spherical_Coordinates = T*R;
end