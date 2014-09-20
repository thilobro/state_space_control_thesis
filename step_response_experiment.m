%This script compares experimental data of a step response to the
%corresponding simulation data in order to judge how good the model fits
%the real system. The .nc files can be downloaded from the git respository
%https://github.com/thilobro/state_space_control_thesis
%ATTENTION: to run this script successfully, carousel_lagrange has to be
%changed to ddelta_motor_sp control

%read in experimental data
ncid1 = netcdf.open('siemensSensorsData6.nc');
timestamp1 = netcdf.getVar(ncid1, 0);
carouselspeed_sp = netcdf.getVar(ncid1, 6);

ncid2 = netcdf.open('lineAngleSensor2Data6.nc');
timestamp2 = netcdf.getVar(ncid2, 0);
azimuth = netcdf.getVar(ncid2, 1);
elevation = netcdf.getVar(ncid2, 2);

ncid3 = netcdf.open('armboneLisaSensorsData6.nc');
timestamp3 = netcdf.getVar(ncid3, 0);
qractualspeed = netcdf.getVar(ncid3, 3);

%resample data
qractualspeed_resamp = interp1(timestamp3, qractualspeed, timestamp1, 'spline');
elevation_resamp = interp1(timestamp2, elevation, timestamp1, 'spline');
azimuth_resamp = interp1(timestamp2, azimuth, timestamp1, 'spline');

%get simulation data
[T, X, U] = carousel_dynamics_SE(carouselspeed_sp, timestamp1);

%shift time axis in order to align simulation and experimental data
Tshift = 20;

%plot data
figure(4);
clf;
ax(1) = subplot(3, 1, 1);
plot(timestamp1, carouselspeed_sp, 'r')
hold on;
plot(timestamp1, qractualspeed_resamp, 'g')
plot(T + Tshift, X(6, :), 'b')
axis([60 100 1.4 1.8]);
xlabel('t [s]')
ylabel('Arm/Motor Speed [rad/s]')
legend('Motor Speed SP','Measured Arm Speed','Simulated Arm Speed')
ax(2) = subplot(3, 1, 2);
plot(T + Tshift, 180/pi*X(3, :), 'b');
hold on;
plot(timestamp1, 180/pi*elevation_resamp, 'g')
axis([60 100 -70 -30]);
legend('Simulated Elevation','Measured Elevation')
xlabel('t [s]')
ylabel('Elevation [deg]')
ax(3) = subplot(3, 1, 3);
plot(T + Tshift, 180/pi*X(4, :), 'b')
hold on;
plot(timestamp1, azimuth_resamp, 'g')
axis([60 100 -3 2]);
legend('Simulated Azimuth', 'Measured Azimuth')
xlabel('t [s]')
ylabel('Azimuth [deg]')

linkaxes(ax, 'x')