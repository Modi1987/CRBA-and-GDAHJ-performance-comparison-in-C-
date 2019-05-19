%% Copyright: Mohammad Safeea, 2019
% Following script is used to plot the results for calculating the inertia
% matrix of serially linked robots using GDAHJ and CRBA methods, results
% plotted here are generated from the file {CRBAvsGDAHJ_print_results.cpp}.

% Tests are carried out on a PC with: Intel(R) Core(TM) i7-8750H, CPU @ 2.20
% GHz, 8 GB RAM

% About software:
% {CRBAvsGDAHJ_print_results.cpp} was compiled with {Visual Studio 2019} on
% Windows 10 (64-bit) Operating System

% Timing results at the output of CRBAvsGDAHJ_print_results.cpp is in
% milliseconds for 500 000 invocations.

% Data format in {res} matrix is:
% [DOF, GDAHJ time, CRBA tima, error]
res=[2,878,116,1.38778e-16;
3,943,209,2.96059e-16;
4,1026,318,6.59195e-16;
5,1127,442,7.36633e-16;
6,1231,593,9.11308e-16;
7,1346,756,1.1139e-15;
8,1461,937,1.52829e-15;
9,1580,1144,1.26177e-15;
10,1694,1357,1.61537e-15;
11,1828,1591,1.39707e-15;
12,1946,1852,1.46758e-15;
13,2082,2111,1.85322e-15;
14,2225,2403,2.59154e-15;
15,2368,2714,3.4062e-15;
16,2512,3053,3.82441e-15;
17,2656,3390,4.43393e-15;
18,2821,3747,4.23161e-15;
19,2982,4123,4.73121e-15;
20,3143,4516,5.50952e-15;
21,3307,4950,5.77033e-15;
22,3486,5381,6.60497e-15;
23,3709,5832,6.94136e-15;
24,3836,6298,7.31441e-15;
25,4026,6780,9.09308e-15;
26,4204,7296,9.63106e-15;
27,4404,7806,1.02536e-14;
28,4589,8357,1.1373e-14;
29,4831,8926,1.21291e-14;
30,4996,9501,1.47543e-14];


plot(res(:,1),res(:,2)/500000);
hold on;
plot(res(:,1),res(:,3)/500000);
xlabel('DOF');
ylabel('Time (milliseconds)')
legend('GDAHJ','CRBA')
title('Time comaprison for CRBA and GDAHJ')

figure
plot(res(:,1),res(:,4))
xlabel('DOF');
ylabel('Error (Kg.m^2)')
title('Error in calculation between CRBA and GDAHJ')