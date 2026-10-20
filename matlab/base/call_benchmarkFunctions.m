function [ y, lb, ub ,x, M  ] = call_benchmarkFunctions( ID )
% Details: Utility function that finds the correct benchmark function from
% input prompt
%
% inputs:
% ID - Input prompt call ID
%
% outputs:
% y - Sample outputs
% lb - Lower bound for sampling
% ub - Higher bound for sampling
% x - Sample inputs
% M - Matlab function to call

%%% 1D Functions
% Single Hump function
if ID ==1
    [ y, lb, ub, x, M ] = Single_hump_1d();
elseif ID == 2
    % Modified Two humps function
    [ y, lb, ub,x, M ] = Two_humps_function_1d();
    % Gramacy-Lee  https://www.sfu.ca/~ssurjano/grlee12.html
elseif ID ==3
    [ y, lb, ub, x, M ] = Gramacy_Lee_1d();
elseif ID == 4
    % Modified Gramacy-Lee
    [ y, lb, ub, x, M ] = Mod_Gramacy_Lee_1d();
elseif ID == 5
    % Perm function 1,0   https://www.sfu.ca/~ssurjano/permdb.html
    [ y, lb, ub, x, M ] = Perm0db_1d();
elseif ID == 6
    % Damped Cosinus function   https://www.sfu.ca/~ssurjano/santetal03dc.html
    [ y, lb, ub, x, M ] = DampedCos_1d();
elseif ID == 7
    % Exploit function
    [ y, lb, ub, x, M ] = Exploit_1d();
elseif ID == 8
    % Forrester function
    [ y, lb, ub, x, M ] = MultiPeak_1d();
    
    
    
    %%% 2D Functions
    %  Michalewicz function https://www.sfu.ca/~ssurjano/michal.html
elseif ID == 9
    [ y, lb, ub, x, M ] = Micha_2d();
    % Drop-Wave  https://www.sfu.ca/~ssurjano/drop.html
elseif ID == 10
    [ y, lb, ub, x, M ] = Drop_wave_2d();
    % Booth function  https://www.sfu.ca/~ssurjano/booth.html
elseif ID == 11
    [ y, lb, ub, x, M ] = Booth_2d();
    % Bohachevsky function https://www.sfu.ca/~ssurjano/boha.html
elseif ID == 12
    [ y, lb, ub, x, M ] = Boha_2d();
    % Branin function https://www.sfu.ca/~ssurjano/branin.html
elseif ID == 13
    [ y, lb, ub, x, M ] = Branin_2d();
    % Franke function https://www.sfu.ca/~ssurjano/franke2d.html
elseif ID == 14
    [ y, lb, ub, x, M ] = Franke_2d();
    % Rosenbrock function https://www.sfu.ca/~ssurjano/rosen.html
elseif ID == 15
    [ y, lb, ub, x, M ] = Rosenbrock_2d();
    % Six Hump Camel function https://www.sfu.ca/~ssurjano/camel6.html
elseif ID == 16
    [ y, lb, ub, x, M ] = SHCamel_2d();
    % Rastrigin function https://www.sfu.ca/~ssurjano/rastr.html
elseif ID == 17
    [ y, lb, ub, x, M ] = Rastrigin_2d();
    % Griewank function https://www.sfu.ca/~ssurjano/griewank.html
elseif ID == 18
    [ y, lb, ub, x, M ] = Griewank_2d();

    % Peaks function
elseif ID == 19
    [ y, lb, ub, x, M ] = Peaks_2d();
elseif ID == 20
    [ y, lb, ub, x, M ] = GoldsteinPrice_2d();

    
    
    %%% 3D Functions
    % Sphere 3D https://www.sfu.ca/~ssurjano/spheref.html
elseif ID == 21
    [ y, lb, ub, x, M ] = Sphere_3d();

elseif ID == 22
    [ y, lb, ub, x, M ] = Hartmann_3d();

elseif ID == 23
    [ y, lb, ub, x, M ] = Ishigami_3d();
elseif ID == 24
    [ y, lb, ub, x, M ] = Ackley_3d();

    %%% 4D Functions
elseif ID == 25
    [ y, lb, ub, x, M ] = Sphere_4d();

elseif ID == 26
    [ y, lb, ub, x, M ] = DixonP_4d();

    %%% 5D Functions
elseif ID == 27
    [ y, lb, ub, x, M ] = Sphere_5d();

    %%% 6D Functions
elseif ID == 28
    [ y, lb, ub ,x, M ] = Hartmann_6d();
else
    error('Error. \nBenchmark ID unknown');
end


end

