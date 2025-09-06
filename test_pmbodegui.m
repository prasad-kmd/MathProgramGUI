% TEST_PMBODEGUI - Test script for the pmbodegui function.
%
% Description:
%   This script creates a sample transfer function and passes it to the
%   pmbodegui function to generate the Bode plot GUI.

% Clear workspace and close all figures
clear;
close all;
clc;

% Define the transfer function
% G(s) = (10s^2 + 20s) / (s^3 + 11s^2 + 110s + 1000)
num = [10 20 0];
den = [1 11 110 1000];
sys = tf(num, den);

% Call the GUI function
pmbodegui(sys);

% Instructions for the user:
% 1. Save both pmbodegui.m and test_pmbodegui.m in the same directory.
% 2. Open MATLAB and navigate to that directory.
% 3. Run the command 'test_pmbodegui' in the MATLAB command window.
% 4. A GUI window should appear with the Bode plots for the sample
%    transfer function.
