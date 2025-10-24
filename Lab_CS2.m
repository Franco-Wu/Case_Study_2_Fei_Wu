%% Lab Case Study 2 Practice
% *ESE 105*
%
% *DUE on Sunday 10/26 11:59pm to Canvas
%
% 20 pts. total: 10 for programming and correct responses, 5 for
% programming style, 5 for presentation
%
% *Name: Minghui Wu (Franco) and Feiyu Ren*
%

%% Instructions
% Run the |Lab CS2.m| script. Follow the instructions inside the "TODO: *****"
% labels below to complete each part: either write code or write your
% response to the question. For example:

% TODO: *************************************************************
% (Replace this comment with code)
% (Add response here)
% *******************************************************************

%%%
%
% *To turn in your assignment:*
%
% * Run the command |publish('Lab_CS2.m','pdf')| in the _Command Window_
% to generate a PDF of your solution |Lab CS2.pdf|. 
% * Submit _both_ the code (|.m| file) and the published output (|.pdf|
% file) to Canvas.
%
clear;
close all;  % uncomment this line if you do not want all figure windows to close when running this code

%% Part 1: Step Equations

%%%
% Implement Equations (8) and (10) in MATLAB to simulate circuit A (Figure
% 1) with R = 1 kΩ and C = 1 µF. Simulate charging the capacitor using a
% step input: set V_C = 0 V at t = 0 and V_in = 1 V for t > 0. Choose a
% suitable h to model the charging process accurately. Plot V_in and V_C
% vs. time to show the charging of the capacitor. You should observe a
% charging curve similar to the Figure 2 in CS2. 
%%%
% Set up simulation parameters:

% TODO: *************************************************************
% Parameters:
R = 1e3; % Resistor 1 kΩ
C = 1e-6; % Capacitor 1 µF 
h = R*C/10000; % interval between updates, how frequent matlab updates the voltages
tEnd = 5e-3; % 5 ms, RC represens how fast the capacitor changes

t = 0:h:tEnd; % time interval from 0s to tEnd in steps of h
vin = ones(size(t)); % input voltage, constant 1V
vC = zeros(size(t)); % voltage across the capacitor
vR = zeros(size(t)); 
vC(1) = 0; % capactiro starts1 uncharged

for k = 1:length(t)-1
    vR(k) = vin(k) - vC(k); % (8)
    vC(k+1) = (1 - h/(R*C))*vC(k) + (h/(R*C))*vin(k); % (10)
end
vR(end) = vin(end) - vC(end); % fills in last resistor value

% Plot to match the example style
figure;
plot(t, vin,'LineWidth', 1.6); 
hold on;
plot(t, vC,'LineWidth', 1.6);
xlabel('time (s)'); ylabel('voltage (V)');
legend('V_{in}','V_C','Location','southeast');
title('Voltage measured across a capacitor in response to a constant 1 V input');
xlim([0 tEnd]); ylim([0 1]); % x an y axis limits
grid on; 
% *******************************************************************

%%%
% Simulate with small h: 
% (please finish the function 'simRCvoltages' in the end and use it here.)

% TODO: *************************************************************
% (Replace this comment with code)
% (Add response here)
% *******************************************************************

%%%
% Plot the figure:

% TODO: *************************************************************
% (Replace this comment with code)
% (Add response here)
% *******************************************************************


%% Part 2: Comparison between sampling intervals
% Now, run several versions of your simulation for various temporal sampling intervals h. As h gets
% larger or smaller, how does your simulation's prediction change? Why is this happening? Does
% the charging behavior of a "real" capacitor change as a function of your choice of h?
% Plot the predicted V_out using 
% 1) an "accurate" choice of h and 
% 2) an "inaccurate" choice of h and
% 3) the theoretical charging curve
% V_C(t) = 1 − exp(−t/RC) [Volts]. 
% Discuss what happens for the "inaccurate" choice of h.
% Note: Be careful to compare the three curves using correct time axes.

%%%
% Simulate with large h:

% TODO: *************************************************************
% (Replace this comment with code)
% (Add response here)
% *******************************************************************

%%%
% Compare charging curves:

% TODO: *************************************************************
% (Replace this comment with code)
% (Add response here)
% *******************************************************************


%% Part 3: Explanation about the sampling intervals
% Relative to the charging curve of the capacitor V_C (t), how do you interpret the meaning of τ =RC,
% called the RC time constant of the circuit?

% TODO: *************************************************************
% (Replace this comment with code)
% (Add response here)
% *******************************************************************

%% Helper functions "simRCvoltages": [V_C, V_R] = simRCvoltages(V_in,V_C0,R,C,h)

% TODO: *************************************************************
% (Replace this comment with code)
% (Add response here)
function [V_C, V_R] = simRCvoltages(V_in,V_C0,R,C,h)
end
% *******************************************************************
