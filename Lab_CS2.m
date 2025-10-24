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

R = 1e3; % Resistor of 1 kilo ohms
C = 1e-6; % capacitor 1 miu farad (coloumb per volt)
h = 0.00001; % interval between updates, how frequent matlab updates the voltage
V_C0 = 0; % capacitor starts uncharged
t_end = 5e-3;
t = 0:h:t_end;
V_in = ones(1,numel(t));

% *******************************************************************

%%%
% Simulate with small h: 
% (please finish the function 'simRCvoltages' in the end and use it here.)
% TODO: *************************************************************
% (Replace this comment with code)
[vC, vR] = simRCvoltages(V_in,V_C0,R,C,h);
% *******************************************************************

%%%
% Plot the figure:

% TODO: *************************************************************
figure
plot(t, vC,'LineWidth', 1.6)
hold on
plot(t,V_in,'LineWidth', 1.6)
xlabel('time (s)')
ylabel('voltage (V)')
xlim([0 t_end])
ylim([0 1.05])
title('Voltage measured across a capacitor in response to a constant 1 V input');
legend('V_{in}','V_C','Location','best')
grid on
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
h = 0.0006;
t_end = 5e-3;
t = 0:h:t_end;
V_in = ones(1,numel(t));
[vC, vR] = simRCvoltages(V_in,V_C0,R,C,h);

figure
plot(t, vC,'LineWidth', 1.6)
hold on
plot(t,V_in,'LineWidth', 1.6)
xlabel('time (s)')
ylabel('voltage (V)')
xlim([0 t_end])
ylim([0 1.05])
title('Voltage measured across a capacitor in response to a constant 1 V input');
legend('V_{in}','V_C','Location','best')
grid on
% How does your simulation's prediction change? 
% If the value of h is large, Matlab seems to be updating the voltage at a
% greater time interval, hence, causing the curve to exhibit relatively
% sharp edges. If the value of h is small, Matlab updates the value of
% voltage at a smaller increment of time, hence, storing more data points,
% allowing the curve to be smoother without sharp edges. 

% Does the charging behavior of a "real" capacitor change as a function of your choice of h?
% No, the charging behavior of a "real" capactiro does not change as a
% function of my decision of h. This is becuse the "real" capacitor's charing
% behavior is only dependent on its Resistance and Capacitance, not on h.
% Changing h only affects how accurately our MATLAB siimulation predicts
% the behavior, it does not change how the physical circuit charges. 
% *******************************************************************

%%%
% Compare charging curves:
% TODO: *************************************************************
h1 = 0.00001; % accurate choice of h
t_end_1 = 5e-3;
t1 = 0:h1:t_end_1;
V_in1 = ones(1,numel(t1));
[vC_1, vR_1] = simRCvoltages(V_in1,V_C0,R,C,h1);

h2 = 0.0006; % inaccurate choice of h
t_end_2 = 5e-3;
t2 = 0:h2:t_end_2;
V_in2 = ones(1,numel(t2));
[vC_2, vR_2] = simRCvoltages(V_in2,V_C0,R,C,h2);

figure
subplot(3,1,1)
plot(t1, vC_1,'LineWidth', 1.6)
hold on
plot(t1,V_in1,'LineWidth', 1.6)
xlabel('time (s)')
ylabel('voltage (V)')
xlim([0 t_end_1])
ylim([0 1.05])
title('Accurate value of h');
legend('V_{in}','V_C','Location','best')
grid on

subplot(3,1,2)
plot(t, vC_2,'LineWidth', 1.6)
hold on
plot(t,V_in2,'LineWidth', 1.6)
xlabel('time (s)')
ylabel('voltage (V)')
xlim([0 t_end_2])
ylim([0 1.05])
title('inaccurate value of h');
legend('V_{in}','V_C','Location','best')
grid on

subplot(3,1,3)
% Compute the continuous theoretical charging curve
t_theory = linspace(0, 5e-3, 1000);          % smooth time vector
V_C_theory = 1 - exp(-t_theory/(R*C));       % theoretical RC charging formula

plot(t_theory, V_C_theory, 'k', 'LineWidth', 1.8); hold on;   % theory (smooth)
xlabel('time (s)')
ylabel('voltage (V)')
xlim([0 5e-3])
ylim([0 1.05])
title('Theoretical Charing Curve')
legend('Theoretical 1 - e^{-t/RC}', 'Accurate h = 0.00001', 'Inaccurate h = 0.0006', 'Location','southeast')
grid on
% inaccurate choice of h will make the curve have sharp edges
% *******************************************************************


%% Part 3: Explanation about the sampling intervals
% Relative to the charging curve of the capacitor V_C (t), how do you interpret the meaning of τ =RC,
% called the RC time constant of the circuit?

% TODO: *************************************************************
% The time constant τ = R*C sets how quickly the capacitor voltage changes.
% It is the time for V_C(t) to reach about 63.2% of its final value after a step input.
% The charging follows V_C(t) = V_in*(1 - exp(-t/(R*C))).
% *******************************************************************

%% Helper functions "simRCvoltages": [V_C, V_R] = simRCvoltages(V_in,V_C0,R,C,h)

% TODO: *************************************************************
% (Replace this comment with code)
% (Add response here)
function [V_C, V_R] = simRCvoltages(V_in,V_C0,R,C,h)
t_end = 5e-3;
t = 0:h:t_end;
V_C = zeros(1,numel(t));
V_C(1) = V_C0;
V_R = zeros(1,numel(t));

    for i = 1:numel(t)-1
        V_R(i) = V_in(i)-V_C(i);
        V_C(i+1) = (1-(h/(R*C)))*V_C(i) + (h/(R*C))*V_in(i);
    end
    V_R(end) = V_in(end) - V_C(end);% fills in last resistor value

end
% *******************************************************************
