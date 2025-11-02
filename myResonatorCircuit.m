%% Case study 3: Circuits as Resonators, Sensors, and Filters
% *ESE 105* 
%
% *Name: Minghui Wu (Franco) and Feiyu Ren*
%
% function myResonatorCircuit(Vin,h) receives a time-series voltage sequence
% sampled with interval h, and returns the output voltage sequence produced
% by a circuit
%
% inputs:
% Vin - time-series vector representing the voltage input to a circuit
% h - scalar representing the sampling interval of the time series in
% seconds
%
% outputs:
% Vout - time-series vector representing the output voltage of a circuit

%{
A tuning fork. Choose your favorite musical tone from the list of piano key frequencies. For
example, symphony orchestras tune to A4 (440 Hz). Tune your RLC oscillator to “ring” at your
favorite frequency after receiving a short voltage/current pulse. Report the circuit component
values you used in your design.
Implement your circuit resonator as the function myResonatorCircuit(Vsound,h) stored in
myResonatorCircuit.m.
Describe the procedure you followed to tune your circuit.
Demonstrate the ringing of your circuit. Comment on the ability of your circuit to function as a
tuning fork. How strong of an input pulse (i.e., what voltage, what duration) is needed to make
it ring for ~5 seconds? 
%}

function Vout = myResonatorCircuit(Vin,h)

% ----- Design choices (tuning-fork) -----
C = 1e-6; % in Farads, Capacitance of 1.0 microF
f0 = 440; %Hz
L = 1/((2*pi*f0)^2*C); % Henries  -> 0.1006 H
R = 0.05; %ohm 

N = numel(Vin); % number of V input discrete time samples
Vout = zeros(N,1); % empty output row vector, for future storate of output when given input

% before applying any input voltage, everything is at "rest", nothing
% happens, so: 
V_C = 0; % Capacitor initially uncharged
i = 0; % current initially 0


% Matrix A vector B 
denom = 1+(h/L)*R;
A = [(1 - (h*h)/(denom*L*C)),  h/(C*denom); 
     -h/(denom*L),              1/denom];

B = [(h*h)/(denom*L*C); h/denom*C];

%{ 
for k = 1:N
    V_C = V_C + (h/C)*i; %eq 5, updates capacitor voltage, more current --> more charge flows --> higher voltage    
    % Eq. (7):  vR = R*i   → output voltage

    Vout(k) = R * i; %eq7, V_R = iR, output voltage = voltage across the resistor  
    % output of resonator circuit, the one that is to be listened and plotted

    i=(i+(h/L)*(Vin(k)-V_C) ) / (1 + (h/L)*R); % eq(16)+KVL with implicit R term for i_{k+1}
end
%}

x = [0; 0]; % x = [V_C; i]
    for k = 1:N
        Vout(k) = R * x(2);%eq7, V_R = iR, output voltage = voltage across the resistor  
        % output of resonator circuit, the one that is to be listened and plotted
        x = A*x + B*Vin(k); % xk+1 = Axk + Buk
    end

% plot output voltage vs time
Fs = 1/h; % sampling frequency
t = (0:N-1)'/Fs;% time vector
figure;
plot(t, Vout, 'LineWidth', 1.4)
xlabel('Time (s)'), ylabel('Output Voltage (V)')
title('RLC Resonator Ringing Response')
grid on
xlim([0, 5]) 

end




