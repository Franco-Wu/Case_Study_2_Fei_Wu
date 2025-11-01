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

R = 100; % Resistance of 100 ohms
L = 100e-3; % Inductance of 100mH
C = 0.1e-6; % Capacitance of 0.11microF

% --- Simulation setup ---
N = numel(Vin); % number of V input discrete time samples
Vout = zeros(N,1); % empty output row vector, for future storate of output when given input

% before applying any input voltage, everything is at "rest", nothing
% happens, so: 
V_C = 0; % Capacitor initially uncharged
i = 0; % current initially 0

Vout(1) = Vin(1); % Initial output voltage is the same as input

for k = 1:N
    V_C = V_C + (h/C)*i; %eq 5, updates capacitor voltage, more current --> more charge flows --> higher voltage    
    % Eq. (7):  vR = R*i   → output voltage

    Vout(k) = R * i; %eq7, V_R = iR, output voltage = voltage across the resistor  
    % output of resonator circuit, the one that is to be listened and plotted

    i=(i+(h/L)*(Vin(k)-V_C) ) / (1 + (h/L)*R); % eq(16)+KVL with implicit R term for i_{k+1}
end

end 