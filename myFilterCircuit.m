%% Case study 3: Circuits as Resonators, Sensors, and Filters
% *ESE 105* 
%
% *Name: Minghui Wu (Franco) and Feiyu Ren*
%
% function myFilterCircuit(Vin,h) receives a time-series voltage sequence
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

function Vout = myFilterCircuit(Vin,h)

N = numel(Vin); % number of V input discrete time samples
Vout = zeros(N,1); % empty output row vector, for future storate of output when given input

% filter out 60Hz
C1 = 1e-6; % F
f0_60 = 60; % Hz
L1 = 1 / ((2*pi*f0_60)^2 * C1); % H (resonates at approximately 60 Hz)
R1 = 8; % ohms (small resistance --> narrow notch)
V_C1 = 0; % capacitor voltage
i1 = 0; % inductor current
vR1 = zeros(N,1);% store resistor voltage
for k = 1:N
    V_C1 = V_C1 + (h/C1)*i1; % update capacitor voltage
    vR1(k) = R1 * i1;% voltage across R (band-pass)
    i1 = ( i1 + (h/L1)*(Vin(k) - V_C1) ) / (1 + (h/L1)*R1); % semi-implicit update
end
y1 = Vin - vR1; % subtract 60Hz noise

% filter out 120 Hz hum
C2 = 1e-6; 
f0_120 = 120;
L2 = 1 / ((2*pi*f0_120)^2 * C2);
R2 = 6; % ohms (small resistance --> narrow notch)

V_C2 = 0;
i2 = 0;
vR2 = zeros(N,1);

for k = 1:N
    V_C2 = V_C2 + (h/C2)*i2;
    vR2(k) = R2 * i2;
    i2 = ( i2 + (h/L2)*(y1(k) - V_C2) ) / (1 + (h/L2)*R2);
end

y2 = y1 - vR2; % subtract 120Hz harmonic noise of 60Hz


% Low-pass filter to remove high-frequency hiss
R3 = 50;% ohms
L3 = 0.0033; % henries
C3 = 3.3e-6;% farads (sets low-pass cutoff near few kHz)
V_C3 = 0;
i3 = 0;

for k = 1:N
    Vout(k) = V_C3;
    V_C3 = V_C3 + (h/C3)*i3; 
    i3 = ( i3 + (h/L3)*(y2(k) - V_C3) ) / (1 + (h/L3)*R3);
end

end