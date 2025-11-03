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

Vin = Vin.';

% denoise 60 Hz
fHum = 60;
bwHum = 6;
C_hum = 0.47e-6;                                        
L_hum = 1 / ((2*pi*fHum)^2 * C_hum);                  
R_hum = 2*pi*L_hum*bwHum;

% A = [1,(h/C_hum);(-h/L_hum),(1-(R_hum*h/L_hum))];
% B = [0;(h/L_hum)];

alpha = 1+(h/L_hum)*R_hum;
A = [(1 - (h*h)/(alpha*L_hum*C_hum)),h/(C_hum*alpha); -h/(alpha*L_hum),1/alpha];
B = [(h*h)/(alpha*L_hum*C_hum); h/(alpha*L_hum)];

N = numel(Vin);
x = zeros(2, N);

for k = 1:N-1
    x(:,k+1) = A*x(:,k) + B*Vin(k);
end

i = x(2,:);
V_hum = i*R_hum;
Vin2 = Vin - V_hum;

% denoise 120 Hz
f2 = 120;
bw2 = 12;
C_2 = 0.47e-6;                                        
L_2 = 1 / ((2*pi*f2)^2 * C_2);                  
R_2 = 2*pi*L_2*bw2;

denom = 1+(h/L_2)*R_2;
A = [(1 - (h*h)/(denom*L_2*C_2)),h/(C_2*denom); -h/(denom*L_2),1/denom];
B = [(h*h)/(denom*L_2*C_2); h/(denom*L_2)];

N = numel(Vin2);
x = zeros(2, N);

for k = 1:N-1
    x(:,k+1) = A*x(:,k) + B*Vin2(k);
end

i = x(2,:);
V_120 = i*R_2;
Vin3 = Vin2 - V_120;

% denoise high frequency
fHigh = 666;
C_high = 0.22e-6;                                        
R_high = 1 / (2*pi*fHigh*C_high); 

V_C = zeros(1,N);

    for k = 1:N-1
        V_C(k+1) = (1-(h/(R_high*C_high)))*V_C(k) + (h/(R_high*C_high))*Vin3(k);
    end

Vout = V_C;
end