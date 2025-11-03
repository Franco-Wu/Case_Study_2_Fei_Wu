%% Case study 3: Circuits as Resonators, Sensors, and Filters
% *ESE 105* 
%
% *Name: Minghui Wu (Franco) and Feiyu Ren*
%
% function mySensorCircuit(Vin,h) receives a time-series voltage sequence
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


function Vout = mySensorCircuit(Vin,h)

f0 = 84;      
bandWidthHz = 12;

C = 0.47e-6;      
L = 1 / ((2*pi*f0)^2 * C);  % set the L accroding to the resonance equation
R = 2*pi*L*bandWidthHz;  % derive R according to the bandwidth

% A = [1,(h/C);(-h/L),(1-(R*h/L))];
% B = [0;(h/L)];

alpha = 1+(h/L)*R;
A = [(1 - (h*h)/(alpha*L*C)),h/(C*alpha); -h/(alpha*L),1/alpha];
B = [(h*h)/(alpha*L*C); h/(alpha*L)];

N = numel(Vin);
x = zeros(2, N);

for k = 1:N-1
    x(:,k+1) = A*x(:,k) + B*Vin(k);
end

i = x(2,:);
Vout = i*R;

end