rng(0)
n = 500; % Number of samples
A3 = [1 .5];
C3 = [1, -0.3, 0.2];
w = sqrt(2) * randn(n + 100, 1);
x = filter(C3, A3, w); % Create the input
A1 = [1, -0.65];
A2 = [1, 0.90, 0.78];
C = 1;
B = [0, 0, 0, 0, 0.4];
e = sqrt(1.5) * randn(n + 100, 1);
y = filter(C, A1, e) + filter(B, A2, x); % Create the output
x = x(101:end); 
y = y(101:end); % Omit initial samples
clear A1 A2 C B e w A3 C3
