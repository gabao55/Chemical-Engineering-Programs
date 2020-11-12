##         Numerical Integration          ##
##  Programado por Gabriel Salateo Rosin  ##


clear
clc
close

## Choose de function you want to integrate:
function f = f(x)
  f = (sin((pi*x)/2)^3.2);
endfunction

## Choose the integration interval:
t = (0:0.1:120);

## Plotting the function you are integrating:
for i = 1:length(t)
  x = t(i);
  F(i) = f(x);
endfor

figure(1)
plot(t,F,'g')
xlabel('x')
ylabel('f(x)')

## Integration algorithm:
h = 2/298;
X = (0:h:2);
n = 2/h;
integral = 0;

for i = 1:2:n
  x0 = X(i);
  x1 = X(i+1);
  x2 = X(i+2);
  f0 = f(x0);
  f1 = f(x1);
  f2 = f(x2);
  I = (h/3)*(f0 + 4*f1 + f2);
  integral = integral + I;
endfor

## Output
integral
