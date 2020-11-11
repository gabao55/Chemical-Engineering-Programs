##           Regressão Linear            ##
## Programado por: Gabriel Salateo Rosin ##

clear
clc
close

#Questão 1:
#Letra a:

x = [5;7.5;10;12.5;15;17.5;20;22.5;25;27.5;30];
y = [3.3;7.5;41.8;51.8;61;101.1;132.9;145.5;171.4;225.8;260.9];

n = length(x);

soma_x = sum(x);
soma_y = sum(y);

soma_x2 = 0;

for i = 1:length(x)
  soma_x2 = soma_x2 + (x(i)^2);
endfor

soma_xy = 0;
for i = 1:length(x)
  soma_xy = soma_xy + (x(i)*y(i));
endfor

a1 = (n*soma_xy - (soma_y*soma_x))/(n*soma_x2 - ((soma_x)^2));
a0 = (soma_y/n) - a1*(soma_x/n);

printf('\nLetra a:\n')
printf('A inclinacao encontrada vale %g e a interseccao com o eixo y vale %g\n\n\n', a1,a0)

#Letra b:

printf('Letra b:\n')

erros_absolutos = [];

for i = 1:n
  erro = abs(a0 + a1*x(i) - y(i));
  erros_absolutos = [erros_absolutos erro];
endfor

printf('\n erros absolutos :\n')

erros_absolutos'

erros_relativos = [];

for i = 1:n
  erro = abs(a0 + a1*x(i) - y(i))/y(i)*100;
  erros_relativos = [erros_relativos erro];
endfor

printf('\n erros percentuais :\n')

erros_relativos'

media_y = sum(y)/n;
desvios = 0;

for i = 1:n;
  desvios = desvios + ((y(i) - media_y)^2);
endfor

printf('\n')

DP = (desvios/(n-1))^0.5

printf('\n')

Sr = 0;

for i = 1:n
  Sr = Sr + ( y(i) - (a0 + a1*x(i)))^2
endfor

erro_padrao = (Sr/(n-2))^0.5

for i = 1:n
  reta(i) = a0 + a1*x(i);
endfor

reta'
plot(x,y,'*')
hold on
plot(x,reta,'r')

#Questão 2:
printf('\n Questão 2: \n')

t = (0:10:100)';
p =[ 0.94; 0.96; 1.0; 1.05; 1.07; 1.09; 1.14; 1.17; 1.21; 1.24; 1.28];

m = length(t);

soma_t = sum(t);
soma_p = sum(p);

soma_t2 = 0;

for i = 1:length(t)
  soma_t2 = soma_t2 + (t(i)^2);
endfor

soma_tp = 0;
for i = 1:length(t)
  soma_tp = soma_tp + (t(i)*p(i));
endfor

A1 = (n*soma_tp - (soma_p*soma_t))/(n*soma_t2 - ((soma_t)^2));
A0 = (soma_p/n) - A1*(soma_t/n);

printf('A inclinacao encontrada vale %g e a interseccao com o eixo y vale %g\n\n\n', A1,A0)

erros_absolutos = [];

for i = 1:n
  erro = abs(A0 + A1*t(i) - p(i));
  erros_absolutos = [erros_absolutos erro];
endfor

printf('\n erros absolutos :\n')

erros_absolutos'

erros_relativos = [];

for i = 1:n
  erro = abs(A0 + A1*t(i) - p(i))/p(i)*100;
  erros_relativos = [erros_relativos erro];
endfor

printf('\n erros percentuais :\n')

erros_relativos'

for i = 1:n
  reta(i) = A0 + A1*t(i);
endfor
