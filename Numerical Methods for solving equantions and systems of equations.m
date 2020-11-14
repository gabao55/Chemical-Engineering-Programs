## Comparação entre Métodos Numéricos para Resolução de Equações ##
## Newton-Raphson, Secante e Newton-Raphson para Sistemas não Lineares ##
##  Programado por Gabriel Salateo Rosin  ##

clc
clear

function saturacao = saturacao(Ta,objetivo)
  Ta = Ta + 273.15;
  saturacao = objetivo - exp( - 139.34411 + (1.575701*(10^5)/Ta) - (6.642308*(10^7)/(Ta^2)) + (1.2438*(10^10)/(Ta^3)) - (8.621949*(10^11)/(Ta^4)));
endfunction

function derivada_saturacao = derivada_saturacao(Ta)
  Ta = Ta + 273.15;
  derivada_saturacao = - (2*6.642308*(10^7)/(Ta^3) + 4*8.621949*(10^11)/(Ta^5) - 1.575701*(10^5)/(Ta^2) - 3*1.2438*(10^10)/(Ta^4))*exp( - 139.34411 + (1.575701*(10^5)/Ta) - (6.642308*(10^7)/(Ta^2)) + (1.2438*(10^10)/(Ta^3)) - (8.621949*(10^11)/(Ta^4)));
endfunction

# Questão 1:
#Método de Newton-Raphson:
for objetivo = 8:2:12
  x0 = 0;
  i = 0;
  imax = 100;
  tol = 0.05;
  while i < imax
    xi = x0 - saturacao(x0,objetivo)/derivada_saturacao(x0);
    i = i + 1;
    xv = x0;
    erro = abs((xi - xv)/(xi));
    if erro < tol
      break
    endif
    x0 = xi;
  endwhile
  iteracao_nh = i;
  printf('\n\nO programa chegou, para uma concentracao de %u mg/L,\na um resultado de temperatura igual a %g °C pelo METODO DE NEWTON-RAPHSON, \ncom um total de %u iteracoes e um erro igual a %f.\n', objetivo, xi, iteracao_nh, erro)
endfor

# Questão 2:
F = 3;
Ca0 = 200;
Cb0 = 200;
V = 40;
k = 0.0045863;
tol = 0.05;
x0 = 200;
y0 = 200;

function f = f(ca,cb,F,Ca0,V,k)
  f = F*Ca0 - F*ca-k*ca*cb*V;
endfunction

function dfca = dfca(cb,F,V,k)
  dfca = -F -k*cb*V;
endfunction

function dfcb = dfcb(ca,V,k)
  dfcb = -k*ca*V;
endfunction

function g = g(ca,cb,F,Cb0,V,k)
  g = F*Cb0 - F*cb-k*ca*cb*V;
endfunction

function dgca = dgca(cb,V,k)
  dgca = -k*cb*V;
endfunction

function dgcb = dgcb(ca,F,V,k)
  dgcb =  -F-k*ca*V;
endfunction

while i < imax
  xi = x0 - (((f(x0,y0,F,Ca0,V,k)*dgcb(x0,F,V,k) - g(x0,y0,F,Cb0,V,k)*dfcb(x0,V,k)))/(dfca(y0,F,V,k)*dgcb(x0,F,V,k) - dfcb(x0,V,k)*dgca(y0,V,k)));
  yi = y0 - (((g(x0,y0,F,Cb0,V,k))*dfca(y0,F,V,k)-f(x0,y0,F,Ca0,V,k)*dgca(y0,V,k))/(dfca(y0,F,V,k)*dgcb(x0,F,V,k)- dfcb(x0,V,k)*dgca(y0,V,k)));
  errox = abs((((f(x0,y0,F,Ca0,V,k)*dgcb(x0,F,V,k) - g(x0,y0,F,Cb0,V,k)*dfcb(x0,V,k)))/(dfca(y0,F,V,k)*dgcb(x0,F,V,k) - dfcb(x0,V,k)*dgca(y0,V,k))));
  erroy = abs((((g(x0,y0,F,Cb0,V,k))*dfca(y0,F,V,k)-f(x0,y0,F,Ca0,V,k)*dgca(y0,V,k))/(dfca(y0,F,V,k)*dgcb(x0,F,V,k)- dfcb(x0,V,k)*dgca(y0,V,k))));
  i = i+ 1;
  if errox < tol && erroy < tol
    break
  endif
  x0 = xi;
  y0 = yi;
endwhile
printf('\n\nO METODO DE NEWTON-RAPHSON PARA SISTEMAS NAO LINEARES nos forneceu um resultado para\na concentracao de A igual a %g e concentracao de B igual a %g,\ncom um total de %u iteracoes e um erro igual a %g.\n',xi,yi,i,erro)
xi;
yi;
i;
errox;
erroy;

#Questão 3:
function fluxo = fluxo(f)
  eD = 0.004;
  Re = 200000;
  fluxo = -1/(f^(0.5)) - 2*log10((eD/3.7) + (2.51/(Re*(f^(0.5)))));
endfunction
tol = 0.01;
x0 = 0.001;
x1 = 0.005;
erro = abs((x1-x0)/(x0));
i = 0;
while erro > tol
  xv = x1;
  xn = x1 - ((fluxo(x1)*(x0 - x1))/(fluxo(x0)-fluxo(x1)));
  erro = abs((xn - xv)/xv);
  x0 = xv;
  x1 = xn;
  i = i+1;
endwhile
iteracao_secante = i;
printf('\n\nO METODO DA SECANTE nos forneceu um resultado para o fluxo igual a %g,\ncom um total de %u iteracoes e um erro relativo igual a %g.\n', xn,iteracao_secante,erro)
