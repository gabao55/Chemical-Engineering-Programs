//      Diagramas Líquido-Vapor de múltiplos compostos  //
//          Programado por: Gabriel Salateo Rosin       //

clc, clear, close

function Psat = Psat(T,Pc,Tc,A,B,C,D)
    //#2: Calculando a Pressão de saturação
    global Pc
    global Tc
    global A
    global B
    global C
    global D
    X = 1 - (T/Tc)
    Psat = Pc * exp((1/(1-X))*(A *X + B*(X^1.5) + C * (X^3) + D * (X^6)))
endfunction

//#1: Declarando o valor das certas constantes para a água:
temperaturas_agua = (280:1:640)
temperaturas_etanol = (273:1:512)
temperaturas_acetona = (273:1:507)
temperaturas_benzeno = (300:1:561)

xx = [220.55 647.1 0.345 -7.76451 1.45838 -2.7758 -1.23303; 61.4 513.9 0.644 -8.51838 0.34163 -5.73683 8.32581; 47 508.1 0.304 -7.45514 1.20200 -2.43926 -3.35590; 48.9 562.2 0.212 -6.98273 1.33213 -2.62863 -3.33399]

Psat_agua = []
Psat_etanol = []
Psat_acetona = []
Psat_benzeno = []

n = length(temperaturas_agua)

Pc = xx(1,1)
Tc = xx(1,2)
A = xx(1,4)
B = xx(1,5)
C = xx(1,6)
D = xx(1,7)

for i=1:n
        T = temperaturas_agua(i)
        doug = Psat(T,Pc,Tc,A,B,C,D)
        Psat_agua = [Psat_agua doug]
end

Psat_agua = Psat_agua

n = length(temperaturas_etanol)

Pc = xx(2,1)
Tc = xx(2,2)
A = xx(2,4)
B = xx(2,5)
C = xx(2,6)
D = xx(2,7)

for i=1:n
        T = temperaturas_etanol(i)
        doug = Psat(T,Pc,Tc,A,B,C,D)
        Psat_etanol = [Psat_etanol doug]
end

Psat_etanol = Psat_etanol

n = length(temperaturas_acetona)

Pc = xx(3,1)
Tc = xx(3,2)
A = xx(3,4)
B = xx(3,5)
C = xx(3,6)
D = xx(3,7)

for i=1:n
        T = temperaturas_acetona(i)
        doug = Psat(T,Pc,Tc,A,B,C,D)
        Psat_acetona = [Psat_acetona doug]
end

Psat_acetona = Psat_acetona

n = length(temperaturas_benzeno)

Pc = xx(4,1)
Tc = xx(4,2)
A = xx(4,4)
B = xx(4,5)
C = xx(4,6)
D = xx(4,7)

for i=1:n
        T = temperaturas_benzeno(i)
        doug = Psat(T,Pc,Tc,A,B,C,D)
        Psat_benzeno = [Psat_benzeno doug]
end

Psat_benzeno = Psat_benzeno

subplot(2,2,1)
plot(temperaturas_agua,Psat_agua, 'c')
xtitle('Psat água','T (K)','Psat(bar)')

subplot(2,2,2)
plot(temperaturas_etanol,Psat_etanol, 'r')
xtitle('Psat etanol', 'T (K)','Psat(bar)')

subplot(2,2,3)
plot(temperaturas_acetona,Psat_acetona, 'm')
xtitle('Psat acetona', 'T (K)','Psat(bar)')

subplot(2,2,4)
plot(temperaturas_benzeno,Psat_benzeno, '-dg')
xtitle('Psat benzeno', 'T(K)','Psat(bar)')
