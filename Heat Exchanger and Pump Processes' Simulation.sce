/*
Programadores:

Gabriel Salatêo Rosin
Daniel Bacellar Souza Vozza 
Lucas Silva Andrade
Márcio Messias de Moraes Neto
Murilo dos Santos Gabriel

*/

clc
clear
close

global R Tc Pc Tref w  
R = 8.314
Tc = 562.05
Pc = 48.94*(10^5)
Tref = 278.15
w = 0.2092
vazao_m = 100 //kg/s
vazao_n = 1280.2458072

function Cp = cp(T)     //função para calcular o valor do cp para uma dada temperatura - Coeficientes retirados do Sandler
    A = -36.193
    B = 0.48444
    C = -31.548*(10^-5)
    D = 77.573*(10^-9)
    Cp =(A + B*T + C*(T^2) + D*(T^3))
endfunction

function integral = integrar(Tinf, Tsup) //função para integrar o cp e calcular o delta H de um GI
    n = 2500
    intervalo = abs(Tsup - Tinf)/n
    valor = 0
    for i = 1:n
        finf = cp(Tinf)
        Tseguinte = Tinf + intervalo
        fsup = cp(Tseguinte)
        Tinf = Tseguinte
        valor = valor + (((fsup + finf)*intervalo)*0.5)
    end
    integral = valor
endfunction

function integral = integrarS(Tinf, Tsup) // Para Cp/T do balanço de entropia 
    n = 2500
    intervalo = abs(Tsup - Tinf)/n
    valor = 0
    for i = 1:n
        finf = cp(Tinf)/Tinf
        Tseguinte = Tinf + intervalo
        fsup = cp(Tseguinte)/Tseguinte
        Tinf = Tseguinte
        valor = valor + (((fsup + finf)*intervalo)*0.5)
    end
    integral = valor
endfunction

//Função de estado SRK. Entrar com os dados na SRK EM K e bar
function v = SRK(T,P)
    P=P*10^5
    teste=0
    k = 0.48508 + 1.55171*w -0.15613*w*w
    alpha = (1+k*(1-((T/Tc)^(0.5))))^2
    b = 0.0867*R*Tc/Pc
    a = 0.42747*(R^2)*(Tc^2)*alpha/Pc
    A = a*P/((R^2)*(T^2))
    B = b*P/(R*T)
    delta = -A*B
    gama = A-B-B^2
    Beta = -1
    q = ((Beta^2) - (3*gama))/9
    r = (2*(Beta)^3 - 9*Beta*gama + 27*delta)/54
    
    if (q^3 - r^2)>= 0 then
        teta = acos(r/(q^1.5))
        Zl = -2*(q^0.5)*(cos(teta/3)) - Beta/3
        Zv = -2*(q^0.5)*(cos((teta + 2*%pi)/3)) - Beta/3
        caso=1
    end 
    if r>0 then
            psi = 1
    else
            psi = -1
    end
    if (q^3 - r^2)< 0 then
        Z = -psi*((((r^2) - (q^3))^0.5 + abs(r))^(1/3) + q/((((r^2) - (q^3))^0.5 + abs(r))^(1/3))) - Beta/3
        Zv = Z
        Zl=Z
        caso=2
    end
       
    
    vl = Zl*R*T/P
    vv = Zv*R*T/P
    fugl = (exp(Zl-1-log(Zl-B)-((A/B)*log(1+(B/Zl)))))
    fugv = (exp(Zv-1-log(Zv-B)-((A/B)*log(1+(B/Zv)))))
    v = [vl vv Zl Zv caso fugl fugv]
endfunction

//Função usada para encontrar a pressão de saturação a uma dada temperatura
function Psat = psat(T)
//Chute inicial feito com a equação de Antoine a fim de fazer um número menor de iterações para achar a pressão de saturação
    P =10^(4.60362-(1701.073/(T+20.806)))
    v_sat = SRK(T,P)
//Psat ocorre quando as fugacidades forem iguais.
    fugl = v_sat(6)
    fugv = v_sat(7)
    while abs(fugl - fugv) > 10^-5
        P = P*(fugl/fugv)
        v_sat = SRK(T, P)
        fugl = v_sat(6)
        fugv = v_sat(7)
    end
    //Psat é dada em bar
    Psat = P
    
endfunction

//------------------------------------------------------ Funções Residuais ----------------------------------------------------------------------------------

// Volume para gás ideal é vgi = R*T/P
// Vale notar que apesar de estarem em maiúsculo, essas propriedades são específicas 

// Energia Interna Residual
function Ur = Uresidual(T,v)
    ac = (0.427480*R^2*Tc^2)/Pc
    K = 0.48508 + 1.55171*w - 0.15613*w^2
    alpha2 = (1+K*(1-(T/Tc)^0.5))^2
    a2 = ac*alpha2
    b2 = 0.086640*R*Tc/Pc
    Ur = (a2/b2)*(((K*T)/((T*Tc*alpha2)^0.5))+ 1)*log(v/(v+b2))
endfunction

// Entalpia Residual
function Hr = Hresidual(T,v,P)
    P=P*10^5
    ac = (0.427480*R^2*Tc^2)/Pc
    K = 0.48508 + 1.55171*w - 0.15613*w^2
    alpha2 = (1+K*(1-(T/Tc)^0.5))^2
    a2 = ac*alpha2
    b2 = 0.086640*R*Tc/Pc
    vgi = R*T/P 
    Hr = (a2/b2)*(((K*T)/((T*Tc*alpha2)^0.5))+ 1)*log(v/(v+b2)) + R*T*((v/vgi) - 1)
endfunction

//Entropia Residual
function Sr = Sresidual(T,v,P)
    P=P*10^5
    ac = (0.427480*R^2*Tc^2)/Pc
    K = 0.48508 + 1.55171*w - 0.15613*w^2
    alpha2 = (1+K*(1-(T/Tc)^0.5))^2
    a2 = ac*alpha2
    b2 = 0.086640*R*Tc/Pc
    vgi = R*T/P
    Sr = ((a2*K)/(b2*(T*Tc*alpha2)^0.5))*log(v/(v+b2)) - R*log(vgi/(v-b2))
endfunction

// Energia Livre de Gibbs Residual
function Gr = Gresidual(T,v,P)
    P=P*10^5
    ac = (0.427480*R^2*Tc^2)/Pc
    K = 0.48508 + 1.55171*w - 0.15613*w^2
    alpha2 = (1+K*(1-(T/Tc)^0.5))^2
    a2 = ac*alpha2
    b2 = 0.086640*R*Tc/Pc
    vgi = R*T/P
    Gr = (a2/b2)*log(v/(v+b2)) + R*T*log(vgi/(v-b2)) + R*T*((v/vgi) - 1)
endfunction

// Energia de Helmholtz Residual
function Ar = Aresidual(T,v,P)
    P=P*10^5
    ac = (0.427480*R^2*Tc^2)/Pc
    K = 0.48508 + 1.55171*w - 0.15613*w^2
    alpha2 = (1+K*(1-(T/Tc)^0.5))^2
    a2 = ac*alpha2
    b2 = 0.086640*R*Tc/Pc
    vgi = R*T/P
    Ar = (a2/b2)*log(v/(v+b2)) + R*T*log(vgi/(v-b2))
endfunction

//----------------------------------------- Análise das Correntes -------------------------------------------

//Referência
//Encontrar a pressão da referência: líquido saturado
Pref = psat(278.5)
v_ref = SRK(Tref,Pref)
//Propriedades na referência
ur_ref=Uresidual(Tref,v_ref(1))
sr_ref = Sresidual(Tref,v_ref(1),Pref)
hr_ref=Hresidual(Tref,v_ref(1),Pref)
gr_ref=Gresidual(Tref,v_ref(1),Pref)
ar_ref=Aresidual(Tref,v_ref(1),Pref)



//Corrente 1
T_1 = 420
P_1 = 6.5
//pressão de saturação na corrente 1
P_1_sat=psat(T_1)
//Verificar se a corrente está líquida ou vapor
v_1 = SRK(T_1,P_1)
    if(P_1_sat>P_1)then
        sr_1 = Sresidual(T_1,v_1(2),P_1)
    else
        sr_1 = Sresidual(T_1,v_1(1),P_1)
    end
//Dalta s de gás ideal
delta_s_gi_1 = integrarS(Tref,T_1)-R*log(P_1/Pref)
//Delta s com a referência
s_1 = sr_1 - sr_ref + delta_s_gi_1

densidade_1 = 1/(v_1(1))

//Corrente 2
// Varificamos que a corrente 1 está líquida. Então como só houve um grande aumento da pressão, sabemos que essa corrente também está líquida. Também verificamos isso no NIST e no site ""Patrick Barrie's program for solving cubic equations of state""

T_2 = 420    //chute inicial: mesma temperatura
P_2 = 25
v_2 = SRK(T_2,P_2)
sr_2 = Sresidual(T_2,v_2(1),P_2)
delta_s_gi_2 = integrarS(Tref,T_2)-R*log(P_2/Pref)
s_2 = sr_2 - sr_ref + delta_s_gi_2

//Balanço na bomba:
//À partir do balanço de entropia na bomba, obtivemos que a variação da entropia é igual a 0
//Assim encontraremos a Temperatura da corrente 2 resolvendo numericamente a equação Sresidual_2 - Sresidual_1 + DeltaS_GI(T1 a T2)) = 0
h = 0.0001 
delta_s_gi_12 = integrarS(T_1, T_2) - R*log(P_2/P_1)
fx = Sresidual(T_2,v_2(1),P_2) - Sresidual(T_1,v_1(1),P_1) + integrarS(T_1, T_2) - R*log(P_2/P_1)
while abs(fx) > h
    derivada = (((Sresidual(T_2 + h,v_2(1),P_2) - Sresidual(T_1,v_1(1),P_1) + integrarS(T_1, T_2 + h) - R*log(P_2/P_1))) - (Sresidual(T_2,v_2(1),P_2) - Sresidual(T_1,v_1(1),P_1) + integrarS(T_1, T_2) - R*log(P_2/P_1)))*10000
    T_2 = T_2 - (fx/derivada)
    fx = Sresidual(T_2,v_2(1),P_2) - Sresidual(T_1,v_1(1),P_1) + integrarS(T_1, T_2) - R*log(P_2/P_1)    
end

//Atualizando os valores da corrente 2 pra nova T_2
v_2 = SRK(T_2,P_2)
sr_2 = Sresidual(T_2,v_2(1),P_2)
delta_s_gi_2 = integrarS(Tref,T_2)-R*log(P_2/Pref)
s_2 = sr_2 - sr_ref + delta_s_gi_2
densidade_2 = 1/(v_2(1))

// Corrente 3 - Pós trocador de calor isobárico
//
T_3 = 515
P_3 = P_2
P_3_sat=psat(T_3)
v_3 = SRK(T_3,P_3)
//Confirmação que a corrente 3 está líquida 
    if(P_3_sat>P_3)then
        sr_3 = Sresidual(T_3,v_3(2),P_3)
    else
        sr_3 = Sresidual(T_3,v_3(1),P_3)
    end
delta_s_gi_3 = integrarS(Tref,T_3)-R*log(P_3/Pref)
s_3 = sr_3 - sr_ref + delta_s_gi_3
densidade_3 = 1/(v_3(2))

//Cálculo das propriedades nas correntes 

//Entalpia
hr_1 = Hresidual(T_1,v_1(1),P_1)
h_1 = hr_1 - hr_ref + integrar(Tref,T_1)
hr_2 = Hresidual(T_2,v_2(1),P_2)
h_2 = hr_2 - hr_ref + integrar(Tref,T_2)
hr_3 = Hresidual(T_3,v_3(2),P_3)
h_3 = hr_3 - hr_ref + integrar(Tref,T_3)

//Energia Interna
u_1 = h_1 - P_1*v_1(1)
u_2 = h_2 - P_2*v_2(1)
u_3 = h_3 - P_3*v_3(2)

//Energia Livre de Helmholtz
a_1 = u_1 - T_1*s_1
a_2 = u_2 - T_2*s_2
a_3 = u_3 - T_3*s_3

//Energia Livre de Gibbs
g_1 = a_1 + P_1*v_1(1)
g_2 = a_2 + P_2*v_2(1)
g_3 = a_3 + P_3*v_3(2)

//Vetores para a tabela

temperatura = [T_1 T_2 T_3]
pressao = [P_1 P_2 P_3]
volume = [v_1(1) v_2(1) v_3(2)]
entropia = [s_1 s_2 s_3]
entalpia = [h_1 h_2 h_3]
energia_interna = [u_1 u_2 u_3]
helmholtz = [a_1 a_2 a_3]
gibbs = [g_1 g_2 g_3]
densidade = [densidade_1 densidade_2 densidade_3]
propriedade = ["Vazão Mássica (kg/s)", "Vazão Molar (mol/s)", "Densidade (mol/m^3)", "Temperatura (K)", "Pressão (bar)", "Volume Específico (m^3/mol)", "Energia Interna (J/mol)", "Entalpia (J/mol)", "Entropia (J/mol.K)", "Energia Livre de Gibbs (J/mol)", "Energia Livre de Helmholtz (J/mol)"]
tabela = [vazao_m vazao_m vazao_m; vazao_n vazao_n vazao_n; densidade; temperatura; pressao; volume; energia_interna; entalpia; entropia; gibbs; helmholtz]


mprintf("A simulação do processo foi realizada com sucesso. Foi gerada uma tabela contendo as informações solicitadas para as três correntes com 5 casas decimais.\n")
mprintf(" Caso deseje uma maior precisão, por favor utilize a aba Navegador de Variáveis.\n")
mprintf("\n\t\t\tPropriedades de cada corrente\n\n\t Corrente 1\t        Corrente 2\t        Corrente 3")
for i=1:11    
    mprintf('\n')
    for j=1:3
        mprintf('\t%10.5f\t',tabela(i,j))
    end
    mprintf('\t%s',propriedade(i))
end

//Por fim mostramos um gráfico PxT que descreve o caminho tomado pelos processos utilizados no projeto:

temperaturas = [400:0.5:514]
for i = 1:length(temperaturas)
    pressoes(i) = psat(temperaturas(i))
end
processo1_T = [T_1 T_2]
processo1_P = [P_1 P_2]
processo2_T = [T_2 T_3]
processo2_P = [P_2 P_3]
plot(temperaturas, pressoes, 'r')
plot(processo1_T, processo1_P, 'b')
plot(processo2_T, processo2_P, 'c')
legend('Curva de saturação (líquido-vapor)','Bomba','Trocador de Calor', 4)
title('Gráfico PxT do processo')
xlabel('Temperatura (K)')
ylabel('Pressão (bar)')
