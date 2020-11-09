//  Cálculo da espessura mínima de anéis de costado de um tanque (API 650 - básico) //
    

clc
clear

// Dados de Entrada://
densidade = input('Insira a densidade do produto a ser armazenado (g/cm³): ')
largura = input('Insira a largura das chapas a serem usadas (ft): ')
C = input('Insira a sobre espessura para corrosão do tanque (pol/ano): ')
diametros = input('Insira os diâmetros em ft em forma de vetor para os quais quer projetar tanques (Exemplo: [100;200;300]): ')
numero_aneis = input('Insira o número de aneis em forma de vetor para o tanque (Exemplo: [5;6;7;8]): ')

alturas = largura*numero_aneis;
densidade_valor = densidade

// Ajuste necessário caso a densidade seja menor que 1 //
if densidade < 1 then
    densidade = 1;
end;


 // Dados de eficiência da junta soldada e tensão admissível //
// Considerando a eficiência da junta soldada com valor de 0.85 e a tensão admissível do material do costado com valor de 21000 psi //

E = 0.85;
S = 21000;



// Cálculo das capacidades dos tanques projetados (conversão: 1 ft³ = 0.178108 bbl)//
for i = 1:length(diametros)
    for j = 1:length(alturas)
    capacidades(j + (i - 1)*length(alturas)) = %pi*((diametros(i)/2)^2)*alturas(j)*0.178108;
    end
end

// Cálculo da espessura mínima de cada anel de costado //
saida = [];
for i = 1:length(diametros)
    t = [];
    for j = 1:length(alturas)
        for k =1:numero_aneis(j)
            t(j,k) = (2.6*diametros(i)*((alturas(j)-largura*(k-1))-1)*densidade/(E*S)) + C;
        end
    end
    saida = [saida;t];
end

// Cálculo da massa aproximada do tanque em kg/ft³ //
densidade_chapa = 18.5239;
for i = 1:length(diametros)
    for j = (1 + (i-1)*length(alturas)):(i*length(alturas))
        espessura_total = 0;
        for k = 1:numero_aneis($)
        espessura_total = espessura_total + saida(j,k);
        end
        espessuras_totais(j) = espessura_total;
    end
end

// Densidade calculada pela tabela de peso de chapas de aço fornecida em aula (18.524 kg/ft³) //
for i = 1:length(diametros)
    for j = 1:length(alturas)
        n = j + (i - 1)*length(alturas);
        massas_tanques(n) = densidade_chapa*%pi*diametros(i)*largura*espessuras_totais(n);
    end
end

// Saída do Programa //

mprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
mprintf('\nProjeto de Tanques de Armazenamento segundo API-650\n')
mprintf('\nCostado construído de chapas de %d ft de largura\n', largura)
mprintf('\nDensidade do produto armazenado  = %.3f g/cm³\n',densidade_valor)
mprintf('\nSobre espessura para corrosão = %.1f pol/ano\n',C)
mprintf('\n\tDiâmetro (ft) \t Altura (ft) \t Capacidade (bbl)\tMassa Aproximada (kg) \t Espessura de cada Anel de Costado (pol)\n')
mprintf('\n\t\t\t\t\t\t\t\t\t\t\t ')

for i = 1:numero_aneis($)
    mprintf('n = %d\t ',i)
end

for i = 1:length(diametros)
    mprintf('\n\t    %d',diametros(i))
    for j = 1:length(alturas)
        n = j + (i - 1)*length(alturas)
        mprintf('\t\t    %d \t\t    %.3f\t\t    %.3f\t', alturas(j), capacidades(n), massas_tanques(n))
        for k = 1:numero_aneis(j)
            mprintf('\t %.3f', saida(n,k))
        end
        mprintf('\n\t')
    end
end
