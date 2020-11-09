##  Programadores: Gabriel Salateo Rosin e Daniel Bascelar Vozza  ##

clc
clear
close

##Banco de dados coletado
anos = [2015 2016 2016 2016 2016 2016 2015 2015 2014 2016 2013 2016 2015 2015 2016 2015 2016 2016 2016 2015 2018 2017 2017 2017 2017 2016 2016 2016 2016 2015 2015 2013 2013 2017 2015 2015 2015 2015 2014 2013 2015 2015 2013 2016 2014 2016 2014 2014 2015 2015 2013 2016 2014 2014 2015 2013 2015 2013 2016 2019 2015 2015 2018 2014 2014 2018 2016 2014 2014 2016 2015 2014 2017 2015 2015 2013 2015 2012 2015 2016 2018 2017 2017 2018 2017];
Km = [28800 15000 29000 33000 31500 37200 80876 60000 64444 34155 69000 36735 34472 66000 99999 27000 32500 30251 24000 37000 39417 21031 35000 37000 42655 32300 33000 34000 50000 30000 48000 74000 81354 64570 33300 56000 60000 61175 51500 77000 102969 61175 81354 34000 110466 32300 65000 51500 102969 61175 81354 34000 65000 51500 69700 40000 72000 77000 61205 3500 45000 30000 12104 51500 82000 14000 50000 110466 81000 32300 79000 55539 34580 55390 46000 45000 69000 98972 69700 49326 32000 33000 13648 44650 42655];
preco = [32700 34900 33799 34990 32900 34970 32990 30990 28900 33490 27900 36890 31980 34990 39913 32990 34500 34500 34900 33000 36900 36490 36990 36990 38730 35890 39890 37990 35990 35990 34990 28990 32990 35990 39990 31990 30990 35990 32890 30990 34900 35990 32990 37990 31900 35890 33990 32890 34900 35990 32990 37990 33990 32890 34900 33990 32200 30990 35990 43800 34990 35990 42500 32890 32490 42990 37990 31900 37120 35890 34990 34890 39900 34900 33900 33990 32990 32990 34900 35500 35900 35900 36990 36900 38730];

##Operações para normalizar
media_anos = sum(anos)/length(anos);
media_Km = sum(Km)/length(Km);
media_preco = sum(preco)/length(preco);
max_anos = max(anos);
min_anos = min(anos);
max_Km = max(Km);
min_Km = min(Km);
max_preco = max(preco);
min_preco = min(preco);
anos_norm = (anos - media_anos)/(max_anos - min_anos);
preco_norm = (preco - media_preco)/(max_preco - min_preco);
Km_norm = (Km - media_Km)/(max_Km - min_Km);
 
anos_norm_transp = anos_norm';
Km_norm_transp = Km_norm';
preco_norm_transp = preco_norm';
Km_norm_quad = Km_norm_transp.^2;
anos_norm_quad = anos_norm_transp.^2;
anos_Km_norm = Km_norm_transp.*anos_norm_transp;
coluna1 = ones(85,1);

##Montagem das matrizes do método
Fi = [coluna1 anos_norm_transp Km_norm_transp anos_norm_quad anos_Km_norm Km_norm_quad];
Fi_transp = Fi';
A = Fi_transp*Fi;
b = Fi_transp*preco_norm_transp;
c = b\A;

Fizao = c*Fi_transp;

x = anos_norm;
y = Km_norm;
z = preco_norm;

##Plot dos gráficos 1 e 2
subplot(2,2,1)
plot3(x,y,z, '.db')
title("Pontos Tabelados")
xlabel("Ano do Veículo Normalizado")
ylabel("Kilometragem Normalizada")
zlabel("Preço Normalizado")

[X,Y] = meshgrid(x,y);
jonson = c(1) + c(2)*X + c(3)*Y + c(4)*X*X + c(5)*X*Y + c(6)*Y*Y;
subplot(2,2,2)
mesh(X,Y,jonson)
title("Superfície Ajustada")
xlabel("Ano do Veículo Normalizado")
ylabel("Kilometragem Normalizada")
zlabel("Preço Ajustado")

##Erro médio quadrático:
for i = 1:length(preco_norm)
  erro_i2 = (preco_norm(i) - Fizao(i))^2;
  erros(1,i) = erro_i2;
endfor
erro_quadratico_medio = sum(erros)/length(preco_norm);

##Laço para calcular o agio (hi = 100 * (pi-p(ai,ki))/p(ai,ki))
H = [];
justo = 0;
barato = 0;
caro = 0;
for i = 1:length(preco)
  hi = 100*((preco_norm(i) - (c(1) + c(2)*anos_norm(i) + c(3)*Km_norm(i) + c(4)*anos_norm_quad(i) + c(5)*anos_norm(i)*Km_norm(i) + c(6)*Km_norm_quad(i)))/(c(1) + c(2)*anos_norm(i) + c(3)*Km_norm(i) + c(4)*anos_norm_quad(i) + c(5)*anos_norm(i)*Km_norm(i) + c(6)*Km_norm_quad(i)));
  if hi <= -10
    barato++;
  elseif hi >= 10
    caro++;
  else
    justo++;
  endif
  H(i) = hi; 
endfor

menor_hi = min(H);
contador = 1;
i = H(contador);
while i!=menor_hi
  contador ++;
  i = H(contador);
endwhile

subplot(2,2,3)
plot(preco,H,'.db')
title("Ágio x Anúncio")
xlabel("Preço do Anúncio (R$)")
ylabel("Ágio Percentual do Anúncio")

ano_celta = (2016 - media_anos)/(max_anos-min_anos);
km_celta = (46000 - media_Km)/(max_Km-min_Km);

preco_celta_norm = c(1) + c(2)*ano_celta + c(3)*km_celta + c(4)*ano_celta^2 + c(5)*ano_celta*km_celta + c(6)* km_celta^2;
format long;
preco_celta = preco_celta_norm*(max_preco-min_preco) + media_preco;
printf("O valor justo para o Celta de Melissa Moira é R$%f\n\n",preco_celta);

##Estudo da depreciação
ano_celta = (2016 - media_anos)/(max_anos-min_anos);
km_celta_8k = (54000 - media_Km)/(max_Km-min_Km);

preco_celta_norm_8k = c(1) + c(2)*ano_celta + c(3)*km_celta_8k + c(4)*ano_celta^2 + c(5)*ano_celta*km_celta_8k + c(6)* km_celta_8k^2;
preco_celta_8k = preco_celta_norm_8k*(max_preco-min_preco) + media_preco;
depreciacao_8k = preco_celta - preco_celta_8k;

ano_celta = (2016 - media_anos)/(max_anos-min_anos);
km_celta_10k = (56000 - media_Km)/(max_Km-min_Km);

preco_celta_norm_10k = c(1) + c(2)*ano_celta + c(3)*km_celta_10k + c(4)*ano_celta^2 + c(5)*ano_celta*km_celta_10k + c(6)* km_celta_10k^2;
preco_celta_10k = preco_celta_norm_10k*(max_preco-min_preco) + media_preco;
depreciacao_10k = preco_celta - preco_celta_10k;

ano_celta_1ano = (2015 - media_anos)/(max_anos-min_anos);
km_celta = (46000 - media_Km)/(max_Km-min_Km);

preco_celta_norm_1ano = c(1) + c(2)*ano_celta_1ano + c(3)*km_celta + c(4)*ano_celta_1ano^2 + c(5)*ano_celta_1ano*km_celta + c(6)* km_celta^2;
preco_celta_1ano = preco_celta_norm_1ano*(max_preco-min_preco) + media_preco;
depreciacao_1ano = preco_celta - preco_celta_1ano;

##Opção para o usuário avaliar seu onix
ano_usuario = input("Caso queira saber o valor justo para o seu Onix, por favor, entre primeiramente com o ano do carro: \n\n");
km_usuario = input("Agora, forneca a kilometragem do veiculo (em km): \n\n");

ano_usuario = (ano_usuario - media_anos)/(max_anos-min_anos);
km_usuario = (km_usuario - media_Km)/(max_Km-min_Km);

preco_usuario = c(1) + c(2)*ano_usuario + c(3)*km_usuario + c(4)*ano_usuario^2 + c(5)*ano_usuario*km_usuario + c(6)*km_usuario^2;
preco_usuario = preco_usuario*(max_preco-min_preco) + media_preco;

printf("O valor justo para o seu Celta é R$%f\n\n",preco_usuario);
