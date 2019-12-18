function total=data2dia_vet(dia,mes,ano)
%=========================================================================%
% fun��o: data2dia_vet retorna um vetor com o total de dias passados at�  %
%         as datas de entrada                                             %
% entrada: dia   = vetor com dias do m�s                                  %
%          m�s   = m�s a que os dias do vetor dia se referem              %
%          ano   = ano em quest�o para corre��es no caso de ano bissexto  %
% sa�da:   total = dia no ano (ex. 02/02 � o 33� dia do ano)              %
%=========================================================================%

% inicializa��o do vetor com os dias do ano
total = zeros(length(dia),1);

% loop que passa cada linha dos vetores de entrada na fun��o data2dia
% que retorna o valor do dia no ano para cada linha.
for i=1:length(dia)
    total(i)=data2dia(dia(i),mes(i),ano(i));
end