function total=data2dia_vet(dia,mes,ano)
%=========================================================================%
% função: data2dia_vet retorna um vetor com o total de dias passados até  %
%         as datas de entrada                                             %
% entrada: dia   = vetor com dias do mês                                  %
%          mês   = mês a que os dias do vetor dia se referem              %
%          ano   = ano em questão para correções no caso de ano bissexto  %
% saída:   total = dia no ano (ex. 02/02 é o 33º dia do ano)              %
%=========================================================================%

% inicialização do vetor com os dias do ano
total = zeros(length(dia),1);

% loop que passa cada linha dos vetores de entrada na função data2dia
% que retorna o valor do dia no ano para cada linha.
for i=1:length(dia)
    total(i)=data2dia(dia(i),mes(i),ano(i));
end