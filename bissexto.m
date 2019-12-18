function bool = bissexto(ano)
%*******************************************************
% Função que retorna verdadeiro se o ano é bissexto e
% falso caso contrário
%*******************************************************

% variavel booleana
bool = false;

% regra geral para anos bissextos
if rem(ano,4)==0
    bool = true;
end
% exceção à regra
if rem(ano,100)==0
    bool = false;
end
% exceção da exceção
if rem(ano,400)==0
    bool = true;
end