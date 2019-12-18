function bool = bissexto(ano)
%*******************************************************
% Fun��o que retorna verdadeiro se o ano � bissexto e
% falso caso contr�rio
%*******************************************************

% variavel booleana
bool = false;

% regra geral para anos bissextos
if rem(ano,4)==0
    bool = true;
end
% exce��o � regra
if rem(ano,100)==0
    bool = false;
end
% exce��o da exce��o
if rem(ano,400)==0
    bool = true;
end