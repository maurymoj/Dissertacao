function total=data2dia(dia,mes,varargin)
%=========================================================================%
% fun��o data2dia: retorna o total de dias passados at� determinada data  %
%                  de determinado ano.                                    %                                      
% entrada: dia   = dia do m�s                                             %
%          m�s   = m�s                                                    %
%          ano   = ano, para casos de ano bissexto                        %
% sa�da:   total = dia no ano                                             %
%=========================================================================%
% inicializa��o da variavel total
total = 0;

% adi��o do valor do dia ao total
total = total + dia;

if isempty(varargin)
    ano = 2007;
else
    ano=varargin{1};
end

% loop de adi��o de dias de acordo com o m�s at� o m�s anterior ao dado
for m=1:(mes-1)
  switch m
    case {1 3 5 7 8 10}     % meses com 31 dias
	  total = total + 31;
    case {4 6 9 11}         % meses com 30 dias
	  total = total + 30;
	case {2}                % no caso de fevereiro verifica se � ano bissexto
      if bissexto(ano)      % e atribui valor de acordo com esse dado
          total = total + 29;
      else
          total = total + 28;
      end
    case {12}               
      
    otherwise  % caso de entrada de n�mero invalido para o m�s
      %error('invalid number');
  end
end