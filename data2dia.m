function total=data2dia(dia,mes,varargin)
%=========================================================================%
% função data2dia: retorna o total de dias passados até determinada data  %
%                  de determinado ano.                                    %                                      
% entrada: dia   = dia do mês                                             %
%          mês   = mês                                                    %
%          ano   = ano, para casos de ano bissexto                        %
% saída:   total = dia no ano                                             %
%=========================================================================%
% inicialização da variavel total
total = 0;

% adição do valor do dia ao total
total = total + dia;

if isempty(varargin)
    ano = 2007;
else
    ano=varargin{1};
end

% loop de adição de dias de acordo com o mês até o mês anterior ao dado
for m=1:(mes-1)
  switch m
    case {1 3 5 7 8 10}     % meses com 31 dias
	  total = total + 31;
    case {4 6 9 11}         % meses com 30 dias
	  total = total + 30;
	case {2}                % no caso de fevereiro verifica se é ano bissexto
      if bissexto(ano)      % e atribui valor de acordo com esse dado
          total = total + 29;
      else
          total = total + 28;
      end
    case {12}               
      
    otherwise  % caso de entrada de número invalido para o mês
      %error('invalid number');
  end
end