function [colunas_nomes, dados] = load_atp_csv(nome_arquivo)
% LOAD_ATP_CSV Le um CSV do ATP com 3 linhas de cabecalho, localizando
% cada coluna por um NOME combinado (linha 2 + linha 3), robusto a
% qualquer ordem de colunas entre arquivos diferentes.
%
% Retorna:
%   colunas_nomes : cell array de strings, ex: 'X0033A@SEG3'
%   dados         : matriz numerica (linhas = amostras, colunas = variaveis)

fid = fopen(nome_arquivo, 'r');
if fid == -1
    error('N�o foi poss�vel abrir o arquivo %s', nome_arquivo);
end

linha1 = fgetl(fid); %#ok<NASGU> % tipo (Voltage/Current/Time) - n�o usada
linha2 = fgetl(fid); % nome principal (PR, SEG3, X0033A, ...)
linha3 = fgetl(fid); % sublabel / origem (PR, SEG2, SEG3, ...)

% Detecta delimitador com base na 1a linha de cabe�alho
if any(linha2 == ';')
    delimitador = ';';
else
    delimitador = ',';
end

nomes_parte = strtrim(strsplit(linha2, delimitador));
subs_parte  = strtrim(strsplit(linha3, delimitador));

n_cols = numel(nomes_parte);
colunas_nomes = cell(1, n_cols);
for i = 1:n_cols
    nome = nomes_parte{i};
    if i <= numel(subs_parte)
        sub = subs_parte{i};
    else
        sub = '';
    end
    if ~isempty(sub)
        colunas_nomes{i} = sprintf('%s@%s', nome, sub);
    else
        colunas_nomes{i} = nome;
    end
end

% Le as linhas de dados
linhas = {};
while true
    linha = fgetl(fid);
    if ~ischar(linha)
        break;
    end
    if isempty(strtrim(linha))
        continue;
    end
    linhas{end+1} = linha; %#ok<AGROW>
end
fclose(fid);

n_lin = numel(linhas);
dados = zeros(n_lin, n_cols);

for i = 1:n_lin
    partes = strtrim(strsplit(linhas{i}, delimitador));
    if delimitador == ';'
        partes = strrep(partes, ',', '.'); % separador decimal = v�rgula neste formato
    end
    for c = 1:n_cols
        dados(i, c) = str2double(partes{c});
    end
end

end
