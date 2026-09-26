function [t_us, V_A, V_B, V_C, I_PR, dt_us] = preparar_dados_seg3(nome_arquivo)
% PREPARAR_DADOS_SEG3 A partir de um CSV do ATP, retorna as 3 tensoes de
% isolador no SEG3 (SEG3 <-> Condutor A, B, C) e a corrente injetada
% (Heidler, PR), com polaridade alinhada, prontas para o metodo DE.
%
% As colunas sao localizadas pelo NOME (nao por indice fixo), pois a
% ordem das colunas pode variar entre arquivos exportados do ATP.

[colunas_nomes, dados] = load_atp_csv(nome_arquivo);

t_s = dados(:, 1);
t_us = t_s * 1e6;

idx_A = buscar_coluna(colunas_nomes, 'X0033A@SEG3');
idx_B = buscar_coluna(colunas_nomes, 'X0033B@SEG3');
idx_C = buscar_coluna(colunas_nomes, 'X0033C@SEG3');

if isempty(idx_A) || isempty(idx_B) || isempty(idx_C)
    error('Colunas SEG3<->condutor n�o encontradas em %s', nome_arquivo);
end

V_A = alinhar_polaridade(dados(:, idx_A));
V_B = alinhar_polaridade(dados(:, idx_B));
V_C = alinhar_polaridade(dados(:, idx_C));

% Busca dinamica da coluna de corrente injetada (nome contem 'XX0036')
idx_corrente = [];
for i = 1:numel(colunas_nomes)
    if ~isempty(strfind(upper(colunas_nomes{i}), 'XX0036')) %#ok<STREMP>
        idx_corrente = i;
        break;
    end
end

if ~isempty(idx_corrente)
    I_PR = alinhar_polaridade(dados(:, idx_corrente));
else
    I_PR = [];
end

dt_us = t_us(2) - t_us(1);

end
