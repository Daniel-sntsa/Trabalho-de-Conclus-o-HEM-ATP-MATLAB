function idx = buscar_coluna(colunas_nomes, alvo)
% BUSCAR_COLUNA Retorna o indice da coluna cujo nome combinado
% (ex: 'X0033A@SEG3') seja igual a 'alvo'. Retorna [] se nao encontrar.

idx = find(strcmp(colunas_nomes, alvo), 1);

end
