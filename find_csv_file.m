function caminho = find_csv_file(nome_arquivo, pasta_dados, script_dir)
% FIND_CSV_FILE Tenta localizar nome_arquivo em varias localizacoes
% plausiveis:
%   1. caminho relativo a pasta_dados
%   2. diretorio do script (script_dir)
%   3. diretorio pai do script (provavel raiz do projeto)
%   4. busca recursiva a partir do diretorio do script
%
% Retorna o caminho absoluto do arquivo se encontrado, ou '' (string
% vazia) caso contrario.

caminho = '';

% 1) caminho direto relativo a pasta_dados
path1 = fullfile(pasta_dados, nome_arquivo);
if exist(path1, 'file') == 2
    caminho = path1;
    return;
end

% 2) mesmo diretorio do script
path2 = fullfile(script_dir, nome_arquivo);
if exist(path2, 'file') == 2
    caminho = path2;
    return;
end

% 3) diretorio pai do script
path3 = fullfile(script_dir, '..', nome_arquivo);
if exist(path3, 'file') == 2
    caminho = path3;
    return;
end

% 4) busca recursiva a partir do diretorio do script
encontrado = busca_recursiva(script_dir, nome_arquivo);
if ~isempty(encontrado)
    caminho = encontrado;
end

end


function resultado = busca_recursiva(pasta, nome_arquivo)
% Subfuncao auxiliar (permitida em arquivos de funcao em qualquer versao
% do MATLAB, inclusive R2016a). Nao usa o campo "folder" de dir(), que
% so existe a partir do R2016b.

resultado = '';

candidato = fullfile(pasta, nome_arquivo);
if exist(candidato, 'file') == 2
    resultado = candidato;
    return;
end

conteudo = dir(pasta);
for i = 1:numel(conteudo)
    item = conteudo(i);
    if item.isdir && ~strcmp(item.name, '.') && ~strcmp(item.name, '..')
        subpasta = fullfile(pasta, item.name);
        r = busca_recursiva(subpasta, nome_arquivo);
        if ~isempty(r)
            resultado = r;
            return;
        end
    end
end

end
