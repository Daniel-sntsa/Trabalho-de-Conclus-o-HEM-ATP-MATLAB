% =========================================================================
% Calculo da Corrente Critica (Ic) via Metodo do Efeito Disruptivo (DE)
% de Hileman (1999), para todos os cenarios de Fase da torre (F1-F4) x
% Resistividade do solo (MB, B, M, A, MA).
%
% Usa as tensoes SEG3 <-> Condutor A/B/C de cada arquivo do ATP como as
% tensoes de isolador, com:
%   - localizacao robusta das colunas por NOME (a ordem varia entre
%     arquivos exportados do ATP);
%   - alinhamento de polaridade (a onda SEG3<->condutor e frequentemente
%     negativo-dominante nesses dados; o metodo DE assume polaridade
%     positiva).
%
% Gera uma tabela com Fase | Zp (Ohm) | Ic (kA) | P(I>=Ic) (%), sinalizando
% cenarios cujo arquivo pareca duplicado/placeholder (pico de tensao
% identico ao de outro arquivo, o que fisicamente nao deveria ocorrer).
%
% A probabilidade de excedencia P(I>=Ic) segue a distribuicao cumulativa
% de corrente de primeira descida de Silveira & Visacro (2020), medida na
% estacao do Morro do Cachimbo:
%       P_I = 1 / (1 + (I/43.3)^3.8)
%
% Requer no mesmo diretorio (compativeis com MATLAB R2016a):
%   load_atp_csv.m, buscar_coluna.m, alinhar_polaridade.m,
%   curve_above.m, calcular_Ic.m, probabilidade_excedencia.m,
%   preparar_dados_seg3.m, find_csv_file.m
% =========================================================================

% -------------------------------------------------------------------
% PARAMETROS DO PROBLEMA
% -------------------------------------------------------------------
CFO = 1200.0; % kV - Tensao Critica de Flashover

FASES = {'F1', 'F2', 'F3', 'F4'};
FASES_DISPLAY = {'Fase I', 'Fase II', 'Fase III', 'Fase IV'};
RESISTIVIDADES = {'MB', 'B', 'M', 'A', 'MA'};

% Valores de Zp (Ohm) para cada Fase e categoria de resistividade,
% conforme as Tabelas 4.1 a 4.4 (colunas: MB, B, M, A, MA)
%              MB        B        M        A        MA
zp_tabela = [ ...
    2.2160,  6.2783, 10.8967, 22.1798, 45.2180;  % Fase I  (Tabela 4.1)
    2.2152,  5.9056,  9.3877, 15.8306, 32.0763;  % Fase II (Tabela 4.2)
    2.2152,  5.8432,  8.9945, 14.2635, 25.1802;  % Fase III (Tabela 4.3)
    2.2152,  5.8364,  8.8824, 13.6816, 21.5681]; % Fase IV (Tabela 4.4)

pasta_dados = pwd; % ajuste para a pasta onde estao os 20 CSVs, se necessario
script_dir = fileparts(mfilename('fullpath'));

tmax_us = 60.0;

% -------------------------------------------------------------------
% LOOP PRINCIPAL: roda para os 20 cenarios (Fase x Resistividade)
% -------------------------------------------------------------------
resultados = struct('Fase', {}, 'Zp_ohm', {}, 'Ic_kA', {}, 'P_Ic_pct', {}, ...
                     'Suspeito', {}, 'Arquivo', {});

% Mapa: pico de tensao (arredondado) -> lista de arquivos com esse pico
% (usado para detectar arquivos duplicados/placeholder)
picos_vistos = containers.Map('KeyType', 'double', 'ValueType', 'any');

idx_resultado = 0;

for fa = 1:numel(FASES)
    for r = 1:numel(RESISTIVIDADES)
        nome_arquivo = sprintf('LTTCC_V1_%s_%s.csv', RESISTIVIDADES{r}, FASES{fa});
        caminho = find_csv_file(nome_arquivo, pasta_dados, script_dir);
        zp_valor = zp_tabela(fa, r);

        idx_resultado = idx_resultado + 1;

        if isempty(caminho)
            fprintf('[AVISO] Arquivo não encontrado, pulando: %s\n', nome_arquivo);
            resultados(idx_resultado).Fase = FASES_DISPLAY{fa};
            resultados(idx_resultado).Zp_ohm = zp_valor;
            resultados(idx_resultado).Ic_kA = NaN;
            resultados(idx_resultado).P_Ic_pct = NaN;
            resultados(idx_resultado).Suspeito = 'ARQUIVO AUSENTE';
            resultados(idx_resultado).Arquivo = nome_arquivo;
            continue;
        end

        [t_us, V_A, V_B, V_C, I_PR, dt_us] = preparar_dados_seg3(caminho);

        if ~isempty(I_PR)
            CurrentPeak = max(I_PR) / 1000; % A -> kA
        else
            CurrentPeak = NaN;
        end

        % Assinatura do pico de tensao (para detectar arquivos duplicados)
        pico_assinatura = round(max(V_A), 2);
        if isKey(picos_vistos, pico_assinatura)
            lista = picos_vistos(pico_assinatura);
            lista{end+1} = nome_arquivo; %#ok<AGROW>
            picos_vistos(pico_assinatura) = lista;
        else
            picos_vistos(pico_assinatura) = {nome_arquivo};
        end

        Ic_A = calcular_Ic(V_A, t_us, tmax_us, CurrentPeak, CFO, dt_us);
        Ic_B = calcular_Ic(V_B, t_us, tmax_us, CurrentPeak, CFO, dt_us);
        Ic_C = calcular_Ic(V_C, t_us, tmax_us, CurrentPeak, CFO, dt_us);

        Ic_vals = [Ic_A, Ic_B, Ic_C];
        if all(isnan(Ic_vals))
            Ic_critico = NaN;
        else
            Ic_critico = min(Ic_vals); % min() ignora NaN por padrao no MATLAB
        end

        % Probabilidade de excedencia P(I >= Ic_critico), segundo a
        % distribuicao cumulativa de corrente de primeira descida de
        % Silveira & Visacro (2020), Eq. (1): P_I = 1 / (1 + (I/43.3)^3.8)
        if isnan(Ic_critico)
            P_Ic = NaN;
        else
            P_Ic = probabilidade_excedencia(Ic_critico) * 100; % em %
        end

        resultados(idx_resultado).Fase = FASES_DISPLAY{fa};
        resultados(idx_resultado).Zp_ohm = zp_valor;
        resultados(idx_resultado).Ic_kA = Ic_critico;
        resultados(idx_resultado).P_Ic_pct = P_Ic;
        resultados(idx_resultado).Suspeito = '';
        resultados(idx_resultado).Arquivo = nome_arquivo;

        fprintf('%s: Zp=%.4f Ohm  ->  Ic_critico=%.3f kA  |  P(I>=Ic)=%.3f %%\n', ...
                nome_arquivo, zp_valor, Ic_critico, P_Ic);
    end
end

% -------------------------------------------------------------------
% Marca como suspeitos os arquivos cujo pico de tensao bateu EXATAMENTE
% com o de outro(s) arquivo(s) - sinal de dados duplicados/placeholder
% -------------------------------------------------------------------
chaves = keys(picos_vistos);
arquivos_suspeitos = {};
for k = 1:numel(chaves)
    lista = picos_vistos(chaves{k});
    if numel(lista) > 1
        arquivos_suspeitos = [arquivos_suspeitos, lista]; %#ok<AGROW>
    end
end

for i = 1:numel(resultados)
    if any(strcmp(resultados(i).Arquivo, arquivos_suspeitos))
        resultados(i).Suspeito = 'PICO DUPLICADO EM OUTRO ARQUIVO - CONFERIR DADOS';
    end
end

% -------------------------------------------------------------------
% Salva a tabela em CSV (fopen/fprintf, sem depender de writetable)
% -------------------------------------------------------------------
caminho_saida = fullfile(pasta_dados, 'Tabela_Corrente_Critica.csv');
fid = fopen(caminho_saida, 'w');
fprintf(fid, 'Fase,Zp_ohm,Ic_kA,P_Ic_pct,Suspeito\n');
for i = 1:numel(resultados)
    if isnan(resultados(i).Ic_kA)
        ic_str = '';
    else
        ic_str = sprintf('%.3f', resultados(i).Ic_kA);
    end
    if isnan(resultados(i).P_Ic_pct)
        p_str = '';
    else
        p_str = sprintf('%.3f', resultados(i).P_Ic_pct);
    end
    fprintf(fid, '%s,%.4f,%s,%s,%s\n', resultados(i).Fase, resultados(i).Zp_ohm, ...
            ic_str, p_str, resultados(i).Suspeito);
end
fclose(fid);

fprintf('\nTabela salva em: %s\n', caminho_saida);
disp('Processamento concluído!');
