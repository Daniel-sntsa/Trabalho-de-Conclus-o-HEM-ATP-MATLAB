function [Ic, VoltageFactor] = calcular_Ic(V_Insulator_V, t_us, tmax_us, ...
    CurrentPeak_kA, CFO_kV, dt_us)
% CALCULAR_IC Determina a corrente critica Ic (kA) para uma dada forma de
% onda de tensao de isolador, usando o metodo do Efeito Disruptivo (DE)
% de Hileman (1999). Busca linear crescente em n (fator de escala da
% tensao), igual a logica do script original em Python.
%
% Entradas:
%   V_Insulator_V : tensao do isolador (V), com polaridade ja alinhada
%   t_us          : vetor de tempo (us)
%   tmax_us       : tempo maximo de integracao (us)
%   CurrentPeak_kA: pico da corrente injetada na simulacao (kA)
%   CFO_kV        : tensao critica de flashover (kV)
%   dt_us         : passo de tempo (us)
%
% Saidas:
%   Ic            : corrente critica (kA), ou NaN se nao convergir
%   VoltageFactor : fator de escala de tensao correspondente, ou NaN

Ic = NaN;
VoltageFactor = NaN;

DE_B = 1.1506 * (CFO_kV ^ 1.36);
V0 = 0.77 * CFO_kV; % kV

% Recorta a serie ate tmax
idx_max = find(t_us <= tmax_us, 1, 'last');
if isempty(idx_max)
    return;
end
V_kV = V_Insulator_V(1:idx_max) / 1000; % V -> kV

VoltagePeak = max(V_kV);
if VoltagePeak <= 0
    return; % onda sem excursao positiva utilizavel
end

n_max = VoltagePeak * 3; % limite de seguranca para evitar loop infinito
n = 1.0;
passo = 0.5;

while true
    VoltageFactor = n / VoltagePeak;
    V_Insul = V_kV * VoltageFactor;

    Vacima = curve_above(V_Insul, V0);
    if ~isempty(Vacima)
        diffv = max(Vacima - V0, 0); % evita valor negativo por erro de ponto flutuante
        Vp = diffv .^ 1.36;
        Area = trapz(Vp) * dt_us; % integracao trapezoidal (grade uniforme)
        if Area >= DE_B
            break;
        end
    end

    n = n + passo;
    if n > n_max
        Ic = NaN;
        VoltageFactor = NaN;
        return; % nao convergiu dentro do limite
    end
end

Ic = VoltageFactor * CurrentPeak_kA;

end
