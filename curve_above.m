function Vacima = curve_above(V, V0)
% CURVE_ABOVE Retorna o trecho de V que fica acima do limiar V0 (da
% primeira a ultima amostra que excede V0). Retorna [] se nada excede.

idx = find(V > V0);
if isempty(idx)
    Vacima = [];
else
    Vacima = V(min(idx):max(idx));
end

end
