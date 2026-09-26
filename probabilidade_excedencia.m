function P = probabilidade_excedencia(Ic_kA)
% PROBABILIDADE_EXCEDENCIA Distribuicao de Silveira/Visacro (Morro do
% Cachimbo) para a probabilidade de a corrente de descarga exceder Ic.

P = 1 ./ (1 + (Ic_kA / 43.3) .^ 3.8);

end
