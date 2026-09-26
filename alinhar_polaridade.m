function Vout = alinhar_polaridade(V)
% ALINHAR_POLARIDADE Inverte o sinal da serie inteira se ela for
% negativo-dominante. O metodo DE (Hileman) assume onda de polaridade
% positiva; inverter o sinal preserva a forma da onda (ao contrario de
% usar abs() ponto a ponto, que distorceria oscilacoes que cruzam zero).

if abs(min(V)) > abs(max(V))
    Vout = -V;
else
    Vout = V;
end

end
