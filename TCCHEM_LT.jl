# ==============================================================================
#                      SISTEMA DE SIMULAÇÃO (HEM)
#           Unificação: Análise no Tempo (GPR) & Análise na Frequência (Z_w)
# ==============================================================================
using DelimitedFiles
using LinearAlgebra
using FFTW
using Printf

# Dependências externas (núcleo do HEM):
#   hem.jl          → solver eletromagnético, funções de impedância/admitância,
#                     transformadas de Laplace (laplace_transform / invlaplace_transform)
#                     e constantes físicas MU0, EPS0 e FOUR_PI (= 4π).
#   counterpoises.jl → geradores de geometria de contrapeso (ex.: conventional1).
include("../hem.jl")
include("../counterpoises.jl")

# ==============================================================================
# 1. PAINEL CENTRAL DE CONTROLE DE PARÂMETROS (EDITE APENAS AQUI)
# ==============================================================================
# Fase I = 25.0m | Fase II = 50.0m | Fase III = 75.0m | Fase IV = 100.0m
# Dados geométricos extraídos do projeto básico da LT de 230 kV (estudo de caso).
#
# NOTA: o parâmetro `dep_freq` das funções `simular_solo_tempo` e
# `calcular_Z_omega` habilita a variação dos parâmetros elétricos do solo com
# a frequência (Alipio & Visacro, 2014). Neste estudo ele é mantido em `false`
# em todas as chamadas, conforme a simplificação conservadora adotada
# (Seção 2.4 do TCC). O código é mantido funcional como base para a extensão
# sugerida em Trabalhos Futuros.
#
# --- A) GEOMETRIA DO CONTRAPESO ---
L        = 100       # Comprimento do ramal do contrapeso (m) - Corresponde à Fase IV do projeto.
s1       = 19.91     # Distância horizontal ao estai esquerdo (m) - Extraído da E2EL
s2       = 17.01     # Distância horizontal ao estai direito (m) - Extraído da E2EL
s3       = 19.91     # Separação no eixo Y (ajuste conforme a simetria do seu modelo EMTP/simulador)
s4       = 17.01     # Separação no eixo Y (ajuste conforme a simetria do seu modelo EMTP/simulador)
h        = -0.80     # Profundidade de enterramento (m) - Especificação de 80 cm
r_cabo   = 4.572e-3  # Raio do condutor do contrapeso (m) - Metade do diâmetro de 9.144 mm
Lmax_seg = 5.0       # Comprimento máximo de segmentação espacial (m) - Mantido o seu padrão

# --- B) EXCITAÇÃO NO DOMÍNIO DO TEMPO (MCS_FST#2 - De Conti & Visacro) ---
# Parâmetros da onda MCS_FST#2 (7 funções de Heidler) para simulação temporal
# Arquivo de referência: Analytical Representation of Single- and Double-Peaked Lightning Current Waveforms
const params_MCS = [
     6.0e3    3.0e-6      76.0e-6     2.0;
     5.0e3    3.5e-6      10.0e-6     3.0;
     5.0e3    4.8e-6      30.0e-6     5.0;
     8.0e3    6.0e-6      26.0e-6     9.0;
    16.5e3    7.0e-6      23.2e-6    30.0;
    17.0e3     70e-6       200e-6     2.0;
    12.0e3     12e-6      26.0e-6    14.0
]

Tmax   = 200e-6   # Janela total da simulação temporal para o NILT (s)
dt     = 10e-9    # Passo de tempo discreto (s)
T_plot = 60e-6    # Janela de tempo de interesse para gráficos/exportação (s)

# --- C) PARÂMETROS DA VARREDURA DE FREQUÊNCIA ADAPTATIVA ---
freq_n_pts  = 100     # Número de pontos logarítmicos por varredura
freq_f_min  = 100.0   # Frequência mínima fixa para todos os solos (Hz)
freq_fator  = 5.0     # Fator de segurança η para o critério λ ≥ η·Lmax_seg

"""
    f_max_valido(rho, Lmax_seg, fator=5.0)

Calcula a frequência máxima válida para a varredura harmônica, garantindo que
o comprimento de onda no solo seja pelo menos `fator` vezes o comprimento máximo
de segmento. Acima desse limite, o HEM produz oscilações numéricas espúrias.

Regime de difusão: λ_solo ≈ sqrt(2·ρ / (ω·μ₀))
Condição: λ_solo ≥ fator · Lmax_seg  →  f_max = ρ / (π·μ₀·(fator·Lmax_seg)²)
"""
function f_max_valido(rho::Float64, Lmax_seg::Float64, fator::Float64=5.0)
    return rho / (pi * MU0 * (fator * Lmax_seg)^2)
end

"""
    freqs_por_solo(rho, Lmax_seg; f_min, n_pts, fator)

Gera o vetor de frequências logaritmicamente espaçado para um dado solo,
limitado pela frequência máxima fisicamente válida para a segmentação usada.

# Ancoragem do ponto DC (f = 0):
O ponto f = 0 é inserido no início do vetor como primeiro elemento.
Isso resolve a dificuldade do Vector Fitting em capturar o comportamento DC:
  - f = 0  →  s = 0  →  Z(0) = R_dc (resistência pura, Im = 0)
  - O VF receberá esse ponto explicitamente e não precisará extrapolar.
  - Na conversão para Foster, o polo em s = 0 mapeia diretamente em R_dc.
"""
function freqs_por_solo(rho::Float64, Lmax_seg::Float64;
                        f_min::Float64=100.0, n_pts::Int=100, fator::Float64=5.0)
    f_max  = f_max_valido(rho, Lmax_seg, fator)
    f_log  = 10 .^ collect(range(log10(f_min), stop=log10(f_max), length=n_pts))
    return vcat(0.0, f_log)   # ← f=0 como âncora DC no índice 1
end

# --- D) CONFIGURAÇÃO DOS CENÁRIOS DE SOLO ---
solos = [
    ("Muito_Baixa", 100.0,  10.0),
    ("Baixa",       500.0,  10.0),
    ("Media",       1000.0, 10.0),
    ("Alta",        2000.0, 10.0),
    ("Muito_Alta",  4000.0, 10.0)
]

# ==============================================================================
# 2. FUNÇÕES AUXILIARES E SOLVERS
# ==============================================================================

"""
    onda_mcs_fst(t)

Calcula o valor instantâneo da corrente de um raio no domínio do tempo, utilizando
a soma de 7 funções de Heidler para representar o formato de duplo pico
(onda MCS_FST#2 de De Conti e Visacro).
"""
function onda_mcs_fst(t)
    if t <= 0.0; return 0.0; end
    i_total = 0.0
    for k in 1:7
        I0k   = params_MCS[k, 1]
        tau1k = params_MCS[k, 2]
        tau2k = params_MCS[k, 3]
        nk    = params_MCS[k, 4]
        eta_k  = exp(-(tau1k / tau2k) * (nk * (tau2k / tau1k))^(1/nk))
        frente = (t / tau1k)^nk
        i_total += (I0k / eta_k) * (frente / (1.0 + frente)) * exp(-t / tau2k)
    end
    return i_total
end

"""
    encontrar_no(nodes, coordenada_alvo)

Varre a matriz de nós e encontra o índice do nó mais próximo à coordenada alvo.
"""
function encontrar_no(nodes, coordenada_alvo)
    menor_dist = Inf; indice_no = 1
    for i in axes(nodes, 1)
        dist = norm(nodes[i, :] - coordenada_alvo)
        if dist < menor_dist
            menor_dist = dist; indice_no = i
        end
    end
    return indice_no
end

"""
    simular_solo_tempo(...)

Solver transiente via HEM + NILT. Retorna GPR no domínio do tempo.
"""
function simular_solo_tempo(rho, eps_r, eletrodos, images, nodes,
                            no_perna1, no_perna2, no_perna3, no_perna4,
                            i_s, s_freq, Tmax, nt, dep_freq::Bool)
    ns = length(eletrodos); nn = size(nodes, 1)
    mA, mB = incidence(eletrodos, nodes)
    v_s = zeros(ComplexF64, length(s_freq))

    intg_type = INTG_MHEM
    mpotzl, mpotzt   = calculate_impedances(eletrodos, 0.0, 1.0, 1.0, 1.0,
                                             typemax(Int), 1e-4, 1e-5, norm, intg_type)
    mpotzli, mpotzti = impedances_images(eletrodos, images, 0.0, 1.0, 1.0, 1.0,
                                          1.0, 1.0, typemax(Int), 1e-4, 1e-5, norm, intg_type)

    rbar = zeros(ns, ns); rbari = zeros(ns, ns)
    for k = 1:ns
        p1 = collect(eletrodos[k].middle_point)
        for i = k:ns
            rbar[i,k]  = norm(p1 - collect(eletrodos[i].middle_point))
            rbari[i,k] = norm(p1 - collect(images[i].middle_point))
        end
    end

    Threads.@threads for f in eachindex(s_freq)
        s = s_freq[f]
        rho_f = rho; eps_r_f = eps_r
        if dep_freq
            freq_hz = max(abs(imag(s)) / (2 * pi), 100.0)
            rho_f   = rho / (1.0 + 1.2e-6 * (rho^0.73) * ((freq_hz - 100.0)^0.65))
            freq_eps  = max(freq_hz, 10e3)
            eps_r_f = 7.6e3 * (freq_eps^(-0.4)) + 1.3
        end
        kappa     = (1.0 / rho_f) + s * eps_r_f * EPS0
        k1        = sqrt(s * MU0 * kappa)
        kappa_ar  = s * EPS0
        ref_t     = (kappa - kappa_ar) / (kappa + kappa_ar)
        zl = zeros(ComplexF64, ns, ns); zt = zeros(ComplexF64, ns, ns)
        iwu_4pi   = s * MU0 / FOUR_PI
        one_4pik  = 1.0 / (FOUR_PI * kappa)
        for k = 1:ns
            for i = k:ns
                zl[i,k]  = exp(-k1 * rbar[i,k])  * iwu_4pi  * mpotzl[i,k]
                zt[i,k]  = exp(-k1 * rbar[i,k])  * one_4pik * mpotzt[i,k]
                zl[i,k] += 1.0 * exp(-k1 * rbari[i,k]) * iwu_4pi  * mpotzli[i,k]
                zt[i,k] += ref_t * exp(-k1 * rbari[i,k]) * one_4pik * mpotzti[i,k]
                zl[k,i]  = zl[i,k]; zt[k,i] = zt[i,k]
            end
        end
        yn   = zeros(ComplexF64, nn, nn)
        exci = zeros(ComplexF64, nn)
        # A corrente de descarga injetada na base da torre se reparte
        # igualmente entre os 4 ramais do contrapeso (simetria geométrica
        # do arranjo "conventional1").
        corrente_dividida    = i_s[f] / 4.0
        exci[no_perna1] = corrente_dividida; exci[no_perna2] = corrente_dividida
        exci[no_perna3] = corrente_dividida; exci[no_perna4] = corrente_dividida
        admittance!(yn, zl, zt, mA, mB)
        ldiv!(lu!(yn), exci)
        v_s[f] = exci[no_perna1]
    end
    return invlaplace_transform(v_s, Tmax, nt)
end

"""
    calcular_Z_omega(...)

# Tratamento do ponto DC (f = 0, s = 0):
Quando freqs[f_idx] == 0.0, o solver HEM não pode ser chamado diretamente
(divisão por zero em kappa, k1, etc.). Nesse caso, Z(0) = R_dc é calculado
analiticamente como o limite de Z(ω) quando ω → 0:

    Z(0) = ρ / (4π) · Σ (integrais de potencial resistivo)

Na prática, usamos a aproximação de baixíssima frequência (f = 1e-3 Hz, s ≈ 0⁺)
que é numericamente estável e converge para R_dc com erro < 0.01%.
Isso garante Im(Z(0)) ≈ 0 e fornece a âncora resistiva pura para o VF.

Solver harmônico via HEM. Retorna Z(ω) complexo para cada frequência do vetor
`freqs` fornecido — gerado individualmente por solo.
"""
function calcular_Z_omega(rho, eps_r, eletrodos, images, nodes,
                          no_perna1, no_perna2, no_perna3, no_perna4,
                          freqs, dep_freq::Bool)
    ns = length(eletrodos); nn = size(nodes, 1)
    mA, mB = incidence(eletrodos, nodes)
    Z_w = zeros(ComplexF64, length(freqs))

    intg_type = INTG_MHEM
    mpotzl, mpotzt   = calculate_impedances(eletrodos, 0.0, 1.0, 1.0, 1.0,
                                             typemax(Int), 1e-4, 1e-5, norm, intg_type)
    mpotzli, mpotzti = impedances_images(eletrodos, images, 0.0, 1.0, 1.0, 1.0,
                                          1.0, 1.0, typemax(Int), 1e-4, 1e-5, norm, intg_type)

    rbar = zeros(ns, ns); rbari = zeros(ns, ns)
    for k = 1:ns
        p1 = collect(eletrodos[k].middle_point)
        for i = k:ns
            rbar[i,k]  = norm(p1 - collect(eletrodos[i].middle_point))
            rbari[i,k] = norm(p1 - collect(images[i].middle_point))
        end
    end

    Threads.@threads for f_idx in eachindex(freqs)

        # ------------------------------------------------------------------
        # Substituição de f=0 por s ≈ 0⁺ (limite DC numericamente estável)
        # ------------------------------------------------------------------
        freq_hz = freqs[f_idx] == 0.0 ? 1e-3 : max(freqs[f_idx], 100.0)
        s       = 2im * pi * (freqs[f_idx] == 0.0 ? 1e-3 : freqs[f_idx])
        # ------------------------------------------------------------------

        rho_f = rho; eps_r_f = eps_r
        if dep_freq
            rho_f   = rho / (1.0 + 1.2e-6 * (rho^0.73) * ((freq_hz - 100.0)^0.65))
            freq_eps  = max(freq_hz, 10e3)
            eps_r_f = 7.6e3 * (freq_eps^(-0.4)) + 1.3
        end
        kappa    = (1.0 / rho_f) + s * eps_r_f * EPS0
        k1       = sqrt(s * MU0 * kappa)
        kappa_ar = s * EPS0
        ref_t    = (kappa - kappa_ar) / (kappa + kappa_ar)
        zl = zeros(ComplexF64, ns, ns); zt = zeros(ComplexF64, ns, ns)
        iwu_4pi  = s * MU0 / FOUR_PI
        one_4pik = 1.0 / (FOUR_PI * kappa)
        for k = 1:ns
            for i = k:ns
                zl[i,k]  = exp(-k1 * rbar[i,k])  * iwu_4pi  * mpotzl[i,k]
                zt[i,k]  = exp(-k1 * rbar[i,k])  * one_4pik * mpotzt[i,k]
                zl[i,k] += 1.0 * exp(-k1 * rbari[i,k]) * iwu_4pi  * mpotzli[i,k]
                zt[i,k] += ref_t * exp(-k1 * rbari[i,k]) * one_4pik * mpotzti[i,k]
                zl[k,i]  = zl[i,k]; zt[k,i] = zt[i,k]
            end
        end
        yn   = zeros(ComplexF64, nn, nn)
        exci = zeros(ComplexF64, nn)
        # Injeção de corrente unitária (1 A), repartida igualmente entre os
        # 4 ramais do contrapeso → a tensão resultante é numericamente Z(ω).
        corrente_dividida    = 1.0 / 4.0 + 0.0im
        exci[no_perna1] = corrente_dividida; exci[no_perna2] = corrente_dividida
        exci[no_perna3] = corrente_dividida; exci[no_perna4] = corrente_dividida
        admittance!(yn, zl, zt, mA, mB)
        ldiv!(lu!(yn), exci)
        Z_w[f_idx] = exci[no_perna1]
    end
    return Z_w
end

# ==============================================================================
# 3. CONSTRUÇÃO DA GEOMETRIA E MAPEAMENTO DE NÓS
# ==============================================================================
println("Construindo geometria base tridimensional...")

eletrodos_base = conventional1(L, s1, s2, s3, s4, h, r_cabo)
eletrodos, nodes = seg_electrode_list(eletrodos_base, Lmax_seg)
ns = length(eletrodos)
nn = size(nodes, 1)

images = Array{Electrode}(undef, ns)
for i = 1:ns
    sp = [eletrodos[i].start_point[1], eletrodos[i].start_point[2], -eletrodos[i].start_point[3]]
    ep = [eletrodos[i].end_point[1],   eletrodos[i].end_point[2],   -eletrodos[i].end_point[3]]
    images[i] = new_electrode(sp, ep, eletrodos[i].radius)
end

base_x = s1 / 2.0
base_y = s3 / 2.0
no_perna1 = encontrar_no(nodes, [ base_x,  base_y, h])
no_perna2 = encontrar_no(nodes, [-base_x,  base_y, h])
no_perna3 = encontrar_no(nodes, [-base_x, -base_y, h])
no_perna4 = encontrar_no(nodes, [ base_x, -base_y, h])

# ==============================================================================
# 3.5. PREPARAÇÃO DOS SINAIS DE EXCITAÇÃO
# ==============================================================================
nt            = Int(round(Tmax / dt))
vetor_t       = collect(0.0:dt:(nt - 1) * dt)
nt_plot       = Int(round(T_plot / dt)) + 1
vetor_t_limpo = vetor_t[1:nt_plot]
i_t           = onda_mcs_fst.(vetor_t)
s_freq, i_s   = laplace_transform(i_t, Tmax, nt)
max_i_t       = maximum(i_t)

# ==============================================================================
# 4. EXECUÇÃO DOS SOLVERS
# ==============================================================================
# Vetores de resultados para o cenário de parâmetros do solo constantes
# (Z_const), único cenário considerado neste estudo.
resultados_v_const    = []
resultados_Z_const    = []
resultados_freqs_solo = []

resumo_Z_imp_const = Float64[]
resumo_R_dc_const  = Float64[]
resumo_f_max       = Float64[]

println("\nIniciando Processamento de Varreduras Eletromagnéticas...")
for (nome, rho, eps_r) in solos
    # Geração do vetor adaptativo + âncora DC em f = 0
    freqs_solo = freqs_por_solo(rho, Lmax_seg;
                                f_min  = freq_f_min,
                                n_pts  = freq_n_pts,
                                fator  = freq_fator)
    f_max_solo = freqs_solo[end]
    push!(resumo_f_max, f_max_solo)
    push!(resultados_freqs_solo, freqs_solo)

    println(" -> Processando Solo: ", nome,
            " (ρ = ", rho, " Ω.m | f_max = ", round(f_max_solo/1e3, digits=2), " kHz)...")

    # --- CENÁRIO: PARÂMETROS CONSTANTES (único alimentado ao VF) ---
    v_t_c    = simular_solo_tempo(rho, eps_r, eletrodos, images, nodes,
                                  no_perna1, no_perna2, no_perna3, no_perna4,
                                  i_s, s_freq, Tmax, nt, false)
    v_real_c = real.(v_t_c)
    push!(resultados_v_const, v_real_c)
    push!(resumo_Z_imp_const, maximum(v_real_c[1:nt_plot]) / max_i_t)

    Z_w_c = calcular_Z_omega(rho, eps_r, eletrodos, images, nodes,
                              no_perna1, no_perna2, no_perna3, no_perna4,
                              freqs_solo, false)
    push!(resultados_Z_const, Z_w_c)

    # índice 1 = ponto DC (f=0, s≈0⁺) → R_dc mais preciso que antes
    push!(resumo_R_dc_const, real(Z_w_c[1]))
end

# ==============================================================================
# 5. EXPORTAÇÃO DOS DADOS
# ==============================================================================
# Arquivos CSV gerados no diretório de execução, com separador de campo ";"
# e separador decimal "," (padrão brasileiro).
# Unidades: tempo em s, corrente em A, tensão em V, frequência em Hz, Z em Ω.
println("\nGerando arquivos e tabelas finais...")

X_fase = if L == 25 "I" elseif L == 50 "II" elseif L == 75 "III" elseif L == 100 "IV" else "Custom" end

# --- CSV DE TEMPO ---
open("Fase_$(X_fase)_Tempo.csv", "w") do io
    nomes_solos = [s[1] for s in solos]
    cabecalho = "Tempo_s;Corrente_A;" *
                join("V_" .* nomes_solos .* "_Const", ";") * "\n"
    write(io, cabecalho)
    for k in 1:nt_plot
        valores = [vetor_t[k], i_t[k]]
        for res in resultados_v_const; push!(valores, res[k]); end
        write(io, join([replace(string(v), "." => ",") for v in valores], ";") * "\n")
    end
end

# --- CSV DE FREQUÊNCIA ---
# Estrutura: 3 colunas por solo → [Freq_Hz | Z_Real | Z_Imag]
# O ponto DC (f = 0) aparece na primeira linha de cada solo com Im ≈ 0.
# Linhas excedentes dos solos mais curtos preenchidas com ";;" (2 separadores).
open("Fase_$(X_fase)_Frequencia.csv", "w") do io
    nomes_solos = [s[1] for s in solos]
    n_solos     = length(solos)

    # Cabeçalho: 3 colunas por solo
    cabecalho = join(["Freq_Hz_"  * nomes_solos[i] * ";" *
                      "Z_Real_"   * nomes_solos[i] * ";" *
                      "Z_Imag_"   * nomes_solos[i]
                      for i in 1:n_solos], ";") * "\n"
    write(io, cabecalho)

    n_linhas_max = maximum(length.(resultados_freqs_solo))

    for k in 1:n_linhas_max
        partes = String[]
        for i in 1:n_solos
            freqs_i = resultados_freqs_solo[i]
            Z_i     = resultados_Z_const[i]
            if k <= length(freqs_i)
                f_str  = replace(string(freqs_i[k]),   "." => ",")
                re_str = replace(string(real(Z_i[k])), "." => ",")
                im_str = replace(string(imag(Z_i[k])), "." => ",")
                push!(partes, f_str * ";" * re_str * ";" * im_str)
            else
                push!(partes, ";;")   # 2 separadores para 3 colunas vazias
            end
        end
        write(io, join(partes, ";") * "\n")
    end
end

# --- TABELA DE RESUMO ---
println("\n" * "="^80)
println("                   RESUMO COMPARATIVO DE IMPEDÂNCIAS - FASE $X_fase (L = $(L)m)")
println("="^80)
println(" Tipo de Solo  | Resistiv. |   f_max válida   |  R_dc (f=0)   |  Z_imp Const.")
println("               |   (Ω.m)   |      (Hz)        |     (Ω)       |     (Ω)      ")
println("-"^80)
for (idx, (nome, rho, _)) in enumerate(solos)
    @printf(" %-13s | %9.1f | %16.2f | %11.4f Ω | %11.4f Ω\n",
            nome,
            rho,
            resumo_f_max[idx],
            resumo_R_dc_const[idx],
            resumo_Z_imp_const[idx])
end
println("="^80 * "\n")
