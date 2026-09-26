# Análise de Risco de Backflashover em Linhas de Transmissão de 230 kV

Códigos e materiais computacionais desenvolvidos para o Trabalho de
Conclusão de Curso *"Risco de Backflashover em Linhas de Transmissão de
230 kV: Influência da Resistividade do Solo e do Dimensionamento do
Aterramento via Modelagem Computacional HEM-ATP"*, apresentado à Escola
Politécnica da UFRJ (Engenharia Elétrica), sob orientação do Prof.
Antonio Carlos Siqueira de Lima.

## Visão geral

O estudo avalia o desempenho de uma linha de transmissão real de 230 kV
frente a descargas atmosféricas, integrando três ferramentas
computacionais em um fluxo único:

1. **Julia** — solução das equações integrais de campo do Modelo
   Eletromagnético Híbrido (*Hybrid Electromagnetic Model* — HEM) para
   obtenção da resposta transitória do sistema de aterramento e extração
   da impedância impulsiva ($Z_P = V_P/I_P$).
2. **ATPDraw** — simulação eletromagnética transitória da torre, dos
   condutores de fase e para-raios e da fonte de corrente de descarga,
   com o aterramento representado por um resistor concentrado equivalente
   a $Z_P$.
3. **MATLAB** — aplicação do Método do Efeito Disruptivo (*Disruptive
   Effect*, DE) para determinação da corrente crítica $I_C$ e da
   probabilidade de excedência associada ao risco de *backflashover*.

## Implementação do HEM utilizada

O núcleo do HEM empregado neste trabalho **não foi desenvolvido pelo autor
deste repositório**. Utiliza-se a implementação em Julia puro do Modelo
Eletromagnético Híbrido desenvolvida por **Pedro Henrique Nascimento
Vieira**, disponível no repositório
[TAGS-julia](https://github.com/pedrohnv/transient-analysis-grounding-systems-julia)
(*Transient Analysis of Grounding Systems*), licenciado sob GPL-3.0.
Segundo o próprio autor, trata-se de um protótipo; uma versão em C, de
alto desempenho, está disponível em
[transient-analysis-grounding-systems](https://github.com/pedrohnv/transient-analysis-grounding-systems)
(DOI: [10.5281/zenodo.2644010](https://doi.org/10.5281/zenodo.2644010)).

Os arquivos `hem.jl` e `counterpoises.jl` pertencem a esse repositório e
**não são redistribuídos aqui**. O script `TCCHEM_LT.jl` deste repositório
os utiliza como dependências, acrescentando a modelagem do contrapeso do
estudo de caso, a excitação pela onda de descarga MCS_FST#2, a varredura
em frequência, a extração de $Z_P$ e a exportação dos resultados.

## Estrutura do repositório

```
├── julia/           # Script de acoplamento ao HEM e extração de Zp (TCCHEM_LT.jl)
├── atpdraw/         # Arquivos de simulação (.acp/.atp) do estudo de caso
├── matlab/          # Rotina de aplicação do Método do Efeito Disruptivo
├── data/            # Parâmetros de entrada (solo, geometria, corrente de descarga)
└── README.md
```

## Requisitos

- **Julia** ≥ 1.6, com os pacotes `DelimitedFiles`, `LinearAlgebra`,
  `FFTW` e `Printf`, além das dependências do TAGS-julia
- **ATPDraw / ATP-EMTP**
- **MATLAB**

## Como executar a etapa em Julia

1. Clone o repositório TAGS-julia:
   ```bash
   git clone https://github.com/pedrohnv/transient-analysis-grounding-systems-julia.git
   ```
2. Copie `julia/TCCHEM_LT.jl` para uma subpasta do repositório clonado
   (por exemplo, `examples/`), já que o script carrega `../hem.jl` e
   `../counterpoises.jl`.
3. Ajuste os parâmetros no bloco **"PAINEL CENTRAL DE CONTROLE DE
   PARÂMETROS"** do script. O comprimento do ramal `L` define a etapa de
   extensão do contrapeso (25, 50, 75 ou 100 m).
4. Execute com múltiplas threads (o solver é paralelizado por frequência):
   ```bash
   julia -t auto TCCHEM_LT.jl
   ```
5. São gerados `Fase_<X>_Tempo.csv` (GPR no tempo para cada solo) e
   `Fase_<X>_Frequencia.csv` (Z(ω) para cada solo), além de uma tabela de
   resumo no terminal com $R_{dc}$ e $Z_P$. Os CSV usam `;` como separador
   de campo e `,` como separador decimal.

## Como citar

Se utilizar este repositório em trabalhos acadêmicos, cite-o:

```bibtex
@Misc{amador2026codigo,
  author       = {Amador, Daniel},
  title        = {Códigos para Análise de Risco de Backflashover em Linhas de Transmissão via Acoplamento HEM-ATP},
  year         = {2026},
  howpublished = {\url{https://github.com/seu-usuario/seu-repositorio}},
}
```

Cite também a implementação do HEM utilizada:

```bibtex
@Misc{vieira_tags,
  author       = {Vieira, Pedro Henrique Nascimento},
  title        = {{TAGS}: Transient Analysis of Grounding Systems},
  howpublished = {\url{https://github.com/pedrohnv/transient-analysis-grounding-systems}},
  doi          = {10.5281/zenodo.2644010},
}
```

## Licença

Este repositório é distribuído sob a licença **GPL-3.0**, em
compatibilidade com a licença do TAGS-julia, do qual depende.
