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
3. **MATLAB** — aplicação do Método do Efeito Disruptivo (Disruptive
   Effect, DE) para determinação da corrente crítica $I_C$ e da
   probabilidade de excedência associada ao risco de *backflashover*.

## Estrutura do repositório

```
├── julia/           # Rotinas de solução do HEM e extração de Zp
├── atpdraw/         # Arquivos de simulação (.acp/.atp) do estudo de caso
├── matlab/          # Rotina de aplicação do Método do Efeito Disruptivo
├── data/            # Parâmetros de entrada (solo, geometria, corrente de descarga)
└── README.md
```

## Requisitos

- Julia ≥ 1.x (pacotes: [listar])
- ATPDraw / ATP-EMTP
- MATLAB ≥ R20XX

## Como citar

Se utilizar este repositório em trabalhos acadêmicos, cite:

```bibtex
@Misc{amador2026codigo,
  author       = {Amador, Daniel},
  title        = {Códigos para Análise de Risco de Backflashover em Linhas de Transmissão via Acoplamento HEM-ATP},
  year         = {2026},
  howpublished = {\url{https://github.com/seu-usuario/seu-repositorio}},
}
```

## Créditos

A implementação original do núcleo do HEM em Julia foi desenvolvida por
[nome, conforme acordado]. Os demais scripts de acoplamento, os arquivos
de simulação do ATPDraw e a rotina do Método do Efeito Disruptivo foram
desenvolvidos pelo autor deste repositório.

## Licença

[MIT / GPL-3.0 / outra — a definir]
