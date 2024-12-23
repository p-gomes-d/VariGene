# VariGene

O VariGene é uma *pipeline* bioinformática desenhada para identificar variantes genéticas a partir de dados de sequenciação genómica. Esta está programada para analisar uma amostra-exemplo - *AMOSTRA_A* - de sequenciação de regiões-alvo segundo a tecnologia Ion AmpliSeq (Thermo Fisher Scientific).

A *pipeline* está configurada para poder ser executada no *Docker* para uma fácil utilização, partilha e reprodutibilidade. Ao ser executada, esta criará um ambiente virtual - *container* - contendo todas os programas e requisitos necessários para realizar a análise. Os resultados serão exportados para este mesmo *container* a partir do qual o utilizador poderá extraí-los conforme a necessidade.

## Como utilizar o VariGene através do *Docker*

Este repositório contém todos os ficheiros necessários para executar a *pipeline*. O utilizador deverá reunir os seguintes ficheiros num único diretório de trabalho:
  1. 
