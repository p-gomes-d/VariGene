# VariGene

O VariGene é uma *pipeline* bioinformática desenhada para identificar variantes genéticas a partir de dados de sequenciação genómica. Esta está programada para analisar uma amostra-exemplo - *AMOSTRA_A* - de sequenciação de regiões-alvo segundo a tecnologia Ion AmpliSeq (Thermo Fisher Scientific).

A *pipeline* está configurada para poder ser executada no *Docker* para uma fácil utilização, partilha e reprodutibilidade. Ao ser executada, esta criará um ambiente virtual - *container* - contendo todas os programas e requisitos necessários para realizar a análise. Os resultados serão exportados para este mesmo *container* a partir do qual o utilizador poderá extraí-los conforme deseja.

## Como utilizar o VariGene através do *Docker*

Este repositório contém todos os ficheiros necessários para executar a *pipeline*. O utilizador deverá reunir os seguintes ficheiros num único diretório de trabalho:

  1. O alinhamento `AMOSTRA_A.bam` que será o ponto de partida da análise;
  2. O *script* `VariGene.sh` que contém todas as instruções necessárias para realizar a análise;
  3. O *script* `report.sh` que no fim da análise produzirá um relatório clínico dos resultados obtidos; e
  4. O ficheiro `Dockerfile` necessário para implementar o VariGene no *Docker*.

NOTA: para conhecer melhor quais os programas utilizados na análise consultar o ficheiro `pseudocodigo.txt`.

NOTA: O sistema-base do VariGene é o Ubuntu (versão 22.04).

AVISO: Antes de inciar a *pipeline* o utilizador deve disponibilizar no mínimo 10 GB de espaço uma vez que o VariGene mantém todos os ficheiros de *output* até terminar, permitindo assim inspecionar algum ficheiro caso seja necessário.

Para implementar o VariGene no *Docker* o utilizador deve executar as seguintes instruções na linha de comandos:

```bash
cd /DIRETÓRIO/COM/OS/FICHEIROS
```
Garantir que os ficheiros necessários se encontram no diretório actual

```bash
docker build -t vari-gene .
```
Implementa o VariGene no *Docker* com a *tag* "vari-gene"

```bash
docker images
```
Garantir que a imagem "vari-gene" se encontra instalada

```bash
docker run vari-gene
```
Inicia a análise e cria um *container* com o resultados

```bash
docker ps -a
```
Listar todos os *containers* que estão em execução ou não

```bash
docker cp CONTAINER-ID:/usr/local/FICHEIRO-DE-INTERESSE .
```
Copiar um resultado do *container* para o diretório actual


## Interpretação dos resultados


