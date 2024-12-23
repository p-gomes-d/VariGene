# VariGene

O VariGene é uma *pipeline* bioinformática desenhada para identificar variantes genéticas a partir de dados de sequenciação genómica. Esta está programada para analisar uma amostra-exemplo - *AMOSTRA_A* - de sequenciação de regiões-alvo segundo a tecnologia Ion AmpliSeq (Thermo Fisher Scientific).

A *pipeline* está configurada para poder ser executada no *Docker* para uma fácil utilização, partilha e reprodutibilidade. Ao ser executada, esta criará um ambiente virtual - *container* - contendo todas os programas e requisitos necessários para realizar a análise. Os resultados serão exportados para este mesmo *container* a partir do qual o utilizador poderá extraí-los conforme desejado.

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

Confirmar que a imagem "vari-gene" se encontra instalada

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

Usando os ficheiros disponibilizados da *AMOSTRA_A*:

  3.1 Apresente métricas de análise que possam ser relevantes em alguns passos da *pipeline*.

  R.: Utilizando os ficheiros disponibilizados, é relevante reunir métricas que permitam averiguar a qualidade dos ficheiros de partida e por conseguinte a qualidade da análise no geral. Analisando o alinhamento `AMOSTRA_A.bam` fornecido, é possível aferir a qualidade do mapeamento.
  Executanto o seguinte comando é possível aferir o número total de reads e o número de reads mapeadas.
  
  ```bash
  samtools flagstat AMOSTRA_A.bam
  ```
  
  Num total de 802979 reads 100% ficaram mapeadas.
  
  Para além da qualidade do mapeamento também é importante averiguar a qualidade da anotação das variantes. Para tal podemos analisar o ficheiro anotado `AMOSTRA_A.anotada.vcf` e executar o seguintes comandos:
  
  ```bash
  bcftools view -H AMOSTRA_A.anotada.vcf | wc -l
  ```
  
  Um total de 83 variantes foram identificadas.

  ```bash
  bcftools view -H AMOSTRA_A.anotada.vcf | grep "CLINVARID" | wc -l
  ```
  
  Das quais 80 ficaram anotadas.
  
  Saber quantos genes foram identificadas variantes e quantas variantes existem por gene pode ser útil e pode ser executado com o seguinte comando:
  
  ```bash
  bcftools query -f "%INFO/GENE\n" AMOSTRA_A.anotada.vcf | grep -v "^.$" | cut -d ":" -f 1 | uniq -c
  ```
  
  Foram encontrados 23 genes com variantes, sendo que o gene com o maior número de variantes identificadas foi o APC, com 11 variantes.

  3.2 Após a filtragem das variantes por qualidade, proponha os passos necessários para a identificação de variantes de interesse clínico.

  R.: As variantes de interesse clínico acrescido são aquelas que estão associadas a uma determinada patologia. Para identificar estas variantes no ficheiro anotado `AMOSTRA_A.anotada.vcf` basta filtrar executando o seguinte commando:

   ```bash
  bcftools view -H AMOSTRA_A.anotada.vcf | grep "CLASSIFICATION=Pathogenic"
  ```
  
  Apenas uma variante patogénica foi identificada, no gene NTHL1.
  
      
