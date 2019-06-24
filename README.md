# BAG2019
**B**ioinformática **A**plicada à **G**enômica (BAG): Introdução à Análise de Dados de Sequenciamento de Nova Geração (Turma SP - 2019)

O desenvolvimento das tecnologias de sequenciamento de nova geração (Next Generation Sequence - NGS) propiciou o sequenciamento global do genoma, exoma ou transcriptoma, revolucionando a área médica. Além disso, trouxe consigo a necessidade do desenvolvimento de ferramentas computacionais para análise do grande volume de dados gerados por estas técnicas. Desta forma, conhecer e saber manipular estes dados com ferramentas de bioinformática tornou-se de grande importância e um diferencial para profissionais da área da saúde atualmente.

Explorando a integração entre as ciências médicas e as ciências computacionais, este curso irá discutir como o NGS tem contribuído para o entendimento das causas genéticas de doenças humanas. Também será abordado, de maneira teórica e prática, o uso das ferramentas de bioinformática necessárias para o processamento dos dados brutos de um sequenciamento, até o seu mapeamento contra um genoma referência e identificação das variantes.

https://www.einstein.br/ensino/atualizacao/bioinformatica_aplicada_a_genomica

# Instrutores

* Dra Andréa Laurato Sertié
* Dra Karina Griesi Oliveira
* Msc Murilo Cervato
* Msc Renato Puga

# Arquivos .FASTQ

* FASTQ [Download Google Drive](https://drive.google.com/open?id=1zAL8GVZuKKmroRe1NKMKtd3L067ZyZBO)

# Instalação

* Mac

```bash
# install homebrew
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# historico completo homebrews
git -C "$(brew --repo homebrew/core)" fetch --unshallow

# install java
brew install java

# install samtools
brew install samtools

# install 
brew install fastqc

# install bwa
brew install bwa

# install wget
brew install wget

# install cutadapt
brew install python
pip3 install --user --upgrade cutadapt

cd /usr/local/bin
sudo ln -s ~/Library/Python/3.7/bin/cutadapt .

# install freebayes
cd
mkdir bioinfo
mkdir bioinfo/app
cd bioinfo/app
git clone --recursive git://github.com/ekg/freebayes.git
make
sudo make install

# install dabases annovar
# NOTA: entreno no site do ANNOVAR com seu e-mail e salve o arquivo no diretorio: ~/bionfo/app/

cd ~/bioinfo/app/
tar -zxvf annovar.latest.tar.gz
cd annovar

# baixar as bases: clinvar e exac03
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/

```

* Windows
	* VirtualBox - Cliente para montar nossa Virtual Machine Linux [VirtualBox](https://www.virtualbox.org)
	* Bio-Linux 8 - [Download .ISO](http://nebc.nerc.ac.uk/downloads/bio-linux-8-latest.iso)

No terminal do Linux, vamos instalar alguns pacotes: (o resto vem instalado).

``` bash
# install freebayes
cd
mkdir bioinfo
mkdir bioinfo/app
cd bioinfo/app
git clone --recursive git://github.com/ekg/freebayes.git
make
sudo make install

# install dabases annovar
# NOTA: entreno no site do ANNOVAR com seu e-mail e salve o arquivo no diretorio: ~/bionfo/app/

cd ~/bioinfo/app/
tar -zxvf annovar.latest.tar.gz
cd annovar

# baixar as bases: clinvar e exac03
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
```


# Bioinformática: Pipelines e Comandos

Olá, bem vindo ao curso de **BIOINFORMÁTICA APLICADA À GENÔMICA: INTRODUÇÃO À ANÁLISE DE DADOS DE SEQUENCIAMENTO DE NOVA GERAÇÃO**, aqui encontrará um guia do curso para executar comandos e pipelines de bioinformática.

# Olá Ambiente de Bioinformática
***Linux:*** O Linux é um sistema operacional, assim como o Windows da Microsoft e o Mac OS da Apple. Ele foi criado pelo finlandês Linus Torvalds, e o nome é a mistura do nome do criador com Unix, um antigo sistema operacional da empresa de mesmo nome. [O que é Linux?](https://www.techtudo.com.br/artigos/noticia/2011/12/o-que-e-linux.html).

***Pipelines:*** Primeiro, o pipeline não é um termo de bioinformática, é na verdade um termo de ciência da computação, envolve encadeamento de processos / threads / funções etc. Em resumo, o resultado de um processo é a entrada de um novo processo na sequência. [A review of bioinformatic pipeline frameworks](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5429012/)

# Qual ambiente vamos utilizar?
Vamos utilizar uma instância na Amazon Cloud com Sistema Operacional Linux Ubuntu 14.04. Os programas para que possamos criar nosso pipeline de análises de sequenciamento de nova geração já foram instalados e serão apenas listados nesse documento. O objetivo é melhorar a experiência de usuários iniciantes em "digitar comandos em um terminal linux" e ser um roteiro para que todos possam se localizar.

# Estrutura de Diretórios
Estrutura de diretórios: sequências, programas e arquivos de referência.

* /bioinfo
	* /bioinfo/app
	* /bioinfo/data
		* /bioinfo/data/fastq
	* /bioinfo/reference


# Programas Instalados
Listamos os programas previamente instalados em nosso ambiente para executar o pipeline de chamada de varinates:

* annovar [Download](http://annovar.openbioinformatics.org/)
* bcftools [Download](https://samtools.github.io/bcftools/bcftools.html)
* bedtools [Download](https://bedtools.readthedocs.io/en/latest/)
* bwa [Download](http://bio-bwa.sourceforge.net/)
* FastQC [Download](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* freebayes [Download](https://github.com/ekg/freebayes)
* samtools [Download](http://samtools.sourceforge.net/)
* cutadapt [Download](https://cutadapt.readthedocs.io/en/stable/guide.html)


# Primeiros Comandos

```bash
# ls para listar arquivos e diretorios
ls

# pwd para informar o diretorio atual
pwd

# cd para entrar e sair de diretorios e voltar para casa
cd resultados/
cd ../
cd
```

Fonte: [Comandos Básicos do Terminal Linux](http://swcarpentry.github.io/shell-novice/)

# Diretório de Resultados

```bash
# volta para casa
cd

# mkdir para criar um diretorio
mkdir resultados

# cd para entrar no diretorio resultados
cd resultados/

# mkdir para criar o diretorio das amostras 003, 017 e 019
mkdir 003 017 019
```

# Download das Referências
Nesta estapa vamos utilizar dois cromossomos humanos (chr13 e chr17), nossos genes de interesse são: BRCA1 e BRCA2. [NCBI BRCA 1 and 2](https://www.ncbi.nlm.nih.gov/books/NBK470239/)

Acessar o site: [Sequence and Annotation Downloads](http://hgdownload.cse.ucsc.edu/downloads.html)


```bash

cd /bioinfo/reference

# wget para fazer download do chr13
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz

# wget para fazer download do chr17 
wget -c http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr17.fa.gz

# cat para concatenar os arquivo do cromossomo 13 e 17 em hg19.fa
zcat chr13.fa.gz chr17.fa.gz > hg19.fa

# rm para deletar os arquivos chr13.fa e chr17.fa
rm chr13.fa chr17.fa
```

# BWA: index reference
Agora, é preciso indexar o arquivo FASTA da referência. Todos os programas de alinhmaneto tem esta etapa que otimiza o processo de alinhamento das nossas sequências na referência.

***NOTA:*** A etapa de index é feita apenas uma vez para cada arquivo de referência:

```bash
# bwa index para gerar o index da referencia hg19.fa
bwa index hg19.fa 
```

# FastQC: Relatório de Controle de Qualidade
Gerar relatório de controle de qualidade com FastQC (Tempo ~10s):

```bash
# fastqc para gerar relatorio de qualidade dos arquivo FASTQ
# opcao (-o ./) diz para salvar o resultado no mesmo arquivo em que o comando esta sendo rodado.

# cd para voltar para a casa
cd

# rodar fastqc e salvar o resultado de cada amostra em seu diretorio
fastqc -o resultados/003/ /bioinfo/data/fastq/003.fastq.gz 
```

FastQC 003 resultado

### Tarefa 01: Repetir o processo para as amostras 017 e 019

~10 min


#cutadapt: localizar e remover adaptadres
O Cutadapt localiza e remove sequências de adaptadores, primers, caudas poly-A e outros tipos de sequência indesejada. Aqui vamos utilizar as funções para "trimar" sequências pequenas e maiores do que o esperado. (Tempo ~3s):

```bash

 -m LEN[:LEN2], --minimum-length=LEN[:LEN2]
                        Discard reads shorter than LEN. Default: 0
 -M LEN[:LEN2], --maximum-length=LEN[:LEN2]
                        Discard reads longer than LEN. Default: no limit
 
# cd para voltar para a casa
cd

# cutadapt para remover sequencias de tamanho menores 100pb e maiores que 220
cutadapt --minimum-length 100 --maximum-length 220 -q 15  -o resultados/003/003.cutadapt.fastq  /bioinfo/data/fastq/003.fastq.gz
```

003.cutadapt.fastq resultado

### Tarefa 02: Repetir o processo para as amostras 017 e 019

# BWA-mem: maximal exact matches
Alinha sequencias de tamanho 70bp-1Mbp com o algoritmo BWA-MEM. Em resumo o algoritmo trabalha com "alinhamento por sementes" com maximal exact matches (MEMs) e então estendendo sementes com o algoritmo Smith-Waterman (SW). Link. Tempo (~60s):

```bash
# cd para voltar para casa
cd

# rodar bwa para alinhar as sequencias contra o genoma de referencia
bwa mem -R '@RG\tID:003\tSM:003_NGSA\tLB:Agilent\tPL:Ion'  ~/references/hg19.fa /bioinfo/data/fastq/003.fastq.gz > resultados/003/003.sam
```

BWA-mem 003 resultado

### Tarefa 03: Repetir o processo para as amostras 017 e 019

# samtools: fixmate, sort e index

## samtools fixmate

Preencha coordenadas de posicionamento, posicione FLAGs relacionadas a partir a alinhamentos classificados por nome. Tempo (~5s):

```bash
samtools fixmate resultados/003/003.sam resultados/003/003.bam
```

SAMTOOLS fixmate 003 resultado

### Tarefa 04: Repetir o processo para as amostras 017 e 019

Tempo: ~10min

## samtools sort

O ***samtools sort*** vai ordenar de nome para ordem de coordenadas. Tempo (5s):

```bash
time samtools sort -O bam -o resultados/003/003_sort.bam -T /tmp/ resultados/003/003.bam 
```

SAMTOOLS sort 003 resultado

### Tarefa 05: Repetir o processo para as amostras 017 e 019

Tempo: ~1min

## samtools index

O ***samtools index*** cria um index (.BAI) do arquivo binário (.BAM):

```bash
samtools index resultados/003/003_sort.bam 
```

SAMTOOLS index 003 resultado

### Tarefa 06: Repetir o processo para as amostras 017 e 019

Tempo: ~1min

# freebayes: chamador de variantes
O FreeBayes é um detector variante genético Bayesiano projetado para encontrar pequenos polimorfismos, especificamente SNPs (polimorfismos de nucleotídeo único), indels (inserções e deleções), MNPs (polimorfismos de múltiplos nucleotídeos) e eventos complexos (eventos compostos de inserção e substituição) menores que os comprimento de um alinhamento de seqüenciamento de leitura curta. Link. Tempo (~6min):

```bash
freebayes -f ~/references/hg19.fa -F 0.01 -C 1 --pooled-continuous resultados/003/003_sort.bam > resultados/003/003.vcf
```

### Parâmetros


```bash
 -C --min-alternate-count N
                   Require at least this count of observations supporting
                   an alternate allele within a single individual in order
                   to evaluate the position.  default: 2

 --pooled-continuous
                   Output all alleles which pass input filters, regardles of
                   genotyping outcome or model.
```


FREEBAYES call variant 003 resultado

### Tarefa 07: Repetir o processo para as amostras 017 e 019

Tempo: ~10min

# annovar: anotador de variantes
ANNOVAR éma ferramenta eficiente para anotar funcionalmente variantes genéticas detectadas a partir de diversos genomas (incluindo o genoma humano hg18, hg19, hg38, bem como mouse, verme, mosca, levedura e muitos outros). [Link]. Aqui vamos converter .VCF para .avinput. Tempo (~5s):


```bash                   
perl /bioinfo/app/annovar/convert2annovar.pl -format vcf4 resultados/003/003.vcf > resultados/003/003.avinput
``` 
 
ANNOVAR convert2annovar 003 resultado

### Tarefa 08: Repetir o processo para as amostras 017 e 019

Anotar as variantes chamadas utilizando algumas bases de dados públicas: Tempo (~5s).

```
perl /bioinfo/app/annovar/table_annovar.pl resultados/003/003.avinput /bioinfo/app/annovar/humandb/ -buildver hg19 -out resultados/003/003 -remove -protocol refGene,exac03,clinvar_20180603 -operation g,f,f -nastring .
```

### Tarefa 09: Repetir o processo para as amostras 017 e 019

ANNOVAR table_annovar 003 resultado

# Documentação Extra


* Criar arquivos BEDs com UCSC table browser [Paste In Identifiers for GENCODE v29](https://genome.ucsc.edu/cgi-bin/hgTables)

* A Importância da Uniformidade de Cobertura sobre a Taxa de Alvos para NGS [Acessar White Paper](https://twistbioscience.com/sites/default/files/resources/2018-10/Twist_NGS_WhitePaper_UniformityonTarget_18Sep18.pdf)

* Referências bibliográficas aula – Bases Moleculares do Sequenciamento de Nova Geração

	* https://www.abmgood.com/marketing/knowledge_base/next_generation_sequencing_introduction.php
	* https://www.youtube.com/watch?time_continue=60&v=jFCD8Q6qSTM
	* https://www.illumina.com/science/technology/next-generation-sequencing.html#
	* https://www.illumina.com/content/dam/illumina-marketing/documents/products/illumina_sequencing_introduction.pdf

* Aplicações do NGS

	* https://www.illumina.com/techniques/sequencing.html

* Bases moleculares do sequenciamento
Preparação de biblioteca e design experimental:
	* https://www.youtube.com/watch?list=PLTt9kKfqE_0FigC-9K6doNq_D83BE9Kzw&v=-kTcFZxP6kM
	* https://www.youtube.com/watch?v=PGAfwSRYv1g&index=2&list=PLTt9kKfqE_0FigC-9K6doNq_D83BE9Kzw
	* https://www.abmgood.com/marketing/knowledge_base/next_generation_sequencing_experimental_design.php
	* https://www.illumina.com/techniques/sequencing/dna-sequencing/targeted-resequencing.html
	* https://www.thermofisher.com/br/en/home/life-science/sequencing/next-generation-sequencing/ion-torrent-next-generation-sequencing-workflow/ion-torrent-next-generation-sequencing-select-targets.html
	* https://www.illumina.com/science/technology/next-generation-sequencing/paired-end-vs-single-read-sequencing.html
	* https://www.illumina.com/science/technology/next-generation-sequencing/multiplex-sequencing.html
	* https://www.illumina.com/science/education/sequencing-coverage.html
	* https://www.illumina.com/science/technology/next-generation-sequencing/deep-sequencing.html
	* https://www.illumina.com/science/technology/next-generation-sequencing/mate-pair-sequencing.html

* Ferramentas para escolha de kits de preparo de biblioteca e de corrida Ion Torrent:
	* https://www.thermofisher.com/br/en/home/life-science/sequencing/next-generation-sequencing/ion-torrent-next-generation-sequencing-workflow/ion-torrent-next-generation-sequencing-construct-library.html
	* https://www.thermofisher.com/br/en/home/life-science/sequencing/next-generation-sequencing/ion-torrent-next-generation-sequencing-workflow/ion-torrent-next-generation-sequencing-run-sequence.html

* Ferramentas para escolha de kits de preparo de biblioteca e de corrida Illumina:
	* https://www.illumina.com/products/selection-tools.html
	* https://www.illumina.com/systems/sequencing-platforms/comparison-tool.html
	* https://www.illumina.com/techniques/sequencing/ngs-library-prep/library-prep-methods.html

* Sequenciamento:
	* Illumina: https://www.youtube.com/watch?v=fCd6B5HRaZ8
	* Ion Torrent: https://www.youtube.com/watch?v=WYBzbxIfuKs

* Análise dos dados
	* https://www.illumina.com/informatics/sequencing-data-analysis.html
	* https://www.abmgood.com/marketing/knowledge_base/next_generation_sequencing_data_analysis.php
	* https://www.youtube.com/watch?v=l4BAfRekohk
