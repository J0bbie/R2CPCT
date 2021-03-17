#!/bin/bash

if [[ "$1" == "-h" || "$1" == "--h" || "$1" == "--help" ]]; then
   echo -e "\e[32mUsage: \e[97m`basename $0` \e[46m[-v /mnt/data/...]\e[49m [-f /home/user/fasta.fasta] [-g gencode.gtf] [-e .vcf.gz] \e[46m[-o /home/user/folder/]\e[49m"
   exit 0
fi

while getopts "v:e:o:f:g:t:" OPTION; do
   case $OPTION in
      v)
      	VCF_PATH=${OPTARG}
      	;;
      f)
        FASTA_PATH=${OPTARG}
        ;;
      g)
        GTF_PATH=${OPTARG}
        ;;
      t)
        THREADS=${OPTARG}
        ;;
      e)
      	EXTENSION=${OPTARG}
	   ;;
      o)
      	OUTPUT_DIR=${OPTARG}
      	;;
      *)
        echo "Unknown option."
	exit 12
      ;;
   esac
done

if [[ -z "${VCF_PATH}" ]]; then
   echo "No path to VCF specified"
   exit 1
fi
if [[ -z "${OUTPUT_DIR}" ]]; then
   echo "No output path specified"
   exit 2
fi
if [[ -z "${FASTA_PATH}" ]]; then
   echo "FASTA needs to be set! Just make a symlink to /mnt/data/ccbc_environment/general/genomes/hsapiens/hg19_HMF/Homo_sapiens.GRCh37.GATK.illumina.fasta in a folder with write permissions."
   exit 3
fi
if [[ -z "${GTF_PATH}" ]]; then
   GTF_PATH="/mnt/data/ccbc_environment/general/annotation/hg19/GENCODE/noChrPrefix_VEP_gencode.v35lift37.annotation.gtf.gz"
fi
if [[ -z "${THREADS}" ]]; then
   THREADS=4
fi
if [[ -z "${EXTENSION}" ]]; then
   EXTENSION=".vcf.gz"
fi

files=`find "$VCF_PATH" -name "*$EXTENSION"`
for vcf in ${files}; do
   output_file=`basename ${vcf/.vcf.gz/.vep.vcf.gz}`
   command="/mnt/data/ccbc_environment/software/general/ensembl-vepv101/vep -i ${vcf} --pick"
   command+=" --ccds --hgvs --symbol --numbers --domains --canonical --protein --biotype --uniprot --tsl --appris"
   command+=" --gene_phenotype --pubmed --variant_class"
   command+=" --fasta ${FASTA_PATH} --vcf_info_field ANN --fork ${THREADS} --cache --offline --no_stats"
   command+=" --gtf ${GTF_PATH}"
   command+=" --dir_plugins /mnt/data/ccbc_environment/software/general/ensembl-vepv101/Plugins/ "
   command+=" --custom /mnt/data/ccbc_environment/software/general/ensembl-vepv101/gnoMAD_2.1.1/gnomad.exomes.r2.1.1.sites.AF.vcf.gz,gnomADe,vcf,exact,0,AF"
   command+=" --custom /mnt/data/ccbc_environment/software/general/ensembl-vepv101/gnoMAD_2.1.1/gnomad.genomes.r2.1.1.sites.AF.vcf.gz,gnomADg,vcf,exact,0,AF"
   command+=" --custom /mnt/data/ccbc_environment/software/general/ensembl-vepv101/cache/Plugins/ClinVar/clinvar_20201003.vcf.gz,ClinVar,vcf,exact,0,CLNDN,CLNHGVS,CLNSIG"
   command+=" --plugin CADD,/mnt/data/ccbc_environment/software/general/ensembl-vepv101/cache/Plugins/CADD/whole_genome_SNVs.tsv.gz,/mnt/data/ccbc_environment/software/general/ensembl-vepv101/cache/Plugins/CADD/InDels.tsv.gz"
   command+=" --plugin FATHMM_MKL,/mnt/data/ccbc_environment/software/general/ensembl-vepv101/cache/Plugins/FATHMM_MKL/fathmm-MKL_Current.tab.gz"
   command+=" --assembly GRCh37 --dir /mnt/data/ccbc_environment/software/general/ensembl-vepv101/cache "
   command+=" --quiet --synonyms /mnt/data/ccbc_environment/software/general/ensembl-vepv101/cache/homo_sapiens/101_GRCh37/chr_synonyms.txt"
   command+=" --format vcf --vcf -o ${OUTPUT_DIR}/${output_file} --compress_output bgzip"
   echo ${command}
done

