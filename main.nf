#!/usr/bin/env nextflow

/*
 *
 *  Pipeline            SnapperRocks
 *  Description         Pipeline for bacterial analysis
 *
 */

params_info = """
===============================================================================
                               SnapperRocks
================================================================================
Input Parameters:
    --run_id      $params.run_id
    --fastq       $params.fastq
    --results     $params.results
    --store_path  $params.store_path
    --db          $params.db
    --scripts     $params.scripts
    --metadata    $params.metadata

Analysis Parameters:
    --trim        $params.trim
    --fastqc      $params.fastqc
    --filterqc    $params.filterqc
    --taxoprofile $params.taxoprofile
    --assembly    $params.assembly
    --qc	  $params.qc
    --quast       $params.quast
    --post        $params.post
    --amr	  $params.amr
    --typing      $params.typing
    --redcap      $params.redcap
    --resprofile  $params.resprofile
    --intrinsic   $params.intrinsic
    --nesoni      $params.nesoni
    --pseudo      $params.pseudo
    --cluster     $params.cluster
    --staph       $params.staph

Execution Parameters:
    --executor    $params.executor
    --container   $params.container
    --njobs       $params.njobs
================================================================================

"""

reads = Channel
  .fromFilePairs("${params.fastq}", flat: true)
	.ifEmpty { exit 1, """ Input read files could not be found.
Have you included the read files in the current directory and do they have the correct naming?
With the parameters specified, reads named ${params.fastq}.
To fix this error either rename your reads to match this formatting or specify the desired format
when initializing SnapperRocks e.g. --fastq *_{1,2}_sequence.fastq.gz"""
}.into{ fastq; strains }

metadata_file = file(params.metadata)

file(params.results).mkdir()
params_file = file(params.results + '/run_params.txt')
params_file.text = params_info

log.info params_info

process Trimmomatic {
  label "general2"
  tag {"$id"}
  publishDir "${params.results}/Trimmed"

  input:
  set id, path(forward), path(reverse) from fastq

  output:
  set id, "${id}_R1.fq.gz", "${id}_R2.fq.gz" into trimmed
  set id, path("${id}_trim_out.log") optional true into qcTrims
  val true into trimDone

  shell:
  if (params.trim) {
    '''
source activate sr_trimmomatic

trimmomatic PE -threads !{task.cpus} "!{forward}" "!{reverse}" \
"!{id}_R1.fq.gz" "!{id}_R1_u.fq.gz" "!{id}_R2.fq.gz" "!{id}_R2_u.fq.gz" \
ILLUMINACLIP:!{params.db}/adapters/trimmomatic.fa:1:30:11 \
LEADING:10 TRAILING:10 MINLEN:50 HEADCROP:10 2> "!{id}_trim_out.log"
    '''
  }
  else {
    '''
mv "!{forward}" "!{id}_R1.fq.gz"
mv "!{reverse}" "!{id}_R2.fq.gz"
    '''
  }
}

process FastQC {
  label "general2"
  tag {"$id"}
  publishDir enabled: false

  input:
  set id, path(forward), path(reverse) from trimmed

  output:
  set id, "${id}_R1.fq.gz", "${id}_R2.fq.gz" optional true into (readsAssembly, readsABRes, readsStaph, readsNesoni, readsST, readsQC)
  env(PASSID) into (ids, krakenIds)
  path("${id}_R1.fq.gz") optional true into krakenForwards
  path("${id}_R2.fq.gz") optional true into krakenReverses
  set id, path("${id}_R1_fastqc.zip"), path("${id}_R2_fastqc.zip") optional true into qcFastQCResult
  set id, path("${id}_R1_fastqc.html"), path("${id}_R2_fastqc.html") optional true into qcFastQCReport
  path("*_fastqc.zip") optional true into fastQCResult
  set id, path("${id}_Failed_QC.csv") optional true into qcFail

  shell:
  if (params.fastqc) {
    '''
source activate sr_fastqc

mv "!{forward}" "R1.fq.gz"
mv "!{reverse}" "R2.fq.gz"

mv "R1.fq.gz" "!{id}_R1.fq.gz"
mv "R2.fq.gz" "!{id}_R2.fq.gz"

fastqc --threads !{task.cpus} --limits "!{params.db}/fastqc_limits.txt" "!{id}_R1.fq.gz" "!{id}_R2.fq.gz"
source deactivate

if [[ "!{params.filterqc}" == "true" ]]; then
  set +eu; source activate sr_fastqcr; set -eu
  cp "!{params.db}/multi-qc-report.Rmd" .
  !{params.scripts}/SR_process_fastqc_output.R "!{id}"
  set +eu; source deactivate; set -eu
  if [ -f "failed_reads.list" ]; then
    while read f; do
      if [ -f "$f" ]; then
        touch failed.qc
      fi
    done < failed_reads.list
  fi

  if [ -f failed.qc ]; then
    rm "!{id}_R1.fq.gz" "!{id}_R2.fq.gz"
    mv fails.csv "!{id}_Failed_QC.csv"
    PASSID="SR_FAILED"
  else
    PASSID="!{id}"
  fi
else
  PASSID="!{id}"
fi
    '''
  }
  else {
    '''
mv "!{forward}" "R1.fq.gz"
mv "!{reverse}" "R2.fq.gz"

mv "R1.fq.gz" "!{id}_R1.fq.gz"
mv "R2.fq.gz" "!{id}_R2.fq.gz"

PASSID="!{id}"
    '''
  }
}

process FastQCFilter {
  label "general2"
  tag {"$id"}
  publishDir "${params.results}/QC"

  input:
  set id, path(failed) from qcFail

  output:
  path("${id}_Failed_FastQC.csv")

  shell:
  '''
mv "!{failed}" "!{id}_Failed_FastQC.csv"
  '''
}

process FastQCSummary {
  label "general"
  publishDir "${params.results}/Summaries"

  input:
  path("*") from fastQCResult.collect()

  output:
  path("${params.run_id}_FastQC_Warnings.csv")
  path("${params.run_id}_FastQC_Failed.csv")
  path("${params.run_id}_FastQC_Report.html")
  val true into fastQCDone

  shell:
  '''
set +eu; source activate sr_fastqcr; set -eu
cp /SnapperRocks/miniconda3/envs/sr_fastqcr/lib/R/library/fastqcr/report_templates/multi-qc-report.Rmd .
!{params.scripts}/SR_process_fastqc_output.R "!{params.run_id}"
set +eu; source deactivate; set -eu

if [ -f warns.csv ]; then
  mv warns.csv "!{params.run_id}_FastQC_Warnings.csv"
fi

if [ -f fails.csv ]; then
  mv fails.csv "!{params.run_id}_FastQC_Failed.csv"
fi

if [ -f report.html ]; then
  mv report.html "!{params.run_id}_FastQC_Report.html"
fi

  '''
}

process BuildStrainList {
  label "general"
  publishDir "${params.results}/Summaries"

  input:
  file ids from ids.collectFile(name: 'strains.txt', newLine: true)

  output:
  path("${params.run_id}.strain.list.txt") into (strainListRedcap, strainListClust, strainListDynamo)
  val true into listDone

  shell:
  '''
sed -i '/SR_FAILED/d' strains.txt
wc -l strains.txt
mv strains.txt "!{params.run_id}.strain.list.txt"
  '''
}

process KrakenProfile {
  label "kraken"
  publishDir "${params.results}/Kraken"

  input:
  file ids from krakenIds.collectFile(name: 'ids.txt', newLine: true)
  path "*" from krakenForwards.collect()
  path "*" from krakenReverses.collect()

  output:
  path("*.report") into (krakenReports, krakenReportsQC)

  when:
  params.taxoprofile

  shell:
  '''
DB="!{params.db}/Kraken2"
USE_RAM="!{params.ramdisk}"
if [[ ${USE_RAM} == 'true' ]]; then
	RAM_PATH="!{params.ramdisk_path}"
	cp -rd "${DB}" "${RAM_PATH}/Kraken2"
	DB="${RAM_PATH}/Kraken2"
fi

source activate sr_kraken2
sed -i '/SR_FAILED/d' ids.txt

while read id;
do
  forward="${id}_R1.fq.gz"
  reverse="${id}_R2.fq.gz"
  if [[ ${USE_RAM} == 'true' ]]; then
    kraken2 --db "${DB}" --report "${id}.report" --gzip-compressed --threads !{task.cpus} --paired "${forward}" "${reverse}" --memory-mapping
  else
    kraken2 --db "${DB}" --report "${id}.report" --gzip-compressed --threads !{task.cpus} --paired "${forward}" "${reverse}"
  fi
done < ids.txt

if [[ ${USE_RAM} == 'true' ]]; then
	rm -rdf "${RAM_PATH}/Kraken2"
fi
  '''
}

process KrakenReport {
  label "general"
  publishDir "${params.results}/Summaries"

  input:
  path "*" from krakenReports.collect()

  output:
  path("${params.run_id}.kraken.summary") into (krakenSummarySRST2, krakenSummaryRedcap, krakenSummaryNesoni, krakenSummaryDynamo, krakenSummaryQC)
  path("*.Species.txt") into species
  path("AllSpecies.txt") into (allSpeciesChrRes, allSpeciesPseudoNesoni, allSpeciesAmrfinder)
  path("Staphylococcus aureus.Species.txt") optional true into (staphListFilter, staphListSummary, staphListReport)
  val true into krakenCluster
  val true into krakenDone

  shell:
  '''
mkdir -p "!{params.store_path}/Kraken"
cp -L *.report "!{params.store_path}/Kraken/"

for F in *.report; do
  "!{params.scripts}/summarize_kraken.sh" "$F" >> tmp.kraken.summary;
done

awk 'BEGIN{print "ID\tReads [%]\tReads [#]\tunclassified\t%\tReads [#]\t#1 Genus\t%\tReads [#]\t#2 Genus\t%\t Reads [#]\t#3 Genus\t%\tT Reads [#]\t#1 Species\t%\tReads [#]\t#2 Species\t%\tReads [#]\t#3 Species"}1' tmp.kraken.summary | cut -f 1-3,5- > "!{params.run_id}.kraken.summary"

rm tmp.kraken.summary

IFS=$'\n'
for SPEC in $(cat "!{params.run_id}.kraken.summary" | awk -F $'\t' 'NR>1 {print $15}' | sort | uniq); do
  echo "Species: ${SPEC}"
  SPECIDS=(`cat "!{params.run_id}.kraken.summary" | awk -F $'\t' -v spec="${SPEC}" '$15==spec {print $1}'`)
  echo "IDS: ${SPECIDS[@]}"
  for ID in ${SPECIDS[@]}; do
    echo ${ID} >> "${SPEC}.Species.txt"
    echo "${SPEC}:${ID}" >> AllSpecies.txt
  done
done
unset IFS
  '''
}

process SpadesAssembly {
  label "general2"
  tag {"$id"}
  publishDir "${params.results}/Assemblies"

  input:
  set id, path(forward), path(reverse) from readsAssembly

  output:
  path("${id}.ass.metrics") into assMetrics
  path("${id}.filtered_ass.metrics") into filteredMetrics
  set id, path("${id}.contigs_filtered.fasta") into (assembliesAbricate, assembliesChrRes, assembliesAmr, assembliesQC)
  path("${id}.contigs_filtered.fasta") into (assembliesMLST, assembliesStaph, assembliesClust)
  val true into spadesDone

  when:
  params.assembly

  shell:
  '''
source activate sr_spades

spades.py -o . -1 "!{forward}" -2 "!{reverse}" -t !{task.cpus} --isolate

source activate sr_perl

COV=10
LEN=100

perl "!{params.scripts}/filterSpadesAss.pl" contigs.fasta ${COV} ${LEN}
perl "!{params.scripts}/fastagrep.pl" -f contigs.fasta.txt contigs.fasta > contigs.fasta.filtered
mv contigs.fasta.filtered "!{id}.contigs_filtered.fasta"

perl "!{params.scripts}/fastasize.pl" contigs.fasta | sed "s/#/!{id}/" > "!{id}.ass.metrics"
perl "!{params.scripts}/fastasize.pl" "!{id}.contigs_filtered.fasta" | sed "s/#/!{id}/" > "!{id}.filtered_ass.metrics"
  '''
}

process SpadesMetrics {
  label "general"
  publishDir "${params.results}/Summaries"
  
  input:
  path("*") from assMetrics.collect()
  path("*") from filteredMetrics.collect()

  output:
  path("${params.run_id}.AssMetrics.stats")
  path("${params.run_id}.FilteredAssMetrics.stats")

  shell:
  '''
cat *.ass.metrics > "!{params.run_id}.AssMetrics.stats"
cat *.filtered_ass.metrics > "!{params.run_id}.FilteredAssMetrics.stats"

sed -i "1iSample,Genomesize(bp),#Contigs,Min,Max,Avg" "!{params.run_id}.AssMetrics.stats"
sed -i "1iSample,Genomesize(bp),#Contigs,Min,Max,Avg" "!{params.run_id}.FilteredAssMetrics.stats"
  '''
}

process QC {
  label "general"
  tag {"$id"}
  publishDir enabled: false

  input:
  each path("${params.run_id}.kraken.summary") from krakenSummaryQC
  each path("*") from krakenReportsQC
  set id, path(qcTrim), path(qcForward), path(qcReverse), path(qcAssembly), path(qcFastQCR1Result), path(qcFastQCR2Result), path(qcFastQCR1Report), path(qcFastQCR2Report) from qcTrims.join(readsQC).join(assembliesQC).join(qcFastQCResult).join(qcFastQCReport)

  output:
  path("${id}/") into qcFiles

  when:
  params.qc

  shell:
  '''
if [[ "!{params.quast}" == "true" ]]; then
    if [ -f "!{qcAssembly}" ]; then
        set +eu
        source activate sr_python
        if REFSPEC="`python !{params.scripts}/SR_refseek.py -q -s !{id} -r !{params.db}/Ref_db.tab -k !{params.run_id}.kraken.summary`"; then
                echo "Found REF: ${REFSPEC}"
        else
                echo "No REF found for species ${REFSPEC}"
                echo "QC output will be limited"
        fi
        set -eu

        if [ ! -z ${REFSPEC} ]; then
            source activate sr_bioconvert
            bioconvert genbank2fasta "!{params.db}/ref/${REFSPEC}.gbk" "${REFSPEC}.fa"
            python "!{params.scripts}/SR_convert_genbank_to_gff3.py" -i "!{params.db}/ref/${REFSPEC}.gbk" -o "${REFSPEC}.gff" --no_fasta
            source deactivate

            set +eu; source activate sr_quast; set -eu
            quast.py "!{qcAssembly}" \
                -r "${REFSPEC}.fa" \
                -g "${REFSPEC}.gff" \
                -1 "!{qcForward}" -2 "!{qcReverse}" \
                -o "!{id}" \
                --conserved-genes-finding \
                --circos \
                --threads !{task.cpus}
            set +eu; source deactivate; set -eu

            rm "${REFSPEC}.fa"
            rm "${REFSPEC}.gff"
        else
            set +eu; source activate sr_quast; set -eu
            quast.py "!{qcAssembly}" \
                -1 "!{qcForward}" -2 "!{qcReverse}" \
                -o "!{id}" \
                --conserved-genes-finding \
                --circos \
                --threads !{task.cpus}
            set +eu; source deactivate; set -eu
        fi

        mv "!{id}/report.html" "!{id}/!{id}_QUAST_Report.html"

        if [ -f "!{id}/circos/circos.png" ]; then
            cp "!{id}/circos/circos.png" "!{id}/!{id}_Circos.png"
        fi

        if [ -f "!{id}/circos/legend.txt" ]; then    
            cp "!{id}/circos/legend.txt" "!{id}/!{id}_Circos_Legend.txt"
        fi
    else
        mkdir -p "!{params.results}/QC"
        echo "!{id}" >> "!{params.results}/QC/no_assemblies.list"
    fi
fi

if [ -f "!{qcTrim}" ]; then
    mkdir -p "!{id}"
    mv "!{qcTrim}" "!{id}/"
fi

if [ -f "!{id}.report" ]; then
    mkdir -p "!{id}/"
    mv "!{id}.report" "!{id}/"
fi

if [ -f "!{qcFastQCR1Result}" ]; then
    mkdir -p "!{id}/"
    mv "!{qcFastQCR1Result}" "!{id}/!{id}_R1_fastqc.zip"
fi

if [ -f "!{qcFastQCR2Result}" ]; then
    mkdir -p "!{id}/"
    mv "!{qcFastQCR2Result}" "!{id}/!{id}_R2_fastqc.zip"
fi

if [ -f "!{qcFastQCR1Report}" ]; then
    mkdir -p "!{id}/"
    mv "!{qcFastQCR1Report}" "!{id}/!{id}_R1_fastqc.html"
fi

if [ -f "!{qcFastQCR2Report}" ]; then
    mkdir -p "!{id}/"
    mv "!{qcFastQCR2Report}" "!{id}/!{id}_R2_fastqc.html"
fi

if [ -d "!{id}" ]; then
    source activate sr_multiqc
    multiqc --outdir "!{id}/" --filename "!{id}_MultiQC_Report.html" --title "SnapperRocks Sample !{id}" --no-megaqc-upload "!{id}/"
    source deactivate
fi
  '''
}

process QCSummary {
  label "general"
  publishDir "${params.results}"

  input:
  path "*" from qcFiles.collect()

  output:
  path("Summaries/${params.run_id}_Summary_QC_Report.html")
  path("QC/")
  val true into qcDone

  shell:
'''
mkdir -p QC
mv */*.html QC/ 2>/dev/null || true
mv */*.png QC/ 2>/dev/null || true
mv */*Legend.txt QC/ 2>/dev/null || true

source activate sr_multiqc
multiqc --outdir "./" --filename "!{params.run_id}_Summary_QC_Report.html" --title "SnapperRocks Run !{params.run_id}" --no-megaqc-upload "."
source deactivate

mkdir -p Summaries
mv "!{params.run_id}_Summary_QC_Report.html" Summaries/
'''
}

process Abricate {
  label "general"
  tag {"$id"}
  publishDir "${params.results}/Abricate"

  input:
  set id, path(assembly) from assembliesAbricate

  output:
  path("${id}.ABRes.tab") into abricateABRes
  path("${id}.Vir.tab") into abricateVir
  path("${id}.Pla.tab") into abricatePla

  when:
  params.post

  shell:
  '''
source activate sr_abricate

abricate "!{assembly}" > "!{id}.ABRes.tab"
abricate --db vfdb "!{assembly}" > "!{id}.Vir.tab"
abricate --db plasmidfinder "!{assembly}" > "!{id}.Pla.tab"
  '''
}

process AbricateSummary {
  label "general"
  publishDir "${params.results}/Summaries"

  input:
  path("*") from abricateABRes.collect()
  path("*") from abricateVir.collect()
  path("*") from abricatePla.collect()

  output:
  path("*Summary.tab") into AbricateSummaries
  val true into abricateDone

  script:
  /
source activate sr_abricate

abricate --summary `ls *.ABRes.tab` > "${params.run_id}.abricate_ResSummary.tab"
sed -i 's@\.contigs_filtered\.fasta@@' "${params.run_id}.abricate_ResSummary.tab"

abricate --summary `ls *.Vir.tab` > "${params.run_id}.abricate_VirSummary.tab"
sed -i 's@\.contigs_filtered\.fasta@@' "${params.run_id}.abricate_VirSummary.tab"

abricate --summary `ls *.Pla.tab` > "${params.run_id}.abricate_PlaSummary.tab"
sed -i 's@\.contigs_filtered\.fasta@@' "${params.run_id}.abricate_PlaSummary.tab"
  /
}

process MLST {
  label "general"
  publishDir "${params.results}/MLST"
  
  input:
  path "*" from assembliesMLST.collect()

  output:
  path("${params.run_id}.MLST.tab") into mlstSummary
  path("${params.run_id}.NOST.txt") into nost

  when:
  params.typing

  script:
  /
source activate sr_mlst

mlst --exclude cronobacter,ecoli_2 --nopath *.fasta > ${params.run_id}.MLST.tab
sed -i 's@[A-Za-z]\t-\t-$@a\t-\t-\t-\t-\t-\t-\t-\t-\t-@' ${params.run_id}.MLST.tab
sed -i 's@\.contigs_filtered\.fasta@@' ${params.run_id}.MLST.tab

set +eu
source activate sr_python
set -eu

python ${params.scripts}\/SR_mlstQC.py -m ${params.run_id}.MLST.tab > ${params.run_id}.NOST.txt

  /
}

process FilterSRST2 {
  label "general"
  tag {"$id"}
  publishDir "${params.results}/NOST"

  input:
  set id, path(forward), path(reverse) from readsST
  each path("${params.run_id}.NOST.txt") from nost

  output:
  set id, path("${id}_R1.fq.gz"), path("${id}_R2.fq.gz") optional true into readsSRST2

  shell:
  '''
echo if grep -Fxq "!{id}" "!{params.run_id}.NOST.txt"
if grep -Fxq "!{id}" "!{params.run_id}.NOST.txt"
then
  echo "!{id} identified for further SRST2 processing"
else
  rm "!{forward}"
  rm "!{reverse}"
fi
  '''
}

process SRST2 {
  label "general"
  tag {"$id"}
  publishDir "${params.results}/SRST2"

  input:
  each path("${params.run_id}.kraken.summary") from krakenSummarySRST2
  set id, path(forward), path(reverse) from readsSRST2

  output:
  path("*_results.txt") optional true into srst2Results
  path("*_no_mlst.txt") optional true

  shell:
  '''
source activate sr_python

STDB_SUM="!{params.db}/MLST/STDB.list"
set +eu; source activate sr_python
if REFSPEC="`python '!{params.scripts}/SR_refseek.py' -q -s !{id} -r '!{params.db}/Ref_db.mlst.tab' -k '!{params.run_id}.kraken.summary'`"; then
        echo "Found REF: ${REFSPEC}"
else
        echo "No REF found for species ${REFSPEC}"
        echo "${REFSPEC}:!{id}" > "!{id}_no_ref.txt"
	exit
fi
set -eu
echo "REFSPEC: ${REFSPEC}"

IFS=$'\n'
for SPEC in $(cat "${STDB_SUM}"); do
  echo "Species:${SPEC}"
  if [ "${SPEC}" == "${REFSPEC}" ]; then
    MLST=${REFSPEC}
    echo "Found REF: ${MLST}"
  fi
done
if [ -z ${MLST} ]; then
  echo "No MLST profile defined for ${REFSPEC}"
  echo "${SPEC}:!{id}" > "!{id}_no_mlst.txt"
  exit
fi
echo "MLST selected:${MLST}"
unset IFS

source activate sr_srst2

MLST_DB="!{params.db}/MLST/${MLST}/${MLST}.fa"
MLST_DEF="!{params.db}/MLST/${MLST}/${MLST}.txt"
echo "${MLST_DB} selected using ${MLST_DEF} definitions"

srst2 --threads !{task.cpus} --input_pe "!{forward}" "!{reverse}" --log --output "./!{id}" --mlst_delimiter '_' --mlst_db "${MLST_DB}" --mlst_definitions "${MLST_DEF}" --forward _R1 --reverse _R2

  '''
}

process SRST2Summary {
  label "general"
  publishDir "${params.results}/Summaries"

  input:
  path "*" from srst2Results.collect()

  output:
  path("${params.run_id}_MLST_srst2.summary") into srst2Summary

  shell:
  '''
source activate sr_srst2

numFiles=$(ls *_results.txt | wc -l)
if (($numFiles > 1)); then 
    srst2 --prev_output *_results.txt --output "!{params.run_id}.MLST"
else 
    mv *_results.txt "!{params.run_id}.MLST__compiledResults.txt"
fi

sed -i 's/_trimmed//' "!{params.run_id}.MLST__compiledResults.txt"
mv "!{params.run_id}.MLST__compiledResults.txt" "!{params.run_id}_MLST_srst2.summary"
  '''
}

process MLSTSummary {
  label "general"
  publishDir "${params.results}/Summaries"

  input:
  file "*" from srst2Summary.ifEmpty { 'EMPTY' }
  path "*" from mlstSummary

  output:
  path("${params.run_id}.MLST.summary.tab") into (stSummaryRedcap, stSummaryStaph, stSummaryDynamo)
  val true into mlstCluster
  val true into stDone

  shell:
  '''
if [ -f "!{params.run_id}_MLST_srst2.summary" ]; then
  source activate sr_python
  python "!{params.scripts}/SR_merge_SRST2_MLST_summaries.py" -s "!{params.run_id}_MLST_srst2.summary" -m "!{params.run_id}.MLST.tab" -o "!{params.run_id}.MLST.summary.tab"
else
  mv "!{params.run_id}.MLST.tab" "!{params.run_id}.MLST.summary.tab"
fi

mkdir -p "!{params.store_path}/MLST"
cp -L "!{params.run_id}.MLST.summary.tab" "!{params.store_path}/MLST/"
  '''
}

process Amrfinderplus {
  label "general"
  tag {"$id"}
  publishDir "${params.results}/Amrfinderplus"

  input:
  set id, path(assembly) from assembliesAmr
  each path("AllSpecies.txt") from allSpeciesAmrfinder

  output:
  path("${id}.amrfinder.tsv") into amrfinder

  when:
  params.amr

  shell:
  '''
source activate sr_python
AMRFlags=`python "!{params.scripts}/amr_filter.py" --id "!{id}"`
source deactivate

source activate sr_amr
export TMPDIR="`pwd`/tmp_`date -u +%s%N`"
mkdir $TMPDIR

amrfinder -n "!{assembly}" --threads !{task.cpus} --plus $AMRFlags --name="!{id}" > "!{id}.amrfinder.tsv"

rm -rdf $TMPDIR
  '''
}

process AmrfinderSummary{
  label "general"
  publishDir "${params.results}/Summaries"

  input:
  path("*") from amrfinder.collect()

  output:
  path("*.csv") into amrfinderSummaries
  val true into amrDone

  shell:
  '''
source activate sr_python

python "!{params.scripts}/amrfinder_summary.py" --output "!{params.run_id}_Amrfinder" *.amrfinder.tsv
  '''
}

process AMRSummary {
  label "general"
  publishDir "${params.results}/Summaries"

  input:
  path "*" from amrfinderSummaries
  path "*" from AbricateSummaries

  output:
  path("*.xlsx")
  val true into amrSDone

  shell:
  '''
source activate sr_python

python "!{params.scripts}/amr_summary.py"
mv 'AMR_Combined.xlsx' "!{params.run_id}_AMR_Combined.xlsx"
  '''
}

process QueryRedcap {
  label "general"
  publishDir "${params.results}/Summaries"

  input:
  path("${params.run_id}.MLST.summary.tab") from stSummaryRedcap
  path("${params.run_id}.kraken.summary") from krakenSummaryRedcap
  path("${params.run_id}.strain.list.txt") from strainListRedcap

  output:
  path("StrainInfo_Redcap.tab")
  path("StrainInfo_Redcap.err") optional true
  path("${params.run_id}.strain.list.txt")
  val true into redcapDone

  when:
  params.redcap

  shell:
  '''
source activate sr_python

python "!{params.scripts}/SR_generate_sample_mlst_report_redcap.py" -s "!{params.run_id}.strain.list.txt" -m "!{params.run_id}.MLST.summary.tab" -k "!{params.run_id}.kraken.summary" -a "!{params.redcap_key}" -o "StrainInfo_Redcap.tab" -e "StrainInfo_Redcap.err"

  '''
}

process ABRes {
  label "general2"
  tag {"$id"}
  publishDir "${params.results}/ABRes"

  input:
  set id, path(forward), path(reverse) from readsABRes

  output:
  path("*__genes__*.txt") into abresResults

  when:
  params.resprofile

  shell:
  '''
source activate sr_srst2

RES_DB="!{params.db}/argannot/argannot0317.v3.clustered.fa"

srst2 --threads !{task.cpus} --input_pe "!{forward}" "!{reverse}" --output "!{id}" --log --gene_db "${RES_DB}" --forward _R1 --reverse _R2
  '''
}

process ABResSummary {
  label "general"
  publishDir "${params.results}/Summaries"

  input:
  path "*" from abresResults.collect()

  output:
  path("${params.run_id}.ResProfile")
  val true into abResDone

  shell:
  /
source activate sr_srst2

srst2 --prev_output *__genes__*.txt  --output "!{params.run_id}.ResProfile"
cat "!{params.run_id}.ResProfile__compiledResults.txt" | sed 's@_trimmed@@' |sed 's@\:[[:graph:]]*@@g' | sed 's@([[:alpha:]]*)@@g' > "!{params.run_id}.ResProfile"
  /
}

process ChrRes {
  label "general"
  tag {"$id"}
  publishDir "${params.results}/ChrRes"

  input:
  set id, path(assembly) from assembliesChrRes
  each path("AllSpecies.txt") from allSpeciesChrRes

  output:
  path("${id}_blastn_*") optional true
  path("${id}_no_ref.txt") optional true
  val true into intrinsicDone

  when:
  params.intrinsic

  shell:
  '''
source activate sr_python

S1="Klebsiella pneumoniae";S2="Escherichia coli";S3="Enterococcus faecium";S4="Enterococcus faecalis";S5="Staphylococcus aureus"

BLASTN_PATH="`which blastn`"
SPEC_SUM=AllSpecies.txt

IFS=$'\n'
for SPEC in $(cat "${SPEC_SUM}" |grep "!{id}" | awk -F $':' '{print $1}' | sort | uniq);do
  echo "Species: ${SPEC}"
  if [ "${SPEC}" == "${S1}" ]; then
    python "!{params.scripts}/pointfinder/PointFinder.py" -i "!{assembly}" -o . -s klebsiella -m blastn -u -r all -p "!{params.db}/pointfinder_db/" -m_p "${BLASTN_PATH}"

  elif [ "${SPEC}" == "${S2}" ]; then
    python "!{params.scripts}/pointfinder/PointFinder.py" -i "!{assembly}" -o . -s escherichia_coli -m blastn -u -r all -p "!{params.db}/pointfinder_db/" -m_p "${BLASTN_PATH}"

  elif [ "${SPEC}" == "${S3}" ] ;then
    python "!{params.scripts}/pointfinder/PointFinder.py" -i "!{assembly}" -o . -s enterococcus_faecium -m blastn -u -r all -p "!{params.db}/pointfinder_db/" -m_p "${BLASTN_PATH}"

  elif [ "${SPEC}" == "${S4}" ];then
    python "!{params.scripts}/pointfinder/PointFinder.py" -i "!{assembly}" -o . -s enterococcus_faecalis -m blastn -u -r all -p "!{params.db}/pointfinder_db/" -m_p "${BLASTN_PATH}"

  elif [ "${SPEC}" == "${S5}" ];then
    python "!{params.scripts}/pointfinder/PointFinder.py" -i "!{assembly}" -o . -s staphylococcus_aureus -m blastn -u -r all -p "!{params.db}/pointfinder_db/" -m_p "${BLASTN_PATH}"

  else
    echo "No reference defined for ${SPEC}"
    echo "${SPEC}:!{id}" > "!{id}_no_ref.txt"
  fi
done
  '''
}

process Nesoni {
  label "nesoni"
  tag "${id}"
  publishDir "${params.results}/Nesoni/${id}"

  input:
  set id, path(forward), path(reverse) from readsNesoni
  each path("${params.run_id}.kraken.summary") from krakenSummaryNesoni

  output:
  set id, path("${id}_report.txt") optional true into nesoniReport
  path("${id}_nesoni/*") optional true
  path("${id}_no_ref.txt") optional true
  path("${id}_done") into nesoniCluster
  val true into nesoniDone

  when:
  params.nesoni

  shell:
  '''
set +eu; source activate sr_python
if REFSPEC="`python !{params.scripts}/SR_refseek.py -q -s !{id} -r !{params.db}/Ref_db.tab -k !{params.run_id}.kraken.summary`"; then
	echo "Found REF: ${REFSPEC}"
else
	echo "No REF found for species ${REFSPEC}"
	echo "${REFSPEC}:!{id}" > "!{id}_no_ref.txt"
	exit
fi
set -eu

echo ${REFSPEC}
REF_GBK="!{params.db}/ref/${REFSPEC}.gbk"
REF_FA="!{params.db}/ref/${REFSPEC}.fa"

set +eu; source activate sr_nesoni; set -eu

nesoni --make-cores !{task.cpus} analyse-sample "!{id}_nesoni" "${REF_GBK}" align: bowtie clip: none pairs: "!{forward}" "!{reverse}"
cp -L "!{id}_nesoni/report.txt" "!{id}_report.txt"

mkdir -p "!{params.store_path}/Nesoni_reports/!{id}"
cp -L "!{id}_nesoni"/*-evidence.txt "!{params.store_path}/Nesoni_reports/!{id}/"
cp -L "!{id}_nesoni/reference/reference.fa" "!{params.store_path}/Nesoni_reports/!{id}/"

touch "!{id}_done"
  '''
}

process PseudoNesoni {
  label "nesoni"
  tag "${id}"
  publishDir "${params.results}/PseudoNesoni/${id}"

  input:
  each path("AllSpecies.txt") from allSpeciesPseudoNesoni
  set id, path("${id}_report.txt") from nesoniReport

  output:
  path("*.fas") optional true
  path("${id}_no_ref.txt") optional true
  val true into pseudoDone

  when:
  params.pseudo

  shell:
  '''
source activate sr_perl

S1="Klebsiella pneumoniae";S2="Escherichia coli";S3="Citrobacter freundii";S4="Enterobacter cloacae";S5="Pseudomonas aeruginosa"
S6="Acinetobacter baumannii";S7="Salmonella enterica";S8="Citrobacter rodentium";S9="Citrobacter koseri";S10="Enterobacter aerogenes"
S11="Enterobacter asburiae";S12="Klebsiella oxytoca";S13="Staphylococcus aureus"
SPEC_SUM=AllSpecies.txt
IFS=$'\n'
for SPEC in $(cat "${SPEC_SUM}" |grep "!{id}" | awk -F $':' '{print $1}' | sort | uniq);do
  echo "Species: ${SPEC}"
  if [ "${SPEC}" == "${S1}" ]; then
    REFSPEC=kpneumoniae2
  elif [ "${SPEC}" == "${S2}" ]; then
    REFSPEC=ecoli2
  elif [ "${SPEC}" == "${S3}" ]; then
    REFSPEC=cfreundii
  elif [ "${SPEC}" == "${S4}" ]; then
    REFSPEC=ecloacae
  elif [ "${SPEC}" == "${S5}" ]; then
    REFSPEC=paeruginosa
  elif [ "${SPEC}" == "${S6}" ]; then
    REFSPEC=abaumannii
  elif [ "${SPEC}" == "${S7}" ]; then
    REFSPEC=senterica
  elif [ "${SPEC}" == "${S8}" ]; then
    REFSPEC=crodentium
  elif [ "${SPEC}" == "${S9}" ]; then
    REFSPEC=ckoseri
  elif [ "${SPEC}" == "${S10}" ]; then
    REFSPEC=eaerogenes
  elif [ "${SPEC}" == "${S11}" ]; then
    REFSPEC=easburiae
  elif [ "${SPEC}" == "${S12}" ]; then
    REFSPEC=koxytoca
  elif [ "${SPEC}" == "${S13}" ]; then
    REFSPEC=Saureus
  else
    echo "No reference defined for ${SPEC}"
    echo "${SPEC}:!{id}" > "!{id}_no_ref.txt"
    exit
  fi
done

unset IFS

REF_FA="!{params.db}/ref/${REFSPEC}.fa"

perl "!{params.scripts}/SNP_mutator.pl" "!{id}_report.txt" "${REF_FA}" "${REFSPEC}"
  '''
}

process ClusterFilter {
  label "nesoni"
  publishDir "${params.results}/Clustering"

  input:
  path "*" from nesoniCluster.collect()
  val mlstFlag from mlstCluster
  val krakenFlag from krakenCluster
  path("${params.run_id}.strain.list.txt") from strainListClust
  file metadata from metadata_file

  output:
  path "*.list" into clusterProcess
  path "metadata.csv" into clusterMetadata
  path "error.txt" optional true
  path "*.snpDists" optional true

  when:
  params.cluster

  shell:
  def metafilter = metadata.name != 'false' ? "--data \"$metadata\"" : ''
  """
#kraken summary
for F in "!{params.store_path}"/Kraken/*.report; do
        "!{params.scripts}/summarize_kraken.sh" \$F >> tmp.kraken.summary
        done
        awk 'BEGIN{print "ID\tReads [%]\tReads [#]\tunclassified\t%\tReads [#]\t#1 Genus\t%\tReads [#]\t#2 Genus\t%\t Reads [#]\t#3 Genus\t%\tT Reads [#]\t#1 Species\t%\tReads [#]\t#2 Species\t%\tReads [#]\t#3 Species"}1' tmp.kraken.summary | cut -f 1-3,5- > "Compiled.Kraken.summary.tab"
        rm tmp.kraken.summary

#MLST summary
cat "!{params.store_path}"/MLST/*.MLST.summary.tab > "Compiled.MLST.summary.tab"

#Generate Strain list
find "!{params.store_path}/Nesoni_reports/" -maxdepth 1 -mindepth 1 -type d | xargs -n 1 basename > "Compiled.strain.list.txt"

#Generate Straininfo report
source activate sr_python

if [[ "!{params.redcap}" == "true" ]]; then
	python "!{params.scripts}/SR_generate_sample_mlst_report_redcap.py" -s "Compiled.strain.list.txt" -m "Compiled.MLST.summary.tab" -k "Compiled.Kraken.summary.tab" $metafilter -a "!{params.redcap_key}" -o StrainInfo.tab -e StrainInfo.err
else
	python "!{params.scripts}/SR_generate_sample_mlst_report_redcap.py" -s "Compiled.strain.list.txt" -m "Compiled.MLST.summary.tab" -k "Compiled.Kraken.summary.tab" $metafilter -o StrainInfo.tab -e StrainInfo.err
fi

#Fix K. michiganensis samples statically
Km_SAMPLES=`cat "!{params.db}/K.michiganensis.samples"`
while read l; do
	Km_RESULT=`echo \$l | grep -q -F "\${Km_SAMPLES}" && echo \${l/Klebsiella oxytoca/Klebsiella michiganensis} || echo \$l`
	echo \${Km_RESULT}
done < "StrainInfo.tab" > "Compiled.StrainInfo.tab"

#Run Clustering
set +eu; source activate sr_nesoni; set -eu
python "!{params.scripts}/SR_process_ST_SNP_tables.py" -e -i "Compiled.StrainInfo.tab" -b "!{params.store_path}/Nesoni_reports/" -s "!{params.run_id}.strain.list.txt"

mv "Compiled.StrainInfo.tab" "metadata.csv"
  """
}

process Cluster {
  label "nesoni"
  publishDir "${params.results}/Clustering"

  input:
  path "*" from clusterProcess.flatten()

  output:
  path "*.snpDists" into clusterSNPDists

  shell:
  '''
set +eu; source activate sr_nesoni; set -eu
python "!{params.scripts}/SR_generate_ST_SNP_tables.py" -i *.list
  '''
}

process NetworkGenerate {
  label "general"
  publishDir "${params.results}/Clustering"

  input:
  path "*" from clusterSNPDists.collect()
  path "*" from clusterMetadata

  output:
  path "*.json"
  val true into clusterDone

  shell:
  '''
source activate sr_python
python "!{params.scripts}/SR_process_snp_dists.py" --input "./" --metadata "./metadata.csv" --output "./"

mkdir -p "!{params.store_path}/CATHAI"
cp -L *.json "!{params.store_path}/CATHAI"
cp -L *.snpDists "!{params.store_path}/CATHAI"
cp -L metadata.csv "!{params.store_path}/CATHAI"
  '''
}

process StaphFilter {
  label "general"
  tag "${id}"
  publishDir enabled: false

  input:
  set id, path(forward), path(reverse), path("Staphylococcus_aureus.Species.txt") from readsStaph.combine(staphListFilter)

  output:
  set id, "${id}_R1.fq.gz", "${id}_R2.fq.gz", "${id}.sa" into readsStaphFiltered optional true

  when:
  params.staph

  shell:
  '''
STRAINLIST="Staphylococcus_aureus.Species.txt"
for SAMPLE in $(cat "${STRAINLIST}"); do
    if [[ "${SAMPLE}" == "!{id}" ]]; then
        touch "!{id}.sa"
    fi
done
  '''
}

process StaphTyperSetup {
  label "general"
  tag "${id}"
  publishDir "${params.results}/StaphTyper"

  input:
  set id, path(forward), path(reverse), path(nix) from readsStaphFiltered

  output:
  path "*_calls.vcf.gz" into staphResultsVCF
  path "*_consensus.fa" into staphResultsConsensus
  path "*.cov" into staphResultsCov

  when:
  params.staph

  shell:
  '''
DBASE="!{params.db}/Staph"
source activate sr_bowtie2

bowtie2 -x "${DBASE}/goi" -1 "!{forward}" -2 "!{reverse}" -S "!{id}.sam"

source activate sr_samtools

samtools view -bS "!{id}.sam" > "!{id}.bam"
samtools sort -@ !{task.cpus} "!{id}.bam" -o "!{id}.sorted.bam"
samtools depth "!{id}.sorted.bam" > "!{id}.cov"

source activate sr_bcftools

bcftools mpileup -Ou -f "${DBASE}/goi.fa" "!{id}.sorted.bam" | bcftools call -mv -Oz -o calls.vcf.gz

tabix "calls.vcf.gz"
cat "${DBASE}/goi.fa" | bcftools consensus calls.vcf.gz > consensus.fa
mv calls.vcf.gz "!{id}_calls.vcf.gz"
mv consensus.fa "!{id}_consensus.fa"
  '''
}

process StaphTyperSummary {
  label "general"
  publishDir "${params.results}/Summaries"

  input:
  path "*" from staphResultsVCF.collect()
  path "*" from staphResultsConsensus.collect()
  path "*" from staphResultsCov.collect()
  path "*" from assembliesStaph.collect()
  path("Staphylococcus_aureus.Species.txt") from staphListSummary

  output:
  path "*.results" into staphSummaries

  shell:
  '''
# summarize gene coverage depth
COV_RESULTS="cov.results"
STRAINLIST="Staphylococcus_aureus.Species.txt"
GOIs=("aroE" "lukF" "lukS" "cna" "sdrE" "pUB110" "pT181")
echo "Summarizing gene coverage depths"
if [ -f "${COV_RESULTS}" ]; then
   echo "Removing old results file.."
   rm "${COV_RESULTS}"
fi

set +eo pipefail

echo "Sample ${GOIs[@]}" > "${COV_RESULTS}"
for SAMPLE in $(cat "${STRAINLIST}"); do
   COV="${SAMPLE}.cov"
   if [[ -f "${COV}" ]]; then
      COV_RESULT=()
      for GOI in "${GOIs[@]}"; do
         COV_RESULT+=("`grep ${GOI} ${COV} | awk '{ sum+=$3 } END { if (NR > 0) print sum/NR; else print "0" }'`")
      done
      echo "${SAMPLE} ${COV_RESULT[@]}" >> ${COV_RESULTS}
   fi
done

set -eo pipefail

# Identify presence of aroE252G
ARO_RESULTS="aro.results"
echo "Identifying presence of aroE252G"
if [ -f "${ARO_RESULTS}" ]; then
   echo "Removing old results file.."
   rm ${ARO_RESULTS}
fi
echo "Sample aroE252G" > ${ARO_RESULTS}
for SAMPLE in $(cat "${STRAINLIST}"); do
   echo "${SAMPLE} `!{params.scripts}/fastaRegexFinder.py -f ${SAMPLE}_consensus.fa -q -r 'atgaaatttgcagttatcggaaatcctatttcacattccttgtcgcccgttatgcatagagcaaattttaattctttaggattagatgatacttatgaagctttaaatattccaattgaagattttcatttaattaaagaaattatttcaaaaaaagaattagatggctttaatatcacaattcctcataaagaacgtatcataccgtatttagattatgttgatgaacaagcgattaatgcaggtgcagttaacactgttttgataaaagatggcaagtggatagggtataatacagatggtatcggttatgtgaaaggattgcacagcgtttatccagatttagaaaatgcatacattttaattttgggcgcaggtggtgcaagtaaaggtattgcttatgaattagcaaaatttgtaaagcccaaattaactgttgcgaatagaacgatggctcgttttgaatcttggaatttaaatataaaccaaatttcattggcagatgctgaaaagtatttagctgaattcgatattgttattaatacaacaccagcgggtatggctggaaataacgaaagtattattaatttaaagcatctttctcccaatactttaatgagtgacattatttatataccgtataaaacacctattttagaggaagcagagcgcaagggaaaccatatttataatggcttagatatgtttgtttaccagggtgcggaaagctttaaaatttggactaataaagatgctgatattaattctatgaaaacagcagttttacaacaattaaaaggagaataa' | awk '{ if ($1=="aroE") print "Present"}END { if (NR==0) print "Absent"; }'`"
done >> ${ARO_RESULTS}

# perform SCCmec typing
SCCMEC_RESULTS="sccmec.results"
GOIs=("aroE" "lukF" "lukS" "cna" "sdrE" "pUB110" "pT181")
echo "Performing SCCmec typing.."
if [ -f "${SCCMEC_RESULTS}" ]; then
   echo "Removing old results file.."
   rm ${SCCMEC_RESULTS}
fi

echo "Sample SCCmec_Type" > "${SCCMEC_RESULTS}"

source activate sr_sccmec

for SAMPLE in $(cat "${STRAINLIST}"); do
   FILE="${SAMPLE}.contigs_filtered.fasta"
   echo "${SAMPLE} `python "!{params.scripts}/SR_sccmec_caller.py" -i "${FILE}" -k 90 -l 60 -o "${SAMPLE}.sccmec" -d "${SAMPLE}" -db_dir "!{params.db}/Staph/sccmecfinder_db" -sc_dir "!{params.scripts}/script_dir" -db_choice extended`"
done >> ${SCCMEC_RESULTS}
  '''
}

process StaphTyperReport {
  label "general"
  publishDir "${params.results}/Summaries"

  input:
  path "*" from staphSummaries
  path "*" from staphListReport
  path("${params.run_id}.MLST.summary.tab") from stSummaryStaph

  output:
  path "*.Styper.*"
  val true into staphDone

  shell:
  '''
STRAINLIST="Staphylococcus aureus.Species.txt"
LOGFILE="!{params.run_id}.Styper.log"

source activate sr_python

python "!{params.scripts}/SR_saureus_report.py" -s "${STRAINLIST}" -t sccmec.results -a aro.results -c cov.results -m "!{params.run_id}.MLST.summary.tab" -l "${LOGFILE}" -o "!{params.run_id}.Styper.csv"
  '''
}

process FinalProcess {
  label "general"
  publishDir "${params.results}"

  input:
  val trimFlag from trimDone.collect().ifEmpty('false')
  val listFlag from listDone.collect().ifEmpty('false')
  val krakenFlag from krakenDone.collect().ifEmpty('false')
  val spadesFlag from spadesDone.collect().ifEmpty('false')
  val abricateFlag from abricateDone.collect().ifEmpty('false')
  val stFlag from stDone.collect().ifEmpty('false')
  val redcapFlag from redcapDone.collect().ifEmpty('false')
  val abResFlag from abResDone.collect().ifEmpty('false')
  val intrinsicFlag from intrinsicDone.collect().ifEmpty('false')
  val nesoniFlag from nesoniDone.collect().ifEmpty('false')
  val pseudoFlag from pseudoDone.collect().ifEmpty('false')
  val clusterFlag from clusterDone.collect().ifEmpty('false')
  val staphFlag from staphDone.collect().ifEmpty('false')
  val amrFlag from amrDone.collect().ifEmpty('false')
  val amrSFlag from amrSDone.collect().ifEmpty('false')
  val qcFlag from qcDone.collect().ifEmpty('false')
  val fastQCFlast from fastQCDone.collect().ifEmpty('false')

  output:
  path "run.success"

  shell:
  '''
touch "run.success"
  '''
}

