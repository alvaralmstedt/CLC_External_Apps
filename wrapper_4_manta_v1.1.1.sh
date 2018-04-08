#!/bin/bash -l

MANTACONFIG=/apps/CLC_ExternalApps/manta_1.1.1/manta-1.1.1.centos5_x86_64/bin/configManta.py
DATE=$(date +"%Y_%m_%d_%H_%M_%S")
TMPOUT=/tmp/manta_$DATE
HELP="Help message"
SAMTOOLS=/apps/bio/apps/samtools/1.3.1/samtools
BAM_FIXER=/apps/CLC_ExternalApps/manta_1.1.1/bam_fixer_2.py
INCOMPLETE_FORMATTER=/apps/CLC_ExternalApps/manta_1.1.1/format_incomplete_manta.sh
echo "Hostname: $HOSTNAME"

while getopts :b:c:f:o:p:q:r:s:t:e:a:u:v:x:n:m:1:2:h opt; do
  case $opt in
    b)
        echo "-b (bam) was input as $OPTARG" >&1
        BAM=$OPTARG
    ;;
    c)
        echo "-c (bam text) was input as $OPTARG" >&1
        TUMOUR_TEXT=$OPTARG
    ;;
    f)
        echo "-f (reference fasta) was input as $OPTARG" >&1
        FASTA=$OPTARG
    ;;
    o)
        echo "-o (output, conactenated) was input as $OPTARG" >&1
        OUTPUT=$OPTARG
    ;;
    p)
        echo "-p (output 2, candidateSV) was input as $OPTARG" >&1
        OUTPUT2=$OPTARG
    ;;
    q)
        echo "-q (output 3, candidateSmallIndels) was input as $OPTARG" >&1
        OUTPUT3=$OPTARG
    ;;
    r)
        echo "-r (output 4, diploidSV) was input as $OPTARG" >&1
        DIPLOID=$OPTARG
    ;;
    s)
        echo "-s (output 5, somaticSV) was input as $OPTARG" >&1
        OUTPUT5=$OPTARG
    ;;
    n)
        echo "-n (tumour + normal) was input as $OPTARG" >&1
        TUMOURNORMAL=$OPTARG
    ;;
    m)
        echo "-m (normal text) was input as $OPTARG" >&1
        NORMAL_TEXT=$OPTARG
    ;;
    t)
        echo "-t (tumour only) was input as $OPTARG" >&1
        TUMOURONLY=$OPTARG
    ;;
    e)
        echo "-e (exome selector) was input as $OPTARG" >&1
        EX=$OPTARG
    ;;
    a)
        echo "-a (Exome manifest bed) was input as $OPTARG" >&1
        EXOME=$OPTARG
    ;;
    u)
        echo "-u (dip_bool) was input as $OPTARG" >&1
        DIP=$OPTARG
    ;;
    v)
        echo "-v (tumour_normal bool) is $OPTARG"
        TN=$OPTARG
    ;;
    x)
        echo "-x (tumour only bool) is $OPTARG"
        TO=$OPTARG
    ;;
    1)
        echo "-1 (incomplete deletions) was input as $OPTARG"
        INCDEL=$OPTARG    
    ;;
    2)
        echo "-2 (incomplete insertions) was input as  $OPTARG"
        INCINS=$OPTARG
    ;;
    h)
        echo "$HELP"
        exit 1
    ;;
    \?)
        echo "Invalid option: -$OPTARG" >&1
        echo "Type $0 -h for usage"
        exit 1
    ;;
  esac
done

THREADS=$(nproc)

if [[ $FASTA == *"hg19"* ]] ; then
    echo "Using default hg19"
    FASTA=/medstore/External_References/hg19/Homo_sapiens_sequence_hg19.fasta
else
    echo "Running samtools faidx on $FASTA"
    $SAMTOOLS faidx $FASTA
fi

if [ ! -z $TUMOUR_TEXT  ] && [ ! -z $NORMAL_TEXT ] ; then 
    BAM=$(cat $TUMOUR_TEXT)
    TUMOURNORMAL=$(cat $NORMAL_TEXT)
    TN_TEXT=True
    echo "Tumour bam (from text): $BAM"
    echo "Normal bam (from text): $TUMOURNORMAL"
elif [ ! -z $TUMOUR_TEXT ] ; then
    BAM=$(cat $TUMOUR_TEXT)
    echo "Tumour only bam (from text): $BAM"
elif [ ! -z $NORMAL_TEXT ] ; then
    TUMOURNORMAL=$(cat $NORMAL_TEXT)
    echo "Normal from text (bam nontext): $TUMOURNORMAL"
#elif [ ! -z $TUMOUR_TEXT ] ; then
#    BAM=$(cat $TUMOUR_TEXT)
fi

#mkdir $TMPOUT
#mkdir $TMPOUT/analysis
#mkdir $OUTPUT

module load samtools/1.3.1

mkdir -p $TMPOUT/bam_fixing
# Move bam out from black box
# if [[ $BAM == *"CLC_External_Tmp"* ]] ; then
echo "Moving $BAM to $TMPOUT/bam_fixing/bam_from_CLC.bam"
cp ${BAM} $TMPOUT/bam_fixing/bam_from_CLC.bam
wait
BAM=$TMPOUT/bam_fixing/bam_from_CLC.bam
# fi

$BAM_FIXER ${BAM} $TMPOUT/bam_fixing/Perfect.sam $TMPOUT/bam_fixing/Secondary.sam $TMPOUT/bam_fixing/bam_fixer_2.log
echo "bam_fixer_2.py completed"
echo "Bam_fixer (main bam file) log:"
cat $TMPOUT/bam_fixing/bam_fixer_2.log
mv $TMPOUT/bam_fixing/merged_sorted.bam $TMPOUT/merged_sorted.bam
mv $TMPOUT/bam_fixing/merged_sorted.bam.bai $TMPOUT/merged_sorted.bam.bai
BAM=$TMPOUT/merged_sorted.bam
rm -rf $TMPOUT/bam_fixing
wait
#$SAMTOOLS sort ${BAM} -@ 25 -m 2G -o ${BAM}_sorted && mv ${BAM}_sorted ${BAM}
#$SAMTOOLS index ${BAM}

#if [[ ! -z $TUMOURNORMAL ]] ; then
if [[ "$TN" = "Enabled" ]] ; then
#    $SAMTOOLS sort ${TUMOURNORMAL} -@ $(($THREADS/4*3)) -m 2G -o ${TUMOURNORMAL}_sorted && mv ${TUMOURNORMAL}_sorted $TUMOURNORMAL
#$SAMTOOLS index ${TUMOURNORMAL}
    echo "Running tumour+normal analysis"
    mkdir $TMPOUT/bam_fixing_normal
    # Move normal bam out from black box
    # if [[ $TUMOURNORMAL == *"CLC_External_Tmp"* ]] ; then
    echo "Moving $TUMOURNORMAL to $TMPOUT/bam_fixing_normal/normal_from_CLC.bam"
    cp ${TUMOURNORMAL} $TMPOUT/bam_fixing_normal/normal_from_CLC.bam
    wait
    TUMOURNORMAL=$TMPOUT/bam_fixing_normal/normal_from_CLC.bam
    # fi
    $BAM_FIXER ${TUMOURNORMAL} $TMPOUT/bam_fixing_normal/Perfect.sam $TMPOUT/bam_fixing_normal/Secondary.sam $TMPOUT/bam_fixing_normal/bam_fixer_2.log
    echo "Bam_fixer (normal bam file) log:"
    cat $TMPOUT/bam_fixing_normal/bam_fixer_2.log
    mv $TMPOUT/bam_fixing_normal/merged_sorted.bam $TMPOUT/merged_sorted_normal.bam
    mv $TMPOUT/bam_fixing_normal/merged_sorted.bam.bai $TMPOUT/merged_sorted_normal.bam.bai
    TUMOURNORMAL=$TMPOUT/merged_sorted_normal.bam
    rm -rf $TMPOUT/bam_fixing_normal
    ${MANTACONFIG} \
--normalBam ${TUMOURNORMAL} \
--tumorBam ${BAM} \
--referenceFasta ${FASTA} \
--runDir $TMPOUT/analysis

#elif [[ ! -z $TUMOURONLY ]] ; then
elif [[ "$TO" = "Enabled" ]] && [[ -z "$EX" ]] ; then
    echo "Running Tumour-only analysis"
    ${MANTACONFIG} \
--tumorBam ${BAM} \
--referenceFasta ${FASTA} \
--runDir $TMPOUT/analysis

#else
elif [[ "$DIP" = "Enabled" ]] && [[ -z "$EX" ]] ; then
#DIP=True
    echo "Running single diploid sample"
    ${MANTACONFIG} \
--bam ${BAM} \
--referenceFasta ${FASTA} \
--runDir $TMPOUT/analysis

elif [[ "$EX" = "Enabled" ]] && [[ "$TO" = "Enabled" ]] ; then
    echo "Preparing to run exome on tumour only sample"
    module load tabix/0.2.6
    cp $EXOME $TMPOUT/exome_manifest.bed
    EXOME=$TMPOUT/exome_manifest.bed
    bgzip $EXOME
    EXOME=${EXOME}.gz
    tabix -0 ${EXOME}
    echo "Exome manifest name before running manta: $EXOME"
    echo "Running exome sample"
    ${MANTACONFIG} \
--tumorBam ${BAM} \
--exome \
--callRegions ${EXOME} \
--referenceFasta ${FASTA} \
--runDir $TMPOUT/analysis

elif [[ "$EX" = "Enabled" ]] && [[ "$DIP" = "Enabled" ]] ; then
    echo "Preparing to run exome on diploid sample"
    module load tabix/0.2.6
    cp $EXOME $TMPOUT/exome_manifest.bed
    EXOME=$TMPOUT/exome_manifest.bed
    bgzip $EXOME
    EXOME=${EXOME}.gz
    tabix -0 ${EXOME}
    echo "Exome manifest name before running manta: $EXOME"
    echo "Running exome sample"
    ${MANTACONFIG} \
--bam ${BAM} \
--exome \
--callRegions ${EXOME} \
--referenceFasta ${FASTA} \
--runDir $TMPOUT/analysis

elif [[ "$EX" = "Enabled" ]] ; then
    echo "Preparing to run exome on diploid sample by default. You forgot to select tumour or diploid"
    module load tabix/0.2.6
    cp $EXOME $TMPOUT/exome_manifest.bed
    EXOME=$TMPOUT/exome_manifest.bed
    bgzip $EXOME
    EXOME=${EXOME}.gz
    tabix -0 ${EXOME}
    echo "Exome manifest name before running manta: $EXOME"
    echo "Running exome sample"
    ${MANTACONFIG} \
--bam ${BAM} \
--exome \
--callRegions ${EXOME} \
--referenceFasta ${FASTA} \
--runDir $TMPOUT/analysis

else
    echo "no analysis type checked"
    exit 1
fi

echo "ls TMPOUT"
ls $TMPOUT

MANTADIR=$TMPOUT/analysis
echo "Starting manta at: `date`"
$MANTADIR/runWorkflow.py -m local -j $THREADS
wait
echo "Manta finished at: `date`"
echo "ls MANTADIR ($MANTADIR)"
ls $MANTADIR
gunzip $MANTADIR/results/variants/*vcf.gz

SWITCH=True
for i in $MANTADIR/results/variants/*.vcf ; do
     if [ $SWITCH = True ] ; then
         grep  "^#" $i > $TMPOUT/concatenated_mantaoutput_final.vcf
         SWITCH=False
     else
         grep -v "^#" $i >> $TMPOUT/concatenated_mantaoutput.vcf
     fi
done

echo "ls MANTADIR/results/variants"
ls $MANTADIR/results/variants

cat $TMPOUT/concatenated_mantaoutput.vcf | sort -u >> $TMPOUT/concatenated_mantaoutput_final.vcf

# Send all the variants that CLC cannot read to plain text files
grep "<INS>" $TMPOUT/concatenated_mantaoutput_final.vcf > $TMPOUT/incomplete_insertions_and_duplications.txt
grep "<DUP:TANDEM>" $TMPOUT/concatenated_mantaoutput_final.vcf >> $TMPOUT/incomplete_insertions_and_duplications.txt
# Put call to convert to csv script here
$INCOMPLETE_FORMATTER $TMPOUT/incomplete_insertions_and_duplications.txt $TMPOUT
cp $TMPOUT/incomplete_insertions_and_duplications.bed $INCINS
grep "<DEL>" $TMPOUT/concatenated_mantaoutput_final.vcf > $TMPOUT/incomplete_deletions_and_breakends.txt
grep "MantaBND" $TMPOUT/concatenated_mantaoutput_final.vcf >> $TMPOUT/incomplete_deletions_and_breakends.txt
# Put call to convert csv script here
$INCOMPLETE_FORMATTER $TMPOUT/incomplete_deletions_and_breakends.txt $TMPOUT
cp $TMPOUT/incomplete_deletions_and_breakends.bed $INCDEL

cp $TMPOUT/*.vcf /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/large_mantaoutput
cp $MANTADIR/results/variants/*.vcf /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/large_mantaoutput

rm $TMPOUT/concatenated_mantaoutput.vcf
echo "Copying $TMPOUT/concatenated_mantaoutput_final.vcf to $OUTPUT"
cp $TMPOUT/concatenated_mantaoutput_final.vcf $OUTPUT

#if [[ ! -z $TUMOURNORMAL ]] ; then
if [[ ! -z $MANTADIR/results/variants/candidateSV.vcf ]] ; then
    echo "Copying $MANTADIR/results/variants/candidateSV.vcf to $OUTPUT2"
    cp $MANTADIR/results/variants/candidateSV.vcf $OUTPUT2
else
    echo "Copying dummy to $OUTPUT2"
    cp /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf $OUTPUT2
fi

if [[ ! -z $MANTADIR/results/variants/candidateSmallIndels.vcf ]] ; then
    echo "Copying $MANTADIR/results/variants/candidateSmallIndels.vcf to $OUTPUT3"
    cp $MANTADIR/results/variants/candidateSmallIndels.vcf $OUTPUT3
else
    echo "Copying dummy to $OUTPUT3"
    cp /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf $OUTPUT3
fi

#mv $MANTADIR/results/variants/diploidSV.vcf $OUTPUT4
#mv $MANTADIR/results/variants/somaticSV.vcf $OUTPUT5
#fi

#if [[ ! -z $TUMOURONLY ]] ; then
if [ "$TO" = "Enabled" ] ; then
    if [[ ! -z "$(grep -v "^#" $MANTADIR/results/variants/tumorSV.vcf)" ]] ; then
        echo "Copying $MANTADIR/results/variants/tumorSV.vcf to $TUMOURONLY"
        cp $MANTADIR/results/variants/tumorSV.vcf $TUMOURONLY
    else
        echo "No variants found in $MANTADIR/results/variants/tumourSV.vcf"
        echo "Copying /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf to $TUMOURONLY"
        cp /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf $TUMOURONLY
    fi
    cp /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf $DIPLOID
    cp /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf $OUTPUT5
elif [ "$DIP" = "Enabled" ] ; then
    if [[ ! -z "$(grep -v "^#" $MANTADIR/results/variants/diploidSV.vcf)" ]] ; then
        echo "Copying $MANTADIR/results/variants/diploidSV.vcf to $DIPLOID"
        cp $MANTADIR/results/variants/diploidSV.vcf $DIPLOID
    else
        echo "No variants found in $MANTADIR/results/variants/diploidSV.vcf"
        echo "Copying /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf to $DIPLOID"
        cp /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf $DIPLOID
    fi
    cp /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf $TUMOURONLY
    cp /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf $OUTPUT5
#elif [[ ! -z $TUMOURNORMAL ]] ; then
elif [ "$TN" = "Enabled" ] ; then
    if [[ ! -z "$(grep -v "^#" $MANTADIR/results/variants/somaticSV.vcf)" ]] ; then
        echo "Copying $MANTADIR/results/variants/somaticSV.vcf to $OUTPUT5"
        cp $MANTADIR/results/variants/somaticSV.vcf $OUTPUT5
    else
        echo "No variants found in $MANTADIR/results/variants/somaticSV.vcf"
        echo "Copying /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf to $OUTPUT5"
        cp /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf $OUTPUT5
    fi
    cp /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf $TUMOURONLY
    cp /medstore/CLC_Import_Export/Alvar_Almstedt/homebrewed/manta_dummy/dummy.vcf $DIPLOID
fi
wait
#rm -rf $TMPOUT
echo "Folder $TMPOUT deleted"
