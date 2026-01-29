process PREPROCESSING {
  tag 'preprocess'

  input:
    val meta

  output:
    val meta

  script:
    """
    export SAMPLE_NAME='${meta.sample}'
    export FASTQ_1='${meta.fastq_1}'
    export FASTQ_2='${meta.fastq_2}'
    export READGROUP_ID='${meta.readgroup_id}'
    export PLATFORM='${meta.platform}'
    export LIBRARY='${meta.library}'
    export LANE='${meta.lane}'
    export CENTER='${meta.center}'
    export INSTRUMENT='${meta.instrument}'
    bash ${projectDir}/bin/02_preprocessing.sh
    """
}
