process HAPPY_DEEPVARIANT {
  tag 'deepvariant-happy'

  input:
    val meta

  output:
    val meta

  script:
    """
    export SAMPLE_NAME='${meta.sample}'
    export READGROUP_ID='${meta.readgroup_id}'
    export PLATFORM='${meta.platform}'
    export LIBRARY='${meta.library}'
    export LANE='${meta.lane}'
    export CENTER='${meta.center}'
    export INSTRUMENT='${meta.instrument}'
    bash ${projectDir}/bin/04_variant_calling_deepvariant_happy.sh
    """
}
