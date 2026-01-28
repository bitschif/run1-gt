process CALL_STRELKA2 {
  tag 'strelka2'

  input:
    val token

  output:
    val token

  script:
    """
    bash ${projectDir}/bin/05_variant_calling_strelka2.sh
    """
}
