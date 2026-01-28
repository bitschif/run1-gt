process CALL_FREEBAYES {
  tag 'freebayes'

  input:
    val token

  output:
    val token

  script:
    """
    bash ${projectDir}/bin/06_variant_calling_freebayes.sh
    """
}
