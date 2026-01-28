process HAPPY_FREEBAYES {
  tag 'freebayes-happy'

  input:
    val token

  output:
    val token

  script:
    """
    bash ${projectDir}/bin/06_variant_calling_freebayes_happy.sh
    """
}
