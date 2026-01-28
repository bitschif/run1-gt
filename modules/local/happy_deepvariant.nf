process HAPPY_DEEPVARIANT {
  tag 'deepvariant-happy'

  input:
    val token

  output:
    val token

  script:
    """
    bash ${projectDir}/bin/04_variant_calling_deepvariant_happy.sh
    """
}
