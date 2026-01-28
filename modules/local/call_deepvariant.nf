process CALL_DEEPVARIANT {
  tag 'deepvariant'

  input:
    val token

  output:
    val token

  script:
    """
    bash ${projectDir}/bin/04_variant_calling_deepvariant.sh
    """
}
