process CALL_GATK {
  tag 'gatk'

  input:
    val token

  output:
    val token

  script:
    """
    bash ${projectDir}/bin/03_variant_calling_gatk.sh
    """
}
