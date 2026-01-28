process HAPPY_GATK {
  tag 'gatk-happy'

  input:
    val token

  output:
    val token

  script:
    """
    bash ${projectDir}/bin/03_variant_calling_gatk_happy.sh
    """
}
