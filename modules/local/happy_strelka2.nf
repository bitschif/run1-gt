process HAPPY_STRELKA2 {
  tag 'strelka2-happy'

  input:
    val token

  output:
    val token

  script:
    """
    bash ${projectDir}/bin/05_variant_calling_strelka2_happy.sh
    """
}
