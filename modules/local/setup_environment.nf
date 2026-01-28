process SETUP_ENVIRONMENT {
  tag 'setup'

  input:
    val token

  output:
    val token

  script:
    """
    bash ${projectDir}/bin/00_setup_environment.sh
    """
}
