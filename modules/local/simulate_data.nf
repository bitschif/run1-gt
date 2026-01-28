process SIMULATE_DATA {
  tag 'simulate'

  input:
    val token

  output:
    val token

  script:
    """
    bash ${projectDir}/bin/01_simulate_data.sh
    """
}
