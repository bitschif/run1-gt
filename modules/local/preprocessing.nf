process PREPROCESSING {
  tag 'preprocess'

  input:
    val token

  output:
    val token

  script:
    """
    bash ${projectDir}/bin/02_preprocessing.sh
    """
}
