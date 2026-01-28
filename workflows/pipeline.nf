include { SETUP_ENVIRONMENT } from '../modules/local/setup_environment'
include { SIMULATE_DATA } from '../modules/local/simulate_data'
include { PREPROCESSING } from '../modules/local/preprocessing'
include { CALL_GATK } from '../modules/local/call_gatk'
include { HAPPY_GATK } from '../modules/local/happy_gatk'
include { CALL_DEEPVARIANT } from '../modules/local/call_deepvariant'
include { HAPPY_DEEPVARIANT } from '../modules/local/happy_deepvariant'
include { CALL_STRELKA2 } from '../modules/local/call_strelka2'
include { HAPPY_STRELKA2 } from '../modules/local/happy_strelka2'
include { CALL_FREEBAYES } from '../modules/local/call_freebayes'
include { HAPPY_FREEBAYES } from '../modules/local/happy_freebayes'

workflow variant_benchmarking {
  def callers_list = params.callers instanceof String ? params.callers.tokenize(',') : params.callers

  def start_ch = Channel.value(1)
  def setup_ch = params.run_setup ? SETUP_ENVIRONMENT(start_ch) : start_ch
  def simulate_ch = params.run_simulation ? SIMULATE_DATA(setup_ch) : setup_ch
  def preprocess_ch = params.run_preprocessing ? PREPROCESSING(simulate_ch) : simulate_ch

  if (params.run_calling && callers_list.contains('gatk')) {
    def gatk_ch = CALL_GATK(preprocess_ch)
    if (params.run_benchmark) {
      HAPPY_GATK(gatk_ch)
    }
  }

  if (params.run_calling && callers_list.contains('deepvariant')) {
    def deepvariant_ch = CALL_DEEPVARIANT(preprocess_ch)
    if (params.run_benchmark) {
      HAPPY_DEEPVARIANT(deepvariant_ch)
    }
  }

  if (params.run_calling && callers_list.contains('strelka2')) {
    def strelka2_ch = CALL_STRELKA2(preprocess_ch)
    if (params.run_benchmark) {
      HAPPY_STRELKA2(strelka2_ch)
    }
  }

  if (params.run_calling && callers_list.contains('freebayes')) {
    def freebayes_ch = CALL_FREEBAYES(preprocess_ch)
    if (params.run_benchmark) {
      HAPPY_FREEBAYES(freebayes_ch)
    }
  }
}
