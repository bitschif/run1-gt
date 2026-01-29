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

def parse_samplesheet(samplesheet_path) {
  Channel.fromPath(samplesheet_path)
    .ifEmpty { error "Samplesheet not found: ${samplesheet_path}" }
    .splitCsv(header: true)
    .map { row ->
      def sample = row.sample ?: row.sample_id ?: row.sample_name
      if (!sample) {
        error "Samplesheet row is missing required 'sample' column: ${row}"
      }
      def fastq1 = row.fastq_1 ?: row.fastq1 ?: row.read1
      if (!fastq1) {
        error "Samplesheet row is missing required 'fastq_1' column: ${row}"
      }
      def fastq2 = row.fastq_2 ?: row.fastq2 ?: row.read2
      if (!fastq2) {
        error "Samplesheet row is missing required 'fastq_2' column: ${row}"
      }
      def meta = [
        sample: sample.toString(),
        fastq_1: fastq1.toString(),
        fastq_2: fastq2.toString(),
        readgroup_id: (row.readgroup_id ?: row.read_group ?: sample).toString(),
        platform: (row.platform ?: params.platform ?: 'ILLUMINA').toString(),
        library: (row.library ?: params.library ?: 'lib1').toString(),
        lane: (row.lane ?: params.lane ?: '1').toString(),
        center: (row.center ?: '').toString(),
        instrument: (row.instrument ?: '').toString()
      ]
      return meta
    }
}

workflow variant_benchmarking {
  def tool_list = params.tools ?: params.callers
  def tools_raw = tool_list instanceof String ? tool_list.tokenize(',') : tool_list
  def tool_map = [
    'haplotypecaller': 'gatk',
    'gatk': 'gatk',
    'deepvariant': 'deepvariant',
    'strelka': 'strelka2',
    'strelka2': 'strelka2',
    'freebayes': 'freebayes'
  ]
  def callers_list = tools_raw.collect { tool ->
    def key = tool.toString().trim().toLowerCase()
    tool_map.get(key, key)
  }.unique()

  def sample_ch = params.input ? parse_samplesheet(params.input) : Channel.value([
    sample: params.sample_name ?: 'SIMULATED_SAMPLE',
    fastq_1: '',
    fastq_2: '',
    readgroup_id: params.readgroup_id ?: params.sample_name ?: 'SIMULATED_SAMPLE',
    platform: params.platform ?: 'ILLUMINA',
    library: params.library ?: 'lib1',
    lane: params.lane ?: '1',
    center: params.center ?: '',
    instrument: params.instrument ?: ''
  ])

  def run_simulation = params.run_simulation && !params.input
  def run_benchmark = params.run_benchmark && (!params.input || params.truth_vcf)

  def setup_ch = params.run_setup ? SETUP_ENVIRONMENT(sample_ch) : sample_ch
  def simulate_ch = run_simulation ? SIMULATE_DATA(setup_ch) : setup_ch
  def preprocess_ch = params.run_preprocessing ? PREPROCESSING(simulate_ch) : simulate_ch

  if (params.run_calling && callers_list.contains('gatk')) {
    def gatk_ch = CALL_GATK(preprocess_ch)
    if (run_benchmark) {
      HAPPY_GATK(gatk_ch)
    }
  }

  if (params.run_calling && callers_list.contains('deepvariant')) {
    def deepvariant_ch = CALL_DEEPVARIANT(preprocess_ch)
    if (run_benchmark) {
      HAPPY_DEEPVARIANT(deepvariant_ch)
    }
  }

  if (params.run_calling && callers_list.contains('strelka2')) {
    def strelka2_ch = CALL_STRELKA2(preprocess_ch)
    if (run_benchmark) {
      HAPPY_STRELKA2(strelka2_ch)
    }
  }

  if (params.run_calling && callers_list.contains('freebayes')) {
    def freebayes_ch = CALL_FREEBAYES(preprocess_ch)
    if (run_benchmark) {
      HAPPY_FREEBAYES(freebayes_ch)
    }
  }
}
