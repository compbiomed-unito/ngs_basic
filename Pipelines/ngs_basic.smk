
# config fields: scratch_root


# SAMPLE SHEET HANDLING #
def make_fastq_table(dirs, pattern=None):
	'''Find fastq files in some directories and assemble them in a read sheet

	dirs -- iterator of pairs of name and path to a directory
	pattern -- a pattern with named group (at least one called SAMPLE)
	
	If pattern is not given it defaults to:
	(?P<SAMPLE>.+)_(?P<INDEX>S[1-9][0-9]*)_(?P<PU>L[0-9]+)_(?P<FRAG>R[12])_001.fastq.gz
	'''
	import os
	import re
	import pandas

	if pattern is None:
		pattern = r'''(?P<SAMPLE>.+)_(?P<INDEX>S[1-9][0-9]*)_(?P<PU>L[0-9]+)_(?P<FRAG>R[12])_001.fastq.gz'''
	pattern = re.compile(pattern)
	fastq_acc = []
	for run_name, dir_path in dirs:
		for file_name in os.listdir(dir_path):
			m = re.match(pattern, file_name)
			if m:
				d = m.groupdict()
				d['RUN'] = run_name
				d['FILENAME'] = file_name
				d['PATH'] = dir_path + '/' + file_name
				fastq_acc.append(d)
	df = pandas.DataFrame(fastq_acc).sort_values('SAMPLE')
	return df.loc[:, ['SAMPLE'] + df.columns.drop('SAMPLE').tolist()] # put SAMPLE as first column without making it the index since SAMPLE is in general not unique


def fastq2workdir(read_sheet, wdir=config['scratch_root']):
	'''Create links to the original fastq in a working directory and produce a sample sheet with a single entry for each sample.
	Files are linked as WDIR/RUN/SAMPLE.PU.FRAG.fastq.gz

	read_sheet -- a read sheet as produced by fastq2sheet
	wdir -- path to working directory, defaults to config['scratch_root']
	
	The sample sheet must contain the columns: SAMPLE, PU, FRAG, RUN
	'''
	import os.path

	read_sheet['WORKING_PATH'] = read_sheet.apply(lambda r: '{}/{RUN}/{SAMPLE}.{FRAG}.{PU}.fastq.gz'.format(wdir, **r), axis=1)
	assert not read_sheet['WORKING_PATH'].duplicated().any()

	# make the links
	linked_count = 0
	for s in read_sheet.itertuples():
		src_path = os.path.abspath(s.PATH)
		dst_path = s.WORKING_PATH
		if not os.path.exists(dst_path):
			os.makedirs(os.path.dirname(dst_path), mode=0o770, exist_ok=True)
			os.symlink(src_path, dst_path)
			print('Created link', dst_path, '->', src_path)
			linked_count += 1
	print('Created', linked_count, 'new links')

	# make the sample sheet
	g = read_sheet.groupby('SAMPLE', as_index=False)
	sample_sheet = g[['RUN']].first()

	# check for samples in more than one run
	for sample, runs in g['RUN']:
		if not all(runs['RUN'] == runs['RUN'].iat[0]):
			raise ValueError('Sample {} belongs to more than one run'.format(sample))

	# add workdir path
	sample_sheet['PATH'] = sample_sheet.apply(lambda r: '{}/{RUN}/{SAMPLE}'.format(wdir, **r), axis=1)
	
	return sample_sheet, read_sheet


if False:
	def sample_sheet_input(ss_path, ct_ext):
		'''create input dict for rules that aggregate samples from sample sheet
		ss_path -- path to the sample sheet tsv

		The sample sheet must contain two columns, SAMPLE and WORKING_PATH
		'''
		import pandas
		ss = pandas.read_csv(ss_path, sep='\t', index_col='SAMPLE', comment='#')
		if ss.index.duplicated().any():
			raise ValueError('Duplicated samples in "{}"'.format(ss_path))
		sample_paths = {s: p + ct_ext for s, p in ss['WORKING_PATH'].items()}
		assert 'SS' not in sample_paths
		sample_paths['SS'] = ss_path
		return sample_paths

###################
# READ PROCESSING #
###################
wildcard_constraints:
	FRAG='R[12]',

def _get_all_lanes(w):
	glob_str = '{SAMPLE}.{FRAG}.{{PU}}.fastq.gz'.format(**w)
	glob_matches = glob_wildcards(glob_str)
	return expand(glob_str, **glob_matches._asdict()) #FIXME this is probably a private method


# create a .fastq pipe, to avoid a useless recompression
rule merge_lanes:
	input:
		_get_all_lanes
	output:
		pipe('{SAMPLE}.{FRAG}.fastq')
	threads: 0.2
	shell:
		'''zcat {input} > {output}'''

# /scratch/shared/ngs_basic/workdir/RUN_180611/TO535POSTIDE.R1.L001.fastq.gz
# /scratch/shared/ngs_basic/workdir/RUN_180611/TO535POSTIDE.R1.L001_fastqc.html
# /scratch/shared/ngs_basic/workdir/RUN_180611/TO535POSTIDE.R1.L001_fastqc.zip

rule fastqc:
	input:
		'/{sample}.fastq.gz'
	output:
		'qc/{sample,}'
	shell:
		'''fastqc --outdir qc {input}'''

############
# CUTADAPT #
############
# Raw read cleaning, especially adapter removal
# Rules for single and pair end
wildcard_constraints:
	a3='[ACGT]*', # adapter sequence
	ml='\d+', # minimum length, integer
	ol='\d+', # minimum overlap, integer

cutadaptpe_stem = '{file}.CAse_a3={a3}_ml={ml}_ol={ol}'
cutadapt_cmd = 'cutadapt {params} --minimum-length={wildcards.ml} --overlap={wildcards.ol} --trim-n --error-rate=0.1'

rule cutadapt_single_end:
	input:
		'{file}.R1.fastq'
	output:
		cutadaptpe_stem + '.R1.fastq.gz', cutadaptpe_stem + '_short.R1.fastq.gz'
	log: cutadaptpe_stem + '.log'
	params: lambda w: '--adapter=' + w.a3 if w.a3 else '' # cutadapt complains if adapter is empty
	#threads: 4 # cutadapt does not support multithreading with --too-short-output...
	shell: cutadapt_cmd + ' --output={output[0]} --too-short-output={output[1]} {input[0]} > {log} 2>&1'

cutadaptpe_stem = cutadaptpe_stem.replace('CAse_', 'CApe_')

rule cutadapt_pair_end:
	input:
		'{file}.R1.fastq',
		'{file}.R2.fastq',
	output:
		[cutadaptpe_stem + suf for suf in ['.R1.fastq.gz', '_short.R1.fastq.gz', '.R2.fastq.gz', '_short.R2.fastq.gz']]
	log:
		cutadaptpe_stem + '.log',
	params:
		lambda w: '--adapter=' + w.a3 if w.a3 else '' # cutadapt complains if adapter is empty
	#threads: 4 # cutadapt does not support multithreading with --too-short-output...
	shell:
		cutadapt_cmd + ''' --pair-filter=any --output={output[0]} --too-short-output={output[1]} --paired-output={output[2]} --too-short-paired-output={output[3]} {input[0]} {input[1]} > {log} 2>&1'''


# SAM/BAM

rule bam_sort:
	input:
		'{file}.bam'
	output:
		'{file}.sorted.bam'
	shell:
		'samtools sort {input[0]} -o {output[0]}'

rule bam_index:
	input:
		'{file}.bam'
	output:
		'{file}.bam.bai'
	shell:
		'samtools index {input[0]} {output[0]} 2> {output[0]}.log'


#######
# BWA #
#######
rule bwa_index:
	input:
		'{file}'
	output:
		directory('{file}.bwa_index')
	shell:
		'''mkdir -p {output[0]} && bwa index -p {output[0]}/bwa_index {input[0]} 2> {output[0]}.log'''

import yaml
with open(config['resources'], 'rt') as f:
	resource_config = yaml.load(f)

def get_resource(refname, restype): # FIXME move to ngs_basic and use this instead of reference_sequences
	if refname not in resource_config:
		msg = 'Reference {} not found. Available references: {}'.format(refname, ', '.join(resource_config))
		raise KeyError(mgs)
	ref_resources = resource_config[refname]
	if restype not in ref_resources:
		msg = 'Resource {} not found for reference {}. Available resources: {}'.format(restype, refname, ', '.join(ref_resources))
		raise KeyError(mgs)
	return ref_resources[restype]

wildcard_constraints:
	ref='|'.join(resource_config.keys())


#reference_sequences = {name: res['sequences'] for name, res in resource_config.items() if 'sequences' in res}

# bwa mem alignment for single end reads, with parameters for handling short alignments in sncRNA pipelines
# sl: seed length, ml: minimum length (score), out: bh (best hit only) or all (all alignments), ref: reference
bwamemse_stem = '{file}.bwamemse_sl={k,\d+}_ml={T,\d+}_out={out,bh|all}_ref={ref}.bam'

rule bwa_mem_single_end:
	input:
		#ref=lambda w: reference_sequences[w.ref] + '.bwa_index',
		ref=lambda w: get_resource(w.ref, 'sequences') + '.bwa_index',
		reads='{file}.R1.fastq.gz'
	output:
		bwamemse_stem
	log:
		bwamemse_stem + '.log'
	benchmark:
		bwamemse_stem + '.benchmark.txt'
	threads: 8
	params:
		all_opt=lambda w: '-a' if w.out == 'all' else ''
	shell: r'''bwa mem -t {threads} -T {wildcards.T} -A 1 -k 10 {params.all_opt} {input.ref}/bwa_index {input.reads} 2> {log} | samtools view -o {output}'''


# bwa mem alignment for pair end reads
bwamempe_stem = '{file}.bwamempe_ref={ref}.bam'

rule bwa_mem_pair_end:
	input:
		ref=lambda w: get_resource(w.ref, 'sequences') + '.bwa_index',
		reads=['{file}.R1.fastq.gz', '{file}.R2.fastq.gz']
	output:
		bwamempe_stem
	log:
		bwamempe_stem + '.log'
	benchmark:
		bwamempe_stem + '.benchmark.txt'
	threads: 8
	shell:
		r'''bwa mem -t {threads} -R '@RG\tID:S1\tSM:S1\tPU:unknown\tPL:ILLUMINA' {input.ref}/bwa_index {input.reads[0]} {input.reads[1]} 2> {log} | samtools view -o {output}'''

# VCF

rule vcf_index:
	input:
		'{file}.vcf.gz'
	output:
		'{file}.vcf.gz.csi'
	shell:
		'''tabix --csi --preset vcf {input} > {output}.log 2>&1'''



