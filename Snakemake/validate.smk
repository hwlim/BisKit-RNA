'''
Configuration & sample sheet validation

'''

## check if reference files exist
if not os.path.isfile(rRNA_fa):
	print( "Error: rRNA reference FASTA does not exist." )
	sys.exit(1)

rRNA_index_check = '/'.join(rRNA_index.split('/')[:-1])
if not os.path.isdir(rRNA_index_check):
	print( "Error: rRNA reference index does not exist. If you don't have the index file ready, run 'hisat-3n-build --base-change C,T your_FASTA_file your_indexFolder_Path/rRNA' and 'hisat-3n-build --base-change C,T --repeat-index your_FASTA_file your_indexFolder_Path/rRNA'. You need to set the path to your index directory and prefix in the Snakefile." )
	sys.exit(1)
if not os.path.isfile(tRNA_fa):
	print( "Error: tRNA reference FASTA does not exist." )
	sys.exit(1)

tRNA_index_check = '/'.join(tRNA_index.split('/')[:-1])
if not os.path.isdir(tRNA_index_check):
	print( "Error: tRNA reference index does not exist. If you don't have the index file ready, run 'hisat-3n-build --base-change C,T your_FASTA_file your_indexFolder_Path/tRNA' and 'hisat-3n-build --base-change C,T --repeat-index your_FASTA_file your_indexFolder_Path/tRNA'. You need to set the path to your index directory and prefix in the Snakefile." )
	sys.exit(1)
if not os.path.isfile(miRNA_hairpin_fa):
	print( "Error: miRNA reference FASTA does not exist." )
	sys.exit(1)

miRNA_index_check = '/'.join(miRNA_hairpin_index.split('/')[:-1])
if not os.path.isdir(miRNA_index_check):
	print( "Error: miRNA reference index does not exist. If you don't have the index file ready, run 'hisat-3n-build --base-change C,T your_FASTA_file your_indexFolder_Path/miRNA' and 'hisat-3n-build --base-change C,T --repeat-index your_FASTA_file your_indexFolder_Path/miRNA'. You need to set the path to your index directory and prefix in the Snakefile." )
	sys.exit(1)

if not os.path.isfile(genome_fa):
	print( "Error: genome reference FASTA does not exist." )
	sys.exit(1)

genome_index_check = '/'.join(genome_index.split('/')[:-1])
if not os.path.isdir(genome_index_check):
	print( "Error: genome reference index does not exist. If you don't have the index file ready, run 'hisat-3n-build --base-change C,T your_FASTA_file your_indexFolder_Path/genome' and 'hisat-3n-build --base-change C,T --repeat-index your_FASTA_file your_indexFolder_Path/genome'. You need to set the path to your index directory and prefix in the Snakefile." )
	sys.exit(1)
if not os.path.isfile(genome_bed):
    print( "Error: genome annotation gtf does not exist." )
    sys.exit(1)
if not os.path.isfile(chrom_size):
    print( "Error: genome chrom.size file doesn't exist." )
    sys.exit(1)
if not os.path.isfile(lookup_table):
    print( "Error: lookup table file doesn't exist." )
    sys.exit(1)

################################################
## Sample information Validation

## Id / Name column must be unique
if not samples.Id.is_unique:
	print( "Error: Id column in sample.tsv is not unique" )
	sys.exit(1)
if not samples.Name.is_unique:
	print( "Error: Name column in sample.tsv is not unique" )
	sys.exit(1)

## Group column must not overlap with Name column
#if not len(set(samples.Group).intersection(set(samples.Name)))==0:
#	print( "Error: Sample Group name must not overlap with Name" )
#	sys.exit(1)

## Only alphanumeric / dash / underbar / dot in sample sheet
## Must start with alphanumeric only
invalid_elem=[]
for col in samples:
    tmp = samples[col].str.count(r'(^[a-zA-Z0-9][a-zA-Z0-9-_\.]+$)')
    index_invalid = (tmp == 0)
    if index_invalid.any():
        invalid_elem = invalid_elem + samples[col][index_invalid].tolist()

if len(invalid_elem) > 0:
    print( "Error: Must be at least two character; Only alphanumeric, dash (-) and underbar (_) are allowed in a sample sheet")
    print( "Invalid values:" )
    for elem in invalid_elem: print( "  - %s" % elem )
    sys.exit(1)

softclipTypes = [True, False]
## check if statistical significance is valid
if soft_clipping not in softclipTypes:
    print("Error: the softClipping parameter in the 'config.bis.yml' file should ONLY be 'True' or 'False', without the quotation marks.")
    sys.exit(1)

sig_types = ["pVal", "FDR"]
## check if statistical significance is valid
if sig_type_call not in sig_types:
    print("Error: the statistical_significance_call variable in the 'config.bis.yml' file should ONLY be 'pVal' or 'FDR'.")
    sys.exit(1)
if sig_type_diff not in sig_types:
    print("Error: the statistical_significance_call variable in the 'config.bis.yml' file should ONLY be 'pVal' or 'FDR'.")
    sys.exit(1)

## check if numerical parameters are actual numbers
if isinstance(cov_thresh, str):
    print("Error: the coverage_threshold variable in the config.bis.yml file needs to be a number.")
    sys.exit(1)
if isinstance(mr_thresh, str):
    print("Error: the methylationRate_threshold variable in the config.bis.yml file needs to be a number.")
    sys.exit(1)
if isinstance(diff_thresh, str):
    print("Error: the differential_threshold variable in the config.bis.yml file needs to be a number.")
    sys.exit(1)
if isinstance(sig_thresh_call, str):
    print("Error: the statistical_threshold_call variable in the config.bis.yml file needs to be a number.")
    sys.exit(1)
if isinstance(sig_thresh_diff, str):
    print("Error: the statistical_threshold_diff variable in the config.bis.yml file needs to be a number.")
    sys.exit(1)