#
source('/athena/masonlab/scratch/users/nai2008/ivanov_functions.R')

### get viral fasta ftp paths

# read in assembly_summary for viral genomes
tt=read.delim('https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt',skip=1, na.strings='')
colnames(tt)[1]='assembly_accession'

tt$organism_name=gsub(" ","_",tt$organism_name)

# drop genomes which were discarded from RefSeq (see ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt, Column 23, for more details)
tt=tt[which(is.na(tt$excluded_from_refseq)),]

path_to_fna=paste0(tt$ftp_path, "/", ss(as.vector(tt$ftp_path),'/',10), '_genomic.fna.gz')

entry=paste0("(\'",tt$organism_name,"\', ","\'",tt$assembly_accession,"\'): ","\'",path_to_fna,"\',")

write.table(entry, file='/athena/masonlab/scratch/users/nai2008/COVID19/viral_fna_files.txt',
	quote = FALSE, row.names = FALSE, col.names = FALSE)

# NOTE: This list includes viruses with the the highest level of assembly of: Chromosomes, Scaffolds, or Contigs
# (see ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt, Column 13, for more details)

# format:
# ('Influenza_A', 'GCF_000865085.1'): 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Influenza_A_virus/reference/GCF_000865085.1_ViralMultiSegProj15622/GCF_000865085.1_ViralMultiSegProj15622_genomic.fna.gz',

### get fasta ftp paths for viruses that are important in human disease

subset=grep("human|Enterovirus_C|Coxsackievirus|rhinovirus|hepatitis_c|Hepatovirus_A|norovirus|yellow_fever|west_nile|dengue|equine_encephalitis|rubella|coronavirus|human_immunodeficiency|Hepatitis_delta_virus|influenza|measles|mumps|respiratory_syncytial_virus|parainfluenza|rabies|ebola|marburg|hantavirus|orthobunyavirus|encephalitis|La_Crosse|rift_valley|Lymphocytic_choriomeningitis|rotavirus|colorado|Human_alphaherpesvirus_1|human_alphaherpesvirus_2|Human_gammaherpesvirus_4|Human_betaherpesvirus_5|Human_alphaherpesvirus_3|betaherpesvirus|gammaherpesvirus|jc_polyomavirus|human_papilloma|Human_parvovirus_B19|human_adenovirus|variola|Cowpox_virus|molluscum|Hepatitis_B_virus\', \'GCF_000861825.2", entry, ignore.case=TRUE)

entry_subset=entry[subset]

write.table(entry_subset, file='/athena/masonlab/scratch/users/nai2008/COVID19/viral_fna_files_subset.txt',
	quote = FALSE, row.names = FALSE, col.names = FALSE)