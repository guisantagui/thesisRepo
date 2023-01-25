__author__ = 'chloeloiseau, modifeyd by MireiaCoscolla'
import os, sys, subprocess,shutil
from subprocess import call
import argparse
from datetime import date
import gzip
import bz2
import itertools
import re
import tempfile
import pysam
import random
import string
from Bio import SeqIO
from collections import defaultdict
import logging
import Modules_for_PipelineTB
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import shutil
import datetime

# Arguments to the CLI
parser = argparse.ArgumentParser(
    description='Whole-genome sequencing pipeline. The pipeline performs trimming, mapping and variant calling. Run: python Pipeline_TB.py [args]')
parser.add_argument('-f', dest='file',
                    help='text file of format: G-number\tpath to the fastq file. Format accepted for the fastq file: *.fastq, *.fastq.gz, *.fastq.bz2. ',
                    required=True)
parser.add_argument('-o', dest='output', help='output directory where output files created by the pipeline will be stored. Should end with /')
parser.add_argument('--read-length',dest='minlen', help='FastQValidator argument: fastQ is considered invalid if a read if < <read-length>. Default value: 20', default=20,type=int)
parser.add_argument('--min-length',dest='MinLength',help='Trimmomatic & SeqPrep argument: drop the read if it is below <min-length>. Post-trimming minimum read length. Default value: 20' , default = 20)
parser.add_argument('-c',dest='contaminant',help='consider contamination if X percent is mapped. Default value: 80.0', type=float,default=80.0)
parser.add_argument('-q',dest='min_mapping_quality', help='Samtools mpileup:Minimum mapping quality. Default value: 20', type=int, default= 20)
parser.add_argument('--min-coverage', dest='min_coverage', help='VARSCAN:min read depth at a position to make  call. Default value: 7', type=int,default=7)
parser.add_argument('--min-reads2', dest='min_reads', help='VARSCAN:min supporting reads at a position to call a variant. Default value: 5', type=int,default=5)
parser.add_argument('--min-avg-qual', dest='avg_qual', help='VARSCAN:min base quality at a position to count a read. Default value: 20', type=int,default=20)
parser.add_argument('--min-var-freq', dest='var_freq', help='VARSCAN:min variant allele frequency threshold. Default value: 0.1', type=float,default=0.1)
parser.add_argument('--min-var-freq-for-hom',dest='var_freq_homo', help='VARSCAN:Minimum frequency to call homozygote. Default value: 0.90', type=float, default= 0.90)
parser.add_argument('--excluded-loci',dest='exclude_positions',help='path to file: Locus_to_exclude_Mtb.txt Default: MTB-pipeline/Pipeline_TB/Locus_to_exclude_Mtb.txt',default='MTB-pipeline/Pipeline_TB/Locus_to_exclude_Mtb.txt')
parser.add_argument('--git',dest='git',help='Specify the path to the MTB-pipeline git repository cloned in your home. Ex: /scicore/home/gagneux/loiseau/. Shoult end with "/".',required=True)
parser.add_argument('-p', dest='threads', help='number of threads to launch trimmomatic and BWA', type = int, default = '4')
parser.add_argument('--overlap-size',dest='overlap',help='SeqPrep argument: minimum overall base pair overlap to merge two reads. default = 15',default=15)
parser.add_argument('--userlog',dest='user',help='append the Pipeline_mapping_log.tsv file',action='store_true')

args = parser.parse_args()

#REF = args.git+'MTB-pipeline/Pipeline_TB/MTB_ref_fasta/MTB_ancestor_reference.fasta'

#REF is B melitensis
#REF = '/scicore/home/gagneux/coscolla/aps/scripts/brucella_pipeline/references/BM.chr1chr2.fasta'
#REF is B melitensis rev1
#REF = '/scicore/home/gagneux/coscolla/aps/scripts/brucella_pipeline/references/rev1_assembly/GCF_000158695.1_ASM15869v1_genomic.fna'

#REF is each chromosomome of B. melitensis Bv1 strainM16
#REF = '/scicore/home/gagneux/coscolla/aps/scripts/brucella_pipeline/references/NC_003317.1.fasta'
#REF = '/scicore/home/gagneux/coscolla/aps/scripts/brucella_pipeline/references/NC_003318.1.fasta'

#REF is MTB_ancestor in garnatxa
REF = '/storage/PGO/store/db_pipe/reference/MTB_ancestor_reference.fasta'



'''create logger: what gets printed to the standard output'''
log = logging.getLogger(__name__)
out_hdlr = logging.StreamHandler(sys.stdout)
out_hdlr.setFormatter(logging.Formatter('%(levelname)s:%(message)s'))
out_hdlr.setLevel(logging.DEBUG)
log.addHandler(out_hdlr)
log.setLevel(logging.DEBUG)

'''Create a dictionary for storing the log'''
LOG_LIST = list()
LOG_DIC = defaultdict(list)

'''Open input_file for reading by the pipeline and throw an exception in case of error (i.e if there is no tab between G number and path to fastq files)'''
fastq_Gnumber_list = []
try:
    with open(args.file, 'r') as input_file:  # open the file with G numbers and paths to fastq files for reading
        fastq_Gnumber_list = [re.split(';|\t',genome.strip()) for genome in input_file if genome != '\n'] # split the lines of the input files by tab and semi colon

        if fastq_Gnumber_list[-1] == ['']: # If the last line of the input file is empty by mistake,  remove it from the list.
            fastq_Gnumber_list = fastq_Gnumber_list[:-1]

except IOError as (errno, strerror):
    log.error("I/O error({0}): {1}".format(errno, strerror))

bad_sample_file = open('/storage/PGO/PGOgit/prueba/genomes.problematic', 'a')  # where should this be stored ? cwd ? chan: se guarda donde ejecutas la pipe

'''Create the files to update LabKey'''
G_number = str
principal_directory = str

'''Pipeline start'''
for fastq in fastq_Gnumber_list:  # fastq[0] is the G number and fastq[1] is the path to the fastq
    print fastq
    G_number = fastq[0]
    principal_directory = args.output+G_number[0:3]+'/'+G_number[3:5]+'/'+G_number[5]+'/' # call this directory the principal_directory

    log.info('Pipeline starting for %s' % G_number)

    LOG_LIST.append(('DAT_OUT',str(date.today())))
    LOG_LIST.append(('GID_IN',G_number))

    str1 = os.path.basename(fastq[1])
    pattern_R1R2 = r'.*_[rR][1-2].*\.fastq'
    pattern_1_2 = r'.*_[1-2]\.fastq'
    match = re.search(pattern_R1R2, str1)
    match2 = re.search(pattern_1_2, str1)
    print match
    print match2

    Forward_pattern_R1 = r'.*_[rR]1.*\.fastq'
    Forward_pattern_1 = r'.*_1\.fastq'
    Reverse_pattern_R1 = r'.*_[rR]2.*\.fastq'
    Reverse_pattern_2 = r'.*_2\.fastq'

    if len(G_number) != 6 or G_number[0] != 'G': # if the G number is not a valid G number

        bad_sample_file.write(G_number+'\tG_number is not a 5 digit number or separator between G_number and FastQ file is not a tab\n')
        log.fatal('G_number is not a 5 digit number or separator between G_number and FastQ file is not a tab')
        continue # go to the next genome on the list

    if os.path.exists(principal_directory): # if the folder in common/mappings already exist
        shutil.rmtree(principal_directory) # Delete the directory
        log.warn('the folder in common_mapping folder already exists for this genome. Deleting it.')

    os.makedirs(principal_directory) # make directory for the G number
    subprocess.call('chmod -R g+w '+ principal_directory,shell=True)
    if len(fastq) == 2 or (len(fastq) == 3 and 'SE' in fastq):  # one run
        if (match2 == None and match == None) or (match == None and match2 != None and 'SE' in fastq) or (match != None and match2 == None and 'SE' in fastq):  # SE
            log.info('Pipeline running mode: single-end, single run')

            path_to_fastq = fastq[1] # second column of the input file
            if not os.path.exists(path_to_fastq):
                bad_sample_file.write(G_number+'\tFastQ file does not exist.\n')
                log.fatal('FastQ file does not exists: "%s"'%path_to_fastq)
                continue

            LOG_LIST.append(('FQ_IN',path_to_fastq))

            basename = Modules_for_PipelineTB.Extract_basename(path=path_to_fastq)[0] # basename
            basename_fastq = Modules_for_PipelineTB.Extract_basename(path=path_to_fastq)[1] # basename without extension

            if basename.endswith('fastq.bz2'):
                subprocess.call(['bzip2', '-dfkq', path_to_fastq]) # Options here are: decompress -overwrite output, quiet and keep input
                path_to_fastq = path_to_fastq[:-4] # new path to the fastq is the old path with without the bz2 extension. This works because it is in the same directory as the bz2. but ideally should be in a temp directory

            '''Check for PHRED encoding. If phred 64, convert to phred33'''
            tf_fastq_phred33 = tempfile.NamedTemporaryFile()
            try: # try opening a normal fastq file if it does not work use the gzip module
                raw_fastq = gzip.open(path_to_fastq,'r')
                first_line =raw_fastq.readline().split(':')
                for line in itertools.islice(raw_fastq,3,None,4): # lopp through the lines starting from the 3rd line and taking a step of 4 (look at the quality string)
                    phred33 = re.search(r'[\d\+\*\(\)\."#\$%&-,:\']',line) # look for the specific phred33 encoding characters in the lines
                    phred64 = re.search(r'[K-Z\[\]^_`a-h]',line)
                    if phred33 != None: #  if the encoding of the sample is
                        LOG_LIST.append(('PHE_OUT','33'))
                        break
                    elif phred64 != None:
                        SeqIO.convert(path_to_fastq,"fastq-illumina", tf_fastq_phred33.name+'.fastq', "fastq-sanger")
                        path_to_fastq = tf_fastq_phred33.name+'.fastq'
                        LOG_LIST.append(('PHE_OUT','64'))
                        break
            except:
                raw_fastq = open(path_to_fastq,'r')
                first_line =raw_fastq.readline().split(':')
                for line in itertools.islice(raw_fastq,3,None,4): # lopp through the first 1000 lines starting from the 3rd line and taking a step of 4
                    phred33 = re.search(r'[\d\+\*\(\)\."#\$%&-,:\']',line) # look for the specific phred33 encoding characters in the lines
                    phred64 = re.search(r'[K-Z\[\]^_`a-h]',line)
                    if phred33 != None: #  if the encoding of the sample is phred33
                        LOG_LIST.append(('PHE_OUT','33'))
                        break

                    elif phred64 != None:
                        SeqIO.convert(path_to_fastq,"fastq-illumina", tf_fastq_phred33.name+'.fastq', "fastq-sanger")
                        path_to_fastq = tf_fastq_phred33.name+'.fastq'
                        LOG_LIST.append(('PHE_OUT','64'))
                        break

            Modules_for_PipelineTB.Run_FastQValidator(fastq=path_to_fastq,minimum_read_length=args.minlen,log=LOG_LIST)
            log.info('FastQValidator: %s'% Modules_for_PipelineTB.stdout_FV)
            if Modules_for_PipelineTB.stdout_FV[Modules_for_PipelineTB.stdout_FV.rfind(':') + 2:-1] == 'FASTQ_SUCCESS':  # if the fastq passes the fastQValidator test. What if only one fastq passes the test?
                LOG_LIST.append(('FV_OUT','FASTQ_SUCCESS'))

                log.info('Starting FastQC on %s'%path_to_fastq)
                Modules_for_PipelineTB.Run_FastQC(fastq=path_to_fastq,output_directory=principal_directory,log=LOG_LIST,gnumber=G_number,basename=basename_fastq, number_of_threads=args.threads)

                if type(Modules_for_PipelineTB.read_length) is str:
                    LOG_LIST.append(('LNG_OUT','NA'))

                elif type(Modules_for_PipelineTB.read_length) is int:
                    LNG_coverage = (float(Modules_for_PipelineTB.total_reads)*Modules_for_PipelineTB.read_length*2)/4411532
                    LOG_LIST.append(('LNG_OUT',str(LNG_coverage)))

                log.info('Starting Trimmomatic')
                Modules_for_PipelineTB.Run_Trimmomatic_SE(number_of_threads=args.threads,fastq=path_to_fastq,minimum_read_length_post_trimming=args.MinLength,log=LOG_LIST)

                '''Parse the fastq header'''
                first_line = [] # where the fastq header will be
                try:  # try opening a normal fastq file if it does not work use the gzip module
                    raw_fastq = gzip.open(path_to_fastq, 'r')
                    first_line = raw_fastq.readline().split(':')
                except:
                    raw_fastq = open(path_to_fastq, 'r')
                    first_line = raw_fastq.readline().split(':')

                flowcell_ID = ''
                flowcell_lane = ''
                if len(first_line) == 10: # if the fastq header is composed of 10 elements
                    flowcell_ID = first_line[2] # then select the second element which is the flowcell ID
                    flowcell_lane = first_line[3] # and select the third element which is the flowcell lane
                elif len(first_line) != 10: # if the fastq header is not composed of 10 elements
                    flowcell_ID = G_number
                    flowcell_lane= 'A'

                #print 'basename is ', basename
                if basename.endswith('fastq.bz2'):
                    print path_to_fastq
                    os.remove(path_to_fastq)
                else:
                    print 'fastq not bz2'
                    print path_to_fastq

                ''' Mapping processed reads to the Mtb Anc genome'''
                log.info('Starting BWA')
                tf_bam = tempfile.NamedTemporaryFile(prefix='mapped_sorted_')
                Modules_for_PipelineTB.BWA_mapping_SE(inreads=Modules_for_PipelineTB.fastq_for_BWA, tempbamfile=tf_bam.name,flowcell=flowcell_ID,lane=flowcell_lane,Gnumber=G_number,number_of_threads=args.threads,log=LOG_LIST,reference_genome=REF)

                log.info('Starting MarkDuplicates')
                Modules_for_PipelineTB.Run_MarkDuplicates(input=tf_bam.name+'.bam', principaldirectory=principal_directory,log=LOG_LIST)

                log.info('Indexing the intermediate BAM file')
                Modules_for_PipelineTB.Index_BAM(input=principal_directory+'*.bam',log=LOG_LIST,name_of_index='IDX1_IN')

                log.info('Starting GATK Realigner Target Creator')
                Modules_for_PipelineTB.GATK_RealignerTargetCreator(log=LOG_LIST,input=principal_directory+'*.bam',number_threads=args.threads, reference_genome=REF)

                tf_bam_postGATK = tempfile.NamedTemporaryFile(prefix='postGATK_bam_')
                
                log.info('Starting GATK IndelRealigner')
                Modules_for_PipelineTB.GATK_IndelRealigner(log=LOG_LIST,input=principal_directory+'*.bam',output=tf_bam_postGATK.name+'.bam',reference_genome=REF,interval_file=Modules_for_PipelineTB.GATK_intervals)

                log.info('Indexing the temp BAM file')
                Modules_for_PipelineTB.Index_BAM(input=tf_bam_postGATK.name+'.bam',log=LOG_LIST,name_of_index='IDX2_IN')

                Modules_for_PipelineTB.Samtools_Flagstat(input_bam=tf_bam_postGATK.name+'.bam',log=LOG_LIST,contaminant_threshold=args.contaminant,gnb=G_number)

                log.info('Filtering BAM with PySam')
                Modules_for_PipelineTB.FilterBAM_AS(infile=tf_bam_postGATK.name+'.bam',gnumber=G_number,principaldirectory=principal_directory)

                log.info('Indexing the final BAM file')
                Modules_for_PipelineTB.Index_BAM(input=principal_directory+G_number+'.bam',log=LOG_LIST,name_of_index='IDX3_IN')

                subprocess.call('rm '+principal_directory+'random_bam_*',shell=True)

                log.info('Starting Qualimap BamQC')
                Modules_for_PipelineTB.Qualimap_BamQC(input_bam=principal_directory+G_number+'.bam',log=LOG_LIST,gnumber=G_number,number_threads=args.threads,output_dir=principal_directory)

                log.info('Starting VarScan')
                Modules_for_PipelineTB.Variant_calling(git=args.git,exclude_positions=args.exclude_positions,min_mapping_quality=args.min_mapping_quality,max_cov=Modules_for_PipelineTB.mean_Coverage_without_duplicates*3,min_cov=args.min_coverage,min_base_qual=args.avg_qual,min_reads2=args.min_reads,hetero_freq=args.var_freq,homo_freq=args.var_freq_homo,reference=REF,Gnumber=G_number,output_directory=principal_directory,log=LOG_LIST)
                Modules_for_PipelineTB.Find_lineage(vcf_file=principal_directory+G_number+'.var.snp.vcf',git=args.git,log=LOG_LIST)

                log.info('Producing whole-genome FASTA')
                Modules_for_PipelineTB.WGSFasta(excluded=args.git+args.exclude_positions,vcf=principal_directory+G_number+'.all.pos.vcf',dir=principal_directory,Gnumber=G_number)

                log.info('Compressing all.pos.vcf')
                with open(principal_directory+G_number+'.all.pos.vcf', 'rb') as f_in, gzip.open(principal_directory+G_number+'.all.pos.vcf.gz', 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                os.remove(principal_directory+G_number+'.all.pos.vcf')

                if args.user:
                    Modules_for_PipelineTB.Mapping_User_log(gnumber=G_number)

                log.info('Pipeline End')

            else:
                log.fatal('Fail at the fastQValidator step')
                LOG_LIST.append(('FV_OUT',Modules_for_PipelineTB.stdout_FV[Modules_for_PipelineTB.stdout_FV.rfind(':') + 2:-1]))
                bad_sample_file.write(G_number + '\t' + fastq[1] +'\tFail at FastQValidator Step\n') # fastq[1] and not path_to_fastq because want to know the original path that causes the problem
                continue  # continue with the next iteration of the loop. ie: the fastq did not pass this step so no point continuing with the other steps of the pipeline

        elif match2 != None or match != None and len(fastq) == 2: # PE
            log.info('Pipeline running mode: paired-end, single run')

            '''Detect if the forward read is given or the reverse and infer the complementary read.'''
            forward_read = ''
            reverse_read = ''

            forward_R1_match = re.search(Forward_pattern_R1, fastq[1])
            forward_1_match = re.search(Forward_pattern_1, fastq[1])
            reverse_R2_match = re.search(Reverse_pattern_R1, fastq[1])
            reverse_2_match = re.search(Reverse_pattern_2, fastq[1])

            if forward_R1_match != None:
                forward_read = fastq[1]
                reverse_read = fastq[1][0:fastq[1].rfind('R1')] + 'R2' + fastq[1][fastq[1].rfind('R1') + 2:]
            elif forward_1_match != None:
                forward_read = fastq[1]
                reverse_read = fastq[1][0:fastq[1].rfind('_1.')] + '_2.'+ fastq[1][fastq[1].rfind('_1.')+3:]
            elif reverse_R2_match != None:
                reverse_read = fastq[1]
                forward_read = fastq[1][0:fastq[1].rfind('R2')] + 'R1' + fastq[1][fastq[1].rfind('R2') + 2:]
            elif reverse_2_match != None:
                reverse_read = fastq[1]
                forward_read = fastq[1][0:fastq[1].rfind('_2.')] + '_1.'+fastq[1][fastq[1].rfind('_2.')+3:]

            LOG_LIST.append(('FQ_IN',forward_read))
            LOG_LIST.append(('FQ_IN',reverse_read))

            basename_forward = Modules_for_PipelineTB.Extract_basename(path=forward_read)[0]
            basename_fastq_forward = Modules_for_PipelineTB.Extract_basename(path=forward_read)[1]

            basename_reverse = Modules_for_PipelineTB.Extract_basename(path=reverse_read)[0]
            basename_fastq_reverse = Modules_for_PipelineTB.Extract_basename(path=reverse_read)[1]

            '''Detect the extension '''
            if basename_forward.endswith('fastq.bz2'):
                subprocess.call(['bzip2', '-dfkq',forward_read])  #  Options here are: decompress -overwrite output, quiet and keep input
                forward_read = forward_read[:-4]  # new path to the fastq is the old path with without the bz2 extension. This works because it is in the same directory as the bz2. but ideally should be in a temp directory
            elif basename_reverse.endswith('fastq.bz2'):
                subprocess.call(['bzip2', '-dfkq',reverse_read])  # Options here are: decompress -overwrite output, quiet and keep input
                reverse_read = reverse_read[:-4]


            '''Check for PHRED encoding. If phred 64, convert to phred33'''
            tf_forward_read_phred33 = tempfile.NamedTemporaryFile()
            tf_reverse_read_phred33 = tempfile.NamedTemporaryFile()

            try: # try opening a normal fasq file if it does not work use the gzip module
                raw_fastq = gzip.open(forward_read,'r')
                first_line =raw_fastq.readline().split(':')
                for line in itertools.islice(raw_fastq,3,None,4): # lopp through the first 1000 lines starting from the 3rd line and taking a step of 4
                    phred33 = re.search(r'[\d\+\*\(\)\."#\$%&-,:\!\']',line) # look for the specific phred33 encoding characters in the lines
                    phred64 = re.search(r'[K-Z\[\]^_`a-h]',line)
                    if phred33 != None: #  if the encoding of the sample is
                        LOG_LIST.append(('PHE_OUT','33'))
                        break

                    elif phred64 != None:
                        SeqIO.convert(forward_read,"fastq-illumina", tf_forward_read_phred33.name+'.fastq', "fastq-sanger")
                        SeqIO.convert(reverse_read,"fastq-illumina", tf_reverse_read_phred33.name+'.fastq', "fastq-sanger")
                        forward_read = tf_forward_read_phred33.name+'.fastq'
                        reverse_read = tf_reverse_read_phred33.name+'.fastq'
                        LOG_LIST.append(('PHE_OUT','64'))
                        break
            except:
                raw_fastq = open(forward_read,'r')
                first_line =raw_fastq.readline().split(':')
                for line in itertools.islice(raw_fastq,3,None,4): # lopp through the first 1000 lines starting from the 3rd line and taking a step of 4
                    phred33 = re.search(r'[\d\+\*\(\)\."#\$%&-,:\!\']',line) # look for the specific phred33 encoding characters in the lines
                    phred64 = re.search(r'[K-Z\[\]^_`a-h]',line)
                    if phred33 != None: #  if the encoding of the sample is phred33
                        LOG_LIST.append(('PHE_OUT','33'))
                        break

                    elif phred64 != None:
                        SeqIO.convert(forward_read,"fastq-illumina", tf_forward_read_phred33.name+'.fastq', "fastq-sanger")
                        SeqIO.convert(reverse_read,"fastq-illumina", tf_reverse_read_phred33.name+'.fastq', "fastq-sanger")

                        forward_read = tf_forward_read_phred33.name+'.fastq'
                        reverse_read = tf_reverse_read_phred33.name+'.fastq'
                        LOG_LIST.append(('PHE_OUT','64'))
                        break
            #print 'stdout is;',Modules_for_PipelineTB.stdout_FV
            Modules_for_PipelineTB.Run_FastQValidator(fastq=forward_read,minimum_read_length=args.minlen,log=LOG_LIST)
            log.info('FastQValidator forward read: %s'% Modules_for_PipelineTB.stdout_FV)

            if Modules_for_PipelineTB.stdout_FV[Modules_for_PipelineTB.stdout_FV.rfind(':') + 2:-1] == 'FASTQ_SUCCESS':  # if the fastq passes the fastQValidator test. What if only one fastq passes the test?
                LOG_LIST.append(('FV_OUT','FASTQ_SUCCESS'))
                Modules_for_PipelineTB.Run_FastQValidator(fastq=reverse_read,minimum_read_length=args.minlen,log=LOG_LIST)
                log.info('FastQValidator reverse read: %s'% Modules_for_PipelineTB.stdout_FV)

                if Modules_for_PipelineTB.stdout_FV[Modules_for_PipelineTB.stdout_FV.rfind(':') + 2:-1] == 'FASTQ_SUCCESS':

                    LOG_LIST.append(('FV_OUT','FASTQ_SUCCESS'))
                    log.info('Starting FastQC on %s'%forward_read)

                    Modules_for_PipelineTB.Run_FastQC(fastq=forward_read,output_directory=principal_directory,log=LOG_LIST,gnumber=G_number,basename=basename_fastq_forward, number_of_threads=args.threads)

                    if type(Modules_for_PipelineTB.read_length) is str:
                        LOG_LIST.append(('LNG_OUT','NA'))

                    elif type(Modules_for_PipelineTB.read_length) is int:
                        LNG_coverage = (float(Modules_for_PipelineTB.total_reads)*Modules_for_PipelineTB.read_length*2)/4411532
                        LOG_LIST.append(('LNG_OUT',str(LNG_coverage)))

                    log.info('Starting Trimmomatic')
                    Modules_for_PipelineTB.Run_Trimmomatic_PE(number_of_threads=args.threads,minimum_read_length_post_trimming=args.MinLength,log=LOG_LIST,forwardRead=forward_read,reverseRead=reverse_read)

                    log.info('Starting SeqPrep')
                    Modules_for_PipelineTB.Run_SeqPrep(forward_paired = Modules_for_PipelineTB.Fastq1P, reverse_paired=Modules_for_PipelineTB.Fastq2P,overlap=args.overlap, minlen=args.MinLength,log=LOG_LIST)


                    '''Parse the fastq header'''
                    first_line = [] # where the fastq header will be
                    try:  # try opening a normal fastq file if it does not work use the gzip module
                        raw_fastq = gzip.open(forward_read, 'r')
                        first_line = raw_fastq.readline().split(':')
                    except:
                        raw_fastq = open(forward_read, 'r')
                        first_line = raw_fastq.readline().split(':')

                    flowcell_ID = ''
                    flowcell_lane = ''
                    if len(first_line) == 10:  # if the fastq header is composed of 10 elements
                        flowcell_ID = first_line[2]  # then select the second element which is the flowcell ID
                        flowcell_lane = first_line[3]  # and select the third element which is the flowcell lane
                    elif len(first_line) != 10:  # if the fastq header is not composed of 10 elements
                        flowcell_ID = G_number
                        flowcell_lane = 'A'

                    #if basename_forward.endswith('fastq.bz2'):
                    #    print "basename forward:", basename_forward
                    #    print forward_read
                        #os.remove(forward_read)
                    #elif basename_reverse.endswith('fastq.bz2'):
                    #    print "basename reverse:", basename_reverse
                    #    print reverse_read
                        #os.remove(reverse_read)


                    tf_bam_se_trimmomatic_1U = tempfile.NamedTemporaryFile(prefix='mapped_sorted_se_trimmomatic_1U')
                    tf_bam_se_trimmomatic_2U = tempfile.NamedTemporaryFile(prefix='mapped_sorted_se_trimmomatic_2U')
                    tf_bam_se_SeqPrep_merge = tempfile.NamedTemporaryFile(prefix='mapped_sorted_se_SeqPrep_merge')
                    tf_bam_pe_SeqPrep_trimmed = tempfile.NamedTemporaryFile(prefix = 'mapped_sorted_pe_SeqPrep_trimmed')
                    tf_bam_pe_SeqPrep_NOTmerged = tempfile.NamedTemporaryFile(prefix = 'mapped_sorted_pe_SeqPrep_NOTmerged')

                    log.info('Start BWA: mapping unpaired forward reads from Trimmomatic')
                    Modules_for_PipelineTB.BWA_mapping_SE(inreads=Modules_for_PipelineTB.Fastq1U,
                                   tempbamfile=tf_bam_se_trimmomatic_1U.name,
                                   flowcell=flowcell_ID,lane=flowcell_lane,Gnumber=G_number,number_of_threads=args.threads,log=LOG_LIST,reference_genome=REF)

                    log.info('Start BWA: mapping unpaired reverse reads from Trimmomatic')
                    Modules_for_PipelineTB.BWA_mapping_SE(inreads=Modules_for_PipelineTB.Fastq2U,
                                   tempbamfile=tf_bam_se_trimmomatic_2U.name,
                                   flowcell=flowcell_ID, lane=flowcell_lane,
                                   Gnumber=G_number,number_of_threads=args.threads,log=LOG_LIST,reference_genome=REF)

                    log.info('Start BWA: mapping merged reads from SeqPrep')
                    Modules_for_PipelineTB.BWA_mapping_SE(inreads=Modules_for_PipelineTB.SeqPrep_Merged_Reads,
                                   tempbamfile=tf_bam_se_SeqPrep_merge.name,
                                   flowcell=flowcell_ID, lane=flowcell_lane,
                                   Gnumber=G_number,number_of_threads=args.threads,log=LOG_LIST,reference_genome=REF)

                    log.info('Start BWA: mapping trimmed reads from SeqPrep')
                    Modules_for_PipelineTB.BWA_mapping_PE(fw_reads=Modules_for_PipelineTB.SeqPrep_Trimmed_Reads_Forward,
                                   rev_reads=Modules_for_PipelineTB.SeqPrep_Trimmed_Reads_Reverse,
                                   tempbamfile=tf_bam_pe_SeqPrep_trimmed.name,flowcell=flowcell_ID,lane=flowcell_lane,Gnumber=G_number,number_of_threads=args.threads,log=LOG_LIST,reference_genome=REF)

                    log.info('Start BWA: mapping non merged reads from Seqprep')
                    Modules_for_PipelineTB.BWA_mapping_PE(fw_reads=Modules_for_PipelineTB.SeqPrep_Non_Merged_Reads_Forward,
                                   rev_reads=Modules_for_PipelineTB.SeqPrep_Non_Merged_Reads_Reverse,
                                   tempbamfile=tf_bam_pe_SeqPrep_NOTmerged.name,flowcell=flowcell_ID,lane=flowcell_lane,Gnumber=G_number,number_of_threads=args.threads,log=LOG_LIST,reference_genome=REF)

                    tf_bam = tempfile.NamedTemporaryFile(prefix = 'mapped_sorted_merged')

                    log.info('Merging all Bam files produced')
                    samtools_merge_command = 'samtools merge {output_bam} {inbam_trimmomatic_1U} {inbam_trimmomatic_2U} {inbam_SeqPrep_merged} {inbam_SeqPrep_trimmed} {inbam_SeqPrep_NOTmerged}'.format(
                    output_bam=tf_bam.name + '.bam',
                    inbam_trimmomatic_1U=tf_bam_se_trimmomatic_1U.name+'.bam',
                    inbam_trimmomatic_2U=tf_bam_se_trimmomatic_2U.name+'.bam',
                    inbam_SeqPrep_merged=tf_bam_se_SeqPrep_merge.name+'.bam',
                    inbam_SeqPrep_trimmed=tf_bam_pe_SeqPrep_trimmed.name+'.bam',
                    inbam_SeqPrep_NOTmerged=tf_bam_pe_SeqPrep_NOTmerged.name+'.bam')

                    subprocess.call(samtools_merge_command,shell=True)
                    LOG_LIST.append(('SM1_IN',samtools_merge_command))

                    log.info('Starting MarkDuplicates')
                    Modules_for_PipelineTB.Run_MarkDuplicates(input=tf_bam.name+'.bam', principaldirectory=principal_directory,log=LOG_LIST)

                    log.info('Indexing the intermediate BAM file')
                    Modules_for_PipelineTB.Index_BAM(input=principal_directory+'*.bam',log=LOG_LIST,name_of_index='IDX1_IN')

                    log.info('Starting GATK Realigner Target Creator')
                    Modules_for_PipelineTB.GATK_RealignerTargetCreator(log=LOG_LIST,input=principal_directory+'*.bam',number_threads=args.threads, reference_genome=REF)

                    tf_bam_postGATK = tempfile.NamedTemporaryFile(prefix='postGATK_bam_')
                    log.info('Starting GATK IndelRealigner')
                    Modules_for_PipelineTB.GATK_IndelRealigner(log=LOG_LIST,input=principal_directory+'*.bam',output=tf_bam_postGATK.name+'.bam',reference_genome=REF,interval_file=Modules_for_PipelineTB.GATK_intervals)

                    log.info('Indexing the temp BAM file')
                    Modules_for_PipelineTB.Index_BAM(input=tf_bam_postGATK.name+'.bam',log=LOG_LIST,name_of_index='IDX2_IN')

                    Modules_for_PipelineTB.Samtools_Flagstat(input_bam=tf_bam_postGATK.name+'.bam',log=LOG_LIST,contaminant_threshold=args.contaminant,gnb=G_number)

                    log.info('Filtering BAM with PySam')
                    Modules_for_PipelineTB.FilterBAM_AS(infile=tf_bam_postGATK.name+'.bam',gnumber=G_number,principaldirectory=principal_directory)

                    log.info('Indexing the final BAM file')
                    Modules_for_PipelineTB.Index_BAM(input=principal_directory+G_number+'.bam',log=LOG_LIST,name_of_index='IDX3_IN')

                    subprocess.call('rm '+principal_directory+'random_bam_*',shell=True)

                    log.info('Starting Qualimap BamQC')
                    Modules_for_PipelineTB.Qualimap_BamQC(input_bam=principal_directory+G_number+'.bam',log=LOG_LIST,gnumber=G_number,number_threads=args.threads,output_dir=principal_directory)

                    Modules_for_PipelineTB.Variant_calling(git=args.git,exclude_positions=args.exclude_positions,
                                                           min_mapping_quality=args.min_mapping_quality,
                                                           max_cov=Modules_for_PipelineTB.mean_Coverage_without_duplicates*3,
                                                           min_cov=args.min_coverage,min_base_qual=args.avg_qual,
                                                           min_reads2=args.min_reads,hetero_freq=args.var_freq,
                                                           homo_freq=args.var_freq_homo,reference=REF,Gnumber=G_number,
                                                           output_directory=principal_directory,log=LOG_LIST)

                    Modules_for_PipelineTB.Find_lineage(vcf_file=principal_directory+G_number+'.var.snp.vcf',git=args.git,log=LOG_LIST)

                    log.info('Producing whole-genome FASTA')
                    Modules_for_PipelineTB.WGSFasta(excluded=args.git+args.exclude_positions,vcf=principal_directory+G_number+'.all.pos.vcf',dir=principal_directory,Gnumber=G_number)

                    log.info('Compressing all.pos.vcf')
                    with open(principal_directory+G_number+'.all.pos.vcf', 'rb') as f_in, gzip.open(principal_directory+G_number+'.all.pos.vcf.gz', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    os.remove(principal_directory+G_number+'.all.pos.vcf')

                    if args.user:
                        Modules_for_PipelineTB.Mapping_User_log(gnumber=G_number)

                    log.info('Pipeline End')

                else:
                    log.fatal('Fail at the fastQValidator step')
                    LOG_LIST.append(('FV_OUT',Modules_for_PipelineTB.stdout_FV[Modules_for_PipelineTB.stdout_FV.rfind(':') + 2:-1]))
                    bad_sample_file.write(G_number + '\t' + fastq[1] +'\tFail at FastQValidator Step\n') # fastq[1] and not path_to_fastq because want to know the original path that causes the problem
                    continue  # continue with the next iteration of the loop. ie: the fastq did not pass this step so no point continuing with the other steps of the pipeline
            else:
                log.fatal('Fail at the fastQValidator step')
                LOG_LIST.append(('FV_OUT',Modules_for_PipelineTB.stdout_FV[Modules_for_PipelineTB.stdout_FV.rfind(':') + 2:-1]))
                bad_sample_file.write(G_number + '\t' + fastq[1] +'\tFail at FastQValidator Step\n') # fastq[1] and not path_to_fastq because want to know the original path that causes the problem
                continue  # continue with the next iteration of the loop. ie: the fastq did not pass this step so no point continuing with the other steps of the pipeline

    elif (len(fastq) > 2 and 'SE' not in fastq) or (len(fastq) > 3 and 'SE' in fastq):  # multiple fastq
        if (match2 == None and match == None) or (match2 != None and 'SE' in fastq) or (match != None and 'SE' in fastq): # SE
            log.info('Pipeline running mode: single-end, multiple runs')
            List_multiple_fastq= []

            if fastq[-1] == 'SE':
                fastq=fastq[:-1]

            for i in range(1,len(fastq)):

                path_to_fastq = fastq[i]

                LOG_LIST.append(('FQ_IN',path_to_fastq))

                basename = Modules_for_PipelineTB.Extract_basename(path=path_to_fastq)[0]
                basename_fastq = Modules_for_PipelineTB.Extract_basename(path=path_to_fastq)[1]

                if basename.endswith('fastq.bz2'):
                    subprocess.call(['bzip2', '-dfkq', path_to_fastq]) # should add the output in a tem directory: ASK Mireia. Options here are: decompress -overwrite output, quiet and keep input
                    path_to_fastq = path_to_fastq[:-4] # new path to the fastq is the old path with without the bz2 extension. This works because it is in the same directory as the bz2. but ideally should be in a temp directory

                List_multiple_fastq += [path_to_fastq]


                tf_fastq_phred33 = tempfile.NamedTemporaryFile()

                try: # try opening a normal fasq file if it does not work use the gzip module
                    raw_fastq = gzip.open(path_to_fastq,'r')
                    first_line =raw_fastq.readline().split(':')
                    for line in itertools.islice(raw_fastq,3,None,4): # lopp through the first 1000 lines starting from the 3rd line and taking a step of 4
                        phred33 = re.search(r'[\d\+\*\(\)\."#\$%&-,:\!\']',line) # look for the specific phred33 encoding characters in the lines
                        phred64 = re.search(r'[K-Z\[\]^_`a-h]',line)
                        if phred33 != None: #  if the encoding of the sample is
                            LOG_LIST.append(('PHE_OUT','33'))
                            break
                        elif phred64 != None:
                            SeqIO.convert(path_to_fastq,"fastq-illumina", tf_fastq_phred33.name+'.fastq', "fastq-sanger")
                            path_to_fastq = tf_fastq_phred33.name+'.fastq'
                            LOG_LIST.append(('PHE_OUT','64'))

                            break
                except:
                    raw_fastq = open(path_to_fastq,'r')
                    first_line =raw_fastq.readline().split(':')
                    for line in itertools.islice(raw_fastq,3,None,4): # lopp through the first 1000 lines starting from the 3rd line and taking a step of 4
                        phred33 = re.search(r'[\d\+\*\(\)\."#\$%&-,:\!\']',line) # look for the specific phred33 encoding characters in the lines
                        phred64 = re.search(r'[K-Z\[\]^_`a-h]',line)
                        if phred33 != None: #  if the encoding of the sample is phred33
                            LOG_LIST.append(('PHE_OUT','33'))
                            break

                        elif phred64 != None:
                            SeqIO.convert(path_to_fastq,"fastq-illumina", tf_fastq_phred33.name+'.fastq', "fastq-sanger")
                            path_to_fastq = tf_fastq_phred33.name+'.fastq'
                            LOG_LIST.append(('PHE_OUT','64'))
                            break

                Modules_for_PipelineTB.Run_FastQValidator(fastq=path_to_fastq,minimum_read_length=args.minlen,log=LOG_LIST)
                log.info('FastQValidator: %s'% Modules_for_PipelineTB.stdout_FV)

                if Modules_for_PipelineTB.stdout_FV[Modules_for_PipelineTB.stdout_FV.rfind(':') + 2:-1] == 'FASTQ_SUCCESS':  # if the fastq passes the fastQValidator test. What if only one fastq passes the test?
                    LOG_LIST.append(('FV_OUT','FASTQ_SUCCESS'))

                    if i < len(fastq) -1: # while there are still more fastq (of the same gnumber) to validate continue
                        continue
                         # continue to see if the other fastq of the multiple run is ok
                    elif i == len(fastq) -1: # if this condition is met it means that fastqValidator has been run on every fastq so now we want to do fastQC to control the LNG and read length
                        sum_LNG = 0
                        count =1
                        for qc in List_multiple_fastq:
                            basename=Modules_for_PipelineTB.Extract_basename(qc)[0]
                            basename_fastq=Modules_for_PipelineTB.Extract_basename(qc)[1]
                            log.info('Starting FastQC on %s'%qc)
                            Modules_for_PipelineTB.Run_FastQC(fastq=qc,output_directory=principal_directory,log=LOG_LIST,gnumber=G_number,basename=basename_fastq, number_of_threads=args.threads)
                            if type(Modules_for_PipelineTB.read_length) is str:
                                sum_LNG = 'NA'
                                if count < len(List_multiple_fastq):  # process the other fastqs
                                    count += 1
                                    continue

                            elif type(Modules_for_PipelineTB.read_length) is int:
                                LNG_coverage = (float(Modules_for_PipelineTB.total_reads)*Modules_for_PipelineTB.read_length*2)/4411532
                                try:
                                    sum_LNG += LNG_coverage # ddo a try except kindof thing
                                except:
                                    sum_LNG = 'NA'
                                    log.info('Could not calculate LNG because read lenght are variable')

                                if count < len(List_multiple_fastq):  # process the other fastqs ## should'nt this be at the same level as the elif???
                                    count += 1
                                    continue

                            LOG_LIST.append(('LNG_OUT',str(sum_LNG)))

                            n = 0
                            for fq in List_multiple_fastq:
                                run_number = str(n+1)
                                log.info('Starting Trimmomatic for run %s' % run_number)
                                Modules_for_PipelineTB.Run_Trimmomatic_SE(number_of_threads=args.threads,fastq=fq,minimum_read_length_post_trimming=args.MinLength,log=LOG_LIST)

                                '''Parse the fastq header'''
                                first_line = [] # where the fastq header will be
                                try:  # try opening a normal fastq file if it does not work use the gzip module
                                    raw_fastq = gzip.open(fq, 'r')
                                    first_line = raw_fastq.readline().split(':')
                                except:
                                    raw_fastq = open(fq, 'r')
                                    first_line = raw_fastq.readline().split(':')

                                letters = list(string.ascii_uppercase) # letters = ['A','B','C','D',...,'Z']

                                flowcell_ID = ''
                                flowcell_lane = ''
                                if len(first_line) == 10:  # if the fastq header is composed of 10 elements
                                    flowcell_ID = first_line[2]  # then select the second element which is the flowcell ID
                                    flowcell_lane = first_line[3]  # and select the third element which is the flowcell lane
                                elif len(first_line) != 10:  # if the fastq header is not composed of 10 elements
                                    flowcell_ID = G_number
                                    flowcell_lane = letters[n] # n starts at 0 with the first run and increments  each time a new run is processed

                                ''' Mapping processed reads to the Mtb Anc genome'''

                                log.info('Starting BWA for run %s' % run_number)
                                tf_bam = tempfile.NamedTemporaryFile(prefix='mapped_sorted_')
                                Modules_for_PipelineTB.BWA_mapping_SE(inreads=Modules_for_PipelineTB.fastq_for_BWA,tempbamfile=tf_bam.name,flowcell=flowcell_ID,lane=flowcell_lane,Gnumber=G_number,number_of_threads=args.threads,log=LOG_LIST,reference_genome=REF)

                                log.info('Starting MarkDuplicates for run %s' % run_number)
                                Modules_for_PipelineTB.Run_MarkDuplicates(input=tf_bam.name+'.bam', principaldirectory=principal_directory,log=LOG_LIST)
                                n+=1
                            print List_multiple_fastq
                            log.info('Merging BAM files from all different runs')
                            tf_bam_merged = tempfile.NamedTemporaryFile(prefix='merged_bam_')
                            command_merge = 'samtools merge {outbam} {inbams}'.format(
                                            outbam=tf_bam_merged.name + '.bam',
                                            inbams=principal_directory + '*.bam')
                            subprocess.call(command_merge, shell=True)
                            LOG_LIST.append(('SM1_IN',command_merge))

                            log.info('Indexing the intermediate BAM file')
                            Modules_for_PipelineTB.Index_BAM(input=tf_bam_merged.name + '.bam',log=LOG_LIST,name_of_index='IDX1_IN')

                            log.info('Starting GATK Realigner Target Creator')
                            Modules_for_PipelineTB.GATK_RealignerTargetCreator(log=LOG_LIST,input=tf_bam_merged.name + '.bam',number_threads=args.threads, reference_genome=REF)

                            tf_bam_postGATK = tempfile.NamedTemporaryFile(prefix='postGATK_bam_')
                            log.info('Starting GATK IndelRealigner')
                            Modules_for_PipelineTB.GATK_IndelRealigner(log=LOG_LIST,input=tf_bam_merged.name + '.bam',output=tf_bam_postGATK.name+'.bam',reference_genome=REF,interval_file=Modules_for_PipelineTB.GATK_intervals)

                            log.info('Indexing the temp BAM file')
                            Modules_for_PipelineTB.Index_BAM(input=tf_bam_postGATK.name+'.bam',log=LOG_LIST,name_of_index='IDX2_IN')

                            Modules_for_PipelineTB.Samtools_Flagstat(input_bam=tf_bam_postGATK.name+'.bam',log=LOG_LIST,contaminant_threshold=args.contaminant,gnb=G_number)

                            log.info('Filtering BAM with PySam')
                            Modules_for_PipelineTB.FilterBAM_AS(infile=tf_bam_postGATK.name+'.bam',gnumber=G_number,principaldirectory=principal_directory)

                            log.info('Indexing the final BAM file')
                            Modules_for_PipelineTB.Index_BAM(input=principal_directory+G_number+'.bam',log=LOG_LIST,name_of_index='IDX3_IN')

                            log.info('Starting Qualimap BamQC')
                            Modules_for_PipelineTB.Qualimap_BamQC(input_bam=principal_directory+G_number+'.bam',log=LOG_LIST,gnumber=G_number,number_threads=args.threads,output_dir=principal_directory)

                            subprocess.call('rm {indivbams}'.format(indivbams=principal_directory + 'random_bam_*'), shell=True)

                            Modules_for_PipelineTB.Variant_calling(git=args.git,exclude_positions=args.exclude_positions,
                                                                   min_mapping_quality=args.min_mapping_quality,
                                                                   max_cov=Modules_for_PipelineTB.mean_Coverage_without_duplicates*3,
                                                                   min_cov=args.min_coverage,min_base_qual=args.avg_qual,min_reads2=args.min_reads,
                                                                   hetero_freq=args.var_freq,homo_freq=args.var_freq_homo,reference=REF,
                                                                   Gnumber=G_number,output_directory=principal_directory,log=LOG_LIST)

                            Modules_for_PipelineTB.Find_lineage(vcf_file=principal_directory+G_number+'.var.snp.vcf',git=args.git,log=LOG_LIST)

                            log.info('Producing whole-genome FASTA')
                            Modules_for_PipelineTB.WGSFasta(excluded=args.git+args.exclude_positions,vcf=principal_directory+G_number+'.all.pos.vcf',dir=principal_directory,Gnumber=G_number)

                            log.info('Compressing all.pos.vcf')
                            with open(principal_directory+G_number+'.all.pos.vcf', 'rb') as f_in, gzip.open(principal_directory+G_number+'.all.pos.vcf.gz', 'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                            os.remove(principal_directory+G_number+'.all.pos.vcf')

                            if args.user:
                                Modules_for_PipelineTB.Mapping_User_log(gnumber=G_number)

                            log.info('Pipeline End')

                else:
                    log.fatal('Fail at the fastQValidator step')
                    LOG_LIST.append(('FV_OUT',Modules_for_PipelineTB.stdout_FV[Modules_for_PipelineTB.stdout_FV.rfind(':') + 2:-1]))
                    bad_sample_file.write(G_number + '\t' + fastq[1] +'\tFail at FastQValidator Step\n') # fastq[1] and not path_to_fastq because want to know the original path that causes the problem
                    continue  # continue with the next iteration of the loop. ie: the fastq did not pass this step so no point continuing with the other steps of the pipeline


        elif match2 != None or match != None and 'SE' not in fastq: # PE
            log.info('Pipeline running mode: paired-end, multiple runs')

            List_multiple_fastq = []
            List_paths_fastq = {}

            for i in range(1, len(fastq)): # loop through the forward fastq paths
                '''Detect if the forward read is given or the reverse and infer the complementary read.'''
                forward_read = ''
                reverse_read = ''
                forward_R1_match = re.search(Forward_pattern_R1, fastq[i])
                forward_1_match = re.search(Forward_pattern_1, fastq[i])
                reverse_R2_match = re.search(Reverse_pattern_R1, fastq[i])
                reverse_2_match = re.search(Reverse_pattern_2, fastq[i])

                if forward_R1_match != None:
                    forward_read = fastq[i]
                    reverse_read = fastq[i][0:fastq[i].rfind('R1')] + 'R2' + fastq[i][fastq[i].rfind('R1') + 2:]
                elif forward_1_match != None:
                    forward_read = fastq[i]
                    reverse_read = fastq[i][0:fastq[i].rfind('_1.')] + '_2.'+ fastq[1][fastq[1].rfind('_1.')+3:]
                elif reverse_R2_match != None:
                    reverse_read = fastq[i]
                    forward_read = fastq[i][0:fastq[i].rfind('R2')] + 'R1' + fastq[i][fastq[i].rfind('R2') + 2:]
                elif reverse_2_match != None:
                    reverse_read = fastq[i]
                    forward_read = fastq[i][0:fastq[i].rfind('_2.')] + '_1.'+ fastq[1][fastq[1].rfind('_2.')+3:]

                LOG_LIST.append(('FQ_IN',forward_read))
                LOG_LIST.append(('FQ_IN',reverse_read))

                basename_forward = Modules_for_PipelineTB.Extract_basename(path=forward_read)[0]
                basename_fastq_forward = Modules_for_PipelineTB.Extract_basename(path=forward_read)[1]

                basename_reverse = Modules_for_PipelineTB.Extract_basename(path=reverse_read)[0]
                basename_fastq_reverse = Modules_for_PipelineTB.Extract_basename(path=reverse_read)[1]

                '''Detect the extension  (gz,bz2,tar, tar.gz,tar.bz2,tgz,zip,Z,) the idea is that the format should be either uncompressed 'tempfile' or gz compressed . WHAT ABOUT THE IDEA GENOMES (BAM)???'''
                if basename_forward.endswith('fastq.bz2'):
                    subprocess.call(['bzip2', '-dfkq',forward_read])  # should add the output in a tem directory: ASK Mireia. Options here are: decompress -overwrite output, quiet and keep input
                    forward_read = forward_read[:-4]  # new path to the fastq is the old path with without the bz2 extension. This works because it is in the same directory as the bz2. but ideally should be in a temp directory
                elif basename_reverse.endswith('fastq.bz2'):
                    subprocess.call(['bzip2', '-dfkq',reverse_read])  # should add the output in a tem directory: ASK Mireia. Options here are: decompress -overwrite output, quiet and keep input
                    reverse_read = reverse_read[:-4]


                List_multiple_fastq += [(forward_read, reverse_read)]

                '''Check for PHRED encoding. If phred 64, convert to phred33'''
                tf_forward_read_phred33 = tempfile.NamedTemporaryFile()
                tf_reverse_read_phred33 = tempfile.NamedTemporaryFile()

                try: # try opening a normal fasq file if it does not work use the gzip module
                    raw_fastq = gzip.open(forward_read,'r')
                    first_line =raw_fastq.readline().split(':')
                    for line in itertools.islice(raw_fastq,3,None,4): # lopp through the first 1000 lines starting from the 3rd line and taking a step of 4
                        phred33 = re.search(r'[\d\+\*\(\)\."#\$%&-,:\!\']',line) # look for the specific phred33 encoding characters in the lines
                        phred64 = re.search(r'[K-Z\[\]^_`a-h]',line)
                        if phred33 != None: #  if the encoding of the sample is 33
                            LOG_LIST.append(('PHE_OUT','33'))
                            break
                        elif phred64 != None:
                            SeqIO.convert(forward_read,"fastq-illumina", tf_forward_read_phred33.name+'.fastq', "fastq-sanger")
                            SeqIO.convert(reverse_read,"fastq-illumina", tf_reverse_read_phred33.name+'.fastq', "fastq-sanger")

                            forward_read = tf_forward_read_phred33.name+'.fastq'
                            reverse_read = tf_reverse_read_phred33.name+'.fastq'
                            LOG_LIST.append(('PHE_OUT','64'))
                            break
                except:
                    raw_fastq = open(forward_read,'r')
                    first_line =raw_fastq.readline().split(':')
                    for line in itertools.islice(raw_fastq,3,None,4): # lopp through the first 1000 lines starting from the 3rd line and taking a step of 4
                        phred33 = re.search(r'[\d\+\*\(\)\."#\$%&-,:\!\']',line) # look for the specific phred33 encoding characters in the lines
                        phred64 = re.search(r'[K-Z\[\]^_`a-h]',line)
                        if phred33 != None: #  if the encoding of the sample is phred33
                            LOG_LIST.append(('PHE_OUT','33'))
                            break

                        elif phred64 != None:
                            SeqIO.convert(forward_read,"fastq-illumina", tf_forward_read_phred33.name+'.fastq', "fastq-sanger")
                            SeqIO.convert(reverse_read,"fastq-illumina", tf_reverse_read_phred33.name+'.fastq', "fastq-sanger")

                            forward_read = tf_forward_read_phred33.name+'.fastq'
                            reverse_read = tf_reverse_read_phred33.name+'.fastq'
                            LOG_LIST.append(('PHE_OUT','64'))
                            break

                Modules_for_PipelineTB.Run_FastQValidator(fastq=forward_read,minimum_read_length=args.minlen,log=LOG_LIST)
                log.info('FastQValidator forward read: %s'% Modules_for_PipelineTB.stdout_FV)
                if Modules_for_PipelineTB.stdout_FV[Modules_for_PipelineTB.stdout_FV.rfind(':') + 2:-1] == 'FASTQ_SUCCESS':  # if the forward fastq passes the fastQValidator test. What if only one fastq passes the test?
                    LOG_LIST.append(('FV_OUT','FASTQ_SUCCESS'))
                    Modules_for_PipelineTB.Run_FastQValidator(fastq=reverse_read,minimum_read_length=args.minlen,log=LOG_LIST)
                    log.info('FastQValidator reverse read: %s'% Modules_for_PipelineTB.stdout_FV)

                    if Modules_for_PipelineTB.stdout_FV[Modules_for_PipelineTB.stdout_FV.rfind(':') + 2:-1] == 'FASTQ_SUCCESS':  # if the forward fastq passes the fastQValidator test.
                        LOG_LIST.append(('FV_OUT','FASTQ_SUCCESS'))

                        if i < len(fastq) - 1:  # while there are still more fastq (of the same gnumber) to validate continue
                                continue

                        elif i == len(fastq) - 1:  # if this condition is met it means that fastqValidator has been run of every fastq so now we want to do fastQC to control the LNG and read length

                            sum_LNG = 0
                            count = 1

                            for qc in List_multiple_fastq: # qc is a tuple (R1, R2). it would be the first run (forward and reverse), then the second run etc...

                                basename_forward=Modules_for_PipelineTB.Extract_basename(qc[0])[0] # qc[0] to select the forward read as qc is a tuple
                                basename_fastq_forward=Modules_for_PipelineTB.Extract_basename(qc[0])[1]

                                log.info('Starting FastQC on %s'%qc[0])

                                Modules_for_PipelineTB.Run_FastQC(fastq=qc[0],output_directory=principal_directory,log=LOG_LIST,gnumber=G_number,basename=basename_fastq_forward, number_of_threads=args.threads)

                                if type(Modules_for_PipelineTB.read_length) is str:  #when a range of read lengths i.e 35-101
                                    sum_LNG = 'NA'
                                    if count < len(List_multiple_fastq):  # process the other fastqs
                                        count += 1
                                        continue

                                elif type(Modules_for_PipelineTB.read_length) is int:
                                    LNG_coverage = (float(Modules_for_PipelineTB.total_reads)*Modules_for_PipelineTB.read_length*2)/4411532
                                    try:
                                        sum_LNG += LNG_coverage # do a try except kindof thing (sometimes not possible because one run has a range of read lengths (str) and the other run has a fixed read leangth (int)
                                    except:
                                        sum_LNG = 'NA'
                                        log.info('Could not calculate LNG because read lengths are variable')

                                    if count < len(List_multiple_fastq):  # process the other fastqs ## should'nt this be at the same level as the elif???
                                        count += 1
                                        continue

                                LOG_LIST.append(('LNG_OUT',str(sum_LNG)))
                                n = 0
                                for fq in List_multiple_fastq:

                                    forward_read = fq[0]
                                    reverse_read = fq[1]
                                    run_number = str(n+1)

                                    log.info('Starting Trimmomatic for run %s' % run_number)
                                    Modules_for_PipelineTB.Run_Trimmomatic_PE(number_of_threads=args.threads,minimum_read_length_post_trimming=args.MinLength,log=LOG_LIST,forwardRead=forward_read,reverseRead=reverse_read)
                            
                                    log.info('Starting SeqPrep for run %s' % run_number)
                                    Modules_for_PipelineTB.Run_SeqPrep(forward_paired = Modules_for_PipelineTB.Fastq1P, reverse_paired=Modules_for_PipelineTB.Fastq2P,overlap=args.overlap, minlen=args.MinLength,log=LOG_LIST)

                                    '''Parse the fastq header'''
                                    first_line = [] # where the fastq header will be
                                    try:  # try opening a normal fastq file if it does not work use the gzip module
                                        raw_fastq = gzip.open(forward_read, 'r')
                                        first_line = raw_fastq.readline().split(':')
                                    except:
                                        raw_fastq = open(forward_read, 'r')
                                        first_line = raw_fastq.readline().split(':')

                                    letters = list(string.ascii_uppercase) # letters = ['A','B','C','D',...,'Z']

                                    flowcell_ID = ''
                                    flowcell_lane = ''
                                    if len(first_line) == 10:  # if the fastq header is composed of 10 elements
                                        flowcell_ID = first_line[2]  # then select the second element which is the flowcell ID
                                        flowcell_lane = first_line[3]  # and select the third element which is the flowcell lane
                                    elif len(first_line) != 10:  # if the fastq header is not composed of 10 elements
                                        flowcell_ID = G_number
                                        flowcell_lane = letters[n] # n starts at 0 with the first run and increments  each time a new run is processed

                                    tf_bam_se_trimmomatic_1U = tempfile.NamedTemporaryFile(prefix='mapped_sorted_se_trimmomatic_1U')
                                    tf_bam_se_trimmomatic_2U = tempfile.NamedTemporaryFile(prefix='mapped_sorted_se_trimmomatic_2U')
                                    tf_bam_se_SeqPrep_merge = tempfile.NamedTemporaryFile(prefix='mapped_sorted_se_SeqPrep_merge')
                                    tf_bam_pe_SeqPrep_trimmed = tempfile.NamedTemporaryFile(prefix = 'mapped_sorted_pe_SeqPrep_trimmed')
                                    tf_bam_pe_SeqPrep_NOTmerged = tempfile.NamedTemporaryFile(prefix = 'mapped_sorted_pe_SeqPrep_NOTmerged')

                                    log.info('Start BWA for run %s: mapping unpaired forward reads from Trimmomatic' % run_number)
                                    Modules_for_PipelineTB.BWA_mapping_SE(inreads=Modules_for_PipelineTB.Fastq1U,
                                                   tempbamfile=tf_bam_se_trimmomatic_1U.name,
                                                   flowcell=flowcell_ID,lane=flowcell_lane,Gnumber=G_number,number_of_threads=args.threads,log=LOG_LIST,reference_genome=REF)

                                    log.info('Start BWA for run %s: mapping unpaired reverse reads from Trimmomatic' % run_number)
                                    Modules_for_PipelineTB.BWA_mapping_SE(inreads=Modules_for_PipelineTB.Fastq2U,
                                                   tempbamfile=tf_bam_se_trimmomatic_2U.name,
                                                   flowcell=flowcell_ID, lane=flowcell_lane,
                                                   Gnumber=G_number,number_of_threads=args.threads,log=LOG_LIST,reference_genome=REF)

                                    log.info('Start BWA for run %s: mapping merged reads from SeqPrep' % run_number)
                                    Modules_for_PipelineTB.BWA_mapping_SE(inreads=Modules_for_PipelineTB.SeqPrep_Merged_Reads,
                                                   tempbamfile=tf_bam_se_SeqPrep_merge.name,
                                                   flowcell=flowcell_ID, lane=flowcell_lane,
                                                   Gnumber=G_number,number_of_threads=args.threads,log=LOG_LIST,reference_genome=REF)

                                    log.info('Start BWA for run %s: mapping trimmed reads from SeqPrep' % run_number)
                                    Modules_for_PipelineTB.BWA_mapping_PE(fw_reads=Modules_for_PipelineTB.SeqPrep_Trimmed_Reads_Forward,
                                                   rev_reads=Modules_for_PipelineTB.SeqPrep_Trimmed_Reads_Reverse,
                                                   tempbamfile=tf_bam_pe_SeqPrep_trimmed.name,flowcell=flowcell_ID,lane=flowcell_lane,Gnumber=G_number,number_of_threads=args.threads,log=LOG_LIST,reference_genome=REF)

                                    log.info('Start BWA for run %s: mapping non merged reads from Seqprep' % run_number)
                                    Modules_for_PipelineTB.BWA_mapping_PE(fw_reads=Modules_for_PipelineTB.SeqPrep_Non_Merged_Reads_Forward,
                                                   rev_reads=Modules_for_PipelineTB.SeqPrep_Non_Merged_Reads_Reverse,
                                                   tempbamfile=tf_bam_pe_SeqPrep_NOTmerged.name,flowcell=flowcell_ID,lane=flowcell_lane,Gnumber=G_number,number_of_threads=args.threads,log=LOG_LIST,reference_genome=REF)

                                    tf_bam = tempfile.NamedTemporaryFile(prefix = 'mapped_sorted_merged')
                                    log.info('Merging all Bam files produced from  run %s' % run_number)
                                    samtools_merge_command = 'samtools merge {output_bam} {inbam_trimmomatic_1U} {inbam_trimmomatic_2U} {inbam_SeqPrep_merged} {inbam_SeqPrep_trimmed} {inbam_SeqPrep_NOTmerged}'.format(
                                    output_bam=tf_bam.name + '.bam',
                                    inbam_trimmomatic_1U=tf_bam_se_trimmomatic_1U.name+'.bam',
                                    inbam_trimmomatic_2U=tf_bam_se_trimmomatic_2U.name+'.bam',
                                    inbam_SeqPrep_merged=tf_bam_se_SeqPrep_merge.name+'.bam',
                                    inbam_SeqPrep_trimmed=tf_bam_pe_SeqPrep_trimmed.name+'.bam',
                                    inbam_SeqPrep_NOTmerged=tf_bam_pe_SeqPrep_NOTmerged.name+'.bam')

                                    subprocess.call(samtools_merge_command,shell=True)
                                    LOG_LIST.append(('SM1_IN',samtools_merge_command))

                                    log.info('Starting MarkDuplicates for run %s' % run_number)
                                    Modules_for_PipelineTB.Run_MarkDuplicates(input=tf_bam.name+'.bam', principaldirectory=principal_directory,log=LOG_LIST)

                                    n+=1
                                log.info('Merging BAM files from all different runs')
                                tf_bam_merged = tempfile.NamedTemporaryFile(prefix='merged_bam_')
                                command_merge = 'samtools merge {outbam} {inbams}'.format(
                                            outbam=tf_bam_merged.name + '.bam',
                                            inbams=principal_directory + '*.bam')
                                subprocess.call(command_merge, shell=True)
                                LOG_LIST.append(('SM2_IN',command_merge))

                                log.info('Indexing the intermediate BAM file')
                                Modules_for_PipelineTB.Index_BAM(input=tf_bam_merged.name + '.bam',log=LOG_LIST,name_of_index='IDX1_IN')

                                subprocess.call('rm {indivbams}'.format(indivbams=principal_directory + 'random_bam_*'), shell=True)

                                log.info('Starting GATK Realigner Target Creator')
                                Modules_for_PipelineTB.GATK_RealignerTargetCreator(log=LOG_LIST,input=tf_bam_merged.name + '.bam',number_threads=args.threads, reference_genome=REF)

                                tf_bam_postGATK = tempfile.NamedTemporaryFile(prefix='postGATK_bam_')
                                log.info('Starting GATK IndelRealigner')
                                Modules_for_PipelineTB.GATK_IndelRealigner(log=LOG_LIST,input=tf_bam_merged.name + '.bam',output=tf_bam_postGATK.name+'.bam',reference_genome=REF,interval_file=Modules_for_PipelineTB.GATK_intervals)

                                log.info('Indexing the temp BAM file')
                                Modules_for_PipelineTB.Index_BAM(input=tf_bam_postGATK.name+'.bam',log=LOG_LIST,name_of_index='IDX2_IN')

                                Modules_for_PipelineTB.Samtools_Flagstat(input_bam=tf_bam_postGATK.name+'.bam',log=LOG_LIST,contaminant_threshold=args.contaminant,gnb=G_number)

                                log.info('Filtering BAM with PySam')
                                Modules_for_PipelineTB.FilterBAM_AS(infile=tf_bam_postGATK.name+'.bam',gnumber=G_number,principaldirectory=principal_directory)

                                log.info('Indexing the final BAM file')
                                Modules_for_PipelineTB.Index_BAM(input=principal_directory+G_number+'.bam',log=LOG_LIST,name_of_index='IDX3_IN')

                                log.info('Starting Qualimap BamQC')
                                Modules_for_PipelineTB.Qualimap_BamQC(input_bam=principal_directory+G_number+'.bam',log=LOG_LIST,gnumber=G_number,number_threads=args.threads,output_dir=principal_directory)


                                Modules_for_PipelineTB.Variant_calling(git=args.git,exclude_positions=args.exclude_positions,
                                                                        min_mapping_quality=args.min_mapping_quality,
                                                                        max_cov=Modules_for_PipelineTB.mean_Coverage_without_duplicates*3,
                                                                        min_cov=args.min_coverage,min_base_qual=args.avg_qual,min_reads2=args.min_reads,
                                                                        hetero_freq=args.var_freq,homo_freq=args.var_freq_homo,reference=REF,
                                                                        Gnumber=G_number,output_directory=principal_directory,log=LOG_LIST)
                                Modules_for_PipelineTB.Find_lineage(vcf_file=principal_directory+G_number+'.var.snp.vcf',git=args.git,log=LOG_LIST)

                                log.info('Producing whole-genome FASTA')
                                Modules_for_PipelineTB.WGSFasta(excluded=args.git+args.exclude_positions,vcf=principal_directory+G_number+'.all.pos.vcf',dir=principal_directory,Gnumber=G_number)
                                log.info('Compressing all.pos.vcf')
                                with open(principal_directory+G_number+'.all.pos.vcf', 'rb') as f_in, gzip.open(principal_directory+G_number+'.all.pos.vcf.gz', 'wb') as f_out:
                                    shutil.copyfileobj(f_in, f_out)
                                os.remove(principal_directory+G_number+'.all.pos.vcf')

                                if args.user:
                                    Modules_for_PipelineTB.Mapping_User_log(gnumber=G_number)


                                log.info('Pipeline End')

                    else:
                        log.fatal('Fail at the fastQValidator step')
                        LOG_LIST.append(('FV_OUT',Modules_for_PipelineTB.stdout_FV[Modules_for_PipelineTB.stdout_FV.rfind(':') + 2:-1]))
                        bad_sample_file.write(G_number + '\t' + fastq[1] +'\tFail at FastQValidator Step\n') # fastq[1] and not path_to_fastq because want to know the original path that causes the problem
                        continue  # continue with the next iteration of the loop. ie: the fastq did not pass this step so no point continuing with the other steps of the pipeline
                else:
                    log.fatal('Fail at the fastQValidator step')
                    LOG_LIST.append(('FV_OUT',Modules_for_PipelineTB.stdout_FV[Modules_for_PipelineTB.stdout_FV.rfind(':') + 2:-1]))
                    bad_sample_file.write(G_number + '\t' + fastq[1] +'\tFail at FastQValidator Step\n') # fastq[1] and not path_to_fastq because want to know the original path that causes the problem
                    continue  # continue with the next iteration of the loop. ie: the fastq did not pass this step so no point continuing with the other steps of the pipeline



'''writing to a log file'''
for ky,vl in LOG_LIST:
    LOG_DIC[ky].append(vl)
ListLogID_PEMR = ['DAT_OUT','GID_IN','FQ_IN','PHE_OUT','FV_IN','FV_OUT','FQC_IN','TRD_OUT','RDL_OUT','LNG_OUT','TRM_IN','TRM_OUT','ASV_OUT','FSV_OUT','RSV_OUT','ADP_OUT','SQP_IN','SIN_OUT','PM_OUT','PA_OUT','PD_OUT','MG_OUT','BWA_IN','SM1_IN','MD_IN','SM2_IN','IDX1_IN','RTC_IN','IDR_IN','BAMF_IN','IDX2_IN','IDX3_IN','SFS_IN','RDC_OUT','DUP_OUT','PP_OUT','MAP_OUT','BQC_IN','COV_OUT','STD_OUT','NC_OUT','N7_OUT','VC_IN','BPF_OUT','VFL_OUT','VAR_OUT','SNP_OUT','HOM_OUT','HET_OUT','HOMO_IND_OUT','ANN_IN','LIN_OUT']
with open(principal_directory+G_number+'.log','w') as outfile:
    for itm in ListLogID_PEMR:
        if LOG_DIC.has_key(itm):
            outfile.write(itm+'\t'+";".join(LOG_DIC[itm])+'\n')
        else:
            outfile.write(itm+'\tNA\n')

bad_sample_file.close()
