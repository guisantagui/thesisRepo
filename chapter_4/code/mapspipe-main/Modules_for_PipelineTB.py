__author__ = 'chloeloiseau'

import subprocess
import pysam
import tempfile
import re
import gzip
import itertools
import random
import string
from Bio import SeqIO
import os
import datetime
import getpass
import csv
import logging


#duda ruta database_update


TRIM_LOC = '/storage/PGO/bin/trimmomatic.jar'
ADAPTOR_FILE_PE = '/storage/PGO/store/db_pipe/adapters/NexteraPE-PE.fa'
ADAPTOR_FILE_SE = '/storage/PGO/store/db_pipe/adapters/TruSeq3-SE.fa'
PICARD = '/storage/PGO/bin/picard.jar'
GATK = '/storage/PGO/bin/GenomeAnalysisTK.jar'
VARSCAN='/storage/PGO/bin/VarScan.v2.4.1.jar'

class AutoVivification(dict):
    """Implementation of autovivification in Python"""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

Data= AutoVivification()

def Extract_basename(path):
    '''From an absolute path, this function will return a tuple where the first component  of the tuple is basename (example:G0001.fastq.gz) and the second component is the basename without extension (example:G0001)'''
    basename = os.path.basename(path)  #  Example: from /Users/chloeloiseau/Desktop/G0001.fastq.gz get G0001.fastq.gz
    index_extension_fastq = basename.find('.') # find the coordinate (index) of the first '.' in the path. Should be the beginning of the extension. What happens if the basename contains a dot.
    basename_without_extension = basename[:index_extension_fastq] # Example: from G0001.fastq.gz get G0001
    return basename, basename_without_extension # returns a tuple
    
def Run_FastQValidator(fastq, minimum_read_length,log):
    global stdout_FV
    logging.warning("Estamos usando Modules_for_PipelineTB adaptado a garnatxa para MTB")
    
    cmd_fastQValidator = 'fastQValidator --noeof --disableSeqIDCheck --minReadLen {minreadlen} --file {path}'.format(
        minreadlen=minimum_read_length, path=fastq)

    output = subprocess.Popen(cmd_fastQValidator, shell=True, stdout=subprocess.PIPE)
    log.append(('FV_IN',cmd_fastQValidator))

    stdout_FV = output.communicate()[0]

def Run_FastQC(fastq, output_directory, log, gnumber, basename,number_of_threads):
    global total_reads
    global read_length
    cmd_fastQC = '/storage/PGO/bin/instalaciones/FastQC/fastqc -t {threads} -q {fastq} --nogroup --noextract -o {outputdirectory}'.format(threads=number_of_threads,fastq=fastq,
                                                                                       outputdirectory=output_directory)
    subprocess.call(cmd_fastQC, shell=True)
    log.append(('FQC_IN', cmd_fastQC))
    subprocess.call('unzip -c ' + output_directory + basename + '_fastqc.zip \*fastqc_data.txt > ' + output_directory + gnumber + '_fastqc_data.txt',shell=True)  # extract the fastqc_data file from the compressed folder without having to extract it.

    ifile_fastqc_data_summary = open(output_directory + gnumber + '_fastqc_data.txt','r')  # open the file created by fastqc and that contains statistics about the fastQ
    table_fastqc_summary = [line.strip().split('\t') for line in ifile_fastqc_data_summary]  # make a list out of the data
    total_reads = int()
    read_length = int()  # sometimes this is not a number. Ex 8-328. If it is a range then do not calculate LNG value.
    for i in table_fastqc_summary:  # loop in the list with all the data in order to find relevant information : total number of reads and read length
        if i[0] == 'Total Sequences':
            total_reads = int(i[1])
        if i[0] == 'Sequence length':
            if '-' in i[1]:
                read_length = str(i[1])  # eg. 8-320
            else:
                read_length = int(i[1])

    ifile_fastqc_data_summary.close()
    log.append(('TRD_OUT', str(total_reads)))
    log.append(('RDL_OUT',str(read_length)))

    # Remove the fastQC outputs:
    subprocess.call('rm ' + output_directory + gnumber + '_fastqc_data.txt', shell=True)
    subprocess.call('rm ' + output_directory + basename + '_fastqc.zip', shell=True)
    subprocess.call('rm ' + output_directory + basename + '_fastqc.html', shell=True)

def Run_Trimmomatic_SE(number_of_threads,fastq,minimum_read_length_post_trimming,log):
    global fastq_for_BWA
    tf_trimmomatic_log = tempfile.NamedTemporaryFile()
    tf_trimmomatic = tempfile.NamedTemporaryFile()
    trimmomatic_cmd = "java -jar {trimmomatic} SE -threads {threads} -phred33 {fastq} {fastqtrimmed} ILLUMINACLIP:{adapter}:2:30:10 SLIDINGWINDOW:5:20 MINLEN:{minlen} 2>> {logfile}".format(
        trimmomatic=TRIM_LOC, threads=number_of_threads, fastq=fastq,fastqtrimmed=tf_trimmomatic.name + '.fastq.gz', adapter=ADAPTOR_FILE_SE, minlen=minimum_read_length_post_trimming,
                                                                                                     logfile=tf_trimmomatic_log.name + '_trimmomatic.log')
    subprocess.call(trimmomatic_cmd, shell=True)
    fastq_for_BWA = tf_trimmomatic.name+'.fastq.gz'

    log.append(('TRM_IN', trimmomatic_cmd))

    trim_log_table = []
    with open(tf_trimmomatic_log.name + '_trimmomatic.log', 'r+') as Trimmomatic_log_file:
        for trim_log_i in Trimmomatic_log_file:
            trim_log_table += [trim_log_i.strip()]

    #logging.warning(trim_log_table)
    TRIM_OUT = trim_log_table[1] + '|' + trim_log_table[2] + '|' + trim_log_table[3]
    log.append(('TRM_OUT', TRIM_OUT))

    trimmomatic_stats = trim_log_table[4].split()

    RIN_OUT = trimmomatic_stats[2]
    log.append(('RIN_OUT', RIN_OUT))

    ASV_OUT = trimmomatic_stats[5][1:-1]
    log.append(('ASV_OUT', ASV_OUT))

    ADP_OUT = trimmomatic_stats[8][1:-1]
    log.append(('ADP_OUT', ADP_OUT))

def Run_Trimmomatic_PE(number_of_threads,minimum_read_length_post_trimming,log,forwardRead,reverseRead):
    global Fastq1P
    global Fastq1U
    global Fastq2P
    global Fastq2U
    tf_trimmomatic_log = tempfile.NamedTemporaryFile()
    tf_1P_trimmomatic = tempfile.NamedTemporaryFile()
    tf_1U_trimmomatic = tempfile.NamedTemporaryFile()
    tf_2P_trimmomatic = tempfile.NamedTemporaryFile()
    tf_2U_trimmomatic = tempfile.NamedTemporaryFile()

    trimmomatic_cmd = "java -jar {trimmomatic} PE -threads {threads} -phred33 {forward} {reverse} {forward_P} {forward_U} {reverse_P} {reverse_U} ILLUMINACLIP:{adapter}:2:30:10 SLIDINGWINDOW:5:20 MINLEN:{minlen} 2>> {logfile}".format(
        trimmomatic=TRIM_LOC, threads=number_of_threads, forward=forwardRead, reverse=reverseRead,
        forward_P=tf_1P_trimmomatic.name + '_trimmomatic_1P.fastq.gz',
        forward_U=tf_1U_trimmomatic.name + '_trimmomatic_1U.fastq.gz',
        reverse_P=tf_2P_trimmomatic.name + '_trimmomatic_2P.fastq.gz',
        reverse_U=tf_2U_trimmomatic.name + '_trimmomatic_2U.fastq.gz', adapter=ADAPTOR_FILE_PE, minlen=minimum_read_length_post_trimming,
        logfile=tf_trimmomatic_log.name + '_trimmomatic.log')
    subprocess.call(trimmomatic_cmd, shell=True)
    Fastq1P = tf_1P_trimmomatic.name + '_trimmomatic_1P.fastq.gz'
    Fastq1U = tf_1U_trimmomatic.name + '_trimmomatic_1U.fastq.gz'
    Fastq2P = tf_2P_trimmomatic.name + '_trimmomatic_2P.fastq.gz'
    Fastq2U = tf_2U_trimmomatic.name + '_trimmomatic_2U.fastq.gz'
    log.append(('TRM_IN', trimmomatic_cmd))

    trim_log_table = []
    with open(tf_trimmomatic_log.name + '_trimmomatic.log', 'r+') as Trimmomatic_log_file:
        for trim_log_i in Trimmomatic_log_file:
            trim_log_table += [trim_log_i.strip()]


    #logging.warning(trim_log_table)
    
    
    TRIM_OUT = trim_log_table[1] + '|' + trim_log_table[2] + '|' + trim_log_table[3] + '|' + trim_log_table[
        4] + '|' + trim_log_table[5] + '|' + trim_log_table[6]
    log.append(('TRM_OUT', TRIM_OUT))

    trimmomatic_stats = trim_log_table[7].split()  # Turn Input Read Pairs: 3712483 Both Surviving: 3496331 (94.18%) Forward Only Surviving: 128190 (3.45%) Reverse Only Surviving: 42805 (1.15%) Dropped: 45157 (1.22%)
    # into: ['Input', 'Read', 'Pairs:', '3712483', 'Both', 'Surviving:', '3496331', '(94.18%)', 'Forward', 'Only', 'Surviving:', '128190', '(3.45%)', 'Reverse', 'Only', 'Surviving:', '42805', '(1.15%)', 'Dropped:', '45157', '(1.22%)']
    
    
    
    RIN_OUT = trimmomatic_stats[3]
    log.append(('RIN_OUT', RIN_OUT))

    ASV_OUT = trimmomatic_stats[7][1:-1]  # 1:-1 to remove the parenthesis
    log.append(('ASV_OUT', ASV_OUT))

    FSV_OUT = trimmomatic_stats[12][1:-1]
    log.append(('FSV_OUT', FSV_OUT))

    RSV_OUT = trimmomatic_stats[17][1:-1]
    log.append(('RSV_OUT', RSV_OUT))

    ADP_OUT = trimmomatic_stats[20][1:-1]
    log.append(('ADP_OUT', ADP_OUT))

def Run_SeqPrep(forward_paired,reverse_paired,overlap,minlen,log):
    '''SeqPrep: merge overlapping PE reads'''
    global SeqPrep_Merged_Reads
    global SeqPrep_Non_Merged_Reads_Forward
    global SeqPrep_Non_Merged_Reads_Reverse
    global SeqPrep_Trimmed_Reads_Forward
    global SeqPrep_Trimmed_Reads_Reverse

    tf_SeqPrep_forward_trimmed = tempfile.NamedTemporaryFile()
    tf_SeqPrep_reverse_trimmed = tempfile.NamedTemporaryFile()
    tf_SeqPrep_forward_notmerged = tempfile.NamedTemporaryFile()
    tf_SeqPrep_reverse_notmerged = tempfile.NamedTemporaryFile()
    tf_SeqPrep_merged = tempfile.NamedTemporaryFile()
    tf_SeqPrep_alignment_log = tempfile.NamedTemporaryFile()
    tf_SeqPrep_log = tempfile.NamedTemporaryFile()  # what should be a good minimum base pair overlap to merge two reads ??????


    #cambiado, anyadida la -S, ennable spinner, porque salia Pairs Processed: 0 siempre
    SeqPrep_cmd = 'SeqPrep -S -f {forward_P} -r {reverse_P} -1 {SeqPrep_forward_trimmed} -2 {SeqPrep_reverse_trimmed} -3 {SeqPrep_forward_notmerged} -4 {SeqPrep_reverse_notmerged} ' \
                  '-A AGATCGGAAGAGCACACGTCT -B AGATCGGAAGAGCGTCGTGTA -L {minlen} -o {OverlapSize} -s {SeqPrep_merged}' \
                  ' -E {SeqPrep_alignment} 2>> {SeqPrep_log}'.format(
        forward_P=forward_paired,
        reverse_P=reverse_paired,
        SeqPrep_forward_trimmed=tf_SeqPrep_forward_trimmed.name + '_SeqPrep_1_trimmed.fastq.gz',
        SeqPrep_reverse_trimmed=tf_SeqPrep_reverse_trimmed.name + '_SeqPrep_2_trimmed.fastq.gz',
        SeqPrep_forward_notmerged=tf_SeqPrep_forward_notmerged.name + '_SeqPrep_1_notmerged.fastq.gz',
        SeqPrep_reverse_notmerged=tf_SeqPrep_reverse_notmerged.name + '_SeqPrep_2_notmerged.fastq.gz',
        minlen=minlen,
        OverlapSize=overlap,
        SeqPrep_merged=tf_SeqPrep_merged.name + '_SeqPrep_merged.fastq.gz',
        SeqPrep_alignment=tf_SeqPrep_alignment_log.name,
        SeqPrep_log=tf_SeqPrep_log.name + '_SeqPrep.log')

    subprocess.call(SeqPrep_cmd, shell=True)

    log.append(('SQP_IN', SeqPrep_cmd))

    SeqPrep_Merged_Reads = tf_SeqPrep_merged.name + '_SeqPrep_merged.fastq.gz'
    SeqPrep_Non_Merged_Reads_Forward = tf_SeqPrep_forward_notmerged.name+'_SeqPrep_1_notmerged.fastq.gz'
    SeqPrep_Non_Merged_Reads_Reverse = tf_SeqPrep_reverse_notmerged.name + '_SeqPrep_2_notmerged.fastq.gz'
    SeqPrep_Trimmed_Reads_Forward = tf_SeqPrep_forward_trimmed.name+ '_SeqPrep_1_trimmed.fastq.gz'
    SeqPrep_Trimmed_Reads_Reverse = tf_SeqPrep_reverse_trimmed.name + '_SeqPrep_2_trimmed.fastq.gz'

    table_seqprep_merging = []  # create an empty list
    SeqPrep_logfile = open(tf_SeqPrep_log.name + '_SeqPrep.log',
                           'r+')  # open the log file for reading
   
    for row in SeqPrep_logfile:

        table_seqprep_merging += [
            row]  # Each element of the list is one line from the log file
    table_seqprep_merging = table_seqprep_merging[1:]  # Discard the first line which corresponds to |/-\|
    pairs_processed = table_seqprep_merging[0]  # Pairs processed stored in a variable : 'Pairs Processed:\t901448\n'
    pairs_merged = table_seqprep_merging[1]  # Pairs Merged stored in a variable
    pairs_with_adapters = table_seqprep_merging[2]
    pairs_discarded = table_seqprep_merging[3]

    regex = r'[0-9]{1,20}'  # regular expression: all digits

    match_merged = re.search(regex, pairs_merged)  # Search for regex in pairs merged
    result_merged = match_merged.group()  # Store the result of the search in a variable

    match_processed = re.search(regex, pairs_processed)  # Search for regex in pairs processed
    result_processed = match_processed.group()  # Store the result of the search in a variable

    match_adapters = re.search(regex, pairs_with_adapters)
    result_adapters = match_adapters.group()

    match_pairs_discarded = re.search(regex, pairs_discarded)
    result_pairs_discarded = match_pairs_discarded.group()

    total_merge = round(float(result_merged) / float(result_processed) * 100,
                        2)  # percentage of reads merged


    SeqPrep_logfile.close()

    log.append(('SIN_OUT', str(result_processed)))
    log.append(('PM_OUT', str(result_merged)))
    log.append(('PA_OUT', str(result_adapters)))
    log.append(('PD_OUT', str(result_pairs_discarded)))
    log.append(('MG_OUT', str(total_merge) + '%'))

def BWA_mapping_SE(inreads, tempbamfile, flowcell, lane, Gnumber, number_of_threads,reference_genome,log):
    #global min_AS
    #min_AS = (0.93*int(RDL))-(int(RDL)*4*0.07)
    bwa_command_se = "bwa mem -t {threads} -v 1 -M -R '@RG\\tID:{flwcl}\\tSM:{samplename}' {reference} {input}".format(
        threads=number_of_threads,flwcl=flowcell + '-' + lane, samplename=Gnumber,
        reference=reference_genome,
        input=inreads)
    samtools_view_command = 'samtools view -Sbhu -'
    bwa_se = subprocess.Popen(bwa_command_se, shell=True,
                              stdout=subprocess.PIPE)

    samtools_view_se = subprocess.Popen(samtools_view_command, shell=True,
                                        stdin=bwa_se.stdout,
                                        stdout=subprocess.PIPE)
    samtools_sort_command_se = 'samtools sort -m 4G - {tempbam}'.format(
        tempbam=tempbamfile)
    samtools_sort_se = subprocess.Popen(samtools_sort_command_se, shell=True,
                                        stdin=samtools_view_se.stdout)

    bwa_se.stdout.close()
    samtools_view_se.stdout.close()
    samtools_sort_se.communicate()

    log.append(('BWA_IN',bwa_command_se + '|' + samtools_view_command + '|' + samtools_sort_command_se))

def BWA_mapping_PE(fw_reads,rev_reads,tempbamfile,flowcell,lane,Gnumber,number_of_threads,reference_genome,log):
    #global min_AS
    #min_AS = (0.93*int(RDL))-(int(RDL)*4*0.07)

    bwa_command_pe = "bwa mem -t {threads} -v 1 -M -R '@RG\\tID:{flwcl}\\tSM:{samplename}' {reference} {R1} {R2}".format(
        threads=number_of_threads,flwcl=flowcell+ '-' + lane, samplename=Gnumber, reference=reference_genome,R1=fw_reads,R2=rev_reads)
    samtools_view_command = 'samtools view -Sbhu -'
    bwa_pe = subprocess.Popen(bwa_command_pe, shell=True, stdout=subprocess.PIPE)
    samtools_view_pe = subprocess.Popen(samtools_view_command, shell=True, stdin=bwa_pe.stdout, stdout=subprocess.PIPE)
    samtools_sort_command_pe = 'samtools sort -m 4G - {tempbam}'.format(
        tempbam=tempbamfile)

    samtools_sort_pe = subprocess.Popen(samtools_sort_command_pe, shell=True, stdin=samtools_view_pe.stdout)
    bwa_pe.stdout.close()
    samtools_view_pe.stdout.close()
    samtools_sort_pe.communicate()


    log.append(('BWA_IN',bwa_command_pe + '|' + samtools_view_command + '|' + samtools_sort_command_pe))

def FilterBAM_AS(infile,gnumber,principaldirectory):
    bam = pysam.AlignmentFile(infile, "rb")
    filt_bam = pysam.AlignmentFile(principaldirectory+gnumber+'.bam', "wb", template=bam)
    bad_bam = pysam.AlignmentFile(principaldirectory+gnumber+'.out.bam', "wb", template=bam)
    for read in bam.fetch():

        if read.has_tag('AS') and read.get_tag('AS') != 0: # Because when 'AS" = 0 (unmapped) the length query return None and therefore I cannot calculate minAS
            AS = read.get_tag('AS')
            length_query = read.infer_query_length()
            minAS=(0.93*length_query)-(length_query*4*0.07)
            if AS > minAS:
                filt_bam.write(read)
            elif AS <= minAS:
                bad_bam.write(read)
        else:
            bad_bam.write(read)



    bam.close()
    bad_bam.close()

def Run_MarkDuplicates(input,principaldirectory,log):
    tf_metrics = tempfile.NamedTemporaryFile(prefix='metrics_')
    MarkDuplicates_command = 'java -jar {picard} MarkDuplicates QUIET=true VERBOSITY=ERROR INPUT={tempbam} OUTPUT={bam} METRICS_FILE={metrics} ASSUME_SORTED=true'.format(picard=PICARD,tempbam=input,
                                                                                                                                                                          bam=principaldirectory+'random_bam_'+''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6))+'.bam',metrics=tf_metrics.name) # TEMPDIR
    subprocess.call(MarkDuplicates_command,shell=True) # Mark duplicates
    log.append(('MD_IN',MarkDuplicates_command))

def Index_BAM(input,log,name_of_index):
    samtools_index_command = 'samtools index {bam}'.format(bam=input)
    subprocess.call(samtools_index_command,shell=True)
    log.append((name_of_index,samtools_index_command))

def GATK_RealignerTargetCreator(log,input,number_threads,reference_genome):
    global GATK_intervals
    tf_GATKintervals = tempfile.NamedTemporaryFile(prefix='intervals_')
    RealignerTargetCreator_command = 'java -jar {GATK} -T RealignerTargetCreator -l ERROR -nt {threads} -R {ref} -o {intervals} -I {input_bam}'.format(GATK = GATK,threads =number_threads,ref=reference_genome,intervals=tf_GATKintervals.name+'.intervals',input_bam=input)
    subprocess.call(RealignerTargetCreator_command,shell=True)
    log.append(('RTC_IN',RealignerTargetCreator_command))
    GATK_intervals = tf_GATKintervals.name+'.intervals'

def GATK_IndelRealigner(log,input,output,reference_genome,interval_file):
    IndelRealigner_command = 'java -jar {GATK} --disable_bam_indexing -T IndelRealigner -l ERROR -R {ref} -targetIntervals {intervals} -I {input_bam} -o {final_bam}'.format(GATK = GATK,ref=reference_genome,intervals=interval_file,input_bam=input,final_bam=output)
    subprocess.call(IndelRealigner_command,shell=True)

    log.append(('IDR_IN',IndelRealigner_command))

def Samtools_Flagstat(input_bam,log,contaminant_threshold,gnb):

    tf_flagstat = tempfile.NamedTemporaryFile()
    samtools_flagstat_command = 'samtools flagstat {bam} > {log}'.format(bam=input_bam,log=tf_flagstat.name+'_flagstat.log')
    subprocess.call(samtools_flagstat_command,shell=True)
    log.append(('SFS_IN',samtools_flagstat_command))

    PostBamReadCount = ''
    PostBamDuplicate = ''
    MappingPercentage = ''
    ProperlyPairedReads = 'NA'

    with open(tf_flagstat.name + '_flagstat.log', 'r') as log_flagstat:
        for stat_line in log_flagstat:
            if 'QC-passed reads' in stat_line:
                PostBamReadCount = stat_line.split()[0]
            if 'duplicates' in stat_line:
                PostBamDuplicate = stat_line.split()[0]
            if 'properly paired' in stat_line:
                tmp = stat_line.split()
                ProperlyPairedReads = tmp[-1].split(':')[0][1:]
            if 'mapped' in stat_line:
                if '%' in stat_line:
                    tmp = stat_line.split()
                    MappingPercentage = tmp[-1].split(':')[0][1:]
                    if MappingPercentage[:-1] <= contaminant_threshold:
                        contaminant_file = open('/storage/PGO/PGOgit/prueba/genomes.contaminants', 'a')
                        contaminant_file.write(gnb + '\t' + MappingPercentage + '\n')
                        contaminant_file.close()

    log.append(('RDC_OUT', PostBamReadCount))
    log.append(('DUP_OUT', PostBamDuplicate))
    log.append(('MAP_OUT', MappingPercentage))
    log.append(('PP_OUT', ProperlyPairedReads))

def Qualimap_BamQC(input_bam,output_dir,gnumber,log,number_threads):
    global mean_Coverage_without_duplicates

    

    #cambiado: le ponemos la ruta al qualimap de garnatxa pq si no no lo detecta
    command_qualimabamqc = '/storage/PGO/bin/instalaciones/qualimap-build-13-01-16/qualimap bamqc -nt {threads} -bam {bamfinal} -sd -sdmode 1 -outdir {outdir} -outfile {outfile}'.format(
        threads =number_threads,bamfinal=input_bam, outdir=output_dir, outfile=gnumber)


    '''command_qualimabamqc = 'qualimap bamqc -nt {threads} -bam {bamfinal} -sd -sdmode 1 -outdir {outdir} -outfile {outfile}'.format(
        threads =number_threads,bamfinal=input_bam, outdir=output_dir, outfile=gnumber)'''
    process_BAMQC =subprocess.call(command_qualimabamqc, shell=True,stdout=subprocess.PIPE)

    
    log.append(('BQC_IN', command_qualimabamqc))

    bamqc_genome_results_file = open(output_dir + 'genome_results.txt', 'r')
    bamqc_genome_results_table = []

    for lines in bamqc_genome_results_file:
        bamqc_genome_results_table += [lines.strip()]

    mean_Coverage_without_duplicates = ''
    std_coverage_without_duplicates = ''
    percentage_coverage_not_covered = ''
    percentage_coverage_below_seven = ''

    for i in range(len(bamqc_genome_results_table)):
        if bamqc_genome_results_table[i] == '>>>>>>> Coverage':
            mean_Coverage_without_duplicates = float(
                bamqc_genome_results_table[i + 2][bamqc_genome_results_table[i + 2].index('=') + 1:-1].replace(
                    ',', '.'))
            std_coverage_without_duplicates = float(
                bamqc_genome_results_table[i + 3][bamqc_genome_results_table[i + 3].index('=') + 1:-1].replace(
                    ',', '.'))
            percentage_coverage_not_covered = 100 - float(
                bamqc_genome_results_table[i + 5][
                bamqc_genome_results_table[i + 5].index('a ') + 1:bamqc_genome_results_table[i + 5].index(
                    '%')].replace(',', '.'))
            percentage_coverage_below_seven = 100 - float(
                bamqc_genome_results_table[i + 10][
                bamqc_genome_results_table[i + 10].index('a ') + 1:bamqc_genome_results_table[i + 10].index(
                    '%')].replace(',', '.'))


    log.append(('COV_OUT', str(mean_Coverage_without_duplicates)))
    log.append(('STD_OUT', str(std_coverage_without_duplicates)))
    log.append(('NC_OUT', str(round(percentage_coverage_not_covered, 2)) + '%'))
    log.append(('N7_OUT', str(round(percentage_coverage_below_seven, 2)) + '%'))

    '''detect the mapping % and if below 80%, record to a file called: potential_contaminations.txt'''

    subprocess.call('rm {bamqcPDF}'.format(bamqcPDF=output_dir + gnumber + '.pdf'), shell=True)
    subprocess.call('rm {bamqcTXT}'.format(bamqcTXT=output_dir + 'genome_results.txt'), shell=True)

def Variant_calling(exclude_positions,git,min_mapping_quality,max_cov,min_cov,min_base_qual,min_reads2,hetero_freq,homo_freq,reference,Gnumber,output_directory,log):
    tf_varscan_log = tempfile.NamedTemporaryFile()
    tf_varscan_vcf = tempfile.NamedTemporaryFile()

    command_samtools_mpileup = "samtools mpileup -ABQ0 -d {maxcoverage} -q {mapqual} -f {reference} {finalbam} ".format(
        maxcoverage=max_cov,
        mapqual=min_mapping_quality, reference=reference,
        finalbam=output_directory + Gnumber + '.bam')
    command_varscan = "java -Xmx5000m -jar {varscan} mpileup2cns --min-coverage {coverage} --min-reads2 {reads2}" \
                      " --min-avg-qual {basequal} --min-var-freq {variantfreqhet} --min-freq-for-hom {variantfreqhomo} " \
                      "--strand-filter {strandfilter} --output-vcf 1 1> {vcf} 2>> {log}".format(
        varscan=VARSCAN,
        coverage=min_cov, reads2=min_reads2,
        basequal=min_base_qual, variantfreqhet=hetero_freq,
        variantfreqhomo=homo_freq, strandfilter=1,
        vcf=tf_varscan_vcf.name+'_VCF', log=tf_varscan_log.name + '_varscan.log')

    p1_samtools_mpileup = subprocess.Popen(command_samtools_mpileup, shell=True,
                                           stdout=subprocess.PIPE)
    p2_varscan = subprocess.Popen(command_varscan, shell=True,
                                  stdin=p1_samtools_mpileup.stdout)
    p1_samtools_mpileup.stdout.close()
    p2_varscan.communicate()
    log.append(('VC_IN', command_samtools_mpileup + '|' + command_varscan))
    bases_pileup_file = str
    bases_failed = str
    number_variants = str

    with open(tf_varscan_log.name + '_varscan.log', 'r') as varscan_log:
        for log_row in varscan_log:
            if 'bases in pileup file' in log_row:
                bases_pileup_file = log_row.split()
            if 'failed by the strand-filter' in log_row:
                bases_failed = log_row.split()
            if 'variant positions reported' in log_row:
                number_variants = log_row.split()

    number_variants = number_variants[0] + ' ' + number_variants[-4] + '+' + number_variants[-2] + ')'
    log.append(('BPF_OUT', bases_pileup_file[0]))
    log.append(('VFL_OUT', bases_failed[0]))
    log.append(('VAR_OUT', number_variants))


    command_snpEFF = 'java -jar /storage/PGO/bin/snpEff.jar ann -c /storage/PGO/store/db_pipe/snpEff_annotation_prueba/snpEff.config -noStats ' \
                                                             '-no-downstream -no-upstream MTB_ANC ' \
                                                             '{inputVCF} > {outputVCF}'.format(inputVCF=tf_varscan_vcf.name+'_VCF', outputVCF=output_directory + Gnumber + '.all.pos.vcf')


    log.append(('ANN_IN',command_snpEFF))
    subprocess.call(command_snpEFF,shell=True)
    command_replace_sample1 = "sed -i 's/FORMAT[[:blank:]]\+Sample1/FORMAT\t{gnumber}/'  {outputVCF}".format(gnumber=Gnumber,outputVCF=output_directory + Gnumber + '.all.pos.vcf')
    subprocess.call(command_replace_sample1,shell=True)

    vcf_file = open(output_directory+ Gnumber + '.all.pos.vcf', 'r')
    vcf_snp = open(output_directory+ Gnumber + '.var.snp.vcf', 'w')
    vcf_indel = open(output_directory+ Gnumber + '.var.homo.indel.vcf', 'w')

    snp_count = 0
    indel_count = 0

    for rows in vcf_file:
        if rows[0] == '#':
            vcf_snp.write(rows)
            vcf_indel.write(rows)
        else:
            if rows.strip().split()[4] != '.' and len(rows.strip().split()[3]) == 1 and len(
                    rows.strip().split()[4]) == 1 and rows.strip().split()[
                6] == 'PASS':  # if 4th columns of vcf is different from a dot and there is only one base (length =1) then it is a SNP
                vcf_snp.write(rows)
                snp_count += 1
            elif (rows.strip().split()[4] != '.' and len(rows.strip().split()[3]) > 1 and rows.strip().split()[6] == 'PASS') or (
                        rows.strip().split()[4] != '.' and len(rows.strip().split()[4]) > 1 and rows.strip().split()[
                6] == 'PASS'):  # indel
                if rows.strip().split()[7].split(';')[3] == 'HOM=1':
                    vcf_indel.write(rows)
                    indel_count += 1
    vcf_snp.close()

    log.append(('SNP_OUT', str(snp_count)))
    log.append(('HOMO_IND_OUT', str(indel_count)))

    snp_file = open(output_directory + Gnumber + '.var.snp.vcf', 'r')
    vcf_homo_snp = open(output_directory + Gnumber + '.var.homo.SNPs.vcf', 'w')
    vcf_het_snp = open(output_directory + Gnumber + '.var.het.SNPs.vcf', 'w')
    '''open the file containing the MTB positions to exclude from the analysis'''
    pos_exclude = []
    with open(git+exclude_positions,'r') as positions_to_exclude_file:
        table = [position.strip().split('\t') for position in positions_to_exclude_file]

    for i in range(1,len(table)):
        StartPosition = int(table[i][1])
        EndPosition = int(table[i][2])
        pos_exclude += [[StartPosition,EndPosition]]

    positions_to_exclude_templist = []
    positions_to_exclude = {}

    for coord in pos_exclude:
        positions_to_exclude_templist += range(coord[0], coord[1])
    for coordinates in positions_to_exclude_templist:
        positions_to_exclude[coordinates] = ''

    vcf_table = [i.split() for i in snp_file]

    homo_count = 0
    het_count = 0

    for line in vcf_table:
        if line[0][0] == '#':
            vcf_homo_snp.write("\t".join(line) + '\n')
            vcf_het_snp.write("\t".join(line) + '\n')
        else:
            if int(line[1]) not in positions_to_exclude:
                if 'HOM=1' in line[7]:
                    homo_count += 1
                    vcf_homo_snp.write("\t".join(line) + '\n')
                elif 'HET=1' in line[7]:
                    het_count += 1
                    vcf_het_snp.write("\t".join(line) + '\n')
    log.append(('HOM_OUT', str(homo_count)))
    log.append(('HET_OUT', str(het_count)))

    snp_file.close()
    vcf_indel.close()
    vcf_het_snp.close()
    vcf_homo_snp.close()
    vcf_file.close()

def lineage(n,vcf,git): # new function
    with open(git+'Lineage_SNPs_with_snps_anim.txt','r') as lineage_file: # open the file where the phylogenetic SNPs are defined
        for line in lineage_file: # Loop into the lines of the file
            line = line.strip()
            line = line.split()

            Lineage = str(line[0]) # Lineage is defined in the first column of the file
            Position = str(line[1]) # Lineage is defined in the second column of the file
            Nucleotide = str(line[2]) # Lineage is defined in the third column of the file

            Data[Lineage][Position] = Nucleotide # Implement the autovivification class

    X = 0
    with open(vcf,'r') as VCF: # open the VCF

        for row in VCF: # Loop into the lines of the VCF

            if row[0] != '#': # Ignore the header of the VCF

                row = row.strip()
                row = row.split()

                POS = str(row[1]) # position is the second column of the VCF
                REF = str(row[3]) # reference is the fourth column of the VCF
                ALT = str(row[4]) # ALT is the fifth column of the VCF

                TUPLE = (POS,ALT) # Define a tuple variable which for each line of the VCF takes the value of (position,alternative base)

                if TUPLE in Data.__getitem__(n).items(): #
                    #print n
                    X+=1

        if X == 3:
            return n, X
        elif X == 2:
            return n, X
        elif X == 1:
            return n, X

def Find_lineage(vcf_file,git,log):
    result_lineage = str
    LineageNames = []
    with open(git+'Lineage_SNPs_with_snps_anim.txt','r') as lineage_file:
        for line in lineage_file:
            line = line.strip()
            line = line.split()
            LineageName = str(line[0])
            LineageNames += [LineageName]
    LineageNames = list(set(LineageNames))
    for name in LineageNames:
        result_lineage = lineage(n=name,vcf=vcf_file,git=git)
        #print result_lineage
        if result_lineage != None:
            LineageInformation = result_lineage[0]
            NumberOfPositions = str(result_lineage[1])
            log.append(('LIN_OUT',LineageInformation+' ('+ NumberOfPositions +' positions)'))

        #I  have to read the first column of the lineage text file into a list, then loop in the list and apply the function on 'i'. result_lineage = Lineage.lineage('L1',principal_directory + G_number + vcf_extension)

def WGSFasta(excluded,vcf,dir,Gnumber):
    '''Make a dictionary of the positions to exclude'''
    pos_exclude = []
    positions_to_exclude_file = open(excluded,'r')

    table = [position.strip().split('\t') for position in positions_to_exclude_file]
    for i in range(1,len(table)):
        pos_exclude += [[int(table[i][1]),int(table[i][2])]]

    positions_to_exclude_file.close()

    positions_to_exclude_templist = []
    positions_to_exclude = {}

    for coord in pos_exclude:
        positions_to_exclude_templist += range(coord[0], coord[1])
    for coordinates in positions_to_exclude_templist:
        positions_to_exclude[coordinates] = ''


    '''Start making the Whole-Genome Fasta file'''
    TableVCF = []
    fasta_list = []
    with open(vcf,'r') as VCF_file:
        for row in VCF_file:
            if row[0] != '#':
                row = row.strip()
                row = row.split()
                TableVCF += [row]
        pos = 1
        for l in range(len(TableVCF)):
            #fasta_list = []
            # For example, l can take the value of: ['MTB_anc', '1', '.', 'T', '.', '.', 'PASS', 'ADP=18;WT=1;HET=0;HOM=0;NC=0', 'GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR', '0/0:32:18:18:18:0:0%:1E0:48:0:7:11:0:0']

            #en TableVCF
            #uno q sale mal ['MTB_anc', '4411525', '.', 'A', '.', '.', 'PASS', 'ADP=1;WT=0;HET=0;HOM=0;NC=1', 'GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR', './.:.:1'], ['MTB_anc', '4411526', '.', 'T', '.', '.', 'PASS', 'ADP=1;WT=0;HET=0;HOM=0;NC=1', 'GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR', './.:.:1']]
            
            #para G47304, que sale mal, se mete en la 5,6,7,8 y 14. El q salia bien entraba siempre en la 2. Y su tablevc tiene ['MTB_anc', '3744', '.', 'G', '.', '.', 'PASS', 'ADP=1;WT=0;HET=0;HOM=0;NC=1', 'GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR', './.:.:1']

            #uno q sale bie ['MTB_anc', '883764', '.', 'T', '.', '.', 'PASS', 'ADP=106;WT=1;HET=0;HOM=0;NC=0', 'GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR', '0/0:182:114:106:105:1:0,94%:5E-1:34:32:54:51:0:1'], ['MTB_anc', '883765', '.', 'A', '.', '.', 'PASS', 'ADP=105;WT=1;HET=0;HOM=0;NC=0', 'GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR', '0/0:198:112:105:105:0:0%:1E0:34:0:56:49:0:0'], 
            Position = TableVCF[l][1]
            RefBase = TableVCF[l][3]
            AltBase = TableVCF[l][4]
            QUAL = TableVCF[l][6]
            Info = TableVCF[l][7].split(';')
            Homo = Info[3] # HOM=0 or HOM=1
            Het = Info[2] # HET=0 or HET=1
            WT = Info[1] # WT=0 or WT=1
            NC = Info[4]
            SAMPLE = TableVCF[l][9]
            a = int(Position)
            b = pos

            if a == b:
                if a in positions_to_exclude or QUAL == 'str10':
                    fasta_list += ['X']
                else:
                    LengthAlt = len(AltBase)
                    LengthRef = len(RefBase)

                    if AltBase != '.' and LengthAlt == 1 and LengthRef > 1: # DELETION
                        #fasta_list += ['-']
                        fasta_list += [AltBase] # which actually is the ANC

                    elif AltBase != '.' and LengthAlt > 1 and LengthRef == 1:  # INSERTION
                        fasta_list += [AltBase[0]]

                    elif LengthRef == 1 and LengthAlt == 1: # Either SNP either ANCESTRAL base

                        if AltBase == '.' and WT == 'WT=1' and './.:.:' not in SAMPLE:  # Ancestral Base
                            fasta_list += [RefBase]
                            #se cumple esta en los que van
                        elif './.:.:' in SAMPLE and NC == 'NC=0': # if the pattern ./.:.: occurs in the vcf line
                            if len(TableVCF[l-1][3]) == 1 and './.:.:' not in TableVCF[l-1][9]: # look at position before and if it is not a deletion
                                fasta_list += ['X']
                            else:
                                fasta_list += ['-']
                                #cuando le bajas mucho el coverage se cumple esta

                        elif AltBase == '.' and NC == 'NC=1':  # Not covered (coverage below threshold)
                            fasta_list += ['X']
                            #se cumple esta en los q no van
                        elif AltBase != '.' and Homo == 'HOM=1':  # Homozygous SNP
                            fasta_list += [AltBase]

                        elif AltBase != '.' and Het == 'HET=1':  # Heterozygous SNP

                            if RefBase == 'A' and AltBase == 'G':
                                fasta_list += ['R']
                            elif RefBase == 'G' and AltBase == 'A':
                                fasta_list += ['R']
                            elif RefBase == 'C' and AltBase == 'T':
                                fasta_list += ['Y']
                            elif RefBase == 'T' and AltBase == 'C':
                                fasta_list += ['Y']
                            elif RefBase == 'G' and AltBase == 'C':
                                fasta_list += ['S']
                            elif RefBase == 'C' and AltBase == 'G':
                                fasta_list += ['S']
                            elif RefBase == 'A' and AltBase == 'T':
                                fasta_list += ['W']
                            elif RefBase == 'T' and AltBase == 'A':
                                fasta_list += ['W']
                            elif RefBase == 'G' and AltBase == 'T':
                                fasta_list += ['K']
                            elif RefBase == 'T' and AltBase == 'G':
                                fasta_list += ['K']
                            elif RefBase == 'C' and AltBase == 'A':
                                fasta_list += ['M']
                            elif RefBase == 'A' and AltBase == 'C':
                                fasta_list += ['M']

                    else:
                        print Position, ':what to do?'

            elif a != b:
                pos = int(a)
                fasta_list += list('-' * (pos - b))
                if a in positions_to_exclude or QUAL == 'str10':
                    fasta_list += ['X']
                else:
                    LengthAlt = len(AltBase)
                    LengthRef = len(RefBase)

                    if AltBase != '.' and LengthAlt == 1 and LengthRef > 1: # DELETION
                        #fasta_list += ['-']
                        fasta_list += [AltBase]
                    elif AltBase != '.' and LengthAlt > 1 and LengthRef == 1:  # INSERTION
                        fasta_list += [AltBase[0]]
                    elif LengthRef == 1 and LengthAlt == 1: # Either SNP either ANCESTRAL base
                        if AltBase == '.' and WT == 'WT=1' and './.:.:' not in SAMPLE:  # Ancestral Base
                            fasta_list += [RefBase]
                        elif './.:.:' in SAMPLE and NC == 'NC=0': # if the pattern ./.:.: occurs in the vcf line

                            if len(TableVCF[l-1][3]) == 1 and './.:.:' not in TableVCF[l-1][9]: # look at position before and if it is not a deletion
                                fasta_list += ['X']
                            else:
                                fasta_list += ['-']

                        elif AltBase == '.' and NC == 'NC=1':  # Not covered (coverage below threshold)
                            fasta_list += ['-']

                        elif AltBase != '.' and Homo == 'HOM=1':  # Homozygous SNP
                            fasta_list += [AltBase]
                            
                        elif AltBase != '.' and Het == 'HET=1':  # Heterozygous SNP

                            if RefBase == 'A' and AltBase == 'G':
                                fasta_list += ['R']
                            elif RefBase == 'G' and AltBase == 'A':
                                fasta_list += ['R']
                            elif RefBase == 'C' and AltBase == 'T':
                                fasta_list += ['Y']
                            elif RefBase == 'T' and AltBase == 'C':
                                fasta_list += ['Y']
                            elif RefBase == 'G' and AltBase == 'C':
                                fasta_list += ['S']
                            elif RefBase == 'C' and AltBase == 'G':
                                fasta_list += ['S']
                            elif RefBase == 'A' and AltBase == 'T':
                                fasta_list += ['W']
                            elif RefBase == 'T' and AltBase == 'A':
                                fasta_list += ['W']
                            elif RefBase == 'G' and AltBase == 'T':
                                fasta_list += ['K']
                            elif RefBase == 'T' and AltBase == 'G':
                                fasta_list += ['K']
                            elif RefBase == 'C' and AltBase == 'A':
                                fasta_list += ['M']
                            elif RefBase == 'A' and AltBase == 'C':
                                fasta_list += ['M']

                    else:
                        print Position, ':what to do?'

            pos += 1

    extension_list = []
    if len(fasta_list) < 2117144:
        extension_list += list('-' * (2117144 - len(fasta_list)))

    output_file = open(dir+Gnumber+'.fasta','w')
    output_file.write(">"+Gnumber+"\n")

    fasta_list.extend(extension_list)

    for base in fasta_list:
            output_file.write(base)
    output_file.close()

def Mapping_User_log(gnumber):
    if os.path.exists("/scicore/home/gagneux/GROUP/tbresearch/genomes/IN_PROGRESS/common_mappings/PipelineTB/database_update/Pipeline_mapping_user_log.tsv"):
        header_exists = True
    else:
        header_exists= False

    headers = 'G_NUMBER\tMAPPED_BY\tDATE\tPIPELINE_VERSION\n'
    with open("/scicore/home/gagneux/GROUP/tbresearch/genomes/IN_PROGRESS/common_mappings/PipelineTB/database_update/Pipeline_mapping_user_log.tsv","a") as tsvfile:
        if not header_exists:
            tsvfile.write(headers)
        tsvfile.write(gnumber+'\t'+getpass.getuser()+'\t'+str(datetime.date.today())+'\tv1\n')









