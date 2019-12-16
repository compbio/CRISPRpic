#!/usr/bin/env python

##############################################
# Author: HoJoon Lee
# Contact: hojoon@stanford.edu
# Date: May 28, 2018
# "Copyright (C) 2018 HoJoon Lee"
##############################################

### LOAD THE NECESSARY PACKAGES ###
import datetime, argparse, re, os, sys, glob, gzip, pandas as pd
from itertools import islice
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main(args):
    input_file = args.input
    window_size = args.window
    seq_file = args.file
    index_size = args.indexsize

    now = datetime.datetime.now()

    sys.stderr.write("Inital Index Size :" + str(index_size) + "\n")

    log_file = open("CRISPRpic_" + input_file.split(".")[0] + ".log","w")
    log_file.write("Date: " + now.strftime("%Y-%m-%d %H:%M") + "\n")
    log_file.write("Command line: python CRISPRpic.py -i " + input_file + " -f " + seq_file + " -w " + str(window_size) +  " -s " + str(index_size) + "\n")
    log_file.close()

    enz_info(input_file,window_size,index_size,log_file)    ## check if gRNAseq is unique in amplicon and determine the lenght of indicies
    seq_freq(seq_file)                  ## count reads by sequences

    df_list = pd.read_table("enzyme_info.txt")
    for i, row in df_list.iterrows():
        locus = row["locus"]
        enz_type = row["enz_type"]
        amp_seq = row["upseq"] + row["downseq"]
        if os.path.exists(locus):
            sys.stderr.write(locus +"\n")
            get_amplicon(locus,amp_seq)                      ## select amplicons by either of two adaptors and both of first and last 8 nucleotides
            make_mutseq(locus,window_size,row["upseq"],row["downseq"],enz_type)            ## make possible sequences with deletions across break point and one nt ins/del/sub
            mut_freq(locus,row["upseq"],row["downseq"],window_size,row["motif_len_up"],row["motif_len_down"],enz_type)               ## count the mutation events by Del_sequences
            freq_summary(locus,row["upseq"],row["downseq"])                      ## summarize the frequencies of mutation types, indels and their size
            plot_freq(locus)                         ## making bar charts
        else:
            sys.stderr.write(locus + "\t" + "amplicon is very low complex" + "\n")

def parse_commandline():
    parser=argparse.ArgumentParser()
    parser.add_argument('-i','--input', help='input enzyme list',required=True)
    parser.add_argument('-f','--file', help='input sequencing file',required=True)
    parser.add_argument('-w','--window', type = int, help='window size',required=True)
    parser.add_argument('-s','--indexsize', type = int, default = 8, help='index initial size',required=False)
    args=parser.parse_args()
    return args

def enz_info(input_file,window_size,index_size,log_file):
    sys.stderr.write("Checking reference sequence and gRNAseqs\n")
    fout = open("enzyme_info.txt","w")
    fout.write('locus' + "\t" + 'enz_type' + "\t" + 'gRNAseq' + "\t" + 'upseq' + "\t" + 'downseq' + "\t" + 'motif_len_up' + "\t" + 'motif_len_down' + "\n")
    with open (input_file) as f:
        for lines in f:
            line = lines.rstrip().split("\t")
            amp_seq = ''
            sys.stderr.write(line[0]+"\n")
            enz_type = int(line[3])
            if enz_type == 3:
                amp_seq = line[1]
            else:
                if re.search(line[2].upper(),line[1]):
                    amp_seq = line[1]
                else:
                    amp_seq = ReverseComplement(line[1])
            split_seq = amp_seq.split(line[2].upper())

            if len(split_seq) > 2:
                sys.stderr.write("gRNAseq site is not unique"+"\n")
            else:
                ## Remove adaptor sequences from amplicon
                enzyme = line[0]                       ### enzyme name
                adpt5 = amp_seq[0:20]
                adpt3 = amp_seq[-20:]

                ## get upstream seq before break & downstream seq after break
                upseq = ''
                downseq = ''
                if enz_type == 1:   # SpCas9
                    upseq = split_seq[0] + line[2].upper()[:-6]
                    downseq = line[2].upper()[-6:] + split_seq[1]
                elif enz_type == 2:  # AsCpf1
                    upseq = split_seq[0] + line[2].upper()[0:25]
                    downseq = line[2].upper()[25:] + split_seq[1]
                elif enz_type == 3: # custom enzyme
                    break_pos = int(line[2])      ## from 5' end of amplicon
                    upseq = amp_seq[0:break_pos]
                    downseq = amp_seq[break_pos:]
                else:
                    sys.stderr.write("Error, please indicate the type of enzyme"+"\n")
                    sys.exit()

                target1 = upseq.split(adpt5)[1]     # upstream seq from break point without 5' adaptor seq
                target2 = downseq.split(adpt3)[0]   # downstream seq from break point without 3' adaptor seq
                wt = target1 + target2
                len1 = len(target1)
                len2 = len(target2)

                ## identify the length of index where all indicies within windowsize are unique in the amplicon staring from 8
                motif_len_up = index_size                       ### stating length is 8
                motif_len_down = index_size
                up_flag = 0
                while up_flag==0:
                    up_check = 0
                    for i in range(0,window_size):
                        firstUpmotif = ''
                        if i ==0:
                            firstUpmotif = target1[-(motif_len_up):]
                        else:
                            firstUpmotif = target1[-(motif_len_up+i):-i]
                        up_test = wt.split(firstUpmotif)
                        if len(up_test) != 2:
                            up_check += 1
                    if up_check == 0:
                        up_flag +=1
                    else:
                        motif_len_up += 1

                down_flag = 0
                while down_flag==0:
                    down_check = 0
                    for i in range(0,window_size):
                        firstDownmotif = target2[i:(motif_len_down+i)]
                        down_test = wt.split(firstDownmotif)
                        if len(down_test) > 2:
                            down_check+=1
                    if down_check == 0:
                        down_flag +=1
                    else:
                        motif_len_down += 1

                if motif_len_up > 15 or motif_len_down > 15:
                    sys.stderr.write(sys.stderr,"very low complex sequences"+"\n")

                if not os.path.exists(enzyme):
                    os.makedirs(enzyme)
                fout.write(enzyme + "\t" + str(enz_type) + "\t" + line[2] + "\t" + upseq + "\t" + downseq + "\t" + str(motif_len_up) + "\t" + str(motif_len_down) + "\n")
    fout.close()

def seq_freq(seq_file):
    sys.stderr.write("Counting reads"+"\n")
    seq_freq = {}
    for fastqs in glob.glob("*" + seq_file):
        n=0
        #with gzip.open(fastqs,'r') as f:
        with open(fastqs,'r') as f:
            while True:
                lines = list(islice(f,4))
                if not lines:
                    break
                seq = lines[1].rstrip()
                if seq in seq_freq:
                    seq_freq[seq] += 1
                else:
                    seq_freq.setdefault(seq,0)
                    seq_freq[seq] += 1
                n += 1

    fout = open("seq_freq.txt","w")
    fout.write("seq\tseq_len\tfreq\n")
    for seqs in seq_freq:
        length = len(seqs)
        if length > 0:
            fout.write(seqs+"\t"+ str(length) + "\t" + str(seq_freq[seqs])+"\n")
    fout.close()

def get_amplicon(locus,amplicon):
    sys.stderr.write("selecting amplicons"+"\n")
    adapt5 = amplicon[0:20]
    first8 = amplicon[20:28]
    adapt3 = amplicon[-20:]
    last8 = amplicon[-28:-20]

    os.mkdir(locus + "/intermediate_files")
    seq_file = open(locus + "/intermediate_files/amplicon.txt","w")
    with open("seq_freq.txt") as f1:
        next(f1)
        for lines1 in f1:
            line1 = lines1.rstrip().split("\t")
            if re.search(first8,line1[0]) and re.search(last8,line1[0]):
                if re.search(adapt5,line1[0]) or re.search(adapt3,line1[0]):
                    amp = line1[0].split(first8)
                    first_loc = amp.pop(0)
                    remain_seq = first8.join(amp)
                    amp1 = remain_seq.split(last8)
                    last_loc = amp1.pop()
                    remain_seq1 = last8.join(amp1)
                    if len(first_loc) < 50 and len(last_loc) < 50:
                        seq_file.write(first8+remain_seq1+last8+"\t"+line1[1]+"\t"+line1[2]+"\t"+line1[0]+"\n")
            else:
                revcom = ReverseComplement(line1[0])
                if re.search(first8,revcom) and re.search(last8,revcom):
                    if re.search(adapt5,revcom) or re.search(adapt3,revcom):
                        amp = revcom.split(first8)
                        first_loc = amp.pop(0)
                        remain_seq = first8.join(amp)
                        amp1 = remain_seq.split(last8)
                        last_loc = amp1.pop()
                        remain_seq1 = last8.join(amp1)
                        if len(first_loc) < 50 and len(last_loc) < 50:
                            seq_file.write(first8+remain_seq1+last8 +"\t"+line1[1]+"\t"+line1[2]+"\t"+line1[0]+"\n")
        seq_file.close()

def make_mutseq(locus,window_size,upstr,downstr,enz_type):
    sys.stderr.write("making possible mutant sequences"+"\n")

    def generate_del_variants(upseq,downseq,dirs,window_size,locus):
        target1 = upseq[20:]
        target2 = downseq[:-20]

        wt =  target1 + target2
        fout = open(locus + "/intermediate_files/" + dirs + "_del_variants.fa","w")
        len1 = len(target1)
        len2 = len(target2)

        aligned_seq = {}
        mut_seq = {}
        if dirs != str(-window_size):
            for i in range(len1-1,-1,-1):
                before = target1[0:i]
                mut = before + target2
                aligned = before + ("-"*(len1-i)) + target2
                mut_type = "5DEL"+str(len1-i-int(dirs))+"_3DEL"+dirs
                if mut in mut_seq:
                    mut_seq[mut].append(mut_type)
                    aligned_seq[mut].add(aligned)
                else:
                    mut_seq.setdefault(mut,[])
                    mut_seq[mut].append(mut_type)
                    aligned_seq.setdefault(mut,set())
                    aligned_seq[mut].add(aligned)

        if dirs != str(window_size):
            for j in range(1,len2+1):
                after = target2[j:]
                mut = target1+after
                aligned = target1 + ("-"*j) + after
                mut_type = "5DEL"+str(-int(dirs))+"_"+"3DEL"+str(j+int(dirs))
                if mut in mut_seq:
                    mut_seq[mut].append(mut_type)
                    aligned_seq[mut].add(aligned)
                else:
                    mut_seq.setdefault(mut,[])
                    mut_seq[mut].append(mut_type)
                    aligned_seq.setdefault(mut,set())
                    aligned_seq[mut].add(aligned)

        for i in range(len1-1,-1,-1):
            before = target1[0:i]
            for j in range(1,len2+1):
                after = target2[j:]
                aligned = before + ("-"*(len1-i)) + ("-"*j) + after
                mut = before + after
                mut_type = "5DEL"+str(len1-i-int(dirs)) + "_3DEL"+str(j+int(dirs))
                if mut in mut_seq:
                    mut_seq[mut].append(mut_type)
                    aligned_seq[mut].add(aligned)
                else:
                    mut_seq.setdefault(mut,[])
                    mut_seq[mut].append(mut_type)
                    aligned_seq.setdefault(mut,set())
                    aligned_seq[mut].add(aligned)

        for seq in mut_seq:
            type_list = ",".join(mut_seq[seq])
            align_list = ",".join(aligned_seq[seq])
            fout.write(seq + "\t" + type_list + "\t" + align_list + "\t" + "Del" + "\n")
        fout.close()
    def onebase_change(upseq,downseq,window_size,locus):
        target1 = upseq[20:]
        target2 = downseq[:-20]
        wt = target1 + target2
        fout = open(locus + "/intermediate_files/onebase_sequences.txt","w")
        fout.write(wt+"\t"+"wt"+ "\t" + wt + "\t" + "WT" + "\n")
        nucls = ['A','C','G','T']
        aligned_seq = {}
        mut_seq = {}
        mut_finaltype = {}
        front_len = len(target1) - window_size -1
        back_len = len(target1) + window_size
        for i in range(0,len(wt)):
            mut = ''
            aligned= ''
            final_type = ''

            for nucl in nucls:
                if i > (front_len) and i < back_len:
                    final_type = "Sub"
                else:
                    final_type = "Unmodified"

                if wt[i:i+1] != nucl:
                    mut = wt[0:i] + nucl + wt[i+1:]
                    aligned = wt[0:i] + nucl.lower() + wt[i+1:]
                    mut_type = "OneSub|"+str(i+1)
                    if mut in mut_seq:
                        mut_seq[mut].append(mut_type)
                        aligned_seq[mut].add(aligned)
                        mut_finaltype[mut].add(final_type)
                    else:
                        mut_seq.setdefault(mut,[])
                        mut_seq[mut].append(mut_type)
                        aligned_seq.setdefault(mut,set())
                        aligned_seq[mut].add(aligned)
                        mut_finaltype.setdefault(mut,set())
                        mut_finaltype[mut].add(final_type)


                mut = wt[0:i] + nucl + wt[i:]
                aligned = wt[0:i] + nucl.lower() + wt[i:]
                mut_type = "OneIns|"+str(i)+"^"+str(i+1)

                if (i-1) > (front_len) and i < back_len:
                    final_type = "Ins"
                else:
                    final_type = "Unmodified"

                if mut in mut_seq:
                    mut_seq[mut].append(mut_type)
                    aligned_seq[mut].add(aligned)
                    mut_finaltype[mut].add(final_type)
                else:
                    mut_seq.setdefault(mut,[])
                    mut_seq[mut].append(mut_type)
                    aligned_seq.setdefault(mut,set())
                    aligned_seq[mut].add(aligned)
                    mut_finaltype.setdefault(mut,set())
                    mut_finaltype[mut].add(final_type)

            if i > (front_len) and i < back_len:
                final_type = "Del"
            else:
                final_type = "Unmodified"
                #mut = wt[0:i-1] + wt[i+1:]
                mut = wt[0:i] + wt[i+1:]
                aligned = wt[0:i] + "-" + wt[i+1:]
                mut_type = "OneDel|"+str(i+1)
                if mut in mut_seq:
                    mut_seq[mut].append(mut_type)
                    aligned_seq[mut].add(aligned)
                    mut_finaltype[mut].add(final_type)
                else:
                    mut_seq.setdefault(mut,[])
                    mut_seq[mut].append(mut_type)
                    aligned_seq.setdefault(mut,set())
                    aligned_seq[mut].add(aligned)
                    mut_finaltype.setdefault(mut,set())
                    mut_finaltype[mut].add(final_type)

        for nucl in nucls:
            mut_type = "OneIns|"+str(len(wt))+"^"
            mut = wt + nucl
            aligned =wt + nucl.lower()
            final_type = "Unmodified"

            if mut in mut_seq:
                mut_seq[mut].append(mut_type)
                aligned_seq[mut].add(aligned)
                mut_finaltype[mut].add(final_type)
            else:
                mut_seq.setdefault(mut,[])
                mut_seq[mut].append(mut_type)
                aligned_seq.setdefault(mut,set())
                aligned_seq[mut].add(aligned)
                mut_finaltype.setdefault(mut,set())
                mut_finaltype[mut].add(final_type)

        for seq in mut_seq:
            type_list = ",".join(mut_seq[seq])
            align_list = ",".join(aligned_seq[seq])
            finaltype_list = ",".join(mut_finaltype[seq])
            fout.write(seq + "\t" + type_list + "\t" + align_list + "\t" + finaltype_list +"\n")
        fout.close()

    if enz_type == 2:
        window_size += 2

    generate_del_variants(upstr,downstr,"0",window_size,locus)
    onebase_change(upstr,downstr,window_size,locus)

    for i in range(1,(window_size+1)):
        upseq = upstr[:-i]
        downseq = upstr[-i:]+downstr
        dirs = str(-i)
        generate_del_variants(upseq,downseq,dirs,window_size,locus)

        upseq = upstr + downstr[:i]
        downseq = downstr[i:]
        dirs = str(i)
        generate_del_variants(upseq,downseq,dirs,window_size,locus)

    del_seq = {}
    del_annot = {}
    for files in glob.glob(locus + "/intermediate_files/*_del_variants.fa"):
        with open(files) as f:
            for lines in f:
                line = lines.rstrip().split("\t")
                if line[3] == "Del":
                    if line[0] in del_seq:
                        aligns = line[2].split(",")
                        for align in aligns:
                            del_seq[line[0]].add(align)
                    else:
                        del_seq.setdefault(line[0],set())
                        aligns = line[2].split(",")
                        for align in aligns:
                            del_seq[line[0]].add(align)

                    if line[0] in del_annot:
                        aligns = line[1].split(",")
                        for align in aligns:
                            del_annot[line[0]].add(align)

                    else:
                        del_annot.setdefault(line[0],set())
                        aligns = line[1].split(",")
                        for align in aligns:
                            del_annot[line[0]].add(align)
    fout = open(locus + "/intermediate_files/Del_sequences.txt","w")
    for seqs in del_seq:
        del_num = len(del_seq[seqs])
        del_list = ",".join(del_seq[seqs])
        annot_list = ",".join(del_annot[seqs])
        fout.write(seqs + "\t" + del_list + "\t" + str(del_num) + "\t" + annot_list +"\n")
    fout.close()

def mut_freq(locus,upstr,downstr,window_size,motif_len_up,motif_len_down,enz_type):
    sys.stderr.write("counting the frequency of mutations"+"\n")
    if enz_type == 2:
        window_size += 2
    target1 = upstr[20:]
    target2 = downstr[:-20]
    wt = target1 + target2
    len_wt = len(wt)
    len1 = len(target1)
    len2 = len(target2)

    ### put possible deletion sequences into hash
    align_seqs = {}
    seq_type = {}
    call_type = {}
    with open(locus + "/intermediate_files/Del_sequences.txt") as f:
        for lines in f:
            line = lines.rstrip().split("\t")
            if line[0] in seq_type:
                seq_type[line[0]].add(line[3])
                align_seqs[line[0]].add(line[1])
                call_type[line[0]].add("Del")
            else:
                seq_type.setdefault(line[0],set())
                seq_type[line[0]].add(line[3])
                align_seqs.setdefault(line[0],set())
                align_seqs[line[0]].add(line[1])
                call_type.setdefault(line[0],set())
                call_type[line[0]].add("Del")
    ### put possible one nucleotide change into hash
    with open(locus + "/intermediate_files/onebase_sequences.txt") as f:
        for lines in f:
            line = lines.rstrip().split("\t")
            if line[0] in seq_type:
                seq_type[line[0]].add(line[1])
                align_seqs[line[0]].add(line[2])
                call_type[line[0]].add(line[3])
            else:
                seq_type.setdefault(line[0],set())
                seq_type[line[0]].add(line[1])
                align_seqs.setdefault(line[0],set())
                align_seqs[line[0]].add(line[2])
                call_type.setdefault(line[0],set())
                call_type[line[0]].add(line[3])

    ### making list of indicies with their uniquness from break point to both direction
    down_indexs = []
    down_index_uniq={}
    up_indexs= []
    up_index_uniq= {}
    for i in range(0,len(target2)-(motif_len_down-1)):
        peptide = target2[i:i+motif_len_down]
        down_indexs.append(peptide)
        if peptide in down_index_uniq:
            down_index_uniq[peptide] += 1
        else:
            down_index_uniq.setdefault(peptide,0)
            down_index_uniq[peptide] += 1

    for i in range(len(target1)-(motif_len_up),-1,-1):
        peptide = target1[i:i+motif_len_up]
        up_indexs.append(peptide)
        if peptide in up_index_uniq:
            up_index_uniq[peptide] += 1
        else:
            up_index_uniq.setdefault(peptide,0)
            up_index_uniq[peptide] += 1

    fout = open(locus + "/intermediate_files/" + locus + "_mut_freq.txt","w")
    fout.write("Aligned_seq" + "\t" + "Seq_len" + "\t" + "Counts" + "\t" + "Type" + "\t" + "Subtyp" + "\t" + "Source" + "\t" + "Original_seq" + "\t" + "Seq_input_fastq" + "\n")
    with open(locus + "/intermediate_files/amplicon.txt") as f:
        for lines in f:
            line = lines.rstrip().split("\t")
            type_list = ''

            if line[0] in seq_type:         ### detect mutations by indexed sequences
                type_list =  ",".join(seq_type[line[0]])
                align_list = ",".join(align_seqs[line[0]])
                call_list = ",".join(call_type[line[0]])
                if re.search("WT",call_list ):
                    fout.write(align_list + "\t" + line[1] + "\t" + line[2] + "\t" + "WT" + "\t" +  type_list + "\t" + "indexed" + "\t" + line[0] + "\t" + line[3] + "\n")
                elif re.search("Del",call_list ):
                    fout.write(align_list + "\t" + line[1] + "\t" + line[2] + "\t" + "Del" + "\t" +  type_list + "\t" + "indexed" + "\t" + line[0] + "\t" + line[3] + "\n")
                elif re.search("Ins",call_list ):
                    fout.write(align_list + "\t" + line[1] + "\t" + line[2] + "\t" + "Ins" + "\t" +  type_list + "\t" + "indexed" + "\t" + line[0] + "\t" + line[3] + "\n")
                elif re.search("Unmodified",call_list ):
                    fout.write(align_list + "\t" + line[1] + "\t" + line[2] + "\t" + "Unmodified" + "\t" +  type_list + "\t" + "indexed" + "\t" + line[0] + "\t" + line[3] + "\n")
                elif re.search("Sub",call_list ):
                    fout.write(align_list + "\t" + line[1] + "\t" + line[2] + "\t" + "Sub" + "\t" +  type_list + "\t" + "indexed" + "\t" + line[0] + "\t" + line[3] + "\n")
                else:
                    sys.stderr.write(call_list+"\n")
                    fout.write(align_list + "\t" + line[1] + "\t" + line[2] + "\t" + "WT_Un" + "\t" +  type_list + "\t" + "indexed" + "\t" + line[0] + "\t" + line[3] + "\n")
            else:
                upflag = 0
                ishift_count = 0
                for i in range(0,len(up_indexs)):
                    if re.search(up_indexs[i],line[0]):
                        split_up = line[0].split(up_indexs[i])
                        up_fraq = len(split_up)
                        if up_index_uniq[up_indexs[i]] == 1:
                            if up_fraq == 2:
                                upflag += 1
                                break
                        else:
                            ishift_count += 1
                if upflag == 0:
                    fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "NA" + "\t" +  type_list + "\t" + "NoUp" + "\t" + "NA" + "\n")
                else:
                    ishift_check=0
                    if ishift_count != 0:
                        if int(i - ishift_count) < window_size:
                            fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "NA" + "\t" +  type_list + "\t" + str(i)+","+str(ishift_count)+",iNA" + "\t" + "NA" + "\t" + line[3] + "\n")
                            ishift_check += 1

                    if ishift_check == 0:
                        down_flag=0
                        jshift_count = 0
                        for j in range(0,len(down_indexs)):
                            if re.search(down_indexs[j],split_up[1]):
                                split_down = split_up[1].split(down_indexs[j])
                                down_fraq = len(split_down)
                                if down_index_uniq[down_indexs[j]] == 1:
                                    if down_fraq == 2:
                                        down_flag += 1
                                        break
                                else:
                                    jshift_count += 1
                        if down_flag == 0:
                            fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "NA" + "\t" +  type_list + "\t" + "NoDown" + "\t" + "NA" + "\t" + line[3] + "\n")
                        else:
                            if i == 0 and j == 0:
                                if split_down[0] == "":
                                    fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "WT" + "\t" +  "wt" + "\t" + "UpDown" + "\t" + line[0] + "\t" + line[3] + "\n")
                                else:
                                    align_seq = split_up[0] + up_indexs[i] + split_down[0].lower() + down_indexs[j] + split_down[1]
                                    fout.write(align_seq + "\t" + line[1] + "\t" + line[2] + "\t" + "Ins" + "\t" +  split_down[0] + "\t" + str(i) + "," + str(j) + "\t" + line[0] + "\t" + line[3] + "\n")
                            else:
                                if split_down[0] == "":
                                    align_seq = split_up[0] + up_indexs[i] + "-"*(i+j) + down_indexs[j] + split_down[1]
                                    fout.write(align_seq + "\t" + line[1] + "\t" + line[2] + "\t" + "Del" + "\t" +  str(i) + "," + str(j) + "\t" + "UpDown" + "\t" + line[0] + "\t" + line[3] + "\n")
                                else:
                                    jshift_check=0
                                    if jshift_count != 0:
                                        if int(j - jshift_count) < window_size:
                                            fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "NA" + "\t" +  type_list + "\t" + str(j)+","+str(jshift_count)+",jNA" + "\t" + "NA" + "\t" + line[3] + "\n")
                                            jshift_check += 1
                                    if jshift_check == 0:
                                        i_index = i - ishift_count
                                        j_index = j - jshift_count
                                        remain_size = len(split_down[0]) - ishift_count - jshift_count
                                        if remain_size > (i_index+j_index):
                                            fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "Ins" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + split_down[0] + "\t" + line[0] + "\t" + line[3] + "\n")
                                        elif remain_size == (i_index+j_index):
                                            if i_index > window_size and j_index > window_size:
                                                ref_frag = target1[-(window_size):] + target2[:window_size]
                                                remain_frag = split_down[0][(i-window_size):(i-window_size+6)]
                                                if ref_frag == remain_frag:
                                                    fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "WT" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + ref_frag + "," + remain_frag + "\t" + line[0] + "\t" + line[3] + "\n")
                                                else:
                                                    fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "Sub" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + ref_frag + "," + remain_frag + "\t" + line[0] + "\t" + line[3] + "\n")
                                            elif i_index == 0:
                                                compare_size = min(remain_size,window_size)
                                                ref_frag = target2[:compare_size]
                                                remain_frag = split_down[0][:compare_size]
                                                if ref_frag == remain_frag:
                                                    fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "WT" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + "UpDown" + "\t" + line[0] + "\t" + line[3] + "\n")
                                                else:
                                                    fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "Sub" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + "UpDown" + "\t" + line[0] + "\t" + line[3] + "\n")
                                            elif j_index == 0:
                                                compare_size = min(remain_size,window_size)
                                                ref_frag = target1[-(compare_size):]
                                                remain_frag = split_down[0][-(compare_size):]
                                                if ref_frag == remain_frag:
                                                    fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "WT" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + ref_frag + "," + remain_frag  + "\t" + line[0] + "\t" + line[3] + "\n")
                                                else:
                                                    fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "Sub" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + ref_frag + "," + remain_frag + "\t" + line[0] + "\t" + line[3] + "\n")
                                            else:
                                                fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "Sub" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + "UpDown" + "\t" + line[0] + "\t" + line[3] + "\n")
                                        else:
                                            if i_index == 0:
                                                compare_size = min(remain_size,window_size)
                                                ref_frag = target2[:compare_size]
                                                remain_frag = split_down[0][:compare_size]
                                                if ref_frag == remain_frag:
                                                    if compare_size == window_size:
                                                        fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "WT" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + "UpDown" + "\t" + line[0] + "\t" + line[3] + "\n")
                                                    else:
                                                        fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "Del" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + "UpDown" + "\t" + line[0] + "\t" + line[3] + "\n")
                                                else:
                                                    fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "Complexs_Del" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + "UpDown" + "\t" + line[0] + "\t" + line[3] + "\n")
                                            elif j_index == 0:
                                                compare_size = min(remain_size,window_size)
                                                ref_frag = target1[-(compare_size):]
                                                remain_frag = split_down[0][-(compare_size):]
                                                if ref_frag == remain_frag:
                                                    if compare_size == window_size:
                                                        fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "WT" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + ref_frag + "," + remain_frag + "\t" + line[0] + "\t" + line[3] + "\n")
                                                    else:
                                                        fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "Del" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + "UpDown" + "\t" + line[0] + "\t" + line[3] + "\n")
                                                else:
                                                    fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "Complexs_Del" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + "UpDown" + "\t" + line[0] + "\t" + line[3] + "\n")
                                            else:
                                                fout.write(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + "Complexs_Del" + "\t" +  str(i) + "," + str(j) +"|"+ str(ishift_count) + "," + str(jshift_count) + "\t" + "UpDown" + "\t" + line[0] + "\t"  + line[3] + "\n")
    fout.close()

def freq_summary(locus,upstr,downstr):
    sys.stderr.write("summarizing frequencies"+"\n")

    total_summary = open("Total_sample_summary.txt","w")
    total_summary.write("locus" + "\t" + "Unmodified" + "\t" + "Del" + "\t" + "Ins" + "\t" + "Sub" + "\t" + "Unmodified" + "\t" + "Del" + "\t" + "Ins" + "\t" + "Sub" + "\t" + "complex_del" + "\t" + "complex_ins" + "\t" + "NA" + "\n")

    uplen = len(upstr) - 19
    downlen = len(downstr) - 19

    freq_table = open(locus + "/" + locus + "_mut_freq.txt","w")
    freq_table.write("Aligned_seq" + "\t" + "Micro_homology" + "\t" + "Seq_len" + "\t" + "Counts" + "\t" + "Type" + "\t" + "Indel_size" + "\n")
    total = 0
    index_count = {}
    groups = {}
    del_index = {}
    del_MH = {}
    with open(locus + "/intermediate_files/Del_sequences.txt") as f:
        for lines in f:
            line = lines.rstrip().split("\t")
            for del_type in line[3].split(","):
                del_index[del_type]=line[3]
                del_MH[del_type]=int(line[2])

    ins_size = {}
    del_size = {}
    del_mut_freq = {}
    with open(locus + "/intermediate_files/" + locus + "_mut_freq.txt") as f:
        next(f)
        for lines in f:
            line = lines.rstrip().split("\t")
            total += int(line[2])

            ### counting deletion / insertion and their size
            ins_bp = 0
            del_bp = 0
            if line[3] == 'Del':
                del_ind = ''
                if line[5] == "indexed":
                    del_ind = line[4].split(",")[0]
                    split = del_ind.split("_")
                    del_bp = 1
                    if len(split) == 2:
                        upDel = int(split[0].replace("5DEL",""))
                        downDel = int(split[1].replace("3DEL",""))
                        del_bp = upDel + downDel
                    else:
                        del_ind = line[4].split(",")[1]

                    if line[4] in del_mut_freq:
                        del_mut_freq[line[4]] += int(line[2])
                    else:
                        del_mut_freq.setdefault(line[4],0)
                        del_mut_freq[line[4]] += int(line[2])
                else:
                    pos_index = line[4].split("|")
                    main = pos_index[0].split(",")
                    index5 = int(main[0])
                    index3 = int(main[1])
                    if len(pos_index)>1:
                        additional = pos_index[1].split(",")
                        index5 -= int(additional[0])
                        index3 -= int(additional[1])
                    del_bp = index5 + index3
                    del_type = "5DEL"+str(index5)+"_"+"3DEL"+str(index3)
                    del_ind = del_type
                    del_seq = del_index[del_type]
                    if del_seq in del_mut_freq:
                        del_mut_freq[del_seq] += int(line[2])
                    else:
                        del_mut_freq.setdefault(del_seq,0)
                        del_mut_freq[del_seq] += int(line[2])
            elif line[3] == 'Ins':
                if line[5] == '0,0':
                    ins_bp = len(line[4])
                    if ins_bp in ins_size:
                        ins_size[ins_bp] += int(line[2])
                    else:
                        ins_size.setdefault(ins_bp,0)
                        ins_size[ins_bp] += int(line[2])
                elif line[5] == 'indexed':
                    ins_bp = 1
                else:
                    pos_index = line[4].split("|")
                    shift = pos_index[1].split(",")
                    ins_bp = len(line[5])-int(shift[0])-int(shift[1])

            ### making new summary table
            align_seqs = line[0].split(",")[0]
            freq_table.write(align_seqs + "\t")
            if line[3] == "Del":
                MH = del_MH[del_ind]-1
                freq_table.write(str(MH))
            else:
                freq_table.write('NA')
            freq_table.write("\t" + line[1] + "\t" + line[2] + "\t")
            if line[3] == "WT":
                freq_table.write("Unmodified")
            else:
                freq_table.write(line[3])
            freq_table.write("\t")

            if line[3] == "Ins":
                freq_table.write(str(ins_bp))
                if ins_bp in ins_size:
                    ins_size[ins_bp] += int(line[2])
                else:
                    ins_size.setdefault(ins_bp,0)
                    ins_size[ins_bp] += int(line[2])
            elif line[3] == "Del":
                freq_table.write(str(del_bp))
                if del_bp in del_size:
                    del_size[del_bp] += int(line[2])
                else:
                    del_size.setdefault(del_bp,0)
                    del_size[del_bp] += int(line[2])
            else:
                freq_table.write('NA')
            freq_table.write("\n")

            #### counting mutation types
            if line[5] == "indexed":
                if line[3] in index_count:
                    index_count[line[3]]+=int(line[2])
                else:
                    index_count.setdefault(line[3],0)
                    index_count[line[3]]+=int(line[2])
            else:
                if line[3] in groups:
                    groups[line[3]] += int(line[2])
                else:
                    groups.setdefault(line[3],0)
                    groups[line[3]] += int(line[2])
    freq_table.close()

    fout = open(locus + "/intermediate_files/ins_size_dist.txt", "w")
    for size in sorted(ins_size):
        fout.write(str(size) + "\t" + str(ins_size[size]) + "\n")
    fout.close()

    del_frame = {}
    fout = open(locus + "/intermediate_files/del_size_dist.txt", "w")
    for size in sorted(del_size):
        delf = size%3
        if delf in del_frame:
            del_frame[delf] += del_size[size]
        else:
            del_frame.setdefault(delf,0)
        fout.write(str(size) + "\t" + str(del_size[size]) + "\n")
    fout.close()

    fout = open(locus + "/intermediate_files/del_frame_count.txt", "w")
    for delf in sorted(del_frame):
        fout.write(str(delf) + "\t" + str(del_frame[delf]) + "\n")
    fout.close()

    fout = open(locus + "/intermediate_files/del_freq.txt", "w")
    for del_seqs in del_mut_freq:
        fout.write(del_seqs + "\t" + str(del_mut_freq[del_seqs]) + "\n")
    fout.close()

    del_pos_freq = {}
    for i in range(0,uplen):
        mut = "5DEL" + str(i)
        del_pos_freq.setdefault(mut,0)
    for i in range(0,downlen):
        mut = "3DEL" + str(i)
        del_pos_freq.setdefault(mut,0)

    if os.path.getsize(locus + "/intermediate_files/del_freq.txt") > 0:
        with open(locus + "/intermediate_files/del_freq.txt") as f:
            next(f)
            for lines in f:
                line = lines.rstrip().split("\t")
                dels = line[0].split(",")
                count = 0
                del_count = {}
                for Del in dels:
                    if not re.search("OneDel",Del):
                        count += 1
                        split = Del.split("_")
                        upDel = int(split[0].replace("5DEL",""))
                        downDel = int(split[1].replace("3DEL",""))
                        if downDel >= 0 and upDel >= 0:
                            for i in range(0,upDel+1):
                                mut = "5DEL"+str(i)
                                if mut in del_count:
                                    del_count[mut] += 1
                                else:
                                    del_count.setdefault(mut,0)
                                    del_count[mut] += 1
                            for j in range(0,downDel+1):
                                mut = "3DEL"+str(j)
                                if mut in del_count:
                                    del_count[mut] += 1
                                else:
                                    del_count.setdefault(mut,0)
                                    del_count[mut] += 1
                        elif downDel < 0:
                            for i in range(-downDel+1,upDel+1):
                                mut = "5DEL"+str(i)
                                if mut in del_count:
                                    del_count[mut] += 1
                                else:
                                    del_count.setdefault(mut,0)
                                    del_count[mut] += 1
                        elif upDel < 0:
                            for i in range(-upDel+1,downDel+1):
                                mut = "3DEL"+str(i)
                                if mut in del_count:
                                    del_count[mut] += 1
                                else:
                                    del_count.setdefault(mut,0)
                                    del_count[mut] += 1
                for del_mut in del_count:
                    ratio = float(del_count[del_mut]) / float(count)
                    freq =  ratio * float(line[1])
                    del_pos_freq[del_mut] += freq
        fout = open(locus + "/intermediate_files/" + locus + "_delPos_freq.txt","w")
        for i in range(uplen-1,0,-1):
            mut = "5DEL" + str(i)
            fout.write(mut + "\t" + str(del_pos_freq[mut]) + "\n")
        for i in range(1,downlen):
            mut = "3DEL" + str(i)
            fout.write(mut + "\t" + str(del_pos_freq[mut]) + "\n")
        fout.close()

    index_wt = 0
    index_unmod = 0
    index_del= 0
    index_ins = 0
    index_sub = 0
    pattern_wt = 0
    pattern_del = 0
    pattern_ins = 0
    pattern_cplx_del = 0
    pattern_sub = 0
    tbd = 0

    if "WT" in index_count:
        index_wt = index_count["WT"]
    if "Unmodified" in index_count:
        index_unmod = index_count["Unmodified"]
    if "Del" in index_count:
        index_del = index_count["Del"]
    if "Ins" in index_count:
        index_ins = index_count["Ins"]
    if "Sub" in index_count:
        index_sub = index_count["Sub"]
    if "WT" in groups:
        pattern_wt = groups["WT"]
    if "Del" in groups:
        pattern_del = groups["Del"]
    if "Ins" in groups:
        pattern_ins = groups["Ins"]
    if "Complexs_Del" in groups:
        pattern_cplx_del = groups["Complexs_Del"]
    if "Sub" in groups:
        pattern_sub = groups["Sub"]
    if "NA" in groups:
        tbd = groups["NA"]

    total_summary.write(locus
               + "\t" + str(index_wt+index_unmod) + "\t" + str(index_del) + "\t" + str(index_ins) + "\t" + str(index_sub)
               + "\t" + str(pattern_wt) + "\t" + str(pattern_del) + "\t" + str(pattern_ins) + "\t" + str(pattern_sub) + "\t" + str(pattern_cplx_del) + "\t"
               + "\t" + str(tbd)
               + "\t" + str(total) + "\n")
    total_summary.close()

    summary = open(locus + "/" + locus + "_freq_table.txt","w")
    summary.write("identified_amplicon" + "\t" + str(total) +"\t" + locus + "\n");
    new_total = total-tbd
    summary.write("classified_amp" + "\t" + str(new_total) + "\n");
    wt_count = index_wt+index_unmod+pattern_wt
    wt_precnt = float(wt_count) / float(total)
    summary.write("WT" + "\t" + str(wt_count) + "\t" + str(wt_precnt) + "\n")
    sub_count = index_sub + pattern_sub
    sub_precnt = float(sub_count)/ float(total)
    summary.write("Sub" + "\t" + str(sub_count) + "\t" + str(sub_precnt) + "\n")
    ins_count = index_ins + pattern_ins
    ins_precnt = float(ins_count) / float(total)
    summary.write("Ins" + "\t" + str(ins_count) + "\t" + str(ins_precnt) + "\n")
    del_count = index_del + pattern_del
    del_precnt = float(del_count) / float(total)
    summary.write("Del" + "\t" + str(del_count) + "\t" + str(del_precnt) + "\n")
    clx_precnt = float(pattern_cplx_del) / float(total)
    summary.write("Complexs_Del" + "\t" + str(pattern_cplx_del) + "\t" + str(clx_precnt) + "\n")
    tbd_precnt = float(tbd) / float(total)
    summary.write("NA" + "\t" + str(tbd) + "\t" + str(tbd_precnt) + "\n")
    summary.close()

def plot_freq(locus):
    sys.stderr.write("plotting figures"+"\n")
    index = []
    label = []
    freq = []
    count = 1
    middle = 0
    if os.path.exists(locus + "/intermediate_files/" + locus + "_delPos_freq.txt"):
        with open(locus + "/intermediate_files/" + locus + "_delPos_freq.txt") as f:
            for lines in f:
                line = lines.rstrip().split("\t")
                if line[0] == '5DEL1':
                    middle = count
                index.append(count)
                label.append(line[0])
                freq.append(float(line[1]))
                count += 1
        plt.bar(index,freq, width = 1, color = "#f08080")
        plt.xlabel('Deletion Position', size=10)
        plt.ylabel('Deletion Count', size=10)
        plt.title(locus)
        plt.axvline(middle+0.5, linestyle='--', lw=0.2)
        plt.savefig(locus + "/" + locus + '_del_dist.pdf', format='pdf')
        plt.close()

    if os.path.exists(locus + "/intermediate_files/del_size_dist.txt"):
        del_freq = {}
        with open(locus + "/intermediate_files/del_size_dist.txt") as f:
            for lines in f:
                line = lines.rstrip().split("\t")
                del_freq[int(line[0])] = int(line[1])
                largest_del = int(line[0])
        index = []
        freq = []
        for i in range(1,largest_del+1):
            index.append(i)
            if i in del_freq:
                freq.append(del_freq[i])
            else:
                freq.append(0)
        plt.bar(index,freq, width = 1, color = "#f08080")
        plt.xlabel('Deletion Size', size=10)
        plt.ylabel('Counts', size=10)
        plt.title(locus)
        plt.savefig(locus + "/" + locus + '_delSize_dist.pdf', format='pdf')
        plt.close()

    mutType_count = {}
    with open(locus + "/" + locus + "_freq_table.txt") as f:
        for lines in f:
            line = lines.rstrip().split("\t")
            mutType_count[line[0]]=int(line[1])

    plt.figure(figsize=(4.5,9))
    p1 = plt.bar(1,mutType_count["Del"], 0.8, color = "#f08080")
    p2 = plt.bar(1,mutType_count["Complexs_Del"], 0.8, bottom=mutType_count["Del"], color = "#ffb6c1")
    p3 = plt.bar(1,mutType_count["Ins"], 0.8, bottom=mutType_count["Del"]+mutType_count["Complexs_Del"], color = "#87cefa")
    p4 = plt.bar(1,mutType_count["Sub"], 0.8, bottom=mutType_count["Del"]+mutType_count["Complexs_Del"]+mutType_count["Ins"], color = "#9acd32")
    p5 = plt.bar(1,mutType_count["WT"], 0.8, bottom=mutType_count["Del"]+mutType_count["Complexs_Del"]+mutType_count["Ins"]+mutType_count["Sub"], color = "#f8f8d4")
    p6 = plt.bar(2,0,1)
    plt.xlabel('MutType Frequency', size=10)
    plt.ylabel('Counts', size=10)
    plt.title(locus)
    plt.xticks([])
    plt.legend((p5[0],p4[0],p3[0],p2[0],p1[0]),('Unmodified', 'Sub', 'Ins', 'Comp_Del', 'Del'), fontsize=10, frameon=False)
    plt.savefig(locus + "/" + locus + '_MutType_dist.pdf', format='pdf')
    plt.close()

    total = 0
    frame_count = {}
    with open(locus + "/intermediate_files/del_frame_count.txt") as f:
        for lines in f:
            line = lines.rstrip().split("\t")
            frame_count[int(line[0])]=int(line[1])
            total += int(line[1])
    plt.figure(figsize=(4,9))
    p1 = plt.bar(1,frame_count[2], 1, color = "#87cefa")
    p2 = plt.bar(1,frame_count[1], 1, bottom=frame_count[2], color = "#ff9846")
    p3 = plt.bar(1,frame_count[0], 1, bottom=frame_count[2]+frame_count[1], color = "#ffdd67")
    p4 = plt.bar(2,0,2)

    plt.xlabel('Deletion Frame', size=10)
    plt.ylabel('Counts', size=10)
    plt.title(locus)
    plt.xticks([])
    plt.legend((p3[0],p2[0],p1[0]),('3n\n' + str(round(float(frame_count[0])/float(total),3)), '3n-1\n' + str(round(float(frame_count[1])/float(total),3)), '3n-2\n' + str(round(float(frame_count[2])/float(total),3))), fontsize=13, frameon=False)
    plt.savefig(locus + "/" + locus + '_DelFrame_dist.pdf', format='pdf')
    plt.close()

    ins_freq = {}
    with open(locus + "/intermediate_files/ins_size_dist.txt") as f:
        for lines in f:
            line = lines.rstrip().split("\t")
            ins_freq[int(line[0])] = int(line[1])
            largest_ins = int(line[0])
    index = []
    freq = []
    for i in range(1,largest_ins+1):
        index.append(i)
        if i in ins_freq:
            freq.append(ins_freq[i])
        else:
            freq.append(0)
    plt.bar(index,freq, width = 1, color = "#87cefa")
    plt.xlabel('Insertion Size', size=10)
    plt.ylabel('Counts', size=10)
    plt.title(locus)
    plt.savefig(locus + "/" + locus + '_InsSize_dist.pdf', format='pdf')
    plt.close()

def ReverseComplement(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	return "".join([seq_dict[base] for base in reversed(seq)])

if __name__ == '__main__':
    args=parse_commandline()
    main(args)
