#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import os

def combine_gfa(input_edges, sorted_depth_file, outputdir):


    myList = []
    for seq_record in SeqIO.parse(input_edges, "fasta"):
        myList.append([seq_record.id, str(seq_record.seq), len(seq_record)])

    myList.sort(key=lambda x: x[2])
    longest = myList[-1][0]
    #print(longest)

    with open(sorted_depth_file) as f:
        first_line = f.readline()
        IR = first_line.split("\t")[0]

        #print(IR)

    #print(IR, longest)


    for seq_record in SeqIO.parse(input_edges,"fasta"):
        # print(seq_record)
        if seq_record.id == longest:
            LSC_seq = seq_record.seq


        if seq_record.id == IR:
            IR_seq = seq_record.seq
            IR_seq_RC = seq_record.seq.reverse_complement()

        if seq_record.id != IR and seq_record.id != longest:
            SSR_seq = seq_record.seq
            SSR2_seq = SSR_seq.reverse_complement()
            #print(SSR)
            #print(SSR2)

    path1 = LSC_seq + IR_seq + SSR_seq + IR_seq_RC
    path2 = LSC_seq + IR_seq + SSR2_seq + IR_seq_RC

    DNA1 = ">path1" + "\n" + path1 + "\n"
    DNA2 = ">path2" + "\n" + path2 + "\n"
    #print(DNA1)
    #print(DNA2)
    #outputdir
    outputfile1 = outputdir + "/path1.fasta"
    outputfile2 = outputdir + "/path2.fasta"

    #saveFasta = open(r'path1.fasta', 'w')
    saveFasta = open(outputfile1, 'w')
    saveFasta.write(str(DNA1))
    saveFasta.close()

    #saveFasta = open(r'path2.fasta', 'w')
    saveFasta = open(outputfile2, 'w')
    saveFasta.write(str(DNA2))
    saveFasta.close()


def main():
    parser = argparse.ArgumentParser(
        description="This script is used to merge the edges from assembly graph.",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-e", "--input_edges_file", dest='input_edges', type=str,
                        help='input edge fasta file')

    parser.add_argument("-d", "--sorted_depth", dest='sorted_depth_file', type=str,
                        help='input depth file')
    
    parser.add_argument("-o", "--outputdir", dest='outputdir', type=str,
                        help='output directory')

    args = parser.parse_args()
    if args.input_edges and args.sorted_depth_file and args.outputdir:
        input_edge_filename = os.path.realpath(args.input_edges)
        input_depth_filename = os.path.realpath(args.sorted_depth_file)
        output_dir = os.path.realpath(args.outputdir)

        combine_gfa(input_edge_filename, input_depth_filename, output_dir)




if __name__ == '__main__':
    main()
