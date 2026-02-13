# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    sub_mat_path = "./substitution_matrices/BLOSUM62.mat"
    nw = NeedlemanWunsch(sub_mat_path, -10, -1)
    list_align = [gg_seq, mm_seq, br_seq, tt_seq]
    list_names = ['Gallus_gallus', 'mus_musculus', 'Balaeniceps_rex_BRD2', 'tursiops_truncatus_BRD2']
    alignment_scores = []
    for seq in list_align:
        nw_alignment_score, seqA_aligned, seqB_aligned = nw.align(hs_seq, seq)
        alignment_scores.append(nw_alignment_score)
    
    # return sorted list (in decreasing order of those most similar to human BRD)
    sorted_list_names = [val for _, val in sorted(zip(alignment_scores, list_names), reverse=True)]
    print(sorted_list_names)

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    for i in range(len(list_align)):
        print(f'Score for {list_names[i]} w/ humans: {alignment_scores[i]}')
    

if __name__ == "__main__":
    main()
