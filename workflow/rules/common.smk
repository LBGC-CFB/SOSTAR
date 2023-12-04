import pandas as pd 

# load config and sample sheet
samples = pd.read_csv(config['samples'], sep='\t', dtype='str').set_index('sample', drop=False)
print(samples)


def filter_ensembl_annot(ensembl_annot, transcripts_list):
    tr_annot_gtf = f"{config['dirs']['outdir']}/ref_transcripts_annotation.gtf"
    tr_annot_bed = f"{config['dirs']['outdir']}/ref_transcripts_annotation.bed"
    tr_list = []
    with open(transcripts_list, "r") as filin:
        for line in filin:
            line = line.split()
            tr_list.append(line[-1])

    with open(ensembl_annot, "r") as filin:
        with open(tr_annot_gtf, "w") as filout_gtf:
            filout_gtf.write("# LoRID: creation of reference transcripts annotation file\n")
            for line in filin:
                if line.startswith("#"):
                    filout.write(line)
                else:
                    lines = line.split()
                    dic_attr = {attr.split()[0]: attr.split()[1].replace('"',"") for attr in line.split("\t")[-1].split(";") if attr.split()}
                    if dic_attr["transcript_id"] in tr_list:
                        filout_gtf.write(line)
