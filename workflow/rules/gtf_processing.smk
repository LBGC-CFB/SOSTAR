rule get_tr_files:
  input:
    ensembl_annot = config['refs']['ensembl_annot'],
    transcripts_list = config['refs']['transcripts_list']
  output:
    tr_annot_gtf = f"{config['dirs']['outdir']}/ref_transcripts_annotation.gtf",
    tr_annot_bed = f"{config['dirs']['outdir']}/ref_transcripts_annotation.bed"
  run:
    tr_list = []
    with open(input[1], "r") as filin:
      for line in filin:
        line = line.split()
        tr_list.append(line[0])
    with open(input[0], "r") as filin:
      with open(output[0], "w") as filout_gtf:
        with open(output[1], "w") as filout_bed:
          filout_gtf.write("# LoRID: creation of reference transcripts annotation file\n")
          for line in filin:
            if line.startswith("#"):
              filout_gtf.write(line)
            else:
              dic_attr = {attr.split()[0]: attr.split()[1].replace('"',"") for attr in line.split("\t")[-1].split(";") if attr.split()}
              if dic_attr["transcript_id"] in tr_list:
                lines = line.split()
                filout_gtf.write(line)
                if lines[2] == "transcript":
                  filout_bed.write(f"{lines[0]}\t{lines[3]}\t{lines[4]}\t{dic_attr['gene_name']}\t755\t{lines[6]}\n")


rule bedtools_intersect:
  input: 
    merged_tr = f"{config['dirs']['outdir']}/assembly/transcripts.merged{{alignment}}.to_inter.gtf",
    bed_file = f"{config['dirs']['outdir']}/ref_transcripts_annotation.bed"
  output: temp(f"{config['dirs']['outdir']}/assembly/transcripts.merged{{alignment}}.inter.gtf")
  container:
    "docker://aucam/lorid:latest"
  shell:
    """
    bedtools intersect \
      -a {input.merged_tr} \
      -wa \
      -wb \
      -header \
      -s \
      -b {input.bed_file} > {output}
    """


rule format_gtf:
  input: f"{config['dirs']['outdir']}/assembly/transcripts.merged{{alignment}}.inter.gtf"
  output: f"{config['dirs']['outdir']}/assembly/transcripts.merged{{alignment}}.filter.gtf"
  run:
    with open(input[0], "r") as filin:
      with open(output[0], "w") as filout:
        filout.write("# LoRID: file cleaned with bedtools intersect results\n")
        for line in filin:
          if line.startswith("#"):
            filout.write(line)
          else:
            lines = line.split()
            gene = lines[-3].split("_")[0]
            dic_attr = {attr.split()[0]: attr.split()[1].replace('"',"") for attr in line.split("\t")[:-6][-1].split(";") if attr.split()}
            dic_attr["gene_name"] = gene
            filout.writelines(str(l) + "\t" for l in lines[:8])
            filout.writelines(f'{key} "{val}"; ' for key, val in dic_attr.items())
            filout.write("\n")


rule add_missing_gene_name:
  input:
    f"{config['dirs']['outdir']}/assembly/transcripts.merged.{{alignment}}.to_correct.gtf"
  output:
    f"{config['dirs']['outdir']}/assembly/transcripts.merged.{{alignment}}.all.gtf"
  run:
    import pandas as pd
    df = pd.read_csv(input[0], sep="\t", names = ['chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'], header=1)
    df_tr = df.loc[df["type"] == "transcript"]
    geneid_gename = {}
    for index, row in df_tr.iterrows():
      attributes = {attr.split()[0]: attr.split()[1].replace('"',"") for attr in row["attributes"].split(";") if attr.split()}
      if "gene_name" in attributes.keys():
        if attributes["gene_id"] not in geneid_gename.keys():
          geneid_gename[attributes["gene_id"]] = attributes["gene_name"]
    with open(input[0], "r") as filin:
      with open(output[0], "w") as filout:
        filout.write("# LoRID: file cleaned up by adding missing gene names\n")
        for line in filin:
          if line.startswith("#"):
            filout.write(line)
          else:
            lines = line.split()
            gene = lines[-3].split("_")[0]
            dic_attr = {attr.split()[0]: attr.split()[1].replace('"',"") for attr in line.split("\t")[-1].split(";") if attr.split()}
            if dic_attr["gene_id"] not in geneid_gename.keys():
              pass
            else:    
              dic_attr["gene_name"] = geneid_gename[dic_attr["gene_id"]] if not ("gene_name" in dic_attr.keys()) else dic_attr["gene_name"]                    
              filout.writelines(str(l) + "\t" for l in lines[:8])
              filout.writelines(f'{key} "{val}"; ' for key, val in dic_attr.items())
              filout.write("\n")


rule gtf_to_bed:
  input:
    f"{config['dirs']['outdir']}/assembly/transcripts.merged.aligned.filter.gtf" if config['options']['bedtools'] else f"{config['dirs']['outdir']}/assembly/transcripts.merged.aligned.all.gtf"
  output:
    f"{config['dirs']['outdir']}/assembly/transcripts.merged.aligned.filter.bed" if config['options']['bedtools'] else f"{config['dirs']['outdir']}/assembly/transcripts.merged.aligned.all.bed"
  run:
    chr = None
    with open(input[0], "r") as filin:
      with open(output[0], "w") as filout:
        filout.write("# LoRID: conversion of the gtf merge file into a bed file\n")
        for line in filin:
          if not line.startswith("#"):
            line = line.split()
            if line[2] == "transcript":
              if chr:
                filout.write(f"{chr}\t{startt}\t{endt}\t{id}\t723\t{strand}\t{startt}\t{endt}\t255,0,0\t{str(exon_count)}\t{str(','.join(exon_length))}\t{str(','.join(exon_start))}\n")
              exon_count = 0
              exon_length = []; exon_start = []
              chr = line[0]; startt = str(int(line[3])-1); endt = line[4]; id = line[11].split('"')[1]; strand = line[6]
            elif line[2] == "exon":
              exon_count += 1
              exon_length.append(str(int(line[4]) - int(line[3])+1))
              exon_start.append(str(int(line[3]) - int(startt)-1))