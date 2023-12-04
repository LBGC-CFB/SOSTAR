rule assembly:
  input: 
    bam = f"{config['dirs']['outdir']}/alignment/{{alignment}}/{{sample}}.{{alignment}}.bam",
    bai = f"{config['dirs']['outdir']}/alignment/{{alignment}}/{{sample}}.{{alignment}}.bam.bai"
  output: f"{config['dirs']['outdir']}/assembly/{{alignment}}/{{sample}}.assembly.{{alignment}}.gtf"
  params: 
    min_isoform = lambda wildcards: '0.01' if wildcards.alignment == "realigned" else '0'
  threads:
    config['threads']
  container:
    "docker://aucam/lorid:latest"
  shell:
    """
    stringtie \
      -L \
      -p {threads} \
      -f {params.min_isoform} \
      -o {output} \
      {input.bam}
    """


rule merge:
  input: expand(f"{config['dirs']['outdir']}/assembly/{{{{alignment}}}}/{{sample}}.assembly.{{{{alignment}}}}.gtf", sample=config['samples'])
  output: 
    merged_gtf = temp(f"{config['dirs']['outdir']}/assembly/transcripts.merged.{{alignment}}.to_inter.gtf") if config['options']['bedtools'] else temp(f"{config['dirs']['outdir']}/assembly/transcripts.merged.{{alignment}}.to_correct.gtf")
  params:
    transcript_guide = config['refs']['ensembl_annot']
  threads:
    config['threads']
  container:
    "docker://aucam/lorid:latest"
  shell:
    """
    stringtie \
      -L \
      --merge \
      -p {threads} \
      -f 0 \
      -o {output.merged_gtf} \
      -G {params.transcript_guide} \
      {input}
    """


rule expression:
  input: 
    bam_files = f"{config['dirs']['outdir']}/alignment/realigned/{{sample}}.realigned.bam",
    transcript_guide =f"{config['dirs']['outdir']}/assembly/transcripts.merged.realigned.filter.gtf" if config['options']['bedtools'] else f"{config['dirs']['outdir']}/assembly/transcripts.merged.realigned.all.gtf",
  output: 
    output_exp = f"{config['dirs']['outdir']}/expression/{{sample}}.expression.gtf",
  params: 
    output_dir = config['dirs']['outdir']
  threads:
    config['threads']
  container:
    "docker://aucam/lorid:latest"
  shell:
    """
    stringtie \
      -L \
      -e \
      -p {threads} \
      -G {input.transcript_guide}  \
      -o {output.output_exp} \
      {input.bam_files}
    """