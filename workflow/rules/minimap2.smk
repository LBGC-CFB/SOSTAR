rule align:
  input: f"{config['dirs']['indir']}/{{sample}}.fastq.gz"
  output: 
    bam = f"{config['dirs']['outdir']}/alignment/aligned/{{sample}}.aligned.bam",
    bai = f"{config['dirs']['outdir']}/alignment/aligned/{{sample}}.aligned.bam.bai",
    sam = temp(f"{config['dirs']['outdir']}/alignment/aligned/{{sample}}.aligned.sam"),
    sorted_bam = temp(f"{config['dirs']['outdir']}/alignment/aligned/{{sample}}_sorted.bam")
  params:
    ref = config['refs']['genome']
  threads: 
    config['threads']
  container:
    "docker://aucam/lorid:latest"
  shell:
    """
    minimap2 \
      -ax splice \
      --MD \
      -t {threads} \
      {params.ref} \
      -o {output.sam} \
      {input}

    samtools view \
      -b \
      -@ {threads} \
      -o  {output.sorted_bam} \
      {output.sam}
    
    samtools sort \
      -@ {threads} \
      -o {output.bam} \
      {output.sorted_bam} && \
      samtools index \
      -@ {threads} \
        {output.bam}
    """


rule realign:
  input: 
    fastq = f"{config['dirs']['indir']}/{{sample}}.fastq.gz",
    bed = f"{config['dirs']['outdir']}/assembly/transcripts.merged.aligned.filter.bed" if config['options']['bedtools'] else f"{config['dirs']['outdir']}/assembly/transcripts.merged.aligned.all.bed"
  output: 
    bam = f"{config['dirs']['outdir']}/alignment/realigned/{{sample}}.realigned.bam",
    bai = f"{config['dirs']['outdir']}/alignment/realigned/{{sample}}.realigned.bam.bai",
    sam = temp(f"{config['dirs']['outdir']}/alignment/realigned/{{sample}}.realigned.sam"),
    sorted_bam = temp(f"{config['dirs']['outdir']}/alignment/realigned/{{sample}}_sorted.realigned.bam")
  params:
    ref = config['refs']['genome']
  threads:
    config['threads']
  container:
    "docker://aucam/lorid:latest"
  shell:
    """
    minimap2 \
      -ax splice \
      --MD \
      -t {threads} \
      --junc-bed {input.bed} \
      {params.ref} \
      -o {output.sam} \
      {input.fastq}

    samtools view \
      -b \
      -@ {threads} \
      -o {output.sorted_bam} \
      {output.sam}
    
    samtools sort \
      -@ {threads} \
      -o {output.bam} \
      {output.sorted_bam} && \
      samtools index \
        -@ {threads} \
        {output.bam}
    """