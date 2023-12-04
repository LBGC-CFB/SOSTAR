rule SOSTAR:
  input:
    exp_file = expand(f"{config['dirs']['outdir']}/expression/{{sample}}.expression.gtf", sample=config['samples']),
    transcripts_ref = f"{config['dirs']['outdir']}/ref_transcripts_annotation.gtf"
  output: 
    output_file = f"{config['dirs']['outdir']}/SOSTAR_annotation_table_results.xlsx",
  params:
    input_folder = f"{config['dirs']['outdir']}/expression",
    output_folder = config['dirs']['outdir']
  shell:
    """
    python3 scripts/SOSTAR.py \
      -I {params.input_folder} \
      -R {input.transcripts_ref} \
      -O {params.output_folder}
    """