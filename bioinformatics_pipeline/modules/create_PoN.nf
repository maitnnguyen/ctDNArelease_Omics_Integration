process create_PoN {

  tag "PoN"
  publishDir "${params.out_dir}/PoN", mode: 'copy'

  input:
    path cov_files  // list of all coverage files in one group!
    path gcmap_file
    path centromere_file

  output:
    path "PoN.rds", emit: pon

  script:
  """
  FILELIST=filelist.txt

  # create header expected by R
  echo -e "Key\\tFile" > \$FILELIST

  # build SampleID + absolute path table
  for f in ${cov_files}; do
    sample=\$(basename \$f _coverage.txt)
    echo -e "\$sample\\t\$(realpath \$f)" >> \$FILELIST
  done

  Rscript ${params.bin_dir}/create_PoN.R \
    --bin_dir ${params.bin_dir} \
    --filelist \$FILELIST \
    --gcmap ${gcmap_file} \
    --centromere ${centromere_file} \
    --outfile PoN \
    --method ${params.pon_method} \
    --pon ${params.pon_group}
  """
}
