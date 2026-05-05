process SCANPY_ANNOTATE {
    label 'process_medium'

    publishDir "${params.outdir}/annotation", mode: params.publish_dir_mode

    input:
    path h5ad

    output:
    path "annotated.h5ad",      emit: h5ad
    path "umap_celltypes.png",  emit: umap
    path "celltypes.csv",       emit: celltypes
    path "versions.yml",        emit: versions

    script:
    """
    annotate.py \\
        --input ${h5ad} \\
        --output annotated.h5ad \\
        --model ${params.celltypist_model}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        celltypist: \$(python3 -c "import celltypist; print(celltypist.__version__)")
        scanpy: \$(python3 -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
