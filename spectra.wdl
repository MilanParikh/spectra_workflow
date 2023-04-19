version 1.0
workflow spectra {
    input {
    	String output_directory
        File anndata_file
        File gene_dict_json_file
        String cell_type_key = 'broad_clusters'
        Boolean use_hvgs = true
        Boolean use_weights = true
        Float lambda = 0.01
        Float delta = 0.001
        Float kappa = 0.00001
        Float rho = 0.00001
        Boolean use_cell_types = true
        Int n_top_vals = 25
        Int num_epochs = 10000
        #general parameters
        Int cpu = 8
        String memory = "64G"
        String docker = "mparikhbroad/spectra:latest"
        Int preemptible = 2
    }

    String output_directory_stripped = sub(output_directory, "/+$", "")

    call run_spectra_model {
        input:
            output_dir = output_directory_stripped,
            anndata_file = anndata_file,
            gene_dict_json_file = gene_dict_json_file,
            cell_type_key = cell_type_key,
            use_hvgs = use_hvgs,
            use_weights = use_weights,
            lambda = lambda,
            delta = delta,
            kappa = kappa,
            rho = rho,
            use_cell_types = use_cell_types,
            n_top_vals = n_top_vals,
            num_epochs = num_epochs,
            cpu=cpu,
            memory=memory,
            docker=docker,
            preemptible=preemptible
    }

    output {
        File spectra_anndata_file = run_spectra_model.spectra_anndata_file
        File spectra_model_file = run_spectra_model.spectra_model_file
    }
}

task run_spectra_model {

    input {
        String output_dir
        File anndata_file
        File gene_dict_json_file
        String cell_type_key
        Boolean use_hvgs
        Boolean use_weights
        Float lambda
        Float delta
        Float kappa
        Float rho
        Boolean use_cell_types
        Int n_top_vals
        Int num_epochs
        String memory
        Int cpu
        String docker
        Int preemptible
    }

    command <<<
        set -e

        mkdir -p outputs

        python <<CODE
        import os
        import scanpy as sc
        import pandas as pd
        import numpy as np
        import json
        import scipy
        import pickle
        from spectra import spectra as spc
        from spectra import spectra_util as util
        from spectra import K_est as kst

        with open("~{gene_dict_json_file}", 'rb') as file:
            annotations = json.load(file)

        adata = sc.read_h5ad("~{anndata_file}")
        cell_type_key = "~{cell_type_key}"
        use_hvgs = ~{true='True' false='False' use_hvgs}
        use_weights = ~{true='True' false='False' use_weights}
        lambda = ~{lambda}
        delta = ~{delta}
        kappa = ~{kappa}
        rho = ~{rho}
        use_cell_types = ~{true='True' false='False' use_cell_types}
        n_top_vals = ~{n_top_vals}
        num_epochs = ~{num_epochs}

        model = spc.est_spectra(adata = adata, gene_set_dictionary = annotations, 
                        use_highly_variable = use_hvgs, cell_type_key = cell_type_key, 
                        use_weights = use_weights, lam = lambda, 
                        delta=delta, kappa = kappa, rho = rho, 
                        use_cell_types = use_cell_types, n_top_vals = n_top_vals, 
                        num_epochs=num_epochs
                       )

        with open('outputs/spectra_model.pickle', 'wb') as f:
            pickle.dump(model, f, pickle.HIGHEST_PROTOCOL)

        adata.write('/outputs/spectra_adata.h5ad')

        CODE

        gsutil -m rsync -r outputs ~{output_dir}
    >>>

    output {
        File spectra_anndata_file = 'outputs/spectra_adata.h5ad'
        File spectra_model_file = 'outputs/spectra_model.pickle'
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(anndata_file, "GB")*2) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}