version 1.0
workflow spectra {
    input {
    	String output_directory
        File anndata_file
        File gene_dict_json_file
        String cell_type_key = 'broad_clusters'
        Boolean use_hvgs = true
        Boolean use_weights = true
        Float lam = 0.01
        Float delta = 0.001
        Float kappa = 0.00001
        Float rho = 0.00001
        Boolean use_cell_types = true
        Int n_top_vals = 25
        Int num_epochs = 10000
        #general parameters
        Boolean use_gpu = false
        String gpuType = "nvidia-tesla-p100"
        Int gpuCount = 0
        Int cpu = 4
        Int memory = 16
        Int extra_disk_space = 0
        String docker = "us.gcr.io/landerlab-atacseq-200218/mehta_spectra:latest"
        Int preemptible = 2
        Array[String] zones = ["us-central1-c"]
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
            lam = lam,
            delta = delta,
            kappa = kappa,
            rho = rho,
            use_cell_types = use_cell_types,
            n_top_vals = n_top_vals,
            num_epochs = num_epochs,
            cpu=cpu,
            memory=memory,
            extra_disk_space = extra_disk_space,
            docker=docker,
            use_gpu = use_gpu,
            gpuCount = gpuCount,
            gpuType = gpuType,
            preemptible=preemptible,
            zones = zones
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
        Float lam
        Float delta
        Float kappa
        Float rho
        Boolean use_cell_types
        Int n_top_vals
        Int num_epochs
        Int memory
        Int extra_disk_space
        Int cpu
        Boolean use_gpu
        Int gpuCount
        String gpuType
        String docker
        Int preemptible
        Array[String] zones
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
        from spectra import spectra_gpu as spc
        from spectra import spectra_util as util
        from spectra import K_est as kst
        import torch

        # check if gpu is detected
        # Check if CUDA is available
        if torch.cuda.is_available():
            print("CUDA is available")

            # Get the number of available GPUs
            num_gpus = torch.cuda.device_count()
            print(f"Number of available GPUs: {num_gpus}")

            # Get the name and memory status of each available GPU
            for i in range(num_gpus):
                gpu_name = torch.cuda.get_device_name(i)
                print(f"GPU {i}: {gpu_name}")

                # Get the memory information
                gpu_memory = torch.cuda.get_device_properties(i).total_memory
                gpu_memory_allocated = torch.cuda.memory_allocated(i)
                gpu_memory_cached = torch.cuda.memory_cached(i)
                gpu_memory_free = gpu_memory - gpu_memory_allocated - gpu_memory_cached

                print(f"\tTotal Memory: {gpu_memory / 1024**3:.2f} GB")
                print(f"\tAllocated Memory: {gpu_memory_allocated / 1024**3:.2f} GB")
                print(f"\tCached Memory: {gpu_memory_cached / 1024**3:.2f} GB")
                print(f"\tFree Memory: {gpu_memory_free / 1024**3:.2f} GB")
        else:
            print("CUDA is not available")

        with open("~{gene_dict_json_file}", 'rb') as file:
            annotations = json.load(file)

        adata = sc.read_h5ad("~{anndata_file}")
        cell_type_key = "~{cell_type_key}"
        use_hvgs = ~{true='True' false='False' use_hvgs}
        use_weights = ~{true='True' false='False' use_weights}
        lam = ~{lam}
        delta = ~{delta}
        kappa = ~{kappa}
        rho = ~{rho}
        use_cell_types = ~{true='True' false='False' use_cell_types}
        n_top_vals = ~{n_top_vals}
        num_epochs = ~{num_epochs}

        model = spc.est_spectra(adata = adata, gene_set_dictionary = annotations, 
                        use_highly_variable = use_hvgs, cell_type_key = cell_type_key, 
                        use_weights = use_weights, lam = lam, 
                        delta=delta, kappa = kappa, rho = rho, 
                        use_cell_types = use_cell_types, n_top_vals = n_top_vals, 
                        num_epochs=num_epochs
                       )

        with open('outputs/spectra_model.pickle', 'wb') as f:
            pickle.dump(model, f, pickle.HIGHEST_PROTOCOL)

        adata.write('outputs/spectra_adata.h5ad')

        CODE

        gsutil -m rsync -r outputs ~{output_dir}
    >>>

    output {
        File spectra_anndata_file = 'outputs/spectra_adata.h5ad'
        File spectra_model_file = 'outputs/spectra_model.pickle'
    }

    runtime {
        docker: docker
        memory: "~{memory}GB"
        bootDiskSizeGb: 12
        disks: "local-disk " + (ceil(size(anndata_file, "GB")*4) + extra_disk_space) + " HDD"
        cpu: cpu
        gpu: use_gpu
        gpuType: gpuType
        gpuCount: gpuCount
        preemptible: preemptible
        zones: zones
    }

}