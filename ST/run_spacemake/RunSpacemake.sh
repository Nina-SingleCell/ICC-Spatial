#!/bin/bash
conda activate spacemake
project_id=ICC_Run
cd /mnt/ExperimentData/SpatialOmics
 
samples=()
for sample_id in "${samples[@]}"; do
    echo "${sample_id}"
    python3.10 /exdisk/SpaceMake/scripts/remove_duplicate_barcodes.py /mnt/ExperimentData/SpatialOmics/tile/${sample_id}/ /mnt/ExperimentData/SpatialOmics/tile/${sample_id} &
done

 
cd /exdisk/SpaceMake
project_id="ICC_Run"
strategy="rRNA:bowtie2->genome:STAR:final"
species_name=human
for sample_id in "${samples[@]}"; do

    data_path="/mnt/ExperimentData/SpatialOmics/fastq_ln/${sample_id}"
    tile_path="/mnt/ExperimentData/SpatialOmics/tile/${sample_id}"

    spacemake projects add_sample --project-id "${project_id}" \
                                  --sample-id "${sample_id}" \
                                  --R1 "${data_path}"/*R1_001.fastq.gz \
                                  --R2 "${data_path}"/*R2_001.fastq.gz \
                                  --species "${species_name}" \
                                  --puck openst \
                                  --run_mode openst \
                                  --barcode_flavor openst \
                                  --puck_barcode_file "${tile_path}"/*.txt.gz \
                                  --map_strategy "${strategy}"

    echo "------------------------------------------"
done
nohup spacemake run --cores 36 --keep-going > logs/report_log_samples.txt &


conda activate OPENSTT 

for sample_id in "${samples[@]}"; do

    source_folder=/exdisk//SpaceMake/projects/${project_id}/processed_data/${sample_id}/illumina/complete_data/dge
    destination_folder=/exdisk//SpaceMake/projects/${project_id}/processed_data/${sample_id}/h5ad

    mkdir -p "$destination_folder"

    count=0
    for file in "$source_folder"/*beads_fc*.h5ad; do
        if [[ -f "$file" ]]; then
            ln -sf "$file" "$destination_folder"/
            ((count++))
        fi
    done
    
    multimodal=/exdisk//SpaceMake/projects/${project_id}/processed_data/${sample_id}/multimodal
    mkdir $multimodal
    
    tileCoord=/exdisk//SpaceMake/puck_data/openst_coordinate_system.csv
    
    openst spatial_stitch --tiles $destination_folder/*.h5ad \
        --tile-id-regex "beads_(.+)\.h5ad$" \
        --tile-coordinates $tileCoord \
        --h5-out $multimodal/spatial_stitched_spots.h5ad
    images=/exdisk//SpaceMake/projects/${project_id}/processed_data/${sample_id}/images
    mkdir -p "$images"
    
    workdir="/exdisk//SpaceMake/projects/${project_id}/processed_data/${sample_id}"
    cd "$workdir"
    openst merge_modalities \
    --h5-in multimodal/spatial_stitched_spots.h5ad \
    --image-in "images/${sample_id}_ICC.tif" \
    --image-key "uns/staining_image"
done

conda activate OPENSTT 
cd  /exdisk/SpaceMake/projects/${project_id}/processed_data/${sample_id}/multimodal
openst manual_pairwise_aligner --h5-in spatial_stitched_spots.h5ad --spatial-key obsm/spatial --image-key uns/staining_image
