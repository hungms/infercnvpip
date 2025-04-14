#' Get gene annotations
#'
#' This function retrieves the gene annotations for a given organism and version.
#'
#' @param org The organism to get annotations for. Can be "human" or "mouse".
#' @param version The version of the annotations to use. Can be "latest" or a specific version number.
#' @return A data frame containing the gene annotations.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' get_annotations(org = "human", version = "latest")
#' }
get_annotations <- function(org = "human", version = "latest"){
    if(version == "latest"){
        annot_files <- list.files(system.file("extdata", package = "infercnvpip"), full.names = T)
        annot_files <- annot_files[grepl(org, annot_files)]
        selected_annot <- annot_files[length(annot_files)]}
    else{
        if(org == "human"){
            selected_annot <- system.file("extdata", paste0("human_", version, "_gene_annotations.txt"), package = "infercnv")}
        else{
            selected_annot <- system.file("extdata", paste0("mouse_", version, "_gene_annotations.txt"), package = "infercnv")}}

    gene_annot <- read.table(selected_annot, sep = "\t", row.names = 1, header = F)
    colnames(gene_annot) <- NULL
    return(gene_annot)}


#' Run infercnv on an individual object
#'
#' This function runs infercnv on an individual object.
#'
#' @param obj The object to run infercnv on.
#' @param ref_obj The reference object to run infercnv on.
#' @param individual_name The name of the individual to run infercnv on.
#' @param annot_column The column name of the annotations to run infercnv on.
#' @param org The organism to run infercnv on.
#' @param version The version of the annotations to use.
#' @param assay The assay to run infercnv on.
#' @param save_dir The directory to save the infercnv results.
#' @param ncores The number of cores to use.
#' @param ... Additional arguments to pass to infercnv::run.
#' @return A list of infercnv objects.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' run_infercnv_individual(obj, ref_obj = NULL, individual_name, annot_column, org = "human", version = "latest", assay = "RNA", save_dir = getwd(), ncores = 1, ...)
#' }
run_infercnv_individual <- function(obj, ref_obj = NULL, individual_name, annot_column, org, version = "latest", assay = "RNA", save_dir = getwd(), ncores = 1, ...){

    stopifnot(assay %in% names(obj@assays))
    stopifnot(annot_column %in% colnames(obj@meta.data))
    stopifnot(org %in% c("human", "mouse"))

    # merge objects
    if(is.null(ref_obj)){
        merged_obj <- obj
        ref_group_names <- NULL
        message("No reference groups are set")
        }
    else{
        message("Merging query and reference objects as a single Seurat object for inferCNV input...")
        stopifnot(assay %in% names(ref_obj@assays))
        stopifnot(annot_column %in% colnames(ref_obj@meta.data))
        merged_obj <- merge(obj, ref_obj)
        merged_obj[[assay]] <- JoinLayers(merged_obj[[assay]])
        ref_group_names <- unique(ref_obj@meta.data[[annot_column]])
        message("Reference groups are set to: ", paste(ref_group_names, collapse = ", "))
        }

    message("Creating infercnv object...")
    # get counts
    counts = as.matrix(merged_obj[[assay]]$counts)

    # get cell annotations
    cell_annot = merged_obj@meta.data %>% dplyr::select(annot_column) # %>% rownames_to_column(var = "cell_id") %>% as.tibble(.)

    # get gene annotations
    gene_annot <- get_annotations(org = org, version = version)

    # create infercnv object
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts,
                                        annotations_file=cell_annot,
                                        delim="\t",
                                        gene_order_file=gene_annot,
                                        ref_group_names=ref_group_names)

    # set output directory
    out_dir <- paste0(save_dir, "/", individual_name)
    dir.create(out_dir, recursive = T)

    message("Running inferCNV on ", individual_name, "...")
    # run infercnv with cluster_by_groups = F
    processed_infercnv_obj = infercnv::run(infercnv_obj,
                                cutoff=0.1,
                                out_dir=paste0(out_dir, "/NO-GROUPING/"),
                                cluster_by_groups=F,
                                denoise=TRUE,
                                HMM=TRUE,
                                analysis_mode = "subclusters",
                                tumor_subcluster_partition_method = c("random_trees"),
                                HMM_type = c("i6"),
                                output_format="pdf",
                                num_threads = ncores,
                                ...)

    # run infercnv with cluster_by_groups = T
    processed_infercnv_obj = infercnv::run(infercnv_obj,
                                cutoff=0.1,
                                out_dir=paste0(out_dir, "/GROUP-BY-", annot_column, "/"),
                                cluster_by_groups=T,
                                denoise=TRUE,
                                HMM=TRUE,
                                analysis_mode = "subclusters",
                                tumor_subcluster_partition_method = c("random_trees"),
                                HMM_type = c("i6"),
                                output_format="pdf",
                                num_threads = ncores,
                                ...)
}



#' Run infercnv on multiple objects
#'
#' This function runs infercnv on multiple objects.
#'
#' @param obj The object to run infercnv on.    
#' @param split.by The column name to split the object by.
#' @param ref_obj The reference object to run infercnv on.
#' @param org The organism to run infercnv on.
#' @param save_dir The directory to save the infercnv results.
#' @param ... Additional arguments to pass to run_infercnv_individual.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' run_infercnv_multi(obj, split.by, ref_obj = NULL, org = "human", annot_column, save_dir = getwd(), ...)
#' }
run_infercnv_multi <- function(obj, split.by, ref_obj = NULL, org, annot_column, save_dir = getwd(), ...){

    stopifnot(split.by %in% colnames(obj@meta.data))
    stopifnot(annot_column %in% colnames(obj@meta.data))
    stopifnot(org %in% c("human", "mouse"))

    # split object
    message("Splitting query object per individual...")
    obj.list <- SplitObject(obj, split.by = split.by)

    if(!is.null(ref_obj)){

        if(split.by %in% colnames(ref_obj@meta.data)){
            message("Splitting reference object per individual...")
            ref_obj.list <- SplitObject(ref_obj, split.by = split.by)
            ref_obj.list <- ref_obj.list[which(names(ref_obj.list) %in% names(obj.list))]
            ref_obj.list <- ref_obj.list[names(obj.list)]
            stopifnot(all(names(obj.list) == names(ref_obj.list)))
            }

        message("Running infercnv individually on ", paste0(names(obj.list), collapse = ", "), "...")

        # run infercnv on each unique value
        for(i in 1:length(obj.list)){
            run_infercnv_individual(
                obj = obj.list[[i]], 
                ref_obj = ref_obj.list[[i]], 
                individual_name = names(obj.list)[i], 
                annot_column = annot_column,
                org = org,
                save_dir = save_dir,
                ...)}
        } else {
        # run infercnv on each unique value
        for(i in 1:length(obj.list)){
            run_infercnv_individual(
                obj = obj.list[[i]], 
                ref_obj = ref_obj, 
                individual_name = names(obj.list)[i], 
                annot_column = annot_column,
                org = org,
                save_dir = save_dir,
                ...)}
    }
}