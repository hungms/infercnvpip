suppressMessages({
library(scUnify)
library(infercnv)
setwd("/nemo/lab/caladod/working/Matthew/project/oscar/OA_SC23423/20250131/output/infercnv")})


gene_annot <- read.table("/flask/scratch/caladod/hungm/reference/alignment/infercnv/GRCm39_gene_annotations.txt", sep = "\t", row.names = 1, header = F)
colnames(gene_annot) <- NULL

SC23423 = qread("/nemo/lab/caladod/working/Matthew/project/oscar/OA_SC23423/seurat/SC23423_merged_rchop_placebo_tumour_integrated.qs")
reference = qread("/nemo/lab/caladod/working/Matthew/project/oscar/OA_SC23423/seurat/20240628_SC23423_B220posGFPneg.qs")
genes <- intersect(rownames(SC23423@assays$RNA$counts), rownames(reference@assays$RNA$counts))

SC23423[["RNA_filter"]] <- CreateAssay5Object(SC23423@assays$RNA$counts[genes,], min.cells = 5)
SC23423@meta.data$group <- gsub("GFPpos", "", SC23423@meta.data$group)
SC23423@meta.data$group	<- gsub("B220pos", "PLACEBO", SC23423@meta.data$group)
SC23423@meta.data$leiden_0.6 = paste0("C", sprintf("%02d", as.numeric(SC23423@meta.data$leiden_0.6)-1))
print(unique(SC23423@meta.data$leiden_0.6))

reference[["RNA_filter"]] <- CreateAssay5Object(reference@assays$RNA$counts[genes,], min.cells = 5)
reference@meta.data$group <- "B220+GFP-"
reference@meta.data$run_id <- "B220+GFP-"
reference@meta.data$leiden_0.6 <- "B220+GFP-"

SC23423_total <- merge(SC23423, reference)
SC23423_total[["RNA_filter"]] <- JoinLayers(SC23423_total[["RNA_filter"]])

SC23423_total@meta.data <- SC23423_total@meta.data %>%
            mutate(mice = ifelse(is.na(`MULTIseq...ID`), as.character(run_id), paste0(as.character(run_id), "_", `MULTIseq...ID`)))

mice = unique(SC23423_total@meta.data$mice)
mice = mice[which(mice != "B220+GFP-")]

rm(SC23423, reference)

for(x in mice){
    obj <- subset(SC23423_total, subset = mice %in% c(mice, "B220+GFP-"))
    rna_counts = obj@assays$RNA_filter$counts
    annotations = obj@meta.data %>%
                    select("mice")
    colnames(annotations) <- NULL

    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=rna_counts,
                                        annotations_file=annotations,
                                        delim="\t",
                                        gene_order_file=gene_annot,
                                        ref_group_names=c("B220+GFP-"))

    outdir <- paste0("/nemo/lab/caladod/working/Matthew/project/oscar/OA_SC23423/20250131/output/infercnv/by_mice/", x)
    dir.create(outdir)
    infercnv_obj = infercnv::run(infercnv_obj,
                                cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir=outdir,
                                cluster_by_groups=F,
                                denoise=TRUE,
                                HMM=TRUE,
				analysis_mode = "subclusters",
				tumor_subcluster_partition_method = c("random_trees"),
HMM_type = c("i6"),
output_format="pdf",
num_threads = 30
    )
}
qsave(infercnv_obj, file = paste0(outdir, "infercnv_obj.qs"))






get_annotations <- function(org = "human", version = "latest"){
    if(version == "latest"){
        annot_files <- list.files(system.file("extdata", package = "infercnv"), full.names = T)
        annot_files <- annot_files[grepl(org, annot_files)]
        annot <- annot_files[length(annot_files)]}
    else{
        if(org == "human"){
            annot <- system.file("extdata", paste0("human_", version, "_gene_annotations.txt"), package = "infercnv")}
        else{
            annot <- system.file("extdata", paste0("mouse_", version, "_gene_annotations.txt"), package = "infercnv")}}

    gene_annot <- read.table(annot, sep = "\t", row.names = 1, header = F)
    colnames(gene_annot) <- NULL
    return(gene_annot)}





run_infercnv_individual <- function(obj, outdir){


list.files(system.file("extdata", package = "strpipe"), full.names = T)