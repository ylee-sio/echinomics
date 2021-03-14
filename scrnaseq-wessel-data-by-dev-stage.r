setwd("~/projects/bioinfo/echinomics/test-datasets/scrna-seq/wessel-lab/")
library(Seurat)
library(patchwork)
library(tidyverse)
library(gridExtra)
library(parallel)
library(plotly)
options(ncpus = 14)

transform_prot_tibble = function(path){
  features_tibble = read_csv(path) %>%
    select(GeneID, Name) %>% 
    unique() %>% 
    drop_na() %>%
    mutate(protein = str_sub(.$Name,4,6) %>% tolower(),
           subfamily = str_sub(.$Name,7,8) %>% tolower(),
           member = str_sub(.$Name,9,) %>% tolower()) %>%
    transmute(GeneID=GeneID,
              protein = protein,
              subfamily = subfamily,
              member = member)
  
  return(feature_tibble)
}
get_gene_locs = function(main_scRNAseq_dataset, features_list){
  
  feature_loc_list = rownames(x = main_scRNAseq_dataset)[rownames(x = main_scRNAseq_dataset) %in% features_list$GeneID]
  return(feature_loc_list)
  
}
get_gene_names = function(features_list, gene_loc_list){
  
  feature_names_list = features_list %>% 
    subset(GeneID %in% gene_loc_list) %>% 
    mutate(named_locus = paste0(.$protein,.$subfamily,"_",.$member)) %>% 
    select(named_locus)
  
  return(feature_names_list)
  
}
get_unique_list = function(features_list, gene_loc_list, by_loc = T, by_name = F){
  
  #loc_unique_list = features_list[which(!(features_list$GeneID %>% duplicated())),]
  loc_unique_list = features_list[which(features_list$GeneID %in% gene_loc_list),] %>% 
    mutate(named_locus = paste0(.$protein,.$subfamily,"_",.$member))
  
  return(loc_unique_list)
  
}
get_dup_list = function(features_list, by_loc = T, by_name = F){
  
  loc_unique_list = features_list[which((features_list$GeneID %>% duplicated())),]%>% 
    mutate(named_locus = paste0(.$protein,.$subfamily,"_",.$member))
  return(loc_unique_list)
  
}
plot_all_stages_specific_protein = function(stage_df, protein_df, stage_name, protein_name) {
  
  # stage_df=cell8
  # protein_df=abcb
  # stage_name="8_cell"
  # protein_name="abcb"

  start = Sys.time()
  print(paste0("Calculating for and plotting: ", stage_name))
  protein_df_geneid = protein_df$GeneID
  protein_df_named_locus = protein_df$named_locus
  
  imax = numeric()
  if (nrow(protein_df) >= 8){
    imax = ceiling(nrow(protein_df)/8)
  }
  
  i=1
  for (i in 1:imax){
    s = (8*i)-7
    t = 8*i
    
    if (t >= nrow(protein_df)){
      t = nrow(protein_df)
    }
    # print("Error1")
    fig1 = VlnPlot(stage_df, features=protein_df_geneid[s:t])
    fig1 = lapply(1:length(protein_df_geneid[s:t]), function(x) {fig1[[x]] + labs(title=protein_df_named_locus[s:t][x])})
    fig1_plotly_list = list()
    
    fig1_plotly_list = lapply(1:length(protein_df_geneid[s:t]), function(z) {fig1_plotly_list[z] = ggplotly(fig1[[z]]) %>% layout(annotations = 
                                                                                                                  list(
                                                                                                                    x = 0.2 , 
                                                                                                                    y = 1.05, 
                                                                                                                    text = protein_df_named_locus[(s-1)+z], 
                                                                                                                    showarrow = F, 
                                                                                                                    xref='paper', 
                                                                                                                    yref='paper', 
                                                                                                                    xanchor="center"
                                                                                                                  )
    ) }
    
    )

    
    # print("Error2")
    a = fig1_plotly_list %>% 
      subplot(nrows = 2, shareX = T, margin = 0.04) %>% 
      layout(title = paste0("SMT: ", protein_name, "; ", "Dev Stage: ", stage_name))
    # print("Error3")

    dir.create(paste0("plotly/", protein_name))
    htmlwidgets::saveWidget(as_widget(a), paste0("plotly/", protein_name, "/", protein_name, "_", stage_name,"_", as.character(i),".html"))
    end = Sys.time()
    elapsed = end-start
    print(paste0(elapsed %>% as.numeric, " seconds/minutes elapsed in creating plot for ", stage_name))
    print(paste0("Plot saved in ", getwd(),"/plotly/", protein_name))
    print("********************************************************************")
  }
 }
temp_efficient_plot = function(protein_df,protein_name){
  
  plot_all_stages_specific_protein(cell8,protein_df,"8cell",as.character(protein_name))
  plot_all_stages_specific_protein(cell64,protein_df,"64cell",as.character(protein_name))
  plot_all_stages_specific_protein(morula,protein_df,"morula",as.character(protein_name))
  plot_all_stages_specific_protein(eb,protein_df,"early_blastula",as.character(protein_name))
  plot_all_stages_specific_protein(hb,protein_df,"hatched_blastula",as.character(protein_name))
  plot_all_stages_specific_protein(mb,protein_df,"mid_blastula",as.character(protein_name))
  plot_all_stages_specific_protein(eg,protein_df,"early_gastrula",as.character(protein_name))
  plot_all_stages_specific_protein(lg,protein_df,"late_gastrula",as.character(protein_name))

}

urchin =readRDS("~/projects/bioinfo/echinomics/test-datasets/scrna-seq/wessel-lab/GSE149221_SpInteg.rds")

features_abc = read_csv("abc_list.csv") %>%
  select(GeneID, protein, subfamily, member)
features_slc = read_csv("slc_list.csv") %>%
  select(GeneID, protein, subfamily, member) 
features_cyps = read_csv("cyps_list.csv") %>%
  select(GeneID, protein, subfamily, member) 
features_nhr = read_csv("nhr_list.csv") %>%
  select(GeneID, Name) %>% 
  drop_na()
features_final = read_csv("final_curation.csv")

urchin = RenameIdents(urchin,
                  '0' = "aboral-ectoderm",
                  '1' = "oral_ectoderm_eb",
                  '2' = "ciliated_cells",
                  '3' = "neural_8_to_eg",
                  '4' = "oral_ectoderm_hb",
                  '5' = "aboral_ectoderm/neural",
                  '6' = "endoderm_hb",
                  '7' = "ciliated_cells_m",
                  '8' = "endoderm_endo-msd",
                  '9' = "neural_AnEctoE",
                  '10' = "transient_8_to_eb",
                  '11' = "nsm_pigment_cells",
                  '12' = "oral_ectoderm_64",
                  '13' = "oral_ectoderm_m",
                  '14' = "endoderm_eb",
                  '15' = "wtf",
                  '16' = "skeleton_64",
                  '17' = "neural_AnEctoL",
                  '18' = "neural_neuro_prog",
                  '19' = "skeleton_hb",
                  '20' = "transient",
                  '21' = "germline")

process_protein_df = function(main_scRNAseq_dataset, features_df, term1, term2){
  
  protein_gene_locs = get_gene_locs(main_scRNAseq_dataset, features_df)
  protein_gene_names = get_gene_names(features_df, protein_gene_locs)
  protein_loc_unique_list = get_unique_list(features_df, protein_gene_locs)
  protein_loc_dup_list = get_dup_list(features_df)
  
  original_list = c(
    rownames(x = urchin)[rownames(x = urchin) %>% str_detect(term1)],
    rownames(x = urchin)[rownames(x = urchin) %>% str_detect(term2)]
  )
  original_list_tibble = tibble(GeneID = original_list, 
                                protein = str_sub(original_list,1,3) %>% tolower(), 
                                subfamily = str_sub(original_list,4,4) %>% tolower(), 
                                member = str_sub(original_list,5) %>% tolower(),
                                named_locus = original_list %>% tolower()
  )
  
  final_combined_df = bind_rows(original_list_tibble, protein_loc_unique_list)
  return(final_combined_df)
}

abcb = process_protein_df(urchin, features_abc, "abc", "ABC") %>% 
  subset(subfamily=="b") 
abcc = process_protein_df(urchin, features_abc, "abc", "ABC") %>% 
  filter(subfamily=="c") 
slc = process_protein_df(urchin, features_slc, "slc", "SLC")
cyps = process_protein_df(urchin, features_cyps, "cyp", "CYP")
nhrs = process_protein_df(urchin, features_nhr, "nr", "nhr") %>% 
  transmute(GeneID = GeneID, named_locus = Name) %>% 
  drop_na()

# 8 cell stage stage *************************************************************************************************************
p1_bc_cell8 = read_tsv("~/projects/bioinfo/echinomics/test-datasets/scrna-seq/wessel-lab/8cell/GSM4494538_Sp1_barcodes.tsv") %>% as_vector()
bc_cell8 = paste0("Sp1_", str_sub(p1_bc_cell8, 1, 16))
cell8 = subset(urchin, cells = bc_cell8)
# 64 cell stage stage *************************************************************************************************************
p1_bc_cell64 = read_tsv("~/projects/bioinfo/echinomics/test-datasets/scrna-seq/wessel-lab/64cell/GSM4494539_Sp2_barcodes.tsv") %>% as_vector()
bc_cell64 = paste0("Sp2_", str_sub(p1_bc_cell64, 1, 16))
cell64 = subset(urchin, cells = bc_cell64)
# morula stage *************************************************************************************************************
p1_bc_morula = read_tsv("~/projects/bioinfo/echinomics//test-datasets/scrna-seq/wessel-lab/morula/GSM4494540_Sp3_barcodes.tsv") %>% as_vector()
bc_morula = paste0("Sp3_", str_sub(p1_bc_morula, 1, 16))
morula = subset(urchin, cells = bc_morula)
# eb stage *************************************************************************************************************
p1_bc_eb = read_tsv("~/projects/bioinfo/echinomics/test-datasets/scrna-seq/wessel-lab/early_blastula//GSM4494541_SpEB_barcodes.tsv") %>% as_vector()
bc_eb = paste0("EB_", str_sub(p1_bc_eb, 1, 16))
eb = subset(urchin, cells = bc_eb)
# hb stage *************************************************************************************************************
p1_bc_hb = read_tsv("~/projects/bioinfo/echinomics/test-datasets/scrna-seq/wessel-lab/hatched_blastula/GSM4494542_AGG1_SpHB_barcodes.tsv") %>% as_vector()
bc_hb = paste0("HB_", str_sub(p1_bc_hb, 1, 18))
hb = subset(urchin, cells = bc_hb)
# mb stage *************************************************************************************************************
p1_bc_mb = read_tsv("~/projects/bioinfo/echinomics/test-datasets/scrna-seq/wessel-lab/mid_blastula/GSM4494543_AGG2_SpMB_barcodes.tsv") %>% as_vector()
bc_mb = paste0("MB_", str_sub(p1_bc_mb, 1, 18))
mb = subset(urchin, cells = bc_mb)
# eg stage ***************************** ********************************************************************************
p1_bc_eg = read_tsv("~/projects/bioinfo/echinomics/test-datasets/scrna-seq/wessel-lab/early_gastrula/GSM4494544_AGG3_SpEG_barcodes.tsv") %>% as_vector()
bc_eg = paste0("EG_", str_sub(p1_bc_eg, 1, 18))
eg = subset(urchin, cells = bc_eg)
# lg stage *************************************************************************************************************
p1_bc_lg = read_tsv("~/projects/bioinfo/echinomics/test-datasets/scrna-seq/wessel-lab/late_gastrula/GSM4494545_AGG4_SpLG_barcodes.tsv") %>% as_vector()
bc_lg = paste0("LG_", str_sub(p1_bc_lg, 1, 18))
lg = subset(urchin, cells = bc_lg)

temp_efficient_plot(abcb, "abcb_x")
temp_efficient_plot(abcc, "abcc_x")
temp_efficient_plot(slc, "slc")
temp_efficient_plot(cyps, "cyp")
temp_efficient_plot(nhrs, "nhr")

