setwd("~/projects/bioinfo/echinomics/test-datasets/scrna-seq/wessel-lab/")
library(Seurat)
library(patchwork)
library(tidyverse)
library(gridExtra)
library(parallel)
library(plotly)
library(cowplot)
library(grid)
library(tidyr)
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
intramutate = function(df, orig_search_list, replacement_list){
  
  for (i in 1:length(replacement_list)){
    
    df = case_when(
      df == orig_search_list[i] ~ replacement_list[i],
      df != orig_search_list[i] ~ df
    )
  }
  
  return(df)
  
}
process_protein_df = function(main_scRNAseq_dataset, features_df, term1, term2){
  
  protein_gene_locs = get_gene_locs(main_scRNAseq_dataset, features_df)
  protein_gene_names = get_gene_names(features_df, protein_gene_locs)
  protein_loc_unique_list = get_unique_list(features_df, protein_gene_locs)
  protein_loc_dup_list = get_dup_list(features_df)
  
  original_list = c(
    rownames(x = main_scRNAseq_dataset)[rownames(x = main_scRNAseq_dataset) %>% str_detect(term1)],
    rownames(x = main_scRNAseq_dataset)[rownames(x = main_scRNAseq_dataset) %>% str_detect(term2)]
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
labelled_feature_plotter = function(scrna_df, feature_df){
  
  plot_list = FeaturePlot(scrna_df, features = feature_df$GeneID)
  
  for (i in 1:nrow(feature_df)){
    plot_list[[i]] = plot_list[[i]] + ggtitle(feature_df$named_locus[i])
  }
  return(plot_list + NoLegend())
}
labelled_dot_plotter = function(scrna_df, feature_df){
  
  dp = DotPlot(scrna_df, features = feature_df$GeneID) + 
    scale_x_discrete(breaks=c(feature_df$GeneID),
                     labels=c(feature_df$named_locus)) +
    RotatedAxis()
  
  return(dp)
}

urchin_0 =readRDS("~/projects/bioinfo/echinomics/test-datasets/scrna-seq/wessel-lab/GSE149221_SpInteg.rds")
urchin_1 = RenameIdents(urchin_0,
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

unique_orig_idents = unique(urchin_1@meta.data$orig.ident)
unique_orig_idents_replacements = c("8 Cell",
                                    "64 Cell",
                                    "Morula",
                                    "Early Blastula",
                                    "Hatched Blastula",
                                    "Mid Blastula",
                                    "Early Gastrula",
                                    "Late Gastrula")


urchin_1@meta.data$orig.ident = intramutate(urchin_1@meta.data$orig.ident, unique_orig_idents, unique_orig_idents_replacements)

cell8 = subset(urchin_1, orig.ident  == "8 Cell")
cell64 = subset(urchin_1, orig.ident == "64 Cell")
morula = subset(urchin_1, orig.ident == "Morula")
eb = subset(urchin_1, orig.ident == "Early Blastula")
hb = subset(urchin_1, orig.ident == "Hatched Blastula")
mb = subset(urchin_1, orig.ident == "Mid Blastula")
eg = subset(urchin_1, orig.ident == "Early Gastrula")
lg = subset(urchin_1, orig.ident == "Late Gastrula")

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

abcb_0 = process_protein_df(urchin_1, features_abc, "abc", "ABC") %>% 
  subset(subfamily=="b") 
abcb = abcb_0[-8,]
abcc_0 = process_protein_df(urchin_1, features_abc, "abc", "ABC") %>% 
  filter(subfamily=="c") 
abcc = abcc_0[-c(6,17,19),]
slc_0 = process_protein_df(urchin_1, features_slc, "slc", "SLC")
slc = slc_0[-c(2,3,4,5,6,10,11,12,13,14,19),]
cyps = process_protein_df(urchin_1, features_cyps, "cyp", "CYP")
nhrs = process_protein_df(urchin_1, features_nhr, "nr", "nhr") %>% 
  transmute(GeneID = GeneID, named_locus = Name) %>% 
  drop_na()

c5_test1 = rownames(urchin_1) %in% c("LOC591982",
                          "LOC752233",
                          "LOC755624",
                          "LOC579421",
                          "LOC590074")
c5_test2 = rownames(urchin_1)[c5_test1]

c5_test3 = tibble(GeneID = c(c5_test2, "ABCC5a", "cell_typeA_score1","ABCC1", "ABCC9a"), named_locus = c("ABCC4a", "ABCC4d","ABCC4c","ABCC4b","ABCC5A_just_added", "ABCC5a","ABCC5a_combined", "ABCC1", "ABCC9a"))

# temp_efficient_plot(abcb, "abcb_x")
# temp_efficient_plot(abcc, "abcc_x")
# temp_efficient_plot(slc, "slc")
# temp_efficient_plot(cyps, "cyp")
# temp_efficient_plot(nhrs, "nhr")

labelled_dot_plotter(mb, slc)
labelled_dot_plotter(eg, abcc[7:10,])
a = labelled_feature_plotter(mb, slc[1:3,]) 
b = DimPlot(mb)
View(abcc_0)
devolist = c(cell8, cell64, morula, eb, hb, mb, eg, lg)

al = map2(.x = devolist, .y = unique_orig_idents_replacements, .f = function(x,y) (DimPlot(x) + ggtitle(y)))
pdf("all_dev_stage_tsne.pdf", onefile = T, width=12, height = 8)
al
dev.off()

bl = map2(.x = devolist, .y = unique_orig_idents_replacements, .f = function(x,y) (labelled_dot_plotter(x,c5_test3) + ggtitle(y)))
pdf("abcc_dotplot_cleaned_test1.pdf", onefile = T, width = 16, height = 9)
bl
dev.off()

cl = map2(.x = devolist, .y = unique_orig_idents_replacements, .f = function(x,y) (labelled_feature_plotter(x,abcc[11:20,]) + plot_annotation(title = y)))
pdf("abcc_featureplot_2.pdf", onefile = T, width = 21, height = 12)
cl
dev.off()

cell_typeA_marker_gene_list <- list(c("ABCC5a", "LOC590074"))
object <- AddModuleScore(object = mb, features = cell_typeA_marker_gene_list, name = "cell_typeA_score")
FeaturePlot(object = object, features = "cell_typeA_score1")

