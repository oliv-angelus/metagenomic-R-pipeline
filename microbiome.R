# ==============================================================================
# ==============================================================================
# SCRIPT DE ANÁLISE DE MICROBIOMA - METAGENÔMICA
# ==============================================================================
# ==============================================================================

# ==============================================================================
# === # 1. INSTALAÇÃO DE PACOTES (Execute apenas na primeira vez) # === #
# ==============================================================================
{
  install.packages("devtools")
  install.packages("BiocManager")
  install.packages("remotes")
  install.packages("tidyverse")
  install.packages("dplyr")
  install.packages("RColorBrewer")
  install.packages("patchwork")
  install.packages("ggplot2")
  install.packages("vegan")
  install.packages("permute")
  install.packages("tibble")
  install.packages("ggrepel")
  install.packages("forcats")
  install.packages("indicspecies")
  install.packages(c("igraph", "ggraph", "Hmisc"))
  
  BiocManager::install("phyloseq")
  BiocManager::install("ape")
  BiocManager::install("DESeq2")
}

# ==============================================================================
# === # 2. CARREGAR BIBLIOTECAS # === #
# ==============================================================================
{
  library(phyloseq)
  library(tidyverse)
  library(ape)
  library(vegan)
  library(RColorBrewer)
  library(patchwork)
  library(ggplot2)
  library(permute)
  library(dplyr)
  library(tibble)
  library(ggrepel)
  library(forcats) 
  library(DESeq2)
  library(indicspecies)
  library(igraph)
  library(ggraph)
  library(Hmisc)
}

# ==============================================================================
# === # 3. DEFINIÇÃO DE VARIÁVEIS E PARÂMETROS (PAINEL DE CONTROLE) # === #
# ==============================================================================

# --- 3.1. Arquivos de Entrada ---
ARQUIVO_ABUNDANCIA <- "abundance.tsv"  # Tabela de abundância (TSV)
ARQUIVO_METADADOS  <- "metadados.tsv"            # Metadados (TSV)

# --- 3.2. Variáveis Experimentais ---
VAR_AGRUPAMENTO  <- "State_of_conservation"      # Coluna dos grupos (Alpha, Beta, Core, DESeq2)
GRUPO_REFERENCIA <- "Native"                     # Grupo de referência para as comparações

# --- 3.3. Parâmetros Gerais ---
SEED_GERAL       <- 123                          # Semente de reprodutibilidade

# --- 3.4. Parâmetros de Controle de Qualidade (QC) ---
COV_CORTE        <- 0.98                         # Corte para Good's Coverage (0.95 = 95%)
RARE_STEP        <- 100                          # Passo da curva de rarefação

# --- 3.5. Parâmetros de Abundância Diferencial (DESeq2) ---
DESEQ_ALPHA      <- 0.01                         # P-valor ajustado de corte (significância)
DESEQ_LFC        <- 7                            # Log2 Fold Change mínimo

# --- 3.6. Parâmetros do Core Microbiome ---
CORE_PREVALENCIA <- 0.90                         # % de amostras que o táxon deve aparecer
CORE_ABUNDANCIA  <- 0.001                        # Abundância relativa mínima
CORE_MICROBIOME_RANK <- "Family"                  # Nível taxonômico ("Phylum", "Family", "Genus" ou "ASV")

# --- 3.7. Parâmetros de Redes de Co-ocorrência ---
NET_TOP_N        <- 25                           # Número de táxons para montar a rede
NET_COR_CUTOFF   <- 0.6                          # Correlação mínima (Spearman rho)
NET_P_CUTOFF     <- 0.05                         # Limiar de significância (apenas p < 0.05)
NETWORK_RANK     <- "Genus"                      # Nível taxonômico para a Rede

# --- 3.8. Parâmetros de Venn e Bioindicadores ---
INDVAL_RANK      <- "Genus"                      # Nível taxonômico para a Análise de Espécies Indicadoras

# ==============================================================================
# ================= FIM DO PAINEL DE CONTROLE / INÍCIO DA EXECUÇÃO =============
# ==============================================================================


# ==============================================================================
# === # 4. IMPORTAÇÃO E CRIAÇÃO DO PHYLOSEQ OBJECT # === #
# ==============================================================================
{
  print("--- Iniciando Importação ---")
  
  # 4.1. Carrega tabelas (check.names = FALSE evita alteração de espaços originais)
  tabela_bruta <- read.delim(ARQUIVO_ABUNDANCIA, header = TRUE, check.names = FALSE) 
  metadados_df <- read.delim(ARQUIVO_METADADOS, header = TRUE, row.names = 1, check.names = FALSE)
  
  # Padronização dos nomes das colunas de metadados
  colnames(metadados_df) <- gsub(" ", "_", colnames(metadados_df))
  
  # 4.2. Limpeza
  tabela_limpa <- tabela_bruta %>% select(-total)
  
  tax_separada <- tabela_limpa %>%
    select(tax) %>% 
    separate(tax, 
             into = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Specie"), 
             sep = ";", 
             fill = "right") %>%
    as.matrix()
  
  # 4.3. Criação da Matriz Numérica
  otu_matriz <- tabela_limpa %>%
    select(-tax) %>%
    as.matrix()
  
  # Criar e atribuir IDs para as linhas
  ids_otas <- paste0("ASV_", 1:nrow(otu_matriz))
  rownames(otu_matriz) <- ids_otas
  rownames(tax_separada) <- ids_otas
  
  # 4.4. Montagem do Phyloseq
  OTU    <- otu_table(otu_matriz, taxa_are_rows = TRUE)
  TAX    <- tax_table(tax_separada)
  SAMPLE <- sample_data(metadados_df)
  
  ps <- phyloseq(OTU, TAX, SAMPLE)
  
  print("Objeto Phyloseq criado com sucesso:")
  print(ps)
}

# ==============================================================================
# === # 5. ANÁLISE DE COBERTURA (GOOD'S COVERAGE) # === #
# ==============================================================================
{
  # Função
  get_goods_coverage <- function(x) {
    n1 <- sum(x == 1)
    N <- sum(x)
    return(1 - (n1 / N))
  }
  
  # Cálculo
  cobertura <- apply(otu_table(ps), 2, get_goods_coverage)
  sample_data(ps)$Goods_Coverage <- cobertura
  
  dados_cov <- data.frame(sample_data(ps)) %>% 
    tibble::rownames_to_column("SampleID") 
  
  # Gráfico
  P_COV <- ggplot(dados_cov, aes(x = reorder(SampleID, Goods_Coverage), y = Goods_Coverage)) +
    geom_col(fill = "#4682B4", color = "black", width = 0.7) +
    geom_hline(yintercept = COV_CORTE, color = "red", linetype = "dashed", size = 0.8) +
    coord_flip() +
    labs(x = NULL, y = NULL) +
    theme(
      text = element_text(family = "serif"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"), 
      plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm"),
      panel.background = element_rect(fill = "gray90"),
      panel.grid.major = element_line(color = "white", size = 0.5),
      panel.grid.minor = element_line(color = "white", size = 0.25)
    )
  
  label_cov <- ggplot() + 
    annotate("text", x = 1, y = 1, label = "Good's Coverage", size = 10,
             fontface = "bold", family = "serif", hjust = 0.5) +
    xlim(0, 2) + ylim(0.5, 1.5) +
    theme_void() +
    theme(panel.background = element_rect(fill = "gray70", color = NA))
  
  grafico_final_cov <- label_cov / P_COV + plot_layout(heights = c(0.1, 1))
  print(grafico_final_cov)
  
  # Exportar TSV antes do gráfico
  write.table(dados_cov, "Tabela_Goods_Coverage.tsv", sep = "\t", quote = FALSE, row.names = FALSE, dec = ".")
  
  ggsave("Goods_Coverage.tiff", plot = grafico_final_cov, device = "tiff",
         width = 20, height = 10, units = "in", dpi = 600, bg = "white", compression = "lzw")
}

# ==============================================================================
# === # 6. CURVA DE RAREFAÇÃO # === #
# ==============================================================================
{
  ps_clean <- prune_samples(sample_sums(ps) > 0, ps)
  
  if (taxa_are_rows(ps_clean)) {
    otu_mat <- t(as(otu_table(ps_clean), "matrix"))
  } else {
    otu_mat <- as(otu_table(ps_clean), "matrix")
  }
  class(otu_mat) <- "matrix"
  
  # Cálculo usando a variável RARE_STEP
  out_rare <- rarecurve(otu_mat, step = RARE_STEP, sample = min(rowSums(otu_mat)), label = FALSE)
  
  names(out_rare) <- rownames(otu_mat)
  rare_df <- lapply(names(out_rare), function(x) {
    df <- data.frame(Reads = attr(out_rare[[x]], "Subsample"),
                     Richness = out_rare[[x]])
    df$SampleID <- x
    return(df)
  }) %>% bind_rows()
  
  dados_labels <- rare_df %>%
    group_by(SampleID) %>%
    summarise(Max_Reads = max(Reads), Max_Richness = max(Richness))
  
  # Gráfico
  P_RARE <- ggplot(rare_df, aes(x = Reads, y = Richness, group = SampleID)) +
    geom_line(color = "black", size = 0.6, alpha = 0.7) +
    geom_text(data = dados_labels, 
              aes(x = Max_Reads, y = Max_Richness, label = SampleID),
              color = "black", hjust = -0.1, size = 3, family = "serif", check_overlap = FALSE) + 
    geom_vline(xintercept = min(sample_sums(ps)), linetype = "dashed", color = "darkred", size = 0.5) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.2))) +
    labs(x = "Sequencing Depth (Reads)", y = "Observed Richness (ASVs)") +
    theme(
      text = element_text(family = "serif"),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 11, face = "bold"),
      legend.position = "none",
      plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm"),
      panel.background = element_rect(fill = "gray90"),
      panel.grid.major = element_line(color = "white", size = 0.5),
      panel.grid.minor = element_line(color = "white", size = 0.25)
    )
  
  label_rare <- ggplot() + 
    annotate("text", x = 1, y = 1, label = "Rarefaction Curve", size = 5,
             fontface = "bold", family = "serif", hjust = 0.5) +
    xlim(0, 2) + ylim(0.5, 1.5) +
    theme_void() +
    theme(panel.background = element_rect(fill = "gray70", color = NA))
  
  grafico_final_rare <- label_rare / P_RARE + plot_layout(heights = c(0.1, 1))
  print(grafico_final_rare)
  
  # Exportar TSV antes do gráfico
  write.table(rare_df, "Tabela_Rarefaction_Curve.tsv", sep = "\t", quote = FALSE, row.names = FALSE, dec = ".")
  
  ggsave("Rarefaction_Curve_BW.tiff", plot = grafico_final_rare, device = "tiff",
         width = 20, height = 10, units = "in", dpi = 600, bg = "white", compression = "lzw")
}

# ==============================================================================
# === # 7. ALPHA DIVERSIDADE # === #
# ==============================================================================
{
  set.seed(SEED_GERAL)
  
  min_reads <- min(sample_sums(ps))
  ps_rar <- rarefy_even_depth(ps, sample.size = min_reads, replace = FALSE, rngseed = SEED_GERAL)
  
  alpha_tab <- estimate_richness(ps_rar, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson"))
  alpha_tab$Pielou <- alpha_tab$Shannon / log(alpha_tab$Observed)
  alpha_tab <- alpha_tab %>% select(!ends_with("se"))
  
  meta <- data.frame(sample_data(ps_rar))
  alpha_final <- cbind(meta, alpha_tab)
  
  # Configuração de Cores
  grupos <- unique(alpha_final[[VAR_AGRUPAMENTO]])
  cores <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(grupos)), grupos)
  
  # Função de Plotagem Personalizada
  plot_alpha <- function(df, metric_col, title_text) {
    formula_stats <- as.formula(paste(metric_col, "~", VAR_AGRUPAMENTO))
    teste <- kruskal.test(formula_stats, data = df)
    p_valor <- teste$p.value
    
    if (p_valor < 0.001) { 
      p_texto <- "p < 0.001" 
    } else { 
      p_texto <- paste("p =", format(round(p_valor, 4), nsmall = 4)) 
    }
    
    titulo_completo <- paste0(title_text, " (", p_texto, ")")
    
    p <- ggplot(df, aes(x = .data[[VAR_AGRUPAMENTO]], y = .data[[metric_col]], fill = .data[[VAR_AGRUPAMENTO]])) +
      geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.8) +
      geom_jitter(width = 0.2, size = 3, shape = 21, fill = "red", color = "black", stroke = 0.3) +
      scale_fill_manual(values = cores) +
      labs(x = NULL, y = NULL) +
      theme(
        text = element_text(family = "serif"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 18, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        legend.position = "none",
        plot.margin = unit(c(0, 0.2, 0.2, 0.2), "cm"),
        panel.background = element_rect(fill = "gray90"),
        panel.grid.major = element_line(color = "white", size = 0.5),
        panel.grid.minor = element_line(color = "white", size = 0.25)
      )
    
    lbl <- ggplot() + 
      annotate("text", x = 1, y = 1, label = titulo_completo, size = 8, fontface = "bold", family = "serif") +
      theme_void() +
      theme(panel.background = element_rect(fill = "gray70", color = NA))
    
    return(lbl / p + plot_layout(heights = c(0.15, 1)))
  }
  
  # Geração e Combinação dos Gráficos
  p1 <- plot_alpha(alpha_final, "Observed", "Observed Richness")
  p2 <- plot_alpha(alpha_final, "Chao1", "Chao1")
  p3 <- plot_alpha(alpha_final, "ACE", "ACE")
  p4 <- plot_alpha(alpha_final, "Shannon", "Shannon Index")
  p5 <- plot_alpha(alpha_final, "Simpson", "Simpson Index (1-D)")
  p6 <- plot_alpha(alpha_final, "Pielou", "Pielou's Evenness (J)")
  
  painel_completo <- (p1 | p2 | p3) / (p4 | p5 | p6) + 
    plot_annotation(title = NULL, theme = theme(plot.title = element_text(family="serif", face="bold", size=16)))
  
  print(painel_completo)
  
  # Salvar Arquivos
  tabela_exportar <- alpha_final %>% tibble::rownames_to_column(var = "SampleID")
  write.table(tabela_exportar, file = "Tabela_Alpha_Kruskal_conservation.tsv", sep = "\t", quote = FALSE, row.names = FALSE, dec = ".")
  
  ggsave("Alpha_Kruskal_conservation.tiff", plot = painel_completo, device = "tiff",
         width = 24, height = 24, units = "in", dpi = 600, compression = "lzw", bg = "white")
}

# ==============================================================================
# === # 8. BETA DIVERSIDADE # === #
# ==============================================================================
{
  METODO_ORD <- "PCoA"
  DISTANCIA  <- "bray"
  
  # Cálculo Estatístico
  dist_matrix <- phyloseq::distance(ps_rar, method = DISTANCIA)
  meta_stats  <- data.frame(sample_data(ps_rar))
  
  set.seed(SEED_GERAL)
  permanova  <- adonis2(dist_matrix ~ meta_stats[[VAR_AGRUPAMENTO]], permutations = 999)
  p_val_beta <- permanova$`Pr(>F)`[1]
  
  if (p_val_beta < 0.001) { 
    p_text_beta <- "p < 0.001" 
  } else { 
    p_text_beta <- paste("p =", format(round(p_val_beta, 4), nsmall = 4)) 
  }
  
  # Ordenação e Preenchimento
  ord_pcoa   <- ordinate(ps_rar, method = METODO_ORD, distance = DISTANCIA)
  p_beta_obj <- plot_ordination(ps_rar, ord_pcoa, color = VAR_AGRUPAMENTO)
  
  df_beta <- p_beta_obj$data
  df_beta$SampleID <- rownames(df_beta)
  eixos_var <- round(ord_pcoa$values$Relative_eig[1:2] * 100, 1)
  
  # Configuração de Cores
  grupos_beta <- unique(df_beta[[VAR_AGRUPAMENTO]])
  cores_beta  <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(grupos_beta)), grupos_beta)
  
  # Plotagem
  P_BETA <- ggplot(df_beta, aes(x = Axis.1, y = Axis.2, color = .data[[VAR_AGRUPAMENTO]], fill = .data[[VAR_AGRUPAMENTO]])) +
    stat_ellipse(geom = "polygon", alpha = 0.2, type = "t", level = 0.95, show.legend = FALSE) +
    geom_point(size = 5, alpha = 0.9, shape = 21, color = "black", stroke = 0.5) +
    geom_text_repel(aes(label = SampleID), size = 5, family = "serif", color = "black", box.padding = 0.5, point.padding = 0.5, max.overlaps = 20, show.legend = FALSE) +
    scale_color_manual(values = cores_beta, name = "Group") +
    scale_fill_manual(values = cores_beta, name = "Group") +
    labs(x = paste0("PCoA 1 (", eixos_var[1], "%)"), y = paste0("PCoA 2 (", eixos_var[2], "%)"), title = NULL) +
    theme(
      text = element_text(family = "serif"),
      axis.title = element_text(size = 16, face = "bold", color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 18, face = "bold"),
      legend.position = "right",
      legend.background = element_rect(fill="gray95", color="gray80"),
      panel.background = element_rect(fill = "gray90"),
      panel.grid.major = element_line(color = "white", size = 0.6),
      panel.grid.minor = element_line(color = "white", size = 0.3),
      plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
    )
  
  # Título e Finalização
  titulo_beta <- paste0("Beta Diversity (Bray-Curtis)  |  PERMANOVA: ", p_text_beta)
  label_beta  <- ggplot() + 
    annotate("text", x = 1, y = 1, label = titulo_beta, size = 10, fontface = "bold", family = "serif", hjust = 0.5) +
    xlim(0, 2) + ylim(0.5, 1.5) +
    theme_void() +
    theme(panel.background = element_rect(fill = "gray70", color = NA))
  
  grafico_final_beta <- label_beta / P_BETA + plot_layout(heights = c(0.1, 1))
  print(grafico_final_beta)
  
  ggsave("Beta_Diversity_ane_PCoA_Bray_Labels-conservation.tiff", plot = grafico_final_beta, device = "tiff", width = 20, height = 18, units = "in", dpi = 600, bg = "white", compression = "lzw")
}

# ==============================================================================
# === # 9. COMPOSIÇÃO (BARRAS DE ABUNDÂNCIA RELATIVA) # === #
# ==============================================================================
{
  gerar_graficos_abundancia_final <- function(ps_obj, out_dir = ".") {
    ranks <- rank_names(ps_obj)
    print(paste("Níveis encontrados:", paste(ranks, collapse = ", ")))
    
    for (i in seq_along(ranks)) {
      nivel_atual <- ranks[i]
      print(paste("Processando nível:", nivel_atual, "..."))
      
      glom <- tax_glom(ps_obj, taxrank = nivel_atual)
      ps_rel <- transform_sample_counts(glom, function(x) x / sum(x) * 100)
      
      df_tax <- psmelt(ps_rel)
      df_tax[[nivel_atual]] <- as.character(df_tax[[nivel_atual]])
      df_tax[[nivel_atual]][is.na(df_tax[[nivel_atual]])] <- "Unclassified"
      
      resumo_tax <- df_tax %>% group_by(.data[[nivel_atual]]) %>% summarise(Mean = mean(Abundance)) %>% arrange(desc(Mean))
      top_20 <- resumo_tax %>% top_n(20, Mean) %>% pull(.data[[nivel_atual]])
      
      df_tax_renomeado <- df_tax %>% mutate(Taxon_Plot = ifelse(.data[[nivel_atual]] %in% top_20, .data[[nivel_atual]], "Others"))
      df_agregado <- df_tax_renomeado %>% group_by(Sample, Taxon_Plot) %>% summarise(Abundance = sum(Abundance), .groups = "drop")
      
      if ("Others" %in% df_agregado$Taxon_Plot) { 
        ordem_fatorial <- c(top_20, "Others") 
      } else { 
        ordem_fatorial <- top_20 
      }
      
      df_agregado$Taxon_Plot <- factor(df_agregado$Taxon_Plot, levels = ordem_fatorial)
      
      n_cores <- length(levels(df_agregado$Taxon_Plot))
      paleta_base <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"))
      cores_finais <- colorRampPalette(paleta_base)(n_cores)
      names(cores_finais) <- levels(df_agregado$Taxon_Plot)
      
      if ("Others" %in% names(cores_finais)) { cores_finais["Others"] <- "black" }
      
      P_ABUND <- ggplot(df_agregado, aes(x = Sample, y = Abundance, fill = Taxon_Plot)) +
        geom_bar(stat = "identity", width = 0.9, color = "black", size = 0.05) +
        scale_fill_manual(values = cores_finais) +
        labs(x = NULL, y = "Relative Abundance (%)", fill = nivel_atual) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 100.1)) +
        theme(
          text = element_text(family = "serif"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 14, face = "bold"),
          legend.text = element_text(size = 14, face = "italic"),
          legend.title = element_text(size = 16, face = "bold"),
          legend.key.size = unit(0.8, "cm"),
          legend.position = "right",
          panel.background = element_rect(fill = "gray90"),
          panel.grid.major.y = element_line(color = "white", size = 0.5),
          panel.grid.major.x = element_blank(),
          plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
        )
      
      titulo_texto <- paste("Relative Abundance:", nivel_atual)
      label_plot <- ggplot() + 
        annotate("text", x = 1, y = 1, label = titulo_texto, size = 8, fontface = "bold", family = "serif", hjust = 0.5) +
        xlim(0, 2) + ylim(0.5, 1.5) + theme_void() + theme(panel.background = element_rect(fill = "gray70", color = NA))
      
      grafico_final <- label_plot / P_ABUND + plot_layout(heights = c(0.1, 1))
      
      nome_arquivo <- paste0("Abundance_ane_", i, "_", nivel_atual, ".tiff")
      caminho_completo <- file.path(out_dir, nome_arquivo)
      ggsave(caminho_completo, plot = grafico_final, device = "tiff", width = 20, height = 14, units = "in", dpi = 600, compression = "lzw", bg = "white")
      
      print(paste("Salvo:", nome_arquivo))
    }
  }
  
  gerar_graficos_abundancia_final(ps)
}

# ==============================================================================
# === # 10. CORE MICROBIOME HEATMAP # === #
# ==============================================================================
{
  print(paste("Analisando Core Microbiome a nível de:", CORE_MICROBIOME_RANK))
  
  # 10.1. Prepara Dados (Abundância Relativa)
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  
  # 10.2. Filtra Core
  core_taxa <- filter_taxa(ps_rel, function(x) sum(x > CORE_ABUNDANCIA) > (CORE_PREVALENCIA * nsamples(ps_rel)), TRUE)
  print(paste("Número de táxons no Core (> ", CORE_PREVALENCIA * 100, "% prev):", ntaxa(core_taxa)))
  
  # 10.3. Agrupamento Taxonômico
  if (CORE_MICROBIOME_RANK == "ASV") {
    ps_core_glom <- core_taxa
  } else {
    ps_core_glom <- tax_glom(core_taxa, taxrank = CORE_MICROBIOME_RANK)
  }
  
  # 10.4. Tranforma para Dataframe
  df_core <- psmelt(ps_core_glom)
  df_core <- df_core %>% arrange(.data[[VAR_AGRUPAMENTO]])
  df_core$Sample <- factor(df_core$Sample, levels = unique(df_core$Sample))
  
  # 10.5. Ajuste dos Nomes
  if (CORE_MICROBIOME_RANK == "ASV") {
    df_core$Taxon_Label <- paste0(df_core$OTU, " (", df_core$Genus, ")")
  } else {
    df_core$Taxon_Label <- as.character(df_core[[CORE_MICROBIOME_RANK]])
    df_core$Taxon_Label[is.na(df_core$Taxon_Label)] <- "Unclassified"
  }
  
  # 10.6. Plotagem
  titulo_core <- paste0("Core Microbiome (", CORE_MICROBIOME_RANK, " > ", CORE_PREVALENCIA * 100, "%)")
  
  P_HEATMAP <- ggplot(df_core, aes(x = Sample, y = Taxon_Label, fill = log10(Abundance * 100 + 0.01))) + 
    geom_tile(color = "white", size = 0.1) +
    scale_fill_distiller(palette = "Spectral", direction = -1, name = "Log10(%)") +
    facet_grid(~ .data[[VAR_AGRUPAMENTO]], scales = "free_x", space = "free_x") +
    labs(x = NULL, y = paste(CORE_MICROBIOME_RANK, "(Core)"), title = NULL) +
    theme(
      text = element_text(family = "serif"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, color = "black"),
      axis.text.y = element_text(size = 10, face = "italic", color = "black"),
      legend.position = "right",
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = 10, face = "bold"),
      panel.background = element_rect(fill = "gray90", color = NA),
      panel.grid = element_blank(),
      plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
    )
  
  label_core <- ggplot() + 
    annotate("text", x = 1, y = 1, label = titulo_core, size = 6, fontface = "bold", family = "serif", hjust = 0.5) +
    xlim(0, 2) + ylim(0.5, 1.5) + theme_void() + theme(panel.background = element_rect(fill = "gray70", color = NA))
  
  grafico_final_core <- label_core / P_HEATMAP + plot_layout(heights = c(0.1, 1))
  print(grafico_final_core)
  
  # 10.7. Salvar Exportações
  tabela_core <- df_core %>% 
    group_by(Taxon_Label, Phylum) %>% 
    summarise(Mean_Abundance_Pct = mean(Abundance * 100), .groups = "drop") %>% 
    arrange(desc(Mean_Abundance_Pct))
  
  write.table(tabela_core, paste0("Tabela_Core_", CORE_MICROBIOME_RANK, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  
  altura_h <- max(6, length(unique(df_core$Taxon_Label)) * 0.3)
  ggsave(paste0("Core_Microbiome_", CORE_MICROBIOME_RANK, ".tiff"), plot = grafico_final_core, device = "tiff", 
         width = 12, height = altura_h, units = "in", dpi = 600, compression = "lzw", bg = "white")
}

# ==============================================================================
# === # 11. REDUNDANCY ANALYSIS (RDA) COM SELEÇÃO AUTOMÁTICA # === #
# ==============================================================================
{
  print("--- Executando Análise de Redundância (RDA) ---")
  
  # 11.1 Extração e Limpeza Direta do Phyloseq
  df_meta <- data.frame(sample_data(ps))
  
  # Converte vírgulas para pontos e transforma todas as colunas em numéricas (se possível)
  df_meta_limpo <- df_meta %>%
    mutate(across(everything(), ~ as.character(.))) %>%
    mutate(across(everything(), ~ str_replace_all(., ",", "."))) %>%
    type.convert(as.is = TRUE)
  
  # 11.2 Transformação da Matriz de Abundância
  abund_matrix <- t(as(otu_table(ps), "matrix"))
  abund_matrix <- abund_matrix[rownames(df_meta_limpo), ]
  abund_hellinger <- decostand(abund_matrix, method = "hellinger")
  
  # 11.3 Seleção Automática de TODAS as Variáveis Numéricas
  variaveis_ambientais <- df_meta_limpo %>%
    select(where(is.numeric)) %>%
    select(where(~ var(., na.rm = TRUE) > 0)) %>% 
    scale() %>% 
    as.data.frame()
  
  # Sincroniza amostras caso existam NAs nas variáveis
  amostras_comuns <- intersect(rownames(abund_hellinger), rownames(drop_na(variaveis_ambientais)))
  abund_hellinger <- abund_hellinger[amostras_comuns, ]
  variaveis_ambientais <- drop_na(variaveis_ambientais)
  
  # 11.4 Seleção de Variáveis (Forward Selection)
  print("Testando todas as variáveis para encontrar as estatisticamente significativas (p < 0.05)...")
  
  modelo_nulo <- rda(abund_hellinger ~ 1, data = variaveis_ambientais)
  modelo_cheio <- rda(abund_hellinger ~ ., data = variaveis_ambientais)
  
  rda_model <- ordistep(modelo_nulo, scope = formula(modelo_cheio), direction = "forward", permutations = 999, trace = FALSE)
  
  if(is.null(rda_model$CCA) || rda_model$CCA$rank == 0) {
    print("AVISO: Nenhuma variável foi significativa. Mostrando modelo completo para fins visuais.")
    rda_model <- rda(abund_hellinger ~ ., data = variaveis_ambientais)
  } else {
    vars_selecionadas <- names(rda_model$CCA$biplot[,1])
    print(paste("--> Sucesso! Variáveis significativas selecionadas:", paste(vars_selecionadas, collapse = ", ")))
  }
  
  set.seed(SEED_GERAL)
  anova_modelo <- anova.cca(rda_model, permutations = 999)
  p_modelo <- anova_modelo$`Pr(>F)`[1]
  
  # Extração de Variância e R2 Ajustado
  r2_adj <- RsquareAdj(rda_model)$adj.r.squared
  sum_rda <- summary(rda_model)
  rda1_var <- round(sum_rda$cont$importance[2, "RDA1"] * 100, 2)
  rda2_var <- round(sum_rda$cont$importance[2, "RDA2"] * 100, 2)
  
  # 11.5 Extração para o Gráfico
  sitios_scores <- as.data.frame(scores(rda_model, display = "sites")) %>%
    mutate(SampleID = rownames(.),
           Group = df_meta_limpo[rownames(.), VAR_AGRUPAMENTO])
  
  setas_scores <- as.data.frame(scores(rda_model, display = "bp")) %>%
    mutate(Variable = rownames(.)) %>%
    mutate(RDA1 = RDA1 * 2.5, RDA2 = RDA2 * 2.5) 
  
  # 11.6 Configuração de Cores Dinâmicas
  grupos_rda <- unique(sitios_scores$Group)
  cores_rda <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(grupos_rda)), grupos_rda)
  
  # 11.7 Geração do Gráfico 
  p_rda <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    
    geom_segment(data = setas_scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                 arrow = arrow(length = unit(0.3, "cm")), color = "blue", linewidth = 0.8) +
    geom_text_repel(data = setas_scores, aes(x = RDA1, y = RDA2, label = Variable),
                    color = "blue", fontface = "bold", family = "serif", size = 5) +
    
    geom_point(data = sitios_scores, aes(x = RDA1, y = RDA2, fill = Group), 
               size = 6, shape = 21, color = "black", stroke = 0.6) +
    
    # === ADICIONADO AQUI: Nomes das amostras puxando linhas com ggrepel ===
    geom_text_repel(data = sitios_scores, aes(x = RDA1, y = RDA2, label = SampleID),
                    size = 4.5, family = "serif", fontface = "italic", color = "black",
                    box.padding = 0.7, point.padding = 0.5, max.overlaps = 50, 
                    min.segment.length = 0, segment.color = "gray30") +
    # ======================================================================
  
  scale_fill_manual(values = cores_rda, name = VAR_AGRUPAMENTO) +
    labs(x = paste0("RDA 1 (", rda1_var, "%)"), y = paste0("RDA 2 (", rda2_var, "%)")) +
    theme(
      text = element_text(family = "serif"),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 12),
      panel.background = element_rect(fill = "gray90"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(color = "white", linewidth = 0.5),
      panel.grid.minor = element_line(color = "white", linewidth = 0.25)
    )
  
  # 11.8 Montagem Final
  titulo_texto <- paste0("Redundancy Analysis (RDA): Key Environmental Drivers\n",
                         "Model p = ", format(round(p_modelo, 4), nsmall = 4), 
                         " | Adjusted R² = ", round(r2_adj, 3))
  
  lbl <- ggplot() + 
    annotate("text", x = 1, y = 1, label = titulo_texto, size = 6, fontface = "bold", family = "serif") +
    theme_void() +
    theme(panel.background = element_rect(fill = "gray70", color = NA))
  
  PAINEL_RDA <- lbl / p_rda + plot_layout(heights = c(0.12, 1))
  print(PAINEL_RDA)
  
  # 11.9 Exportações
  write.table(sitios_scores, paste0("Tabela_RDA_Sites_Scores_", VAR_AGRUPAMENTO, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(setas_scores, paste0("Tabela_RDA_Arrow_Scores_", VAR_AGRUPAMENTO, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  ggsave(paste0("Figure_RDA_Selected_Drivers_", VAR_AGRUPAMENTO, ".tiff"), plot = PAINEL_RDA, 
         width = 12, height = 10, dpi = 600, compression = "lzw", bg = "white")
}

# ==============================================================================
# === # 12. ABUNDÂNCIA DIFERENCIAL (DESeq2) # === #
# ==============================================================================
{
  ps_deseq <- ps
  d <- sample_data(ps_deseq)
  
  if (any(is.na(d[[VAR_AGRUPAMENTO]]))) {
    print(paste("Aviso: Removendo amostras com NA em", VAR_AGRUPAMENTO))
    ps_deseq <- subset_samples(ps_deseq, !is.na(get(VAR_AGRUPAMENTO)))
  }
  
  # Define Referência Base
  sample_data(ps_deseq)[[VAR_AGRUPAMENTO]] <- factor(sample_data(ps_deseq)[[VAR_AGRUPAMENTO]])
  sample_data(ps_deseq)[[VAR_AGRUPAMENTO]] <- relevel(sample_data(ps_deseq)[[VAR_AGRUPAMENTO]], ref = GRUPO_REFERENCIA)
  
  otu_table(ps_deseq) <- otu_table(ps_deseq) + 1
  formula_deseq <- as.formula(paste("~", VAR_AGRUPAMENTO))
  diagdds <- phyloseq_to_deseq2(ps_deseq, formula_deseq)
  
  print(paste("Rodando DESeq2 para:", VAR_AGRUPAMENTO, "..."))
  diagdds <- DESeq(diagdds, test = "Wald", fitType = "parametric")
  results_deseq <- results(diagdds, cooksCutoff = FALSE)
  
  tax_table_df <- as.data.frame(tax_table(ps_deseq))
  res_tab <- cbind(as.data.frame(results_deseq), tax_table_df[rownames(results_deseq), ])
  
  # Exportar a Tabela Estatística COMPLETA (Antes dos Filtros)
  tabela_deseq_full <- res_tab %>% tibble::rownames_to_column(var = "ASV_ID") %>% arrange(padj)
  write.table(tabela_deseq_full, paste0("Tabela_DESeq2_FULL_", VAR_AGRUPAMENTO, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE, dec = ".")
  
  # Filtragem Estatística (Significativos)
  sigtab <- res_tab %>%
    filter(padj < DESEQ_ALPHA) %>%
    filter(abs(log2FoldChange) > DESEQ_LFC) %>%
    drop_na(Genus)
  
  if (nrow(sigtab) == 0) {
    stop("Nenhum biomarcador encontrado com esses cortes estatísticos.")
  } else {
    print(paste("Biomarcadores encontrados:", nrow(sigtab)))
  }
  
  niveis <- levels(sample_data(ps_deseq)[[VAR_AGRUPAMENTO]])
  NIVEL_TESTE <- setdiff(niveis, GRUPO_REFERENCIA)
  
  sigtab <- sigtab %>% mutate(Enriched_In = ifelse(log2FoldChange > 0, NIVEL_TESTE, GRUPO_REFERENCIA))
  sigtab$Genus <- factor(as.character(sigtab$Genus), levels = unique(as.character(sigtab$Genus[order(sigtab$log2FoldChange)])))
  
  # Configuração de Cores
  grupos_no_grafico <- unique(sigtab$Enriched_In)
  n_cores_necessarias <- length(grupos_no_grafico)
  paleta_cores <- brewer.pal(max(3, n_cores_necessarias), "Dark2")
  cores_grupos <- setNames(paleta_cores[1:n_cores_necessarias], grupos_no_grafico)
  
  # === ALTERAÇÃO: Rank explícito no Título ===
  titulo_texto <- paste("Differential Abundance (Genus):", VAR_AGRUPAMENTO)
  subtitulo_texto <- paste("Significant Biomarkers (p <", DESEQ_ALPHA, "| LFC >", DESEQ_LFC, ")")
  texto_completo_faixa <- paste0(titulo_texto, "\n", subtitulo_texto)
  # ==========================================
  
  # Plotagem
  P_DESEQ <- ggplot(sigtab, aes(x = Genus, y = log2FoldChange, color = Enriched_In)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.2, color = "black", alpha = 0.6) +
    geom_point(size = 4) +
    scale_color_manual(values = cores_grupos, name = "More abundant in:") +
    labs(y = "Log2 Fold Change", x = NULL) +
    coord_flip() +
    theme(
      text = element_text(family = "serif"),
      axis.text.y = element_text(size = 11, face = "italic", color = "black"),
      axis.title.x = element_text(face = "bold", size = 12, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      legend.position = "top", 
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 11),
      legend.background = element_rect(fill = "gray90", color = NA),
      legend.key = element_rect(fill = "gray90", color = NA),
      panel.background = element_rect(fill = "gray90", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "white", size = 0.6),
      panel.grid.minor = element_line(color = "white", size = 0.3),
      plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")
    )
  
  label_deseq <- ggplot() + 
    annotate("text", x = 1, y = 1, label = texto_completo_faixa, size = 5, fontface = "bold", family = "serif", hjust = 0.5, lineheight = 0.9) +
    xlim(0, 2) + ylim(0.5, 1.5) + theme_void() + theme(panel.background = element_rect(fill = "gray70", color = NA))
  
  grafico_final_deseq <- label_deseq / P_DESEQ + plot_layout(heights = c(0.15, 1))
  print(grafico_final_deseq)
  
  # Salvar Exportações Filtradas e Gráfico
  tabela_exportar <- sigtab %>% tibble::rownames_to_column(var = "ASV_ID") %>% arrange(desc(abs(log2FoldChange)))
  nome_tabela <- paste0("Tabela_DESeq2_Significant_log7_", VAR_AGRUPAMENTO, ".tsv")
  write.table(tabela_exportar, file = nome_tabela, sep = "\t", quote = FALSE, row.names = FALSE, dec = ".")
  
  altura_dinamica <- max(6, nrow(sigtab) * 0.25)
  ggsave(paste0("DiffAbund_DESeq2_log7_", VAR_AGRUPAMENTO, "_Final.tiff"), plot = grafico_final_deseq, device = "tiff", width = 14, height = altura_dinamica, units = "in", dpi = 600, compression = "lzw", bg = "white", limitsize = FALSE)
}

# ==============================================================================
# === # 13. REDES DE CO-OCORRÊNCIA # === #
# ==============================================================================
{
  grupos_rede <- unique(as.character(sample_data(ps)[[VAR_AGRUPAMENTO]]))
  grupos_rede <- grupos_rede[!is.na(grupos_rede)]
  
  for (grupo in grupos_rede) {
    
    print(paste("--- Gerando rede para o grupo:", grupo, "---"))
    
    # 13.1. Filtragem de Amostras
    ps_sub <- subset_samples(ps, get(VAR_AGRUPAMENTO) == grupo)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    
    # 13.2. Aglomeração
    ps_glom <- tax_glom(ps_sub, taxrank = NETWORK_RANK)
    top_taxa <- names(sort(taxa_sums(ps_glom), decreasing = TRUE)[1:min(NET_TOP_N, ntaxa(ps_glom))])
    ps_net <- prune_taxa(top_taxa, ps_glom)
    
    # 13.3. Matriz Numérica
    otu_net <- as.matrix(otu_table(ps_net))
    if (!taxa_are_rows(ps_net)) { otu_net <- t(otu_net) }
    
    # 13.4. Correlação (Spearman)
    cor_result <- Hmisc::rcorr(t(otu_net), type = "spearman")
    cor_matrix <- cor_result$r
    p_matrix <- cor_result$P
    
    # 13.5. Limiares
    cor_matrix[abs(cor_matrix) < NET_COR_CUTOFF] <- 0
    cor_matrix[p_matrix > NET_P_CUTOFF] <- 0
    diag(cor_matrix) <- 0
    
    if (sum(cor_matrix != 0) == 0) {
      print(paste("Aviso: Nenhuma correlação significativa atingiu os critérios para", grupo))
      next
    }
    
    # 13.6. Objeto igraph e Metadados
    net_obj <- graph_from_adjacency_matrix(cor_matrix, mode = "undirected", weighted = TRUE)
    
    tax_info <- as.data.frame(tax_table(ps_net))
    total_reads_rede <- sum(taxa_sums(ps_net))
    
    V(net_obj)$Abundance <- (taxa_sums(ps_net)[V(net_obj)$name] / total_reads_rede) * 100
    V(net_obj)$Phylum    <- tax_info[V(net_obj)$name, "Phylum"]
    V(net_obj)$TaxonName <- tax_info[V(net_obj)$name, NETWORK_RANK]
    
    E(net_obj)$Interaction <- E(net_obj)$weight
    E(net_obj)$Sign <- ifelse(E(net_obj)$Interaction > 0, "Positive", "Negative")
    E(net_obj)$weight <- abs(E(net_obj)$weight)
    
    # =======================================================================
    # CORREÇÃO: EXTRAÇÃO DIRETA DA TABELA DE ARESTAS (À PROVA DE ERROS)
    # =======================================================================
    nome_limpo <- gsub(" ", "_", grupo) 
    
    edge_df_raw <- igraph::as_data_frame(net_obj, what = "edges")
    edge_df_clean <- data.frame(
      Taxon_A = edge_df_raw$from,
      Taxon_B = edge_df_raw$to,
      Spearman_Rho = E(net_obj)$Interaction,
      Interaction_Type = E(net_obj)$Sign
    )
    write.table(edge_df_clean, paste0("Tabela_Network_Edges_", nome_limpo, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Tabela de Nós (Métricas Topológicas de cada Táxon)
    node_df <- igraph::as_data_frame(net_obj, what = "vertices") %>%
      mutate(Degree = degree(net_obj),                  
             Betweenness = betweenness(net_obj),        
             Closeness = closeness(net_obj)) %>%        
      arrange(desc(Degree))                             
    
    write.table(node_df, paste0("Tabela_Network_Nodes_Stats_", nome_limpo, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    # =======================================================================
    
    # 13.7. Configuração Visual
    filos_unicos <- unique(V(net_obj)$Phylum)
    cores_phylum <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(filos_unicos)), filos_unicos)
    
    titulo_texto <- paste("Co-occurrence Network:", grupo)
    subtitulo_texto <- paste("Top", NET_TOP_N, NETWORK_RANK, "| Spearman \u03C1 >", NET_COR_CUTOFF, "| p <", NET_P_CUTOFF)
    texto_completo_faixa <- paste0(titulo_texto, "\n", subtitulo_texto)
    
    # 13.8. Plotagem
    set.seed(SEED_GERAL) 
    
    P_NET <- ggraph(net_obj, layout = "fr") +
      geom_edge_link(aes(edge_width = abs(Interaction), edge_alpha = abs(Interaction), color = Sign)) +
      scale_edge_color_manual(values = c("Positive" = "green4", "Negative" = "red3"), name = "Correlation") +
      scale_edge_width(range = c(0.5, 2.5), guide = "none") +
      scale_edge_alpha(range = c(0.4, 0.8), guide = "none") +
      geom_node_point(aes(size = Abundance, fill = Phylum), shape = 21, color = "black", stroke = 0.5) +
      geom_node_text(aes(label = TaxonName), repel = TRUE, size = 3.5, fontface = "italic", family = "serif", color = "black") +
      scale_fill_manual(values = cores_phylum, name = "Phylum") +
      scale_size_continuous(range = c(4, 14), name = "Relative Abundance (%)") +
      guides(fill = guide_legend(override.aes = list(size = 8))) +
      theme_void() +
      theme(
        text = element_text(family = "serif"),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 11, family = "serif", color = "black"),
        legend.text = element_text(size = 10, family = "serif", color = "black"),
        legend.background = element_rect(fill = "gray90", color = NA),
        legend.key = element_rect(fill = "gray90", color = NA),
        legend.spacing.y = unit(0.3, "cm"),
        panel.background = element_rect(fill = "gray90", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        plot.margin = unit(c(0.2, 0.5, 0.5, 0.5), "cm")
      )
    
    label_net <- ggplot() + 
      annotate("text", x = 1, y = 1, label = texto_completo_faixa, size = 5, fontface = "bold", family = "serif", hjust = 0.5, lineheight = 0.9) +
      xlim(0, 2) + ylim(0.5, 1.5) + theme_void() + theme(panel.background = element_rect(fill = "gray70", color = NA))
    
    grafico_final_net <- label_net / P_NET + plot_layout(heights = c(0.12, 1))
    print(grafico_final_net)
    
    # 13.9. Salvar Gráfico Final
    nome_arquivo <- paste0("Network_", nome_limpo, ".tiff")
    ggsave(nome_arquivo, plot = grafico_final_net, device = "tiff", 
           width = 12, height = 10, units = "in", dpi = 600, compression = "lzw", bg = "white")
  }
}

# ==============================================================================
# === # 14. INDICATOR SPECIES ANALYSIS (IndVal) # === #
# ==============================================================================
{
  print(paste("--- Executando Análise de Espécies Indicadoras (IndVal) a nível de:", INDVAL_RANK, "---"))
  
  # 14.1. Aglomeração Taxonômica
  ps_indval <- tax_glom(ps, taxrank = INDVAL_RANK)
  
  # 14.2. Extração de Abundância e Metadados
  # O indicspecies requer que os táxons sejam colunas e as amostras linhas
  abund_matrix <- t(as(otu_table(ps_indval), "matrix"))
  grupos_indval <- as.character(sample_data(ps_indval)[[VAR_AGRUPAMENTO]])
  
  # 14.3. Execução do Modelo Multipatt
  # duleg = TRUE força o algoritmo a procurar bioindicadores de apenas UM grupo específico (sem combinações complexas)
  set.seed(SEED_GERAL)
  inv_model <- multipatt(abund_matrix, grupos_indval, func = "IndVal.g", duleg = TRUE, control = how(nperm = 999))
  
  # 14.4. Limpeza e Extração dos Resultados Estatísticos
  df_inv <- inv_model$sign %>%
    rownames_to_column("ASV_ID") %>%
    filter(p.value <= 0.05) %>%           # Filtra apenas os bioindicadores significativos
    arrange(desc(stat))                   # Ordena pelo maior Valor Indicador (stat)
  
  if(nrow(df_inv) == 0) {
    print("Aviso: Nenhuma espécie indicadora significativa encontrada com p <= 0.05.")
  } else {
    
    # Cruza o ID da tabela com o nome real da taxonomia
    tax_info <- as.data.frame(tax_table(ps_indval))
    df_inv$TaxonName <- tax_info[df_inv$ASV_ID, INDVAL_RANK]
    
    # Mapeia qual grupo aquela bactéria está indicando
    # O inv_model$comb guarda os nomes dos grupos. O index guarda qual grupo venceu.
    nomes_dos_grupos <- colnames(inv_model$comb)
    df_inv$Indicated_Group <- nomes_dos_grupos[df_inv$index]
    
    # Prepara a tabela de exportação removendo colunas residuais do algoritmo
    tabela_exportar <- df_inv %>%
      select(TaxonName, ASV_ID, Indicated_Group, Indicator_Value = stat, p_value = p.value)
    
    write.table(tabela_exportar, paste0("Tabela_IndVal_Bioindicators_", INDVAL_RANK, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE, dec = ".")
    
    # 14.5. Configuração de Cores
    grupos_presentes <- unique(df_inv$Indicated_Group)
    cores_indval <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(grupos_presentes)), grupos_presentes)
    
    # Limita aos Top 30 maiores bioindicadores para o gráfico não virar uma bagunça
    df_inv_plot <- head(df_inv, 30)
    
    # Ordena os fatores para o ggplot construir o gráfico do maior para o menor
    df_inv_plot$TaxonName <- factor(df_inv_plot$TaxonName, levels = rev(unique(df_inv_plot$TaxonName)))
    
    # 14.6. Plotagem (Lollipop Chart elegante)
    P_INDVAL <- ggplot(df_inv_plot, aes(x = stat, y = TaxonName, color = Indicated_Group)) +
      geom_segment(aes(x = 0, xend = stat, y = TaxonName, yend = TaxonName), color = "gray60", linewidth = 1) +
      geom_point(size = 5) +
      scale_color_manual(values = cores_indval, name = "Bioindicator of:") +
      labs(x = "Indicator Value (IndVal Stat)", y = NULL) +
      theme(
        text = element_text(family = "serif"),
        axis.text.y = element_text(size = 11, face = "italic", color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x = element_text(face = "bold", size = 12),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11),
        panel.background = element_rect(fill = "gray90"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.major.x = element_line(color = "white", linewidth = 0.6),
        panel.grid.major.y = element_blank(), # Remove linha Y para destacar o "Lollipop"
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.2, 0.5, 0.5, 0.5), "cm")
      )
    
    # 14.7. Montagem Final Estilo Patchwork
    titulo_texto <- paste("Indicator Species Analysis (IndVal):", INDVAL_RANK)
    subtitulo_texto <- "Top Significant Bioindicators (p < 0.05)"
    texto_completo_faixa <- paste0(titulo_texto, "\n", subtitulo_texto)
    
    label_indval <- ggplot() + 
      annotate("text", x = 1, y = 1, label = texto_completo_faixa, size = 5, fontface = "bold", family = "serif", hjust = 0.5, lineheight = 0.9) +
      xlim(0, 2) + ylim(0.5, 1.5) + theme_void() + theme(panel.background = element_rect(fill = "gray70", color = NA))
    
    grafico_final_indval <- label_indval / P_INDVAL + plot_layout(heights = c(0.12, 1))
    print(grafico_final_indval)
    
    # Altura dinâmica para acomodar a quantidade de táxons sem achatar o gráfico
    altura_dinamica <- max(6, nrow(df_inv_plot) * 0.25)
    ggsave(paste0("Indicator_Species_IndVal_", INDVAL_RANK, ".tiff"), plot = grafico_final_indval, device = "tiff", 
           width = 12, height = altura_dinamica, units = "in", dpi = 600, compression = "lzw", bg = "white")
  }
}
