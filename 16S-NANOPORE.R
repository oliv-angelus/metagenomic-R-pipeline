# ==============================================================================
# SCRIPT MESTRE: ANÁLISE DE METAGENÔMICA 16S (PHYLOSEQ)
# Autor: 
# Data: 2024-06
# ==============================================================================

# === # 1. INSTALAÇÃO DE PACOTES (Execute apenas na primeira vez) # === #

{install.packages("devtools")
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
  install.packages(c("igraph", "ggraph", "Hmisc"))
  #
  BiocManager::install("phyloseq")
  BiocManager::install("ape")
  BiocManager::install("DESeq2")}

# === # 2. CARREGAR BIBLIOTECAS # === #

{library(phyloseq)
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
  library(igraph)
  library(ggraph)
  library(Hmisc)}

# ==============================================================================
# === # 3. DEFINIÇÃO DE VARIÁVEIS E PARÂMETROS (PAINEL DE CONTROLE) # === #
# ==============================================================================
# ALTERE AQUI PARA MUDAR A ANÁLISE INTEIRA SEM MEXER NO CÓDIGO ABAIXO

# --- 3.1. Arquivos de Entrada ---
ARQUIVO_ABUNDANCIA <- "ABUND.tsv"  # Nome do arquivo da tabela de abundância (TSV)
ARQUIVO_METADADOS  <- "DATA.tsv"         # Nome do arquivo de metadados (TSV)

# --- 3.2. Variáveis Experimentais ---
# Nome exato da coluna no metadado que define seus grupos (ex: "Anthropization", "Local", "Tratamento")
VAR_AGRUPAMENTO    <- "Beach_type"       # Usado em: Alpha, Beta, Core, DESeq2
# Nome do grupo que serve de CONTROLE ou REFERÊNCIA (ex: "Native", "Control", "Zero")
GRUPO_REFERENCIA   <- "island"               # Usado em: DESeq2 (para calcular o FoldChange contra ele)

# --- 3.3. Parâmetros Gerais ---
SEED_GERAL         <- 123                    # Semente para reprodutibilidade (garante mesmos resultados sempre)

# --- 3.4. Parâmetros de Controle de Qualidade (QC) ---
COV_CORTE          <- 0.95                   # Corte para Good's Coverage (0.95 = 95%)
RARE_STEP          <- 100                    # Passo da curva de rarefação (diminuir se tiver poucos reads, aumentar se tiver muitos)

# --- 3.5. Parâmetros de Abundância Diferencial (DESeq2) ---
DESEQ_ALPHA        <- 0.05                   # P-valor ajustado de corte (significância)
DESEQ_LFC          <- 5                      # Log2 Fold Change mínimo (2 = 4x mudança; 5 = 32x mudança - mais rigoroso)

# --- 3.6. Parâmetros do Core Microbiome ---
CORE_PREVALENCIA   <- 0.90                   # % de amostras que o táxon deve aparecer (0.75 = 75%)
CORE_ABUNDANCIA    <- 0.001                  # Abundância relativa mínima para ser considerado "presente" (0.001 = 0.1%)
CORE_RANK          <- "Genus"                # Nível taxonômico: "Phylum", "Family", "Genus" ou "ASV"

# --- 3.7. Parâmetros de Redes de Co-ocorrência ---
NET_TOP_N          <- 25                     # Número de gêneros mais abundantes para montar a rede (evita "hairball")
NET_COR_CUTOFF     <- 0.6                    # Correlação mínima (Spearman rho) para criar uma conexão
NET_P_CUTOFF       <- 0.05                   # P-valor máximo para a correlação ser aceita

# ==============================================================================
# FIM DA DEFINIÇÃO DE VARIÁVEIS - INÍCIO DO PROCESSAMENTO
# ==============================================================================

# === # 4. e 5. IMPORTAÇÃO, LIMPEZA E CRIAÇÃO DO PHYLOSEQ OBJECT # === #

{
  
  print("--- Iniciando Importação ---")
  
  # 4.1. Carrega tabelas usando as variáveis definidas
  tabela_bruta <- read.delim(ARQUIVO_ABUNDANCIA, header = TRUE, check.names = FALSE) 
  metadados_df <- read.delim(ARQUIVO_METADADOS, header = TRUE, row.names = 1)
  
  # 4.2. Limpeza
  # Remove a coluna total
  tabela_limpa <- tabela_bruta %>% select(-total)
  
  # Separa a coluna de taxonomia
  tax_separada <- tabela_limpa %>%
    select(tax) %>% 
    separate(tax, 
             into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Specie"), 
             sep = ";", 
             fill = "right") %>%
    as.matrix()
  
  # 4.3. Criação da Matriz Numérica
  otu_matriz <- tabela_limpa %>%
    select(-tax) %>%
    as.matrix()
  
  # Criar IDs para as linhas e atribuir
  ids_otas <- paste0("ASV_", 1:nrow(otu_matriz))
  rownames(otu_matriz) <- ids_otas
  rownames(tax_separada) <- ids_otas
  
  # === # MONTAGEM DO OBJETO PHYLOSEQ # === #
  
  OTU <- otu_table(otu_matriz, taxa_are_rows = TRUE)
  TAX <- tax_table(tax_separada)
  SAMPLE <- sample_data(metadados_df)
  
  ps <- phyloseq(OTU, TAX, SAMPLE)
  
  print("Objeto Phyloseq criado:")
  print(ps)
  }

# === # 6. ANÁLISE DE COBERTURA (GOOD'S COVERAGE) # === #

{# Função
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
    geom_hline(yintercept = COV_CORTE, color = "red", linetype = "dashed", size = 0.8) + # Usa variável COV_CORTE
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
  
  ggsave("Goods_Coverage.tiff", plot = grafico_final_cov, device = "tiff",
         width = 20, height = 10, units = "in", dpi = 600, bg = "white", compression = "lzw")
}

# === # 7. CURVA DE RAREFAÇÃO # === #

{ps_clean <- prune_samples(sample_sums(ps) > 0, ps)
  
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
  
  ggsave("Rarefaction_Curve_BW.tiff", plot = grafico_final_rare, device = "tiff",
         width = 20, height = 10, units = "in", dpi = 600, bg = "white", compression = "lzw")
}

# === # 8. ALPHA DIVERSIDADE # === #

{set.seed(SEED_GERAL) # Usa a variável SEED_GERAL
  min_reads <- min(sample_sums(ps))
  ps_rar <- rarefy_even_depth(ps, sample.size = min_reads, replace = FALSE, rngseed = SEED_GERAL)
  
  alpha_tab <- estimate_richness(ps_rar, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson"))
  alpha_tab$Pielou <- alpha_tab$Shannon / log(alpha_tab$Observed)
  alpha_tab <- alpha_tab %>% select(!ends_with("se"))
  
  meta <- data.frame(sample_data(ps_rar))
  alpha_final <- cbind(meta, alpha_tab)
  
  # Configuração de Cores usando VAR_AGRUPAMENTO
  grupos <- unique(alpha_final[[VAR_AGRUPAMENTO]])
  cores <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(grupos)), grupos)
  
  # Função Plot
  plot_alpha <- function(df, metric_col, title_text) {
    # Usa VAR_AGRUPAMENTO globalmente definido
    formula_stats <- as.formula(paste(metric_col, "~", VAR_AGRUPAMENTO))
    teste <- kruskal.test(formula_stats, data = df)
    p_valor <- teste$p.value
    
    if (p_valor < 0.001) { p_texto <- "p < 0.001" } else { p_texto <- paste("p =", format(round(p_valor, 4), nsmall = 4)) }
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
  
  # Geração dos Gráficos
  p1 <- plot_alpha(alpha_final, "Observed", "Observed Richness")
  p2 <- plot_alpha(alpha_final, "Chao1", "Chao1")
  p3 <- plot_alpha(alpha_final, "ACE", "ACE")
  p4 <- plot_alpha(alpha_final, "Shannon", "Shannon Index")
  p5 <- plot_alpha(alpha_final, "Simpson", "Simpson Index (1-D)")
  p6 <- plot_alpha(alpha_final, "Pielou", "Pielou's Evenness (J)")
  
  painel_completo <- (p1 | p2 | p3) / (p4 | p5 | p6) + 
    plot_annotation(title = NULL, theme = theme(plot.title = element_text(family="serif", face="bold", size=16)))
  
  print(painel_completo)
  
  ggsave("Alpha_Diversity_Selected_Kruskal.tiff", plot = painel_completo, device = "tiff",
         width = 24, height = 24, units = "in", dpi = 600, compression = "lzw", bg = "white")
  
  # Exportar TSV
  tabela_exportar <- alpha_final %>% tibble::rownames_to_column(var = "SampleID")
  write.table(tabela_exportar, file = "Tabela_Alpha_Diversidade.tsv", sep = "\t", quote = FALSE, row.names = FALSE, dec = ".")
}

# === # 9. BETA DIVERSIDADE # === #

{METODO_ORD <- "PCoA"
  DISTANCIA <- "bray"
  
  # Estatística
  dist_matrix <- phyloseq::distance(ps_rar, method = DISTANCIA)
  meta_stats <- data.frame(sample_data(ps_rar))
  
  set.seed(SEED_GERAL)
  # Usa VAR_AGRUPAMENTO
  permanova <- adonis2(dist_matrix ~ meta_stats[[VAR_AGRUPAMENTO]], permutations = 999)
  p_val_beta <- permanova$`Pr(>F)`[1]
  
  if (p_val_beta < 0.001) { p_text_beta <- "p < 0.001" } else { p_text_beta <- paste("p =", format(round(p_val_beta, 4), nsmall = 4)) }
  
  # Ordenação
  ord_pcoa <- ordinate(ps_rar, method = METODO_ORD, distance = DISTANCIA)
  p_beta_obj <- plot_ordination(ps_rar, ord_pcoa, color = VAR_AGRUPAMENTO)
  df_beta <- p_beta_obj$data
  df_beta$SampleID <- rownames(df_beta)
  eixos_var <- round(ord_pcoa$values$Relative_eig[1:2] * 100, 1)
  
  # Cores
  grupos_beta <- unique(df_beta[[VAR_AGRUPAMENTO]])
  cores_beta <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(grupos_beta)), grupos_beta)
  
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
  
  titulo_beta <- paste0("Beta Diversity (Bray-Curtis)  |  PERMANOVA: ", p_text_beta)
  label_beta <- ggplot() + 
    annotate("text", x = 1, y = 1, label = titulo_beta, size = 10, fontface = "bold", family = "serif", hjust = 0.5) +
    xlim(0, 2) + ylim(0.5, 1.5) +
    theme_void() +
    theme(panel.background = element_rect(fill = "gray70", color = NA))
  
  grafico_final_beta <- label_beta / P_BETA + plot_layout(heights = c(0.1, 1))
  print(grafico_final_beta)
  ggsave("Beta_Diversity_PCoA_Bray_Labels.tiff", plot = grafico_final_beta, device = "tiff", width = 20, height = 18, units = "in", dpi = 600, bg = "white", compression = "lzw")
}

# === # 10. COMPOSIÇÃO (BARRAS DE ABUNDÂNCIA) # === #

{gerar_graficos_abundancia_final <- function(ps_obj, out_dir = ".") {
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
    
    if ("Others" %in% df_agregado$Taxon_Plot) { ordem_fatorial <- c(top_20, "Others") } else { ordem_fatorial <- top_20 }
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
    
    nome_arquivo <- paste0("Abundance_", i, "_", nivel_atual, ".tiff")
    caminho_completo <- file.path(out_dir, nome_arquivo)
    ggsave(caminho_completo, plot = grafico_final, device = "tiff", width = 20, height = 14, units = "in", dpi = 600, compression = "lzw", bg = "white")
    print(paste("Salvo:", nome_arquivo))
  }
}
  gerar_graficos_abundancia_final(ps)
}

# === # 11. ABUNDÂNCIA DIFERENCIAL (DESeq2) # === #

{ps_deseq <- ps
  d <- sample_data(ps_deseq)
  
  # Usa VAR_AGRUPAMENTO (que foi definido como VAR_INTERESSE no seu pedido anterior)
  if (any(is.na(d[[VAR_AGRUPAMENTO]]))) {
    print(paste("Aviso: Removendo amostras com NA em", VAR_AGRUPAMENTO))
    ps_deseq <- subset_samples(ps_deseq, !is.na(get(VAR_AGRUPAMENTO)))
  }
  
  # Define Referência (GRUPO_REFERENCIA)
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
  
  # Filtra usando variáveis DESEQ_ALPHA e DESEQ_LFC
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
  
  grupos_no_grafico <- unique(sigtab$Enriched_In)
  n_cores_necessarias <- length(grupos_no_grafico)
  paleta_cores <- brewer.pal(max(3, n_cores_necessarias), "Dark2")
  cores_grupos <- setNames(paleta_cores[1:n_cores_necessarias], grupos_no_grafico)
  
  titulo_texto <- paste("Differential Abundance:", VAR_AGRUPAMENTO)
  subtitulo_texto <- paste("Significant Biomarkers (p <", DESEQ_ALPHA, "| LFC >", DESEQ_LFC, ")")
  texto_completo_faixa <- paste0(titulo_texto, "\n", subtitulo_texto)
  
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
  
  altura_dinamica <- max(6, nrow(sigtab) * 0.25)
  ggsave(paste0("DiffAbund_DESeq2_", VAR_AGRUPAMENTO, "_Final.tiff"), plot = grafico_final_deseq, device = "tiff", width = 14, height = altura_dinamica, units = "in", dpi = 600, compression = "lzw", bg = "white")
  
  # Exportar TSV
  tabela_exportar <- sigtab %>% tibble::rownames_to_column(var = "ASV_ID") %>% arrange(desc(abs(log2FoldChange)))
  nome_tabela <- paste0("Tabela_Biomarcadores_", VAR_AGRUPAMENTO, ".tsv")
  write.table(tabela_exportar, file = nome_tabela, sep = "\t", quote = FALSE, row.names = FALSE, dec = ".")
}

# === # 12. CORE MICROBIOME HEATMAP # === #

{print(paste("Analisando Core Microbiome a nível de:", CORE_RANK))
  
  # 1. Preparar dados (Abundância Relativa)
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  
  # 2. Filtrar Core (Prevalência e Abundância)
  core_taxa <- filter_taxa(ps_rel, function(x) sum(x > CORE_ABUNDANCIA) > (CORE_PREVALENCIA * nsamples(ps_rel)), TRUE)
  
  print(paste("Número de táxons no Core (> ", CORE_PREVALENCIA*100, "% prev):", ntaxa(core_taxa)))
  
  # 3. Agrupamento Taxonômico (Lógica Inteligente)
  if (CORE_RANK == "ASV") {
    # Se for ASV, não fazemos glom (agrupamento), usamos direto
    ps_core_glom <- core_taxa
  } else {
    # Se for Genus, Family, etc, agrupamos
    ps_core_glom <- tax_glom(core_taxa, taxrank = CORE_RANK)
  }
  
  # 4. Transformar para Dataframe
  df_core <- psmelt(ps_core_glom)
  
  # Ordenar amostras
  df_core <- df_core %>% arrange(.data[[VAR_AGRUPAMENTO]])
  df_core$Sample <- factor(df_core$Sample, levels = unique(df_core$Sample))
  
  # 5. Tratamento de Nomes para o Eixo Y
  # Se for ASV, queremos ver o nome da ASV. Se for Gênero, o nome do Gênero.
  if (CORE_RANK == "ASV") {
    # Cria um nome composto: ASV_1 (Genus)
    df_core$Taxon_Label <- paste0(df_core$OTU, " (", df_core$Genus, ")")
  } else {
    # Usa o nível escolhido
    df_core$Taxon_Label <- as.character(df_core[[CORE_RANK]])
    df_core$Taxon_Label[is.na(df_core$Taxon_Label)] <- "Unclassified"
  }
  
  # 6. Gráfico
  titulo_core <- paste0("Core Microbiome (", CORE_RANK, " > ", CORE_PREVALENCIA*100, "%)")
  
  P_HEATMAP <- ggplot(df_core, aes(x = Sample, y = Taxon_Label, fill = log10(Abundance * 100 + 0.01))) + 
    geom_tile(color = "white", size = 0.1) +
    scale_fill_distiller(palette = "Spectral", direction = -1, name = "Log10(%)") +
    facet_grid(~ .data[[VAR_AGRUPAMENTO]], scales = "free_x", space = "free_x") +
    labs(x = NULL, y = paste(CORE_RANK, "(Core)"), title = NULL) +
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
  
  # 7. Salvar
  altura_h <- max(6, length(unique(df_core$Taxon_Label)) * 0.3)
  ggsave(paste0("Core_Microbiome_", CORE_RANK, ".tiff"), plot = grafico_final_core, device = "tiff", 
         width = 12, height = altura_h, units = "in", dpi = 600, compression = "lzw", bg = "white")
  
  # 8. Tabela
  tabela_core <- df_core %>% 
    group_by(Taxon_Label, Phylum) %>% 
    summarise(Mean_Abundance_Pct = mean(Abundance * 100), .groups = "drop") %>% 
    arrange(desc(Mean_Abundance_Pct))
  
  write.table(tabela_core, paste0("Tabela_Core_", CORE_RANK, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)}

# === # 13. REDES DE CO-OCORRÊNCIA # === #

{# Usa variáveis NET_TOP_N, NET_COR_CUTOFF, NET_P_CUTOFF
  ps_genus <- tax_glom(ps, "Genus")
  ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))
  top_taxa <- names(sort(taxa_sums(ps_rel), decreasing = TRUE)[1:NET_TOP_N])
  ps_net <- prune_taxa(top_taxa, ps_rel)
  
  otu_mat <- t(as(otu_table(ps_net), "matrix"))
  res_cor <- rcorr(otu_mat, type = "spearman")
  corr_r <- res_cor$r
  corr_p <- res_cor$P
  
  flatten_corr_matrix <- function(cormat, pmat) {
    ut <- upper.tri(cormat)
    data.frame(row = rownames(cormat)[row(cormat)[ut]], column = rownames(cormat)[col(cormat)[ut]], cor = cormat[ut], p = pmat[ut])
  }
  
  edge_list <- flatten_corr_matrix(corr_r, corr_p)
  edge_list_final <- edge_list %>% filter(p < NET_P_CUTOFF) %>% filter(abs(cor) > NET_COR_CUTOFF) %>% mutate(Interaction = ifelse(cor > 0, "Positive", "Negative"))
  
  print(paste("Conexões mantidas:", nrow(edge_list_final)))
  
  net <- graph_from_data_frame(edge_list_final, directed = FALSE)
  net <- delete.vertices(net, which(degree(net) == 0))
  
  tax_info <- as.data.frame(tax_table(ps_net))
  tax_info$Genus[is.na(tax_info$Genus)] <- "Unclassified"
  mean_abund <- taxa_sums(ps_net) / nsamples(ps_net)
  
  V(net)$Phylum <- tax_info[V(net)$name, "Phylum"]
  V(net)$Genus  <- tax_info[V(net)$name, "Genus"]
  V(net)$Size   <- mean_abund[V(net)$name]
  V(net)$Degree <- degree(net)
  
  filos_presentes <- unique(V(net)$Phylum)
  cores_filos <- setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(filos_presentes)), filos_presentes)
  cores_arestas <- c("Positive" = "green", "Negative" = "red")
  
  titulo_net <- paste("Microbial Co-occurrence Network (Top", NET_TOP_N, "Genera)")
  subtitulo_net <- paste("Spearman rho >", NET_COR_CUTOFF, "| p <", NET_P_CUTOFF)
  
  P_NET <- ggraph(net, layout = "fr") + 
    geom_edge_link(aes(color = Interaction, width = abs(cor)), alpha = 0.6) +
    scale_edge_width(range = c(0.5, 2), guide = "none") + 
    scale_edge_color_manual(values = cores_arestas) +
    geom_node_point(aes(size = Size, fill = Phylum), shape = 21, color = "black", stroke = 0.5) +
    scale_size_continuous(range = c(3, 10), name = "Abundance") +
    scale_fill_manual(values = cores_filos) +
    guides(fill = guide_legend(override.aes = list(size = 8), title = "Phylum")) +
    geom_node_text(aes(label = Genus), repel = TRUE, size = 5, family = "serif", fontface = "italic") +
    labs(x = NULL, y = NULL) +
    theme_void() + 
    theme(
      text = element_text(family = "serif"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 15),
      legend.text = element_text(size = 12),
      panel.background = element_rect(fill = "gray90", color = NA),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  
  label_net <- ggplot() + 
    annotate("text", x = 1, y = 1, label = paste(titulo_net, "\n", subtitulo_net), size = 5, fontface = "bold", family = "serif", hjust = 0.5, lineheight = 0.9) +
    xlim(0, 2) + ylim(0.5, 1.5) + theme_void() + theme(panel.background = element_rect(fill = "gray70", color = NA))
  
  grafico_final_net <- label_net / P_NET + plot_layout(heights = c(0.15, 1))
  print(grafico_final_net)
  
  ggsave("Network_CoOccurrence.tiff", plot = grafico_final_net, device = "tiff", width = 16, height = 14, units = "in", dpi = 600, compression = "lzw", bg = "white")
  
  hubs_df <- data.frame(Genus = V(net)$Genus, Phylum = V(net)$Phylum, Degree = V(net)$Degree, Abundance = V(net)$Size) %>% arrange(desc(Degree))
  write.table(hubs_df, "Tabela_Network_Hubs.tsv", sep="\t", quote=FALSE, row.names=FALSE)
  print("Análise Finalizada com Sucesso!")}

                                    # --- 1. IMPORTAÇÃO E LIMPEZA ---
df_meta <- read.delim("DATA.tsv", sep = "\t", header = TRUE)
rownames(df_meta) <- df_meta$Sample

df_meta_limpo <- df_meta %>%
  mutate(across(pH:alcalinity, ~ as.numeric(str_replace_all(as.character(.), ",", ".")))) %>%
  mutate(Beach_type = str_to_title(Beach_type))

abund_matrix <- t(as(otu_table(ps), "matrix"))
abund_matrix <- abund_matrix[rownames(df_meta_limpo), ]

# --- 2. TRANSFORMAÇÃO E SELEÇÃO A PRIORI ---
abund_hellinger <- decostand(abund_matrix, method = "hellinger")

# SOLUÇÃO: Selecionamos apenas as 5 variáveis que contam a sua história ecológica!
variaveis_ambientais <- df_meta_limpo %>%
  select(pH, temperature, salinity, dissolved_oxygen, turbidity) %>%
  scale() %>% 
  as.data.frame()

# --- 3. EXECUÇÃO DA RDA E TESTE GLOBAL ---
rda_model <- rda(abund_hellinger ~ ., data = variaveis_ambientais)

# Teste global do modelo
set.seed(123)
anova_modelo <- anova.cca(rda_model, permutations = 999)
p_modelo <- anova_modelo$`Pr(>F)`[1]

# Extração de variância e R2 Ajustado
r2_adj <- RsquareAdj(rda_model)$adj.r.squared
sum_rda <- summary(rda_model)
rda1_var <- round(sum_rda$cont$importance[2, "RDA1"] * 100, 2)
rda2_var <- round(sum_rda$cont$importance[2, "RDA2"] * 100, 2)

# --- 4. EXTRAÇÃO PARA O GRÁFICO ---
sitios_scores <- as.data.frame(scores(rda_model, display = "sites")) %>%
  mutate(SampleID = rownames(.),
         Group = df_meta_limpo$Beach_type)

setas_scores <- as.data.frame(scores(rda_model, display = "bp")) %>%
  mutate(Variable = rownames(.)) %>%
  # Retiramos o filtro! Agora as 5 setas ecológicas vão aparecer.
  mutate(RDA1 = RDA1 * 2.5, RDA2 = RDA2 * 2.5) 

# --- 5. GERAÇÃO DO GRÁFICO (Estilo Alfa) ---
cores_projeto <- c("Urban" = "#d95f02", "Island" = "#1b9e77")

p_rda <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  
  geom_segment(data = setas_scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), color = "blue", linewidth = 0.8) +
  geom_text_repel(data = setas_scores, aes(x = RDA1, y = RDA2, label = Variable),
                  color = "blue", fontface = "bold", family = "serif", size = 5) +
  
  geom_point(data = sitios_scores, aes(x = RDA1, y = RDA2, fill = Group), 
             size = 6, shape = 21, color = "black", stroke = 0.6) +
  
  scale_fill_manual(values = cores_projeto) +
  labs(x = paste0("RDA 1 (", rda1_var, "%)"), y = paste0("RDA 2 (", rda2_var, "%)")) +
  theme(
    text = element_text(family = "serif"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.background = element_rect(fill = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "white", linewidth = 0.5),
    panel.grid.minor = element_line(color = "white", linewidth = 0.25)
  )

# --- 6. MONTAGEM FINAL ---
titulo_texto <- paste0("Redundancy Analysis (RDA): Key Environmental Drivers\n",
                       "Model p = ", format(round(p_modelo, 4), nsmall = 4), 
                       " | Adjusted R² = ", round(r2_adj, 3))

lbl <- ggplot() + 
  annotate("text", x = 1, y = 1, label = titulo_texto, size = 6, fontface = "bold", family = "serif") +
  theme_void() +
  theme(panel.background = element_rect(fill = "gray70", color = NA))

PAINEL_RDA <- lbl / p_rda + plot_layout(heights = c(0.12, 1))

print(PAINEL_RDA)

ggsave("Figure_RDA_Selected_Drivers.tiff", plot = PAINEL_RDA, 
       width = 12, height = 10, dpi = 600, compression = "lzw", bg = "white")
