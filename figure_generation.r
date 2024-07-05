library(tidyverse)
library(ggsignif)
library(forcats)
library(networkD3)
library(RColorBrewer)
library(viridis)

full_cell_metadata <- read.csv("B:/TFM/AMP_PD_cluster/04-meta_data/full_cell_metadata.tsv", row.names=1)
corrected_full_cell_metadata <- read.csv("B:/TFM/AMP_PD_cluster/04-meta_data/corrected_full_cell_metadata.tsv", row.names=1)
corrected_full_cell_metadata$brain_region <- as.factor(corrected_full_cell_metadata$brain_region)

patient_metadata <- corrected_full_cell_metadata %>% 
  select(participant_id, cohort, case_control, diagnosis_latest, age_at_baseline, sex, ethnicity, race, 
         education_level_years, hoehn_yahr_stage, path_braak_lb, 
         path_braak_nft) %>% 
  filter(!duplicated(.))


full_cell_metadata <- full_cell_metadata %>% 
  mutate(participant_id = fct_reorder(participant_id, 
                                      as.numeric(as.factor(case_control_other_latest)))) %>% 
  mutate(brain_region = as.factor(brain_region)) 

#### STATS #####
# Contar número de casos, controles, other
rownames(patient_metadata) <- patient_metadata$participant_id
patient_metadata[1] <- NULL
table(patient_metadata$diagnosis_latest)  

# Wilcox EDAD
patient_metadata$participant_id <- rownames(patient_metadata)
patient_metadata <- patient_metadata %>%
  full_join(total_genes_counts, by = "participant_id")
wilcox.test(patient_metadata$age_at_baseline[patient_metadata$case_control=="Control"], 
            patient_metadata$age_at_baseline[!(patient_metadata$case_control == "Control")])

write.csv(
  patient_metadata,
  "04-meta_data/patient_metadata.tsv",
  sep = "\t",
  row.names = T,
  col.names = T,
)

# IQR y media counts por muestra
(IQR_mean_counts <- full_cell_metadata %>%
  summarise(
    Q1_nCount = quantile(nCount_RNA, 0.25, na.rm = TRUE),
    Q3_nCount = quantile(nCount_RNA, 0.75, na.rm = TRUE),
    mean_nCount = mean(nCount_RNA, na.rm = TRUE)
  ))

# IQR y media counts por muestra (por caso o control)
# Recodificar "Other" como "Case"
full_cell_metadata <- full_cell_metadata %>%
  mutate(only_case_control = ifelse(case_control_other_latest == "Other", "PD", case_control_other_latest)) %>% 
  mutate(only_case_control = ifelse(case_control_other_latest == "Case", "PD", only_case_control))
(IQR_mean_counts_cases_control <- full_cell_metadata %>%
    group_by(only_case_control) %>% 
    summarise(
      Q1_nCount = quantile(nCount_RNA, 0.25, na.rm = TRUE),
      Q3_nCount = quantile(nCount_RNA, 0.75, na.rm = TRUE),
      mean_nCount = mean(nCount_RNA, na.rm = TRUE)
    ))
# IQR y media counts por muestra (por región)
(IQR_mean_counts_region <- full_cell_metadata %>%
    group_by(brain_region) %>% 
    summarise(
      Q1_nCount = quantile(nCount_RNA, 0.25, na.rm = TRUE),
      Q3_nCount = quantile(nCount_RNA, 0.75, na.rm = TRUE),
      mean_nCount = mean(nCount_RNA, na.rm = TRUE)
    ))

# IQR y media genes por muestra (por región)
(IQR_mean_genes_region <- full_cell_metadata %>%
    group_by(brain_region) %>% 
    summarise(
      Q1_nFeature_RNA = quantile(nFeature_RNA, 0.25, na.rm = TRUE),
      Q3_nFeature_RNA = quantile(nFeature_RNA, 0.75, na.rm = TRUE),
      mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE)
    ))


# IQR y media genes por muestra
(IQR_mean_genes <- full_cell_metadata %>%
    summarise(
      Q1_nFeature_RNA = quantile(nFeature_RNA, 0.25, na.rm = TRUE),
      Q3_nFeature_RNA = quantile(nFeature_RNA, 0.75, na.rm = TRUE),
      mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE)
    ))

# IQR y media genes por muestra (por caso o control)
(IQR_mean_genes_cases_control <- full_cell_metadata %>%
    group_by(only_case_control) %>% 
    summarise(
      Q1_nFeature_RNA = quantile(nFeature_RNA, 0.25, na.rm = TRUE),
      Q3_nFeature_RNA = quantile(nFeature_RNA, 0.75, na.rm = TRUE),
      mean_nFeature_RNA= mean(nFeature_RNA, na.rm = TRUE)
    ))


# Calcular IQR y núcleos medios por muestra (por caso o control)
(mean_nuclei_cases_control <- full_cell_metadata %>%
  group_by(sample_id) %>%
  summarise(total_nuclei = n()) %>%
  summarise(
    Q1_nuclei_per_sample = quantile(total_nuclei, 0.25, na.rm = TRUE),
    Q3_nuclei_per_sample = quantile(total_nuclei, 0.75, na.rm = TRUE),
    mean_nuclei_per_sample = mean(total_nuclei, na.rm = TRUE)
    
  ))

# Calcular núcleos totales por muestra (por caso o control)
(total_nuclei_cases_control <- full_cell_metadata %>%
  group_by(only_case_control) %>%
  summarise(total_nuclei = n()))

# Agrupar los datos por sample_id y case_control_other_latest para obtener el número de núcleos por muestra
(nuclei_per_participant_id <- full_cell_metadata %>%
  group_by(participant_id, only_case_control) %>%
  summarise(total_nuclei = n()) %>%
  ungroup())

(mean_nuclei_per_participant_id <- nuclei_per_participant_id %>%
    summarise(mean_nuclei = mean(total_nuclei),
              Q1_nuclei_per_sample = quantile(total_nuclei, 0.25, na.rm = TRUE),
              Q3_nuclei_per_sample = quantile(total_nuclei, 0.75, na.rm = TRUE)))

# Calcular la media del número de núcleos por condición
(mean_nuclei_cases_control <- nuclei_per_participant_id %>%
  group_by(only_case_control) %>%
  summarise(mean_nuclei = mean(total_nuclei),
            Q1_nuclei_per_sample = quantile(total_nuclei, 0.25, na.rm = TRUE),
            Q3_nuclei_per_sample = quantile(total_nuclei, 0.75, na.rm = TRUE)))



# Calcular núcleos totales por región
(total_nuclei_region <- full_cell_metadata %>%
  group_by(brain_region) %>%
  summarise(total_nuclei = n()))

  
  
  
  
  
  
#### GRÁFICOS ####
# Parámetros por defecto
beauty_colors = c("#C40834", "#055D7A", "#048B67", "#E6A205","#022836")
beauty_lighter_colors = c("#EF476F", "#118AB2", "#06D6A0", "#FFD166", "#073B4C")

path = "07-final_plots_for_figures/"

my_theme <-  function() {
  theme_minimal(base_size = 26) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 34),
      plot.subtitle = element_text(hjust = 0.5, size = 24),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 30),
      legend.text = element_text(size = 24),
      axis.title.x = element_text(size = 28),
      axis.title.y = element_text(size = 28),
      panel.background = element_rect(fill = "white", colour = "white"),  # Fondo blanco
      panel.grid.minor = element_blank(),  # Eliminar la cuadrícula menor
      panel.grid.major.y = element_line(size = 0.5)  # Ajustar el tamaño de la cuadrícula principal
    )
}

# Función para guardar el gráfico con configuraciones estándar
save_plot <- function(plot, filename, path, width = 12, height = 8, dpi = 500, bg = "white") {
  ggsave(
    filename = filename,
    plot = plot,
    path = path,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg
  )
}



# Boxplots #
data = full_cell_metadata
# Agrupar los participantes por grupo para mejorar la claridad
set.seed(123)  # Para reproducibilidad
random = sample(nrow(full_cell_metadata), size = 10000)
test = full_cell_metadata[random,]
#### Boxplot 1: participant vs UMIs ####
p <- ggplot(data = data) +
  geom_boxplot(mapping = aes(x = participant_id, y = nCount_RNA, fill = case_control_other_latest, colour = case_control_other_latest), 
               outliers = FALSE) +  # Outliers con baja opacidad y tamaño pequeño
  scale_fill_manual(values = beauty_lighter_colors) +  # Colores personalizados para el interior de los boxplots
  scale_color_manual(values = beauty_colors) +  # Paleta de colores personalizada para los bordes
  scale_y_continuous(
    limits = c(0, 60000), 
    breaks = seq(0, 60000, by = 5000), 
    labels = function(x) x / 1000,
    expand = c(0,0)
  ) +  # Ajustar el eje Y
  labs(
    x = "Participant",
    y = bquote('Number of UMIs per Cell (x ' * 10^3 * ')'),
    colour = "Group",
    fill = "Group"
  ) +
  my_theme() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_line(size = 0.25),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Guardar el gráfico con dimensiones específicas y resolución ajustada, asegurando un fondo blanco
save_plot(
  filename = "boxplot_UMI_per_cell_by_participant.png",
  plot = p,
  path = path
)


# Conteos y genes totales por paciente
total_genes_counts <- full_cell_metadata %>%
  group_by(participant_id) %>%
  summarise(
    total_nCount_RNA = sum(nCount_RNA),
    total_nFeature_RNA = sum(nFeature_RNA)
  )
#### Boxplot 2: participant vs genes ####
p <- ggplot(data = data) +
  geom_boxplot(mapping = aes(x = participant_id, y = nFeature_RNA, fill = case_control_other_latest, colour = case_control_other_latest), 
               outliers = FALSE) +  # Outliers con baja opacidad y tamaño pequeño
  scale_fill_manual(values = beauty_lighter_colors) +  # Colores personalizados para el interior de los boxplots
  scale_color_manual(values = beauty_colors) +  # Paleta de colores personalizada para los bordes
  scale_y_continuous(
    limits = c(0, 15000), 
    breaks = seq(0, 15000, by = 2500), 
    labels = function(x) x / 1000,
    expand = c(0,0)
  ) +  # Ajustar el eje Y
  labs(
    x = "Participant",
    y = bquote('Number of Genes per Cell (x ' * 10^3 * ')'),
    colour = "Group",
    fill = "Group"
  ) +
  my_theme() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_line(size = 0.25),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 

# Guardar el gráfico con dimensiones específicas y resolución ajustada, asegurando un fondo blanco
save_plot(
  filename = "boxplot_genes_per_cell_by_participant.png",
  plot = p,
  path = path
)


#### Boxplot 3: brain regions vs UMIs ####
p <- ggplot(data = data) +
  geom_boxplot(mapping = aes(x = brain_region, y = nCount_RNA, fill = brain_region, colour = brain_region), 
               outliers = FALSE) +  # Outliers con baja opacidad y tamaño pequeño
  scale_fill_manual(values = beauty_lighter_colors) +  # Colores personalizados para el interior de los boxplots
  scale_color_manual(values = beauty_colors) +  # Paleta de colores personalizada para los bordes
  scale_y_continuous(
    limits = c(0, 45000), 
    breaks = seq(0, 45000, by = 5000), 
    labels = function(x) x / 1000,
    expand = c(0,0)
  ) +  # Ajustar el eje Y
  labs(
    x = "Brain Region",
    y = bquote('Number of UMIs per Cell (x ' * 10^3 * ')'),
    colour = "Brain Region",
    fill = "Brain Region"
  ) +
  my_theme() +
  theme(legend.position = "none")

save_plot(
  filename = "boxplot_UMI_per_cell_by_brain_region.png",
  plot = p,
  path = path
)

#### Boxplot 4: brain regions vs genes ####
p <- ggplot(data = data) +
  geom_boxplot(mapping = aes(x = brain_region, y = nFeature_RNA, fill = brain_region, colour = brain_region), 
               outliers = FALSE) +  # Outliers con baja opacidad y tamaño pequeño
  scale_fill_manual(values = beauty_lighter_colors) +  # Colores personalizados para el interior de los boxplots
  scale_color_manual(values = beauty_colors) +  # Paleta de colores personalizada para los bordes
  scale_y_continuous(
    limits = c(0, 12500), 
    breaks = seq(0, 12500, by = 2500), 
    labels = function(x) x / 1000,
    expand = c(0,0)
  ) +  # Ajustar el eje Y
  labs(
    x = "Brain Region",
    y = bquote('Number of Genes per Cell (x ' * 10^3 * ')'),
    colour = "Brain Region",
    fill = "Brain Region"
  ) +
  my_theme() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank())

save_plot(
  filename = "boxplot_genes_per_cell_by_brain_region.png",
  plot = p,
  path = path
)

# Barplots #
#### Barplot1: patient vs nuclei ####

total_nuclei_per_patient <- full_cell_metadata %>%
  group_by(participant_id, case_control_other_latest) %>% # Incluimos case_control_other_latest en el group_by
  summarise(total_nuclei = n()) %>%
  ungroup()

p <- ggplot(data = total_nuclei_per_patient) +
  geom_bar(mapping = aes(x = participant_id, y = total_nuclei, fill = case_control_other_latest)
           , stat = "identity") +  
  scale_fill_manual(values = beauty_lighter_colors) +
  scale_y_continuous(
    limits = c(0, 50000),
    breaks = seq(0, 50000, by = 10000), 
    labels = function(x) x / 1000,
    expand = c(0,0)
  ) +
  labs(
    x = "Participant",
    y = bquote('Number of Nuclei (x ' * 10^3 * ')'),
    fill = "Condition"
  ) +
  my_theme() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_line(size = 0.25),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

save_plot(
  filename = "barplot_nuclei_by_participant.png",
  plot = p,
  path = path
)

#### Barplot2: region vs nuclei ####

total_nuclei_per_region <- full_cell_metadata %>%
  group_by(brain_region) %>% # Incluimos case_control_other_latest en el group_by
  summarise(total_nuclei = n()) %>%
  ungroup()

p <- ggplot(data = total_nuclei_per_region) +
  geom_bar(mapping = aes(x = brain_region, y = total_nuclei, fill = brain_region)
           , stat = "identity") +  
  scale_fill_manual(values = beauty_lighter_colors) +
  scale_y_continuous(
    limits = c(0, 600000),
    breaks = seq(0, 600000, by = 100000), 
    labels = function(x) x / 1000,
    expand = c(0,0)
  ) +
  labs(
    x = "Brain Region",
    y = bquote('Number of Nuclei (x ' * 10^3 * ')'),
    fill = "Condition"
  ) +
  my_theme() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_line(size = 0.25)
        )

save_plot(
  filename = "barplot_nuclei_by_region.png",
  plot = p,
  path = path
)








# Cargo las anotaciones finales y elimino los unknown.
annotations <- read.delim("B:/TFM/AMP_PD_cluster/06-results/02-outputs/annotations.tsv")
(table_annotations <- table(annotations$annotations))
annotations <- annotations[!(annotations$annotations == "Unknown"),]



# Extender la paleta si es necesario
beauty_colors_ext = c(
           "#EF476F",
           "#118AB2",
           "#06D6A0",
           "#FFD166",
           "#073B4C",
           "#FF5DF3",
           "#ff8e5d",
           "#8eba00",
           "#79d7ff",
           "#D33322",
           "#964d51",
           "#C46904",
           "#5d8e75",
           "#565656",
           "#790000",
           "#A40052",
           "#008200",
           "#975EFE")

# Ordenar los factores por abundancia
annotations$annotations <- factor(annotations$annotations, levels = names(sort(table(annotations$annotations), decreasing = TRUE)))
annotations$supercluster_term <- factor(annotations$supercluster_term, levels = names(sort(table(annotations$supercluster_term), decreasing = TRUE)))


paleta_28_colores <- c(
       brewer.pal(8, "Set1"),   # 12 colores
       brewer.pal(8, "Dark2"),      # 8 colores
       brewer.pal(12, "Set3"),       # 9 colores
       brewer.pal(8, "Accent")
)
# 1. Barplot: X = "leiden_0_1", stacks = "supercluster_term"
p <- ggplot(annotations, aes(x = as.factor(leiden_0_10_corrected), fill = supercluster_term)) +
  geom_bar(position = "fill") +
  labs(x = "Leiden Clusters", y = "Proportion", fill = "Predicted Cell Types") +
  scale_y_continuous(
    expand = c(0,0)
  ) +
  my_theme() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 34),
        plot.subtitle = element_text(hjust = 0.5, size = 24),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 20),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.2, "lines"),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(20, face = "bold"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 16),
        panel.background = element_rect(fill = "white", colour = "white"),  # Fondo blanco
        panel.grid.minor = element_blank(),  # Eliminar la cuadrícula menor
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(10, 10, 10, 10)
        )+
  scale_fill_manual(values = paleta_28_colores, guide = guide_legend(ncol = 1))
  

save_plot(
  filename = "barplot_predicted_cell_types_by_leiden.png",
  plot = p,
  path = path
)

# 2. Barplot: X = "leiden_0_1", stacks = "annotations"
p <- ggplot(annotations, aes(x = as.factor(leiden_0_10_corrected), fill = annotations)) +
  geom_bar(position = "fill") +
  labs(x = "Leiden Clusters", y = "Proportion", fill = "Final Cell Types") +
  scale_y_continuous(
    expand = c(0,0)
  ) +
  my_theme() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 34),
        plot.subtitle = element_text(hjust = 0.5, size = 24),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 20),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.2, "lines"),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(20, face = "bold"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 16),
        panel.background = element_rect(fill = "white", colour = "white"),  # Fondo blanco
        panel.grid.minor = element_blank(),  # Eliminar la cuadrícula menor
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(10, 10, 10, 10)
  )+
  scale_fill_manual(values = beauty_colors_ext, guide = guide_legend(ncol = 1))


save_plot(
  filename = "barplot_cell_types_by_leiden.png",
  plot = p,
  path = path
)

# 3. Barplot: X = "supercluster_term", stacks = "annotations"
p <- ggplot(annotations, aes(x = as.factor(supercluster_term), fill = annotations)) +
  geom_bar(position = "fill") +
  labs(x = "Predicted Cell Types", y = "Proportion", fill = "Final Cell Types") +
  scale_y_continuous(
    expand = c(0,0)
  ) +
  my_theme() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 34),
        plot.subtitle = element_text(hjust = 0.5, size = 24),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 20),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.2, "lines"),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 16),
        panel.background = element_rect(fill = "white", colour = "white"),  # Fondo blanco
        panel.grid.minor = element_blank(),  # Eliminar la cuadrícula menor
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(10, 10, 10, 10)
  )+
  scale_fill_manual(values = beauty_colors_ext, guide = guide_legend(ncol = 1))


save_plot(
  filename = "barplot_cell_types_by_predicted_cell_types.png",
  plot = p,
  path = path
)

# 4. Barplot: X = "annotations", stacks = "supercluster_term"
p <- ggplot(annotations, aes(x = as.factor(annotations), fill = supercluster_term )) +
  geom_bar(position = "fill") +
  scale_y_continuous(
    expand = c(0,0)
  ) +
  labs(x = "Final Cell Types", y = "Proportion", fill = "Predicted Cell Types") +
  my_theme() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 34),
        plot.subtitle = element_text(hjust = 0.5, size = 24),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "lines"),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 16),
        panel.background = element_rect(fill = "white", colour = "white"),  # Fondo blanco
        panel.grid.minor = element_blank(),  # Eliminar la cuadrícula menor
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(10, 10, 10, 10)
  )+
  scale_fill_manual(values = paleta_28_colores, guide = guide_legend(ncol = 1))

save_plot(
  filename = "barplot_predicted_cell_types_by_cell_types.png",
  plot = p,
  path = path
)


# 5. Piechart
# Calcular las proporciones de cada categoría en annotations
annotation_counts <- annotations %>%
  count(annotations) %>%
  mutate(proportion = n / sum(n),
         percentage = paste0(round(proportion * 100, 1), "%"))
# Crear el pie chart
p <- ggplot(annotation_counts, aes(x = "", y = proportion, fill = annotations)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(x = NULL, y = NULL, fill = "Cell Types") +
  my_theme() +
  scale_fill_manual(values = beauty_colors_ext) +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),  
        axis.line = element_blank())  

save_plot(
  filename = "piechart_cell_types.png",
  plot = p,
  path = path
)

# 6. Barplot: X = "annotations"
p <- ggplot(annotation_counts, aes(x = annotations, y = proportion, fill = annotations)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(
    limits = c(0, 0.6),
    breaks = seq(0, 0.6, by = 0.1), 
    expand = c(0,0)
  ) +
  labs(x = "Cell Types", y = "Proportion", fill = "Cell Types") +
  my_theme() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 34),
        plot.subtitle = element_text(hjust = 0.5, size = 24),
        legend.position = "none",
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 16),
        panel.background = element_rect(fill = "white", colour = "white"),  # Fondo blanco
        panel.grid.minor = element_blank(),  # Eliminar la cuadrícula menor
        panel.grid.major.x = element_blank(),
        plot.margin = margin(10, 10, 10, 10)
  ) +
  geom_text(aes(label = percentage), vjust = -0.5, size = 5, face = "bold") +
  scale_fill_manual(values = beauty_colors_ext, guide = guide_legend(ncol = 1))

save_plot(
  filename = "barplot_cell_types.png",
  plot = p,
  path = path
)


# 6. Barplot de tipos celulares en función de condición
# Calcular las proporciones de cada categoría en annotations
colnames(annotations)[1] <- "barcodekey"
combined_data <- full_cell_metadata %>%
  right_join(annotations, by = "barcodekey")
combined_data <- combined_data[!(combined_data$case_control_other_latest == "Other"),]
# Calcular las proporciones de cada categoría en annotations
annotation_counts <- combined_data %>%
  group_by(annotations, case_control_other_latest) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(annotations) %>%
  mutate(total = sum(count),
         proportion = count / total * 100)

# Calcular el porcentaje de núcleos que pertenecen a muestras de casos (PD)
total_cases <- combined_data %>%
  filter(case_control_other_latest == "Case") %>%
  nrow()

total_cells <- combined_data %>%
  nrow()

case_percentage <- total_cases / total_cells * 100

# Crear el gráfico de barras horizontal apilado
p <- ggplot(annotation_counts, aes(x = annotations, y = proportion, fill = case_control_other_latest)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 5, fontface = "bold") + # Añadir números de células
  geom_hline(yintercept = 100-case_percentage, linetype = "dashed", color = "black", size = 1) +
  scale_y_continuous(
    expand = c(0,0)
  ) +
  coord_flip() +
  labs(x = NULL, y = "Proportion (%)", fill = "Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        plot.subtitle = element_text(hjust = 0.5, size = 15),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(size = 0.5)) +
  scale_fill_manual(values = beauty_colors_ext)

save_plot(
  filename = "barplot_cell_types_case_control.png",
  plot = p,
  path = path
)

# Barplot: X = "leiden_0_1", stacks = "supercluster_term" (Para la referencia SIN filtrar por disección!)
combined_data$supercluster_term.x <- factor(combined_data$supercluster_term.x, levels = names(sort(table(combined_data$supercluster_term.x), decreasing = TRUE)))
p <- ggplot(combined_data, aes(x = as.factor(leiden_0_10_corrected), fill = supercluster_term.x)) +
  geom_bar(position = "fill") +
  labs(x = "Leiden Clusters", y = "Proportion", fill = "Predicted Cell Types") +
  scale_y_continuous(
    expand = c(0,0)
  ) +
  my_theme() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 34),
        plot.subtitle = element_text(hjust = 0.5, size = 24),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "lines"),
        axis.title.x = element_text(size = 24),
        axis.text.x = element_text(20, face = "bold"),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 16),
        panel.background = element_rect(fill = "white", colour = "white"),  # Fondo blanco
        panel.grid.minor = element_blank(),  # Eliminar la cuadrícula menor
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.major.x = element_blank(),
        plot.margin = margin(10, 10, 10, 10)
  )+
  scale_fill_manual(values = paleta_28_colores, guide = guide_legend(ncol = 1))


save_plot(
  filename = "barplot_predicted_cell_types_by_leiden_all_dissections.png",
  plot = p,
  path = path
)


table(combined_data$supercluster_term.x)
table(annotations$supercluster_term)

head(combined_data)
head(annotations)

# Heatmap diseccion vs original
# Normalizar los conteos por fila
cross_tab_norm <- cross_tab %>%
  mutate(across(-supercluster_term.x, ~ . / sum(.)))

# Identificar términos comunes y únicos
terms_x <- unique(cross_tab_norm$supercluster_term.x)
terms_y <- colnames(cross_tab_norm)[-1]
common_terms <- intersect(terms_x, terms_y)
unique_x <- setdiff(terms_x, common_terms)
unique_y <- setdiff(terms_y, common_terms)

# Identificar términos específicos que deben ir al final
specific_terms <- c("Thalamic excitatory", "Hippocampal CA4", "Hippocampal dentate gyrus", "Mammillary body")
remaining_unique_x <- unique_x[!unique_x %in% specific_terms]
unique_x_ordered <- c(remaining_unique_x, specific_terms)

# Ordenar términos: comunes primero y luego únicos
ordered_terms_x <- c(common_terms, unique_x_ordered)
ordered_terms_y <- c(common_terms, unique_y)

# Reordenar filas y columnas en cross_tab_norm
cross_tab_norm <- cross_tab_norm %>%
  arrange(factor(supercluster_term.x, levels = ordered_terms_x))

# Convertir la tabla de conteo normalizada en formato largo adecuado para ggplot2
long_data_norm <- cross_tab_norm %>%
  gather(key = "supercluster_term.y", value = "proportion", -supercluster_term.x) %>%
  mutate(supercluster_term.x = factor(supercluster_term.x, levels = ordered_terms_x),
         supercluster_term.y = factor(supercluster_term.y, levels = ordered_terms_y))
# Crear el heatmap con la paleta continua magma
heatmap <- ggplot(long_data_norm, aes(x = supercluster_term.x, y = supercluster_term.y, fill = proportion)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0, 1), option = "magma") +
  labs(x = "Predicted Cell Types (all-dissections)", y = "Predicted Cell Types (specific-dissections)", fill = "Proportion") +
  my_theme() +
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 8),
        legend.key.size = unit(1, "lines"),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 10))
# Guardar el heatmap
save_plot(
  filename = "heatmap_supercluster_term_comparison.png",
  plot = heatmap,
  path = path,
  dpi = 500
)



# Barplot region cerebral y tipo celular final
# Agrupar los datos por brain_region y annotation
total_nuclei_per_region_annotation <- full_cell_metadata %>%
  right_join(annotations, by = "barcodekey") %>%
  group_by(brain_region, annotations) %>%
  summarise(total_nuclei = n()) %>%
  ungroup()

# Crear el gráfico de barras apilado
p <- ggplot(data = total_nuclei_per_region_annotation) +
  geom_bar(mapping = aes(x = brain_region, y = total_nuclei, fill = annotations), stat = "identity", position = "stack") +
  scale_fill_manual(values = beauty_colors_ext) +
  scale_y_continuous(
    limits = c(0, 600000),
    breaks = seq(0, 600000, by = 100000), 
    labels = function(x) x / 1000,
    expand = c(0,0)
  ) +
  labs(
    x = "Brain Region",
    y = bquote('Number of Nuclei (x ' * 10^3 * ')'),
    fill = "Annotation"
  ) +
  my_theme() +
  theme(legend.position = "right",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_line(size = 0.25)
  )

# Guardar el gráfico
save_plot(
  filename = "barplot_nuclei_by_region_by_cell_types.png",
  plot = p,
  path = path
)








# Agrupar los datos por brain_region y annotation
total_nuclei_per_region_annotation <- full_cell_metadata %>%
  right_join(annotations, by = "barcodekey") %>%
  group_by(brain_region, annotations) %>%
  summarise(total_nuclei = n()) %>%
  ungroup()

# Calcular el total de núcleos por brain_region
total_nuclei_per_region <- total_nuclei_per_region_annotation %>%
  group_by(brain_region) %>%
  summarise(total_nuclei_region = sum(total_nuclei)) %>%
  ungroup()

# Calcular el porcentaje relativo de cada tipo celular en cada región cerebral
relative_nuclei_per_region_annotation <- total_nuclei_per_region_annotation %>%
  left_join(total_nuclei_per_region, by = "brain_region") %>%
  mutate(relative_nuclei = (total_nuclei / total_nuclei_region) * 100) %>%
  select(brain_region, annotations, relative_nuclei)

# Mostrar la tabla resultante
print(relative_nuclei_per_region_annotation, n=100)
write.csv(relative_nuclei_per_region_annotation, "06-results/02-outputs/relative_nuclei_per_region_annotation.csv", row.names = FALSE)
