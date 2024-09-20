getwd()
#===============================Introduction===================================#

#Publication :
#https://www.nature.com/articles/s41597-020-00617-9

#Figure : 
#https://www.nature.com/articles/s41597-020-00617-9/figures/1

#Données requises
#Un arbre phylogénétique au format newick/tree.
#4 tables de dimension 2*(d'au moins le nombre de feuilles/génome de l'arbre) avec chaques 
#avec chaques données associé à un génome.
#Pour l'annotation il y a besoin des noeuds des taxa.

#Lien des données :
#https://figshare.com/articles/online_resource/910_metagenome-assembled_genomes_from_the_phytobiomes_of_three_urban-farmed_leafy_Asian_greens/12472673
#"glv_mags_qual_tax_summary.tsv" et "glv_mag_de_novo_unrooted.tree" ont été ici
#utilisés et renommés, et l'arbre a été réenraciné au préalable avec iTol.

#Tuto utiles :
#https://guangchuangyu.github.io/ggtree-book/chapter-treeio.html#introduction
#https://yulab-smu.top/treedata-book/chapter10.html
#https://www.datanovia.com/en/blog/ggplot-legend-title-position-and-labels/
#https://viz-ggplot2.rsquaredacademy.com/ggplot2-modify-legend
#https://ggplot2-book.org/scales-guides.html
#https://www.datanovia.com/en/fr/blog/magnifique-diagramme-de-venn-ggplot-avec-r/

#=================================Packages=====================================#

install.packages("ggtree")
install.packages("ggtreeExtra") 
install.packages("ggnewscale")
install.packages("ggVenndiagram")
install.packages("ggimage")
install.packages("ape")
install.packages("readr")

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggtreeExtra))
suppressPackageStartupMessages(library(ggtree))
library(ggVennDiagram)
library(ggnewscale)
library(ggimage)
library(readr)
library(ape)

#==================================Données=====================================#

#Données de la heatmap et du bar-plot
data        <- read_tsv("glv_mags_data.tsv",
                        show_col_types = FALSE)

#Données de l'arbre phylogénétique
Tree        <- read.tree("glv_mags_rooted.tree")
labels      <- Tree$tip.label 
n           <- length(labels) 
taxa_nod    <- c(517, 450, 438, 434, 432, 429, 428, 392, 
                 386, 373, 372, 325, 322, 321, 312, 294
                 )
#Données pour l'annotation de l' arbre
info_ann_v  <- c(438, 434, 432, 429, 428, 386, 373, 322, 321, 312, 294)
info_ann_v  <- sort(info_ann_v)
info_ann_h  <- c("GLV1604141", "GLV1604246",
                 "GLV1604651", "GLV1604689"
                 )
taxa_lab_v  <- c("Patescibacteria", "Chloroflexora", "Cyanobacteria", 
                 "Firmucutes", "Verrucomicrobiota", "Gemmatimonadota", 
                 "Nitrospirota", "Binatota", "Bdellovibrionota", 
                 "Myxococcota", "Acidobacteria"
                 )
taxa_lab_h  <- c("Actinobacteriota", "Bacteroidota",
                 "Gammeproteobacteria", "Alphaproteobacteria"
                 )

#données des noms de génomes
mag_id      <- data$mag_id

#Données pour les heatmap
redundancy  <- data$redundancy
completness <- data$completeness

#Génération de données pour le paramètre détection
set.seed(123)
sample_det  <- c("a", "b", "c", "d", "e", "f", "g", NA)
proba_det   <- c(0.45, 0.1, 0.06, 0.06, 0.06, 0.06, 0.06, 0.15)
detection   <- sample(x       = sample_det,
                      prob    = proba_det,
                      replace = TRUE, 
                      size    = n
                      )

#Données Bar-plot
size        <- data$genome_size_mbp

#============================Formatage des données=============================#

data_det            <- data.frame(labels,
                                  detection
                                  )
data_red            <- data.frame(mag_id,
                                  redundancy
                                  )
data_comp           <- data.frame(mag_id,
                                  completness
                                  )
data_bar            <- data.frame(mag_id, 
                                  size
                                  )
data_lab            <- data.frame(c("GLV1604497", "GLV1604497", 
                                    "GLV1604497", "GLV1604497"
                                    ), 
                                  c("Genome size", "1 Mpb", 
                                    "5 Mbp", "10 Mpb"
                                    )
                                  )
data_taxa_v         <- data.frame(info_ann_v, 
                                  taxa_lab_v
                                  )
data_taxa_h         <- data.frame(info_ann_h, 
                                  taxa_lab_h
                                  )
data_Venn           <- data.frame(sample(x       = "GLV1604141", 
                                         size    = 5, 
                                         replace = TRUE), 
                                  c("Detection", "Choy sum", 
                                    "Bayam", "Kai lan", 
                                    "Venn_detection_legend.png"
                                    )
                                  )

colnames(data_det)     <- c("ID", "Detection")
colnames(data_red)     <- c("ID", "Redundancy")
colnames(data_comp)    <- c("ID", "Completness")
colnames(data_bar)     <- c("ID", "Genome_size")
colnames(data_lab)     <- c("ID", "Label")
colnames(data_taxa_v)  <- c("Node", "Label")
colnames(data_taxa_h)  <- c("ID", "Label")
colnames(data_Venn)    <- c("ID", "Label")

#============================Tracer de la figure===============================#

#paramètres de positionnement
size_text    <- 2.5
size_text2   <- 3.2
width        <- 0.15
max_tree_len <- 1.46
offset_tree  <- 0.11  #Offset = écart entre deux graphiques
offset_hm    <- 0.045
offset_bar   <- 0.040
offset_taxa  <- 0.08
ratio        <- c(-0.9, 0.4, 0.5)
text_angle   <- 21
taxa_angle   <- c(-48, -123, 148, 65)

#Fonction pour annoter les taxas
Annotaxa <- function(i, off){
            geom_fruit(data    = data_taxa_h[i,],
                       geom    = geom_text,
                       mapping = aes(y     = ID,
                                     label = Label
                                     ),
                       offset  = off,
                       size    = size_text2,
                       angle   = taxa_angle[i]
                       )
            }

#Fonction pour annoter le bar_plot
Annoplot <- function(i){
            geom_fruit(data    = data_lab[i+1,],
                       geom    = geom_text,
                       mapping = aes(y     = ID,
                                     label = Label
                                     ),
                       offset  = 0 + ratio[i]*width, 
                       size    = size_text2,
                       angle   = text_angle,
                       hjust   = -0.22,
                       vjust   = 0,
                       grid.params = list(linetype   = 2,
                                          color      = "black",
                                          )
                       )
            }

#Fonction pour annoter le diagramme de Venn
Annovenn <- function(i, h, v){
            geom_fruit(data    = data_Venn[i,],
                       geom    = geom_text,
                       mapping = aes(y     = ID,
                                     label = Label
                                     ),
                       size    = 3,
                       angle   = 0,
                       hjust   = h,
                       vjust   = v,
                       )
            }

#Création d'une image d'un diagramme de Venn
x = list(A=c(1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15),
         B=c(1, 2, 6, 7, 8, 9, 10, 21, 22, 23, 24, 25, 26),
         C=c(1, 3, 4, 5, 7, 8, 9, 10, 31, 32, 33, 34, 35, 36, 37)
         )
v <- ggVennDiagram(x,
                   category.names = NA,
                   label = "none"
                   ) + 
     scale_fill_gradientn(colors = c("#f2efe5", "#ecbe97", "#d0c59b", "#86979d", 
                                     "#efd5a5", "#ec9098", "#79a9ab"
                                     ),
                          guide  = FALSE
                          )
ggsave("Venn_detection_legend.png", v, height = 8, width = 9)

#Arbre phylogénétique circulaire
p <- ggtree(Tree,
            right  = TRUE,
            layout = "circular",
            size   = 0.2
            ) +
     geom_hilight(node      = taxa_nod, 
                  extendto  = max_tree_len + 0.05,
                  to.bottom = TRUE,
                  alpha     = 1,
                  fill      = c("#e1dacd", "#e1ccb3", "#fdf697", "#7fb39b",
                                "#7ed4cb", "#de8d96", "#97e1f1", "#f8b1a2",
                                "#8ebadd", "#f57fc4", "#f57fc4", "#9bdba3",
                                "#c8cacb", "#c5e29f", "#fed7a0", "#f68a8f"
                                ),
                  ) +
     geom_treescale(x        = 0.44, 
                    y        = 52, 
                    width    = 0.1, 
                    linesize = 0.3,
                    fontsize = 2.5,
                    offset   = -2
                    ) 

     #Ouverture et rotation de l'arbre
p <- suppressMessages(open_tree(p, 11))
p <- suppressMessages(rotate_tree(p, 112))

     #Annotation des clades
p <- p +
     #Clade radial à l'arbre
     geom_cladelab(data        = data_taxa_v,
                   mapping     = aes(node  = Node,
                                     label = Label
                                     ),
                   align       = TRUE,
                   offset      = -0.28,
                   barsize     = NA ,
                   fontsize    = size_text,
                   angle       = "auto",
                   horizontal  = TRUE
                   ) +
     #Clade tangent à l'arbre
     Annotaxa(i   = 1,
              off = offset_taxa
              ) +
     Annotaxa(i   = 2:4,
              off = 0
              ) + 
    
     #Heatmap detection
     new_scale_fill() + #Permet d'associer une autre palette de couleurs 
     geom_fruit(data    = data_det,
                geom    = geom_tile,
                mapping = aes(y    = ID,
                              fill = Detection,
                ),
                color    = "grey15",
                width  = 0.05,
                offset  = offset_tree,
                size    = 0.001,
                ) + 
     scale_fill_manual(values   = c("#86979d", "#f2efe5", "#ec9098" ,"#efd5a5",
                                    "#79a9ab", "#ecbe97", "#d0c59b"
                                     ),
                       guide    = FALSE,
                       na.value = "#ffffff",
                       ) +
     #Étiquette heatmap détection
     geom_fruit(data    = data_lab[1,], 
                geom    = geom_text,
                mapping = aes(y     = ID,
                              label = "Detection"
                ),
                offset  = 0,
                size    = size_text,
                angle   = text_angle,
                hjust   = -0.14,
                vjust   = 0.30,
                ) +
     
     #Heatmap completness
     new_scale_fill() +
     geom_fruit(data        = data_comp,
                geom        = geom_tile,
                mapping     = aes(y    = ID,
                                  fill = Completness,
                ),
                color       = "grey15",
                width       = 0.05,
                offset      = offset_hm,
                size        = 0.001,
                ) +
     scale_fill_gradientn(colors = c("#c6deed", "#2c6aad"),
                          limits = c(50, 100),
                          breaks = c(50, 100),
                          #Format de la légende
                          name   = "Completness (%)",
                          guide  = guide_colorbar(ticks = FALSE, 
                                                  title.position = "top",
                                                  title.vjust    = -1,
                                                  label.position = "top",
                                                  barwidth       = 5.7,
                                                  barheight      = 1,
                                                  frame.colour = "black",
                                                  )
                          )+
     #Étiquette heatmap completness
     geom_fruit(data    = data_lab[1,],
                geom    = geom_text,
                mapping = aes(y     = ID,
                              label = "Completness"
                              ),
                offset  = 0,
                size    = size_text,
                angle   = text_angle,
                hjust   = -0.1,
                vjust   = 0.30,
                ) +
    
     #Heatmap redundancy
     new_scale_fill() + 
     geom_fruit(data    = data_red, 
                geom    = geom_tile,
                mapping = aes(y     = ID,
                              fill  = Redundancy
                              ),
                color   = "grey15",
                width   = 0.05,
                offset  = offset_hm,
                size    = 0.001
                ) +
     scale_fill_gradientn(colors = c("#fbd5c2", "#d46350"),
                          limits = c(0, 10),
                          breaks = c(0, 10),
                          #Format de la légende
                          name   = "Redundancy (%)",
                          guide  = guide_colorbar(ticks          = FALSE,
                                                  title.position = "top",
                                                  title.vjust    = -1,
                                                  label.position = "top",
                                                  barwidth       = 5.55,
                                                  barheight      = 1,
                                                  frame.colour   = "black")
                          ) +
     #Étiquette heatmap redundancy
     geom_fruit(data    = data_lab[1,], 
                geom    = geom_text,
                mapping = aes(y     = ID,
                              label = "Redundancy"
                              ),
                offset  = 0,
                size    = size_text,
                angle   = text_angle,
                hjust   = -0.1,
                vjust   = 0.30,
                ) +
  
     #Bar_plot taille du génome 
     new_scale_fill() +
     geom_fruit(data    = data_bar, 
                geom    = geom_col,
                mapping = aes(y    = ID, 
                              x    = Genome_size,
                              fill = "Mpb"
                              ),
                color   = "#d3d3d3",
                pwidth  = width,
                offset  = offset_bar,
                size    = 0.8
                ) +
    scale_fill_manual(values = "#d3d3d3",
                      guide  = FALSE
                      ) +
    #Étiquette bar-plot
    Annoplot(1)+
    Annoplot(2)+
    Annoplot(3)+
    geom_fruit(data    = data_lab[1,],
               geom    = geom_text,
               mapping = aes(y     = ID,
                             label = Label
                             ),
               pwidth  = width,
               offset  = -0.5*width,
               size    = size_text2,
               angle   = 31,
               hjust   = 1.3,
               vjust   = -1
               ) +
    #Étiquette bacteria au centre
    geom_label(mapping    = aes(x     = 0,
                                y     = 0,
                                label = "Bacteria"),
               vjust      = 1.3,
               size       = size_text2,
               label.size = 0
               ) +
    #Légende heatmap detection
    geom_fruit(data    = data_Venn[5,], #ajoute l'image du diagramme de Venn
               geom    = geom_image,
               mapping = aes(y     = ID,
                             image = Label
                             ),
               size    = 0.09,
               angle   = 180,
               nudge_x = 0.67,
               nudge_y = 2.7
               ) +
    geom_fruit_list(Annovenn(i = 1,     #geom_fruit_list pour éviter d'avoir un
                             h = -0.7,  #offset entre chaque Annovenn
                             v = -16.56
                             ),
                    Annovenn(i = 2,
                             h = -0.64,
                             v = -15.1
                             ),
                    Annovenn(i = 3,
                             h = -0.39,
                             v = -6
                             ),
                    Annovenn(i = 4,
                             h = -2.05,
                             v = -6
                             )
                    ) +
    #position de la légende
    theme(legend.position      = c(1, 1.01),
          legend.direction     = "horizontal",
          legend.box           = "horizontal",
          legend.justification = c("right", "top"),
          )           

#Génère un pdf de la figure
pdf(file = "Arbre circulaire Djibril.pdf", height = 10, width = 10)
p
dev.off()