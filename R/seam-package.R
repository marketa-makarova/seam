#' seam: Sparse Eigengene Analysis of Modules
#'
#' Tools for extended analysis of weighted gene co-expression networks
#' using sparse principal component analysis. Includes
#' functions for group comparison, visualizing of drivers,
#' and clinical trait correlation analysis.
#' @name seam
NULL
## usethis namespace: start
#' @importFrom dplyr across
#' @importFrom dplyr all_of
#' @importFrom dplyr arrange
#' @importFrom dplyr bind_cols
#' @importFrom dplyr bind_rows
#' @importFrom dplyr c_across
#' @importFrom dplyr case_when
#' @importFrom dplyr desc
#' @importFrom dplyr everything
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr if_any
#' @importFrom dplyr if_else
#' @importFrom dplyr join_by
#' @importFrom dplyr left_join
#' @importFrom dplyr matches
#' @importFrom dplyr mutate
#' @importFrom dplyr n
#' @importFrom dplyr pull
#' @importFrom dplyr rename
#' @importFrom dplyr rename_with
#' @importFrom dplyr right_join
#' @importFrom dplyr rowwise
#' @importFrom dplyr select
#' @importFrom dplyr slice_head
#' @importFrom dplyr starts_with
#' @importFrom dplyr summarize
#' @importFrom dplyr ungroup
#' @importFrom dplyr where
#' @importFrom elasticnet spca
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 theme_minimal
#' @importFrom graphics abline
#' @importFrom graphics boxplot
#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom graphics text
#' @importFrom impute impute.knn
#' @importFrom rlang .data
#' @importFrom rlang :=
#' @importFrom stats as.dist
#' @importFrom stats cor
#' @importFrom stats hclust
#' @importFrom stats mad
#' @importFrom stats na.omit
#' @importFrom stats reorder
#' @importFrom stats sd
#' @importFrom stats setNames
#' @importFrom stats var
#' @importFrom tidyr pivot_longer
#' @importFrom utils head
#' @importFrom utils tail
#' @importFrom utils write.csv
#' @importFrom WGCNA blueWhiteRed
#' @importFrom WGCNA labeledHeatmap
## usethis namespace: end
NULL
