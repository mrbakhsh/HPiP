    #' FreqInteractors
    #' @title Plot the Pathogen Proteins' frequency of Interactions with Host
    #' Proteins
    #' @param ppi A data.frame containing pathogen proteins in the first column
    #' and host proteins in the second column.
    #' @param cex.size Text size.
    #' @return A frequency plot.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @importFrom ggplot2 element_rect
    #' @importFrom ggplot2 coord_flip
    #' @importFrom ggplot2 element_rect
    #' @importFrom stats aggregate
    #' @description This function plots the pathogen proteins' Frequency of
    #' interactions with host proteins
    #' @export
    #' @examples
    #' ppi <- data.frame(
    #'     node1 = c("A", "A", "A", "B", "B", "B", "B"),
    #'     node2 = c("C", "E", "D", "F", "G", "H", "I")
    #' )
    #' FreqInteractors(ppi)
    FreqInteractors <- function(ppi, cex.size = 12) {
      Group.1 <- NULL
      sizeInt <- NULL

      agg.data <-
        aggregate(list(ppi[, 2]),
          by = list(ppi[, 1]),
          FUN = function(x) {
            paste(unique(x),
              collapse = ";"
            )
          }
        )
      colnames(agg.data)[2] <- "interactors"
      agg.data$interactors <-
        vapply(lapply(strsplit(agg.data$interactors, ";"), unique),
          paste, character(1L),
          collapse = ";"
        )

      agg.data$sizeInt <- lengths(strsplit(agg.data$interactors, ";"))
      agg.data <- agg.data[, c(1, 3)]

      print(ggplot2::ggplot(agg.data, aes(Group.1, y = sizeInt)) +
        geom_point(size = 3) +
        coord_flip() +
        ggthemes::theme_hc() +
        theme(axis.text.y = element_text(
          color = "black", # label color and size
          size = cex.size, lineheight = 0.9, colour = "black"
        )) +
        theme(axis.text.x = element_text(
          color = "black", # label color and size
          size = cex.size, lineheight = 0.9, colour = "black"
        )) +
        theme(axis.line = element_line(
          colour = "black",
          size = 0.5, linetype = "solid"
        )) +
        theme(axis.ticks = element_line(
          colour = "black",
          size = 0.5, linetype = "solid"
        )) +
        labs(title = "", x = "", y = "Number of Interactors") +
        theme(
          legend.title = element_blank(),
          legend.spacing.y = unit(0, "mm"),
          panel.border = element_rect(colour = "black", fill = NA),
          aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
          legend.position = "top",
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black")
        ) +
        theme(text = element_text(size = cex.size, colour = "black")) +
        theme(axis.ticks.length = unit(.25, "cm")))
    }
