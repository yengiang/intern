# Load libraries

require(remotes)
library(Matrix)
library(SeuratObject)
library(Seurat)
library(future)
library(rgeos)
library(ggplot2)
library(rlang)
library(patchwork)
library(qrcode)
library(devtools)
library(InSituType)
library(readr)
library(RColorBrewer)
library(stringr)
library(ggrepel)

set.seed(0)


# Define function for visualizing cells

function(object, fov, boundaries) {
  if (!is.list(x = boundaries)) {
    if (is.null(x = names(x = boundaries))) {
      boundaries <- rep_len(x = list(boundaries), length.out = length(x = fov))
      names(x = boundaries) <- fov
    } else {
      boundaries <- .AsList(x = boundaries)
    }
  }
  if (any(!nchar(x = names(x = boundaries)))) {
    missing <- setdiff(x = fov, y = names(x = boundaries))
    idx <- which(x = !nchar(x = names(x = boundaries)))
    boundaries <- c(
      boundaries[intersect(x = names(x = boundaries), y = fov)],
      rep_len(x = boundaries[idx], length.out = length(x = missing))
    )
    names(x = boundaries)[!nchar(x = names(x = boundaries))] <- missing
  }
  if (any(!fov %in% names(x = boundaries))) {
    for (i in setdiff(x = fov, y = names(x = boundaries))) {
      boundaries[[i]] <- Boundaries(object = object[[i]])[1L]
    }
  }
  fov <- union(x = fov, y = names(x = boundaries))
  if (length(x = boundaries) != length(x = fov)) {
    fov <- intersect(x = fov, y = names(x = boundaries))
  }
  boundaries <- boundaries[fov]
  for (i in fov) {
    boundaries[[i]] <- Filter(
      f = function(x) {
        return(x %in% Boundaries(object = object[[i]]) || is_na(x = x))
      },
      x = boundaries[[i]]
    )
  }
  boundaries <- Filter(f = length, x = boundaries)
  return(boundaries)
}



.MolsByFOV <- function(object, fov, molecules) {
  keys <- Key(object = object)[fov]
  keyed.mols <- sapply(
    X = names(x = keys),
    FUN = function(img) {
      if (is.null(x = Molecules(object = object[[img]]))) {
        return(NULL)
      }
      key <- keys[img]
      mols <- grep(pattern = paste0('^', key), x = molecules, value = TRUE)
      names(x = mols) <- mols
      mols <- gsub(pattern = paste0('^', key), replacement = '', x = mols)
      keyed <- sapply(
        X = SeuratObject::Keys(object = object[[img]]),
        FUN = function(x) {
          return(grep(pattern = paste0('^', x), x = mols, value = TRUE))
        }
      )
      keyed <- unlist(x = keyed)
      names(x = keyed) <- gsub(
        pattern = '^.*\\.',
        replacement = '',
        x = names(x = keyed)
      )
      missing <- mols[!mols %in% keyed]
      missing <- missing[missing %in% Features(x = object[[img]])]
      if (length(x = missing)) {
        # TODO: replace with default molecules
        default <- Molecules(object = object[[img]])[1L]
        mn <- names(x = missing)
        missing <- paste0(
          SeuratObject::Key(object = object[[img]][[default]]),
          missing
        )
        names(x = missing) <- mn
      }
      return(c(missing, keyed))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  found <- names(x = unlist(x = keyed.mols))
  found <- gsub(pattern = '^.*\\.', replacement = '', x = found)
  missing <- setdiff(x = molecules, y = found)
  names(x = missing) <- missing
  for (img in fov) {
    imissing <- missing
    for (i in seq_along(along.with = imissing)) {
      for (lkey in Keys(object = object[[img]])) {
        imissing[[i]] <- gsub(
          pattern = paste0('^', lkey),
          replacement = '',
          x = imissing[[i]]
        )
      }
    }
    imissing <- names(x = imissing[imissing %in% Features(x = object[[img]])])
    keyed.mols[[img]] <- c(keyed.mols[[img]], imissing)
  }
  keyed.mols <- Filter(f = length, x = keyed.mols)
  keyed.mols <- sapply(X = keyed.mols, FUN = unname, simplify = FALSE)
  return(keyed.mols)
}



imagecells<-function (object, fov = NULL, boundaries = NULL, group.by = NULL, 
                      split.by = NULL, cols = NULL, shuffle.cols = FALSE, size = 0.5, 
                      molecules = NULL, mols.size = 0.1, mols.cols = NULL, mols.alpha = 1, 
                      nmols = 100000, alpha = 1, border.color = "white", border.size = NULL, 
                      na.value = "grey50", dark.background = TRUE, crop = FALSE, 
                      cells = NULL, overlap = FALSE, axes = FALSE, combine = TRUE, 
                      coord.fixed = TRUE) 
{
  cells <- cells %||% Cells(x = object)
  fov <- fov %||% DefaultFOV(object = object)
  fov <- Filter(f = function(x) {
    return(x %in% Images(object = object) && inherits(x = object[[x]], 
                                                      what = "FOV"))
  }, x = fov)
  if (!length(x = fov)) {
    stop("No compatible spatial coordinates present")
  }
  boundaries <- boundaries %||% sapply(X = fov, FUN = function(x) {
    return(DefaultBoundary(object = object[[x]]))
  }, simplify = FALSE, USE.NAMES = TRUE)
  boundaries <- .BoundariesByImage(object = object, fov = fov, 
                                   boundaries = boundaries)
  fov <- names(x = boundaries)
  overlap <- rep_len(x = overlap, length.out = length(x = fov))
  crop <- rep_len(x = crop, length.out = length(x = fov))
  names(x = crop) <- fov
  group.by <- boundaries %!NA% group.by %||% "ident"
  vars <- c(group.by, split.by)
  md <- if (!is_na(x = vars)) {
    FetchData(object = object, vars = vars[!is.na(x = vars)], 
              cells = cells, slot="counts")
  }
  else {
    NULL
  }
  pnames <- unlist(x = lapply(X = seq_along(along.with = fov), 
                              FUN = function(i) {
                                return(if (isTRUE(x = overlap[i])) {
                                  fov[i]
                                } else {
                                  paste(fov[i], boundaries[[i]], sep = "_")
                                })
                              }))
  pdata <- vector(mode = "list", length = length(x = pnames))
  names(x = pdata) <- pnames
  for (i in names(x = pdata)) {
    ul <- unlist(x = strsplit(x = i, split = "_"))
    img <- paste(ul[1:length(ul) - 1], collapse = "_")
    lyr <- ul[length(ul)]
    if (is.na(x = lyr)) {
      lyr <- boundaries[[img]]
    }
    pdata[[i]] <- lapply(X = lyr, FUN = function(l) {
      if (l == "NA") {
        return(NA)
      }
      df <- fortify(model = object[[img]][[l]])
      df <- df[df$cell %in% cells, , drop = FALSE]
      if (!is.null(x = md)) {
        df <- merge(x = df, y = md, by.x = "cell", by.y = 0, 
                    all.x = TRUE)
      }
      df$cell <- paste(l, df$cell, sep = "_")
      df$boundary <- l
      return(df)
    })
    pdata[[i]] <- if (!is_na(x = pdata[[i]])) {
      do.call(what = "rbind", args = pdata[[i]])
    }
    else {
      unlist(x = pdata[[i]])
    }
  }
  if (!is.null(x = molecules)) {
    molecules <- .MolsByFOV(object = object, fov = fov, molecules = molecules)
    mdata <- vector(mode = "list", length = length(x = fov))
    names(x = mdata) <- fov
    for (img in names(x = mdata)) {
      idata <- object[[img]]
      if (!img %in% names(x = molecules)) {
        mdata[[img]] <- NULL
        next
      }
      if (isTRUE(x = crop[img])) {
        idata <- Overlay(x = idata, y = idata)
      }
      imols <- gsub(pattern = paste0("^", Key(object = idata)), 
                    replacement = "", x = molecules[[img]])
      mdata[[img]] <- FetchData(object = idata, vars = imols, 
                                nmols = nmols)
    }
  }
  else {
    mdata <- NULL
  }
  plots <- vector(mode = "list", length = length(x = pdata) * 
                    ifelse(test = length(x = group.by), yes = length(x = group.by), 
                           no = 1L))
  idx <- 1L
  for (group in group.by) {
    for (i in seq_along(along.with = pdata)) {
      img <- unlist(x = strsplit(x = names(x = pdata)[i], 
                                 split = "_"))[1L]
      p <- SingleImagePlot(data = pdata[[i]], col.by = pdata[[i]] %!NA% 
                             group, molecules = mdata[[img]], cols = cols, 
                           shuffle.cols = shuffle.cols, size = size, alpha = alpha, 
                           mols.size = mols.size, mols.cols = mols.cols, 
                           mols.alpha = mols.alpha, border.color = border.color, 
                           border.size = border.size, na.value = na.value, 
                           dark.background = dark.background)
      if (!is.null(x = split.by)) {
        p <- p + facet_wrap(facets = vars(!!sym(x = split.by)))
      }
      if (!isTRUE(x = axes)) {
        p <- p + NoAxes(panel.background = element_blank())
      }
      if (!anyDuplicated(x = pdata[[i]]$cell)) {
        p <- p + guides(fill = guide_legend(override.aes = list(size = 4L, 
                                                                alpha = 1)))
      }
      if (isTRUE(coord.fixed)) {
        p <- p + coord_fixed()
      }
      plots[[idx]] <- p
      idx <- idx + 1L
    }
  }
  if (isTRUE(x = combine)) {
    plots <- wrap_plots(plots)
  }
  return(plots)
}

<bytecode: 0x000002a632875778>
  <environment: namespace:Seurat>
  
  

# Define functions for counting the cells


.BoundariesByImage <- function(object, fov, boundaries) {
  if (!is.list(x = boundaries)) {
    if (is.null(x = names(x = boundaries))) {
      boundaries <- rep_len(x = list(boundaries), length.out = length(x = fov))
      names(x = boundaries) <- fov
    } else {
      boundaries <- .AsList(x = boundaries)
    }
  }
  if (any(!nchar(x = names(x = boundaries)))) {
    missing <- setdiff(x = fov, y = names(x = boundaries))
    idx <- which(x = !nchar(x = names(x = boundaries)))
    boundaries <- c(
      boundaries[intersect(x = names(x = boundaries), y = fov)],
      rep_len(x = boundaries[idx], length.out = length(x = missing))
    )
    names(x = boundaries)[!nchar(x = names(x = boundaries))] <- missing
  }
  if (any(!fov %in% names(x = boundaries))) {
    for (i in setdiff(x = fov, y = names(x = boundaries))) {
      boundaries[[i]] <- Boundaries(object = object[[i]])[1L]
    }
  }
  fov <- union(x = fov, y = names(x = boundaries))
  if (length(x = boundaries) != length(x = fov)) {
    fov <- intersect(x = fov, y = names(x = boundaries))
  }
  boundaries <- boundaries[fov]
  for (i in fov) {
    boundaries[[i]] <- Filter(
      f = function(x) {
        return(x %in% Boundaries(object = object[[i]]) || is_na(x = x))
      },
      x = boundaries[[i]]
    )
  }
  boundaries <- Filter(f = length, x = boundaries)
  return(boundaries)
}













#Here we modify imagefeaturecounts to plot actual count data


imagecounts<-function (object, features, fov = NULL, boundaries = NULL, cols = if (isTRUE(x = blend)) {
  c("lightgrey", "#ff0000", "#00ff00")
} else {
  c("lightgrey", "firebrick1")
}, size = 0.5, min.cutoff = NA, max.cutoff = NA, split.by = NULL, 
molecules = NULL, mols.size = 0.1, mols.cols = NULL, nmols = 1000, 
alpha = 1, border.color = "white", border.size = NULL, dark.background = TRUE, 
blend = FALSE, blend.threshold = 0.5, crop = FALSE, cells = NULL, 
scale = c("feature", "all", "none"), overlap = FALSE, axes = FALSE, 
combine = TRUE, coord.fixed = TRUE) 
{
  cells <- cells %||% Cells(x = object)
  scale <- scale[[1L]]
  scale <- match.arg(arg = scale)
  no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
                    axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold", 
                                                                                           size = 14, margin = margin(r = 7)))
  fov <- fov %||% DefaultFOV(object = object)
  fov <- Filter(f = function(x) {
    return(x %in% Images(object = object) && inherits(x = object[[x]], 
                                                      what = "FOV"))
  }, x = fov)
  if (!length(x = fov)) {
    stop("No compatible spatial coordinates present")
  }
  boundaries <- boundaries %||% sapply(X = fov, FUN = function(x) {
    return(DefaultBoundary(object = object[[x]]))
  }, simplify = FALSE, USE.NAMES = TRUE)
  boundaries <- .BoundariesByImage(object = object, fov = fov, 
                                   boundaries = boundaries)
  fov <- names(x = boundaries)
  if (isTRUE(x = blend) || !is.null(x = split.by)) {
    type <- ifelse(test = isTRUE(x = "blend"), yes = "Blended", 
                   no = "Split")
    if (length(x = fov) != 1L) {
      fov <- fov[1L]
      warning(type, " image feature plots can only be done on a single image, using \"", 
              fov, "\"", call. = FALSE, immediate. = TRUE)
    }
    if (any(!overlap) && length(x = boundaries[[fov]]) > 
        1L) {
      warning(type, " image feature plots require overlapped segmentations", 
              call. = FALSE, immediate. = TRUE)
    }
    overlap <- TRUE
  }
  overlap <- rep_len(x = overlap, length.out = length(x = fov))
  crop <- rep_len(x = crop, length.out = length(x = fov))
  names(x = crop) <- names(x = overlap) <- fov
  if (isTRUE(x = blend)) {
    if (length(x = features) != 2L) {
      stop("Blended feature plots only works with two features")
    }
    default.colors <- eval(expr = formals(fun = ImageFeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)), 
                   `0` = {
                     warning("No colors provided, using default colors", 
                             immediate. = TRUE)
                     default.colors
                   }, `1` = {
                     warning("Only one color provided, assuming specified is double-negative and augmenting with default colors", 
                             immediate. = TRUE)
                     c(cols, default.colors[2:3])
                   }, `2` = {
                     warning("Only two colors provided, assuming specified are for features and augmenting with '", 
                             default.colors[1], "' for double-negatives", 
                             immediate. = TRUE)
                     c(default.colors[1], cols)
                   }, `3` = cols, {
                     warning("More than three colors provided, using only first three", 
                             immediate. = TRUE)
                     cols[1:3]
                   })
  }
  md <- FetchData(object = object, vars = c(features, split.by[1L]), 
                  cells = cells, slot= "counts")
  split.by <- intersect(x = split.by, y = colnames(x = md))
  if (!length(x = split.by)) {
    split.by <- NULL
  }
  imax <- ifelse(test = is.null(x = split.by), yes = ncol(x = md), 
                 no = ncol(x = md) - length(x = split.by))
  features <- colnames(x = md)[1:imax]
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(md[[feature]]), 
                  no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(md[[feature]]), 
                  no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(features, min.cutoff, 
                                              max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, 
  ]$maxcolors, no = length(x = cols))
  for (i in seq_along(along.with = features)) {
    f <- features[[i]]
    data.feature <- md[[f]]
    min.use <- SetQuantile(cutoff = min.cutoff[i], data = data.feature)
    max.use <- SetQuantile(cutoff = max.cutoff[i], data = data.feature)
    data.feature[data.feature < min.use] <- min.use
    data.feature[data.feature > max.use] <- max.use
    if (brewer.gran != 2) {
      data.feature <- if (all(data.feature == 0)) {
        rep_len(x = 0, length.out = length(x = data.feature))
      }
      else {
        as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), 
                                         breaks = brewer.gran)))
      }
    }
    md[[f]] <- data.feature
  }
  if (is.null(x = split.by)) {
    split.by <- RandomName()
    md[[split.by]] <- factor(x = split.by)
  }
  if (!is.factor(x = md[[split.by]])) {
    md[[split.by]] <- factor(x = md[[split.by]])
  }
  if (isTRUE(x = blend)) {
    md <- lapply(X = levels(x = md[[split.by]]), FUN = function(x) {
      df <- md[as.character(x = md[[split.by]]) == x, , 
               drop = FALSE]
      no.expression <- features[colMeans(x = df[, features]) == 
                                  0]
      if (length(x = no.expression)) {
        stop("The following features have no value: ", 
             paste(no.expression, collapse = ", "))
      }
      return(cbind(df[, split.by, drop = FALSE], BlendExpression(data = df[, 
                                                                           features])))
    })
    md <- do.call(what = "rbind", args = md)
    features <- setdiff(x = colnames(x = md), y = split.by)
  }
  pnames <- unlist(x = lapply(X = seq_along(along.with = fov), 
                              FUN = function(i) {
                                return(if (isTRUE(x = overlap[i])) {
                                  fov[i]
                                } else {
                                  paste(fov[i], boundaries[[i]], sep = "_")
                                })
                              }))
  pdata <- vector(mode = "list", length = length(x = pnames))
  names(x = pdata) <- pnames
  for (i in names(x = pdata)) {
    ul <- unlist(x = strsplit(x = i, split = "_"))
    img <- paste(ul[1:length(ul) - 1], collapse = "_")
    lyr <- ul[length(ul)]
    if (is.na(x = lyr)) {
      lyr <- boundaries[[img]]
    }
    pdata[[i]] <- lapply(X = lyr, FUN = function(l) {
      df <- fortify(model = object[[img]][[l]])
      df <- df[df$cell %in% cells, , drop = FALSE]
      if (!is.null(x = md)) {
        df <- merge(x = df, y = md, by.x = "cell", by.y = 0, 
                    all.x = TRUE)
      }
      df$cell <- paste(l, df$cell, sep = "_")
      df$boundary <- l
      return(df)
    })
    pdata[[i]] <- if (!is_na(x = pdata[[i]])) {
      do.call(what = "rbind", args = pdata[[i]])
    }
    else {
      unlist(x = pdata[[i]])
    }
  }
  if (!is.null(x = molecules)) {
    molecules <- .MolsByFOV(object = object, fov = fov, molecules = molecules)
    mdata <- vector(mode = "list", length = length(x = fov))
    names(x = mdata) <- fov
    for (img in names(x = mdata)) {
      idata <- object[[img]]
      if (!img %in% names(x = molecules)) {
        mdata[[img]] <- NULL
        next
      }
      if (isTRUE(x = crop[img])) {
        idata <- Overlay(x = idata, y = idata)
      }
      imols <- gsub(pattern = paste0("^", Key(object = idata)), 
                    replacement = "", x = molecules[[img]])
      mdata[[img]] <- FetchData(object = idata, vars = imols, 
                                nmols = nmols)
    }
  }
  else {
    mdata <- NULL
  }
  if (isTRUE(x = blend)) {
    ncol <- 4
    color.matrix <- BlendMatrix(two.colors = cols[2:3], col.threshold = blend.threshold, 
                                negative.color = cols[1])
    cols <- cols[2:3]
    colors <- list(color.matrix[, 1], color.matrix[1, ], 
                   as.vector(x = color.matrix))
    blend.legend <- BlendMap(color.matrix = color.matrix)
  }
  limits <- switch(EXPR = scale, all = range(unlist(x = md[, 
                                                           features])), NULL)
  plots <- vector(mode = "list", length = length(x = levels(x = md[[split.by]])))
  names(x = plots) <- levels(x = md[[split.by]])
  for (i in seq_along(along.with = levels(x = md[[split.by]]))) {
    ident <- levels(x = md[[split.by]])[i]
    plots[[ident]] <- vector(mode = "list", length = length(x = pdata))
    names(x = plots[[ident]]) <- names(x = pdata)
    if (isTRUE(x = blend)) {
      blend.key <- suppressMessages(expr = blend.legend + 
                                      scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = md[[split.by]])) > 
                                                                                             1, yes = ident, no = "")), expand = c(0, 0)) + 
                                      labs(x = features[1L], y = features[2L], title = if (i == 
                                                                                           1L) {
                                        paste("Color threshold:", blend.threshold)
                                      }
                                      else {
                                        NULL
                                      }) + no.right)
    }
    for (j in seq_along(along.with = pdata)) {
      key <- names(x = pdata)[j]
      img <- unlist(x = strsplit(x = key, split = "_"))[1L]
      plots[[ident]][[key]] <- vector(mode = "list", length = length(x = features) + 
                                        ifelse(test = isTRUE(x = blend), yes = 1L, no = 0L))
      data.plot <- pdata[[j]][as.character(x = pdata[[j]][[split.by]]) == 
                                ident, , drop = FALSE]
      for (y in seq_along(along.with = features)) {
        feature <- features[y]
        cols.use <- if (isTRUE(x = blend)) {
          cc <- as.numeric(x = as.character(x = data.plot[, 
                                                          feature])) + 1
          colors[[y]][sort(unique(x = cc))]
        }
        else {
          NULL
        }
        colnames(data.plot) <- gsub("-", "_", colnames(data.plot))
        p <- SingleImagePlot(data = data.plot, col.by = gsub("-", 
                                                             "_", feature), size = size, col.factor = blend, 
                             cols = cols.use, molecules = mdata[[img]], 
                             mols.size = mols.size, mols.cols = mols.cols, 
                             alpha = alpha, border.color = border.color, 
                             border.size = border.size, dark.background = dark.background) + 
          CenterTitle() + labs(fill = feature)
        if (isTRUE(x = blend)) {
          p <- p + guides(fill = "none")
        }
        if (isTRUE(coord.fixed)) {
          p <- p + coord_fixed()
        }
        if (!isTRUE(x = axes)) {
          p <- p + NoAxes(panel.background = element_blank())
        }
        else if (isTRUE(x = blend) || length(x = levels(x = md[[split.by]])) > 
                 1L) {
          if (y != 1L) {
            p <- p + theme(axis.line.y = element_blank(), 
                           axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
                           axis.title.y.left = element_blank())
          }
          if (i != length(x = levels(x = md[[split.by]]))) {
            p <- p + theme(axis.line.x = element_blank(), 
                           axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
                           axis.title.x = element_blank())
          }
        }
        if (!isTRUE(x = blend)) {
          if (length(x = cols) == 1L) {
            p <- p + scale_fill_brewer(palette = cols)
          }
          else {
            cols.grad <- cols
            fexp <- data.plot[data.plot[[split.by]] == 
                                ident, feature, drop = TRUE]
            fexp <- unique(x = fexp)
            if (length(x = fexp) == 1L) {
              warning("All cells have the same value (", 
                      fexp, ") of ", feature, call. = FALSE, 
                      immediate. = TRUE)
              if (fexp == 0) {
                cols.grad <- cols.grad[1L]
              }
            }
            if (scale == "feature") {
              limits <- range(pdata[[j]][[feature]])
            }
            p <- p + ggplot2::scale_fill_gradientn(colors = cols.grad, 
                                                   guide = "colorbar", limits = limits)
          }
        }
        p <- p + if (i == 1L) {
          ggplot2::labs(title = feature)
        }
        else {
          ggplot2::labs(title = NULL)
        }
        plots[[ident]][[key]][[y]] <- p
      }
      if (isTRUE(x = blend)) {
        plots[[ident]][[key]][[length(x = plots[[ident]][[key]])]] <- blend.key
      }
      else if (length(x = levels(x = md[[split.by]])) > 
               1L) {
        plots[[ident]][[key]][[y]] <- suppressMessages(expr = plots[[ident]][[key]][[y]] + 
                                                         scale_y_continuous(sec.axis = dup_axis(name = ident)) + 
                                                         no.right)
      }
    }
    plots[[ident]] <- unlist(x = plots[[ident]], recursive = FALSE, 
                             use.names = FALSE)
  }
  plots <- unlist(x = plots, recursive = FALSE, use.names = FALSE)
  if (isTRUE(x = combine)) {
    if (isTRUE(x = blend) || length(x = levels(x = md[[split.by]])) > 
        1L) {
      plots <- wrap_plots(plots, ncol = ifelse(test = isTRUE(x = blend), 
                                               yes = 4L, no = length(x = features)), nrow = length(x = levels(x = md[[split.by]])), 
                          guides = "collect")
    }
    else {
      plots <- wrap_plots(plots)
    }
  }
  return(plots)
}
<bytecode: 0x000002a632680470>
  <environment: namespace:Seurat>