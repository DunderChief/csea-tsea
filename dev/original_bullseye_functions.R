# The lab's original bullseye plotting code (David R. O'Brien), as used by the CSEA/TSEA web tools.
# Kept verbatim for reference: app/R/bullseye.R is a port of draw.bullseyes() and add.heatmap() to
# plain base graphics (polygon() instead of plotrix::draw.circle(nv = 6), a base-graphics legend
# instead of gridBase). The tree functions (create.igraph, draw.tree) built the dendrogram whose node
# positions were placed by hand in tkplot; the app ships those trees as PNG overlays instead.
# Not sourced by the app or the tests.

options(stringsAsFactors=FALSE)

create.tkplot <- function(igraph){
  tkplot(igraph, 1600, 900)
  return(igraph)
}

##  Find only the nodes that are on the edge of graph
outsideNodes <- function(igraph){
  isSample <- c()
  for(i in 1:vcount(igraph)){
    samp <- length(neighbors(igraph, i)) == 1
    isSample <- c(isSample, samp)
  }
  return(isSample)
}

##  Plot an empty graph
create.blank <- function(coords, xmargins=40, ymargins=50, asp=NULL){
  plot(coords, axes=FALSE, xlab='',ylab='', type='n',
       xlim=c(min(coords$x)-xmargins, max(coords$x)+xmargins),
       ylim=c(min(coords$y)-ymargins, max(coords$y)+ymargins),
       asp=asp)
}

##  Create dendrogram and convert to network graph
create.igraph <- function(exprs_dat, top=3000){
  require(ape)
  require(igraph)
  d <- dist(t(exprs_dat[1:top, ])) # Euclidean distances
  phylo <- as.phylo(hclust(d)) # Create a phylo tree
  g <- as.igraph(phylo, directed=FALSE, )
  labels <- V(g)$name
  edges <- get.edgelist(g)
  return(g)
}

##  Creates a line segment with some text going over the top of it
makeLineText <- function(x0, y0, x1, y1, text, lwd=1, text_space=8, text_size=1.2,
                         col='black'){
  segments(x0, y0, x1, y1, lwd=lwd, col=col)
  mid <- c(mean(c(x1,x0)), mean(c(y1, y0)))
  slope <- (y1-y0) / (x1-x0)
  angle <- atan(slope)
  angle <- angle * (180/pi) #Radians to degrees
  x <- mid[1] - (text_space * (abs(slope)/slope))
  y <- mid[2] + text_space
  text(x, y, labels=text, srt=angle, cex=text_size)
}

##  Draw a tree with line segments
draw.tree <- function(coords, igraph, draw.text=FALSE, lwd=1, col='black'){
  require(igraph)
  labels <- V(igraph)$name
  for(i in 1:length(labels)){
    x0 <- coords[i, 1] ## Start of segment is where node starts
    y0 <- coords[i, 2]
    for(j in neighbors(igraph, i)){ ## Now draw lines to it's neighbors
      if(length(neighbors(igraph, i))==1) {
        if(draw.text==TRUE){
          makeLineText(x0, y0, coords[j, 1], coords[j, 2], 
                       text=labels[i], lwd=lwd, col=col)
        } else {
          segments(x0, y0, coords[j, 1], coords[j, 2], lwd=lwd, col=col)
        }      
      } else {
        segments(x0, y0, coords[j, 1], coords[j, 2], lwd=lwd, col=col)
      }
    }
  }
}



draw.bullseyes <- function(coords,
                           pSI,
                           pvals,
                           min_pval=1e-5,
                           sigcutoff=.1,
                           circle_size=4.5,
                           lwd=3,
                           border='black',
                           bottom_label=NULL,
                           top_label=NULL,
                           nonsig.col='lightgrey',
                           hmcols=colorRampPalette(c('yellow2','#FF6700', '#D11500',
                                                     'darkred'), 
                                                   bias=.15, space='Lab')(128),
                           shader_width=2,
                           add.text=TRUE,
                           text.only=FALSE,
                           text_spacer=5,
                           text_size=1){
  require(plotrix)
  require(RColorBrewer)
  coords <- coords[order(rownames(coords)), ]
  ##  Find number of genes per pSI cutoff
  gene_num <- data.frame(colSums(pSI < .05, na.rm=TRUE),
                         colSums(pSI < .01, na.rm=TRUE),
                         colSums(pSI < .001, na.rm=TRUE),
                         colSums(pSI < .0001, na.rm=TRUE))
  names(gene_num) <- c('05', '01', '001', '0001')
  
  
  pvals[pvals <= min_pval] <- min_pval
  sig <- as.data.frame(which(pvals<sigcutoff, arr.ind=TRUE), 
                       row.names=1:length(pvals<sigcutoff))
  sig$pval <- pvals[pvals < sigcutoff]
  bins <- cut(c(min_pval, sigcutoff, sig$pval), 
              breaks=length(hmcols), 
              labels=hmcols[length(hmcols):1])
  if(length(bins) > 2){ 
    bins <- bins[3:length(bins)]
  } else {
    bins <- c()
  }
  
  sig$color <- bins #as.data.frame(matrix(bins, nrow=nrow(sig)), row.names=rownames(sig))
  bins_empty <- as.data.frame(matrix(rep(nonsig.col, length(unlist(pvals))), 
                                     nrow=nrow(pvals)),
                              row.names=rownames(pvals))
  bins_empty[, 1] <- 'white'
  if(nrow(sig) > 0){
    for(i in 1:nrow(sig)){
      bins_empty[sig[i,'row'], sig[i,'col']] <- as.character(sig[i, 'color'])
    }
  }
  bins <- bins_empty
  ##  Draw circles
  gene_num <- gene_num[order(rownames(gene_num)), ]
  
  bins <- bins[order(rownames(bins)), ]
  circle_params <- data.frame(coords, gene_num, bins)
  circle_params$alt_txt1 <- top_label
  circle_params$alt_txt2 <- bottom_label
  
  ##  Draw Circles
  if(!text.only){
    for(i in 1:nrow(circle_params)){
      xx <- circle_params
      x <- xx[i,1]
      y <- xx[i,2]
      r <- sqrt(xx[i, 3]/pi)*circle_size
      line_col='white'
      draw.circle(x, y, radius=r, col=xx[i, 7], nv=6, 
                  border='lightgray', lwd=lwd)
      x1 <- x - cos(60 * pi/180) * (r -4)
      y1 <- y - sin(60 * pi/180) * (r -4)
      x2 <- x + cos(60 * pi/180) * (r -4)
      y2 <- y + sin(60 * pi/180) * (r -4)
      segments(x1, y1, x2, y2, col=line_col, lwd=shader_width)
      side_len <- cos(60 * pi/180) * r * 2
      line_spacer = 7
      for(j in 1:floor((side_len-1)/line_spacer)){
        d <- j*line_spacer
        newx1 <- x1 - d * sin(30 * pi/180)
        newy1 <- y1 + d * cos(30 * pi/180)
        newx2 <- x2 - d 
        newy2 <- y2
        segments(newx1, newy1, newx2, newy2, col=line_col, lwd=shader_width)
        newx1 <- x1 + d 
        newy1 <- y1 
        newx2 <- x2 + d * sin(30 * pi/180)
        newy2 <- y2 - d * cos(30 * pi/180)
        segments(newx1, newy1, newx2, newy2, col=line_col, lwd=shader_width)
      }
      draw.circle(x, y, radius=sqrt(xx[i, 4]/pi)*circle_size, col=xx[i, 8], nv=6, 
                  border=border, lwd=lwd)
      draw.circle(x, y, radius=sqrt(xx[i, 5]/pi)*circle_size, col=xx[i, 9], nv=6, 
                  border=border, lwd=lwd)
      draw.circle(x, y, radius=sqrt(xx[i, 6]/pi)*circle_size, col=xx[i, 10], nv=6, 
                  border=border, lwd=lwd)
    }
  }
  
  ##  Add Text
  if(add.text){
    for(i in 1:nrow(circle_params)){
      xx <- circle_params
      x <- xx[i,1]
      y <- xx[i,2]
      #ifelse(labels$n_lines[i]==1, up_spacer <- 8, up_spacer <- 22)
      #text(x, y + sqrt(xx[i, 3]/pi)*circle_size + up_spacer, labels=xx$alt_txt1[i])
      text(x, y - sqrt(xx[i, 3]/pi)*circle_size - text_spacer, 
           labels=xx$alt_txt2[i], cex=text_size)
    }
  }
}

##  Add heatmap legend to plot
add.heatmap <- function(sigcutoff=0.1, min_pval=1e-5, n_ticks=4,
                        hmcols=colorRampPalette(c('yellow2','#FF6700', '#D11500', 
                                                  'darkred'), 
                                                bias=.15, space='Lab')(128),
                        title_height=.06,
                        tick_height=.025,
                        tick_lab_height=.04
) {
  require(gridBase)
  vp <- baseViewports()
  pushViewport(vp$inner)
  grid.raster(t(hmcols),
              y=unit(0, 'npc'),
              width=unit(.5,"npc"),
              height=unit(.05,"npc"))
  grid.text('p-values', y=unit(title_height,'npc'), gp=gpar(cex=1))
  
  ##  Draw the number of ticks based on variable n_ticks
  for(i in n_ticks:0){
    x <- unit(0.25 + (1-(i/n_ticks)) * .5, 'npc')
    grid.text(as.character(sigcutoff * (i/n_ticks)),
              x=x,
              y=unit(tick_lab_height, 'npc'))
    grid.text('|', 
              x=x,
              y=unit(tick_height, 'npc'))
  }  
  popViewport(1)
}



