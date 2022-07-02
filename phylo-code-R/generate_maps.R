
set.seed(123)

library(ape)
library(picante)
library(ggplot2)
#library(ggtree)
library(ggrepel)
library(rentrez)
library(tidyverse)
library(magrittr)

setwd("~/Github/betacov_bats_mapping")

tree <- read.nexus("bCoVmapping_240bp_RdRp_it4.nex")
evo1 <- evol.distinct(tree, type = 'fair.proportion') # 'equal.splits' # 'fair.proportion'

evo1 %>% 
  as_tibble() %>%
  mutate(Species = str_replace(Species, "NC_", "NC")) %>%
  mutate(Species = word(Species, start = 1, end = 1, sep = "_")) %>%
  mutate(Species = str_replace(Species, "NC", "NC_")) -> evo1

evo1 %>% pull(Species) %>% unique() %>% str_replace("\\'", "") -> all_cov_ids

# Grab the host information from the NCBI accession information

get_metadata <- function(x){
  
  query_index <- split(seq(1,length(x)), ceiling(seq_along(seq(1,length(x)))/300))
  Seq_result <- vector("list", length(query_index))
  
  for (i in 1:length(query_index)) {
    Seq_result[[i]] <- entrez_summary(db = "nuccore",id = x[unlist(query_index[i])])
    Sys.sleep(1)
  }
  
  if(length(x) == 1){
    return(Seq_result)
  } else {
    return(Seq_result %>% 
             purrr::flatten() %>% 
             unname() %>% 
             structure(class = c("list", "esummary_list"))) # coerce to esummary list
  }
}

all_cov_summary <- get_metadata(all_cov_ids)

all_accession <- extract_from_esummary(all_cov_summary, elements = c("accessionversion"))
all_subname_list <- extract_from_esummary(all_cov_summary, elements = c("subname", "subtype"), 
                                          simplify=FALSE)
all_subname_list <- all_subname_list %>%
  lapply(function(x) {
    subname_row <- x["subname"] %>%
      stringr::str_split("\\|") %>%
      unlist() %>%
      matrix(nrow = 1, byrow = FALSE) %>%
      as.data.frame()
    subtype_names <- x["subtype"] %>%
      stringr::str_split("\\|") %>%
      unlist()
    set_colnames(subname_row, subtype_names)
  })

all_metadata_df <-  data.frame(accession = all_accession,
                               suppressMessages(bind_rows(all_subname_list)))
str(all_metadata_df[,1:21]) 

all_metadata_df %>% dplyr::select(accession, host) -> host_df

# Clean the data frame

host_df %>% 
  as_tibble() %>%
  filter(!str_detect(host, ".sp")) %>%
  filter(!is.na(host)) -> host_df

host_df %<>%
  mutate(host = recode(host, !!!c("Pipistrellys abramus" = "Pipistrellus abramus",
                                  "Chaerephon plicata" = "Chaerephon plicatus",
                                  "Rhinilophus ferrumequinum" = "Rhinolophus ferrumequinum",
                                  "Hipposideroa pratti" = "Hipposideros pratti")))

host_df

raw.hosts <- host_df$host

# Pull NCBI taxonomy

library(taxize)

mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}
hdict <- function(names) { 
  names.orig <- names
  names <- str_replace(names, " sp\\.","")
  names <- str_replace(names, " gen\\.","")
  u <- get_uid(names, rank_filter = c("subspecies", "species", "genus", "family", "order", "class"), 
               division_filter = "vertebrates", ask = FALSE)
  c <- classification(u)
  n <- !is.na(u)
  attributes(u) <- NULL
  s <- unlist(lapply(c, function(x){tryCatch(x$name[[which(x$rank=="species")]], error = function(e) {NA})}), use.names = FALSE)
  g <- unlist(lapply(c, function(x){tryCatch(x$name[[which(x$rank=="genus")]], error = function(e) {NA})}), use.names = FALSE)
  f <- unlist(lapply(c, function(x){tryCatch(x$name[[which(x$rank=="family")]], error = function(e) {NA})}), use.names = FALSE)
  o <- unlist(lapply(c, function(x){tryCatch(x$name[[which(x$rank=="order")]], error = function(e) {NA})}), use.names = FALSE)
  c2 <- unlist(lapply(c, function(x){tryCatch(x$name[[which(x$rank=="class")]], error = function(e) {NA})}), use.names = FALSE)
  
  levels <- c("species", "genus", "family", "order", "class")
  u <- unlist(lapply(c, function(x){tryCatch(last(na.omit(x[x$rank %in% levels,'id'])), 
                                             error = function(e) {NA})}), use.names = FALSE)
  
  data.frame(HostOriginal = names.orig,
             HostTaxID = u,
             HostNCBIResolved = n, 
             Host = s,
             HostGenus = g,
             HostFamily = f,
             HostOrder = o, 
             HostClass = c2) %>% mutate_cond(HostNCBIResolved == FALSE, Host = HostOriginal) %>% return()
}

matches <- hdict(unique(word(raw.hosts, 1, 2)))

matches %>% filter(HostNCBIResolved == 'TRUE' | HostOriginal %in% c('Neoromicia capensis')) %>% 
  rename(host = 'Host') -> matches

evo1 %>%
  mutate(Species = str_replace(Species, "\\'","")) %>%
  left_join(host_df, by = c("Species" = "accession")) %>%
  left_join(matches) -> hell.df

hell.df %>% dplyr::select(host, w) %>%
  rename(Host = 'host',
         EvolDist = 'w') %>%
  na.omit() %>%
  unique() %>%
  group_by(Host) %>%
  summarize(MaxEvolDist = max(EvolDist),
            MeanEvolDist = mean(EvolDist)) -> hell.df

# Grab the IUCN object which you know you simply love to work with 

library(fasterize)
library(sf)
library(maps)

iucn <- st_read(dsn = 'C:/Users/cjcar/Dropbox/HowManyHelminths2019',
                 layer = 'TERRESTRIAL_MAMMALS')

iucn <- iucn[iucn$order_ == 'CHIROPTERA',]

iucn.all <- iucn[iucn$binomial %in% unique(hell.df$Host),]

iucn.all %>% rename(Host = binomial) %>% left_join(hell.df) -> iucn.all

r <- raster::getData("worldclim",var="alt",res=5) # Make a blank raster
map.num <- fasterize(iucn.all, r, field = NULL, fun = 'count')
map.1 <- fasterize(iucn.all, r, field = 'MaxEvolDist', fun = 'sum')
map.2 <- fasterize(iucn.all, r, field = 'MeanEvolDist', fun = 'sum')

par(mfrow = c(2,1))
par(mar = c(0,0,0,0))

plot(map.num, axes = F, box = F, main = '', col = rev(RColorBrewer::brewer.pal(11,"Spectral")))
maps::map('world', interior = F, add = T)
title("(A)", adj = 0.05, line = -1)
plot(map.1, axes = F, box = F, main = '', col = rev(RColorBrewer::brewer.pal(11,"Spectral")))
maps::map('world', interior = F, add = T)
title("(B)", adj = 0.05, line = -1)

library(raster)
raster::writeRaster(map.1, "VirusDistinctiveness.tif", overwrite = TRUE)

#plot(map.2, main = 'Mean evolutionary distance of BCoV')

###########################

dist.df <- cophenetic(tree)

library(vegan)


# nmds <- metaMDS(as.matrix(dist.df), distance = "bray")
# plot(nmds)

# pcnm <- pcnm(as.matrix(dist.df))
# pca <- princomp(pcnm$vectors)
# plot(pca$scores[,2] ~ pca$scores[,1])

hell <- pcoa(as.matrix(dist.df))

hell.species <- rownames(hell$vectors[hell$vectors[,1] < -5,]) %>% unique()

good <- pcoa(as.matrix(dist.df[!(rownames(dist.df) %in% hell.species),
                               !(colnames(dist.df) %in% hell.species)]))

plot(good$vectors[,2] ~ good$vectors[,1])

# Bind the virus vectors to the hosts 

virus.all <- good$vectors[,1:2]

virus.all %>% as_tibble() %>% 
  mutate(Virus = rownames(good$vectors)) %>%
  dplyr::select(Virus, Axis.1, Axis.2) %>%
  rename(PCoA1 = 'Axis.1',
         PCoA2 = 'Axis.2') %>% 
  mutate(Species = str_replace(Virus, "NC_", "NC")) %>%
  mutate(Species = word(Species, start = 1, end = 1, sep = "_")) %>%
  mutate(Species = str_replace(Species, "NC", "NC_")) %>%
  mutate(Species = str_replace(Species, "\\'", "")) %>%
  dplyr::select(Species, PCoA1, PCoA2) %>%
  group_by(Species) %>%
  summarize(PCoA1 = mean(PCoA1), 
            PCoA2 = mean(PCoA2)) %>%
  left_join(host_df, by = c("Species" = "accession")) %>%
  left_join(matches) %>%
  dplyr::select(host, PCoA1, PCoA2) %>%
  na.omit() %>%
  group_by(host) %>%
  summarize(PCoA1 = mean(PCoA1),
            PCoA2 = mean(PCoA2)) -> pc.df# %>%
 # rename(Host = 'match') 

iucn.all <- iucn[iucn$binomial %in% unique(pc.df$host),]

iucn.all %>% rename(host = binomial) %>% left_join(pc.df) -> pc.all

pc.1 <- fasterize(pc.all, r, field = 'PCoA1', fun = 'sum')
pc.2 <- fasterize(pc.all, r, field = 'PCoA2', fun = 'sum')

par(mfrow = c(3,1))
plot(map.num, main = 'Number of bat species with BCoV in GenBank')
plot(pc.1/map.num, main = 'PCoA1')
writeRaster(pc.1/map.num, 'VirusPC1.tif', overwrite = TRUE)
plot(pc.2/map.num, main = 'PCoA2')
writeRaster(pc.2/map.num, 'VirusPC2.tif', overwrite = TRUE)

############################# BIVARIATE CODE

library(classInt)
library(raster)
library(rgdal)
library(dismo)
library(XML)
library(maps)
library(sp)

colmat<-function(nquantiles=10, upperleft=rgb(255,255,0, maxColorValue=255), upperright=rgb(50,205,50, maxColorValue=255), bottomleft="red", bottomright="blue", xlab="x label", ylab="y label"){
  my.data<-seq(0,1,.01)
  my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
  my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
  my.pal.2<-findColours(my.class,c(upperright, bottomright))
  col.matrix<-matrix(nrow = 101, ncol = 101, NA)
  for(i in 1:101){
    my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
    col.matrix[102-i,]<-findColours(my.class,my.col)}
  plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, xlab=xlab, ylab=ylab,cex.lab=1.3)
  for(i in 1:101){
    col.temp<-col.matrix[i-1,]
    points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
  seqs<-seq(0,100,(100/nquantiles))
  seqs[1]<-1
  col.matrix<-col.matrix[c(seqs), c(seqs)]}

col.matrix<-colmat(nquantiles=10,
                   xlab="PCoA1", ylab="PCoA2")

bivariate.map<-function(rasterx, rastery, colormatrix=col.matrix, nquantiles=10){
  quanmean<-getValues(rasterx)
  temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
  brks<-with(temp, quantile(unique(temp),na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r1<-within(temp, quantile <- cut(quanmean, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr<-data.frame(r1[,2]) 
  quanvar<-getValues(rastery)
  temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
  brks<-with(temp, quantile(unique(temp),na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r2<-within(temp, quantile <- cut(quanvar, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr2<-data.frame(r2[,2])
  as.numeric.factor<-function(x) {as.numeric(levels(x))[x]}
  col.matrix2<-colormatrix
  cn<-unique(colormatrix)
  for(i in 1:length(col.matrix2)){
    ifelse(is.na(col.matrix2[i]),col.matrix2[i]<-1,col.matrix2[i]<-which(col.matrix2[i]==cn)[1])}
  cols<-numeric(length(quantr[,1]))
  for(i in 1:length(quantr[,1])){
    a<-as.numeric.factor(quantr[i,1])
    b<-as.numeric.factor(quantr2[i,1])
    cols[i]<-as.numeric(col.matrix2[b,a])}
  r<-rasterx
  r[1:length(r)]<-cols
  return(r)}
### STOP COPYING AND PASTE INTO R ###

bivmap<-bivariate.map(pc.1/map.num,pc.2/map.num, colormatrix=col.matrix, nquantiles=10)
plot(bivmap,frame.plot=F,axes=F,box=F,add=F,legend=F,col=as.vector(col.matrix))
map(interior=F,add=T)

############################################################
############################################################
############################################################
############################################################
############################################################

library(sf)

iucn <- st_read(dsn = 'C:/Users/cjcar/Dropbox/HowManyHelminths2019',
                layer = 'TERRESTRIAL_MAMMALS')

iucn <- iucn[iucn$order_ == 'CHIROPTERA',]
batnames <- unique(iucn$binomial)

# Same map but bats

tree2 <- read.nexus("~/Github/clover/clover/clover_0.1_mammalviruses/phylogenies/upham_tree_666.nex")
dist.df <- cophenetic(tree2)
dist.df <- dist.df[rownames(dist.df) %in% gsub(" ", "_", batnames), 
        colnames(dist.df) %in% gsub(" ", "_", batnames)]

good <- pcoa(as.matrix(dist.df))
plot(good$vectors[,2] ~ good$vectors[,1])

# Bind the virus vectors to the hosts 

iucn.all <- iucn[iucn$binomial %in% gsub("_"," ", rownames(good$vectors)),]
vec.df <- data.frame(binomial = rownames(good$vectors),
                     PC1 = good$vectors[,1],
                     PC2 = good$vectors[,2])

iucn.all %>% left_join(vec.df %>% mutate(binomial = gsub("_"," ", binomial))) -> pc.all

map.num <- fasterize(iucn.all, r, field = NULL, fun = 'count')
pc.1 <- fasterize(pc.all, r, field = 'PC1', fun = 'sum')
pc.2 <- fasterize(pc.all, r, field = 'PC2', fun = 'sum')

par(mfrow = c(2,1))
plot(map.num, main = 'Number of bat species with BCoV in GenBank')
plot(pc.1/map.num, main = 'PCoA1')
writeRaster(pc.1/map.num, 'BatPC1.tif', overwrite = TRUE)
plot(pc.2/map.num, main = 'PCoA2')
writeRaster(pc.2/map.num, 'BatPC2.tif', overwrite = TRUE)

############################# BIVARIATE CODE

library(classInt)
library(raster)
library(rgdal)
library(dismo)
library(XML)
library(maps)
library(sp)

colmat<-function(nquantiles=10, upperleft=rgb(255,255,0, maxColorValue=255), upperright=rgb(50,205,50, maxColorValue=255), bottomleft="red", bottomright="blue", xlab="x label", ylab="y label"){
  my.data<-seq(0,1,.01)
  my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
  my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
  my.pal.2<-findColours(my.class,c(upperright, bottomright))
  col.matrix<-matrix(nrow = 101, ncol = 101, NA)
  for(i in 1:101){
    my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
    col.matrix[102-i,]<-findColours(my.class,my.col)}
  plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, xlab=xlab, ylab=ylab,cex.lab=1.3)
  for(i in 1:101){
    col.temp<-col.matrix[i-1,]
    points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
  seqs<-seq(0,100,(100/nquantiles))
  seqs[1]<-1
  col.matrix<-col.matrix[c(seqs), c(seqs)]}

col.matrix<-colmat(nquantiles=10,
                   xlab="PCoA1", ylab="PCoA2")

bivariate.map<-function(rasterx, rastery, colormatrix=col.matrix, nquantiles=10){
  quanmean<-getValues(rasterx)
  temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
  brks<-with(temp, quantile(unique(temp),na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r1<-within(temp, quantile <- cut(quanmean, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr<-data.frame(r1[,2]) 
  quanvar<-getValues(rastery)
  temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
  brks<-with(temp, quantile(unique(temp),na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r2<-within(temp, quantile <- cut(quanvar, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr2<-data.frame(r2[,2])
  as.numeric.factor<-function(x) {as.numeric(levels(x))[x]}
  col.matrix2<-colormatrix
  cn<-unique(colormatrix)
  for(i in 1:length(col.matrix2)){
    ifelse(is.na(col.matrix2[i]),col.matrix2[i]<-1,col.matrix2[i]<-which(col.matrix2[i]==cn)[1])}
  cols<-numeric(length(quantr[,1]))
  for(i in 1:length(quantr[,1])){
    a<-as.numeric.factor(quantr[i,1])
    b<-as.numeric.factor(quantr2[i,1])
    cols[i]<-as.numeric(col.matrix2[b,a])}
  r<-rasterx
  r[1:length(r)]<-cols
  return(r)}
### STOP COPYING AND PASTE INTO R ###

bivmap2<-bivariate.map(pc.1/map.num,pc.2/map.num, colormatrix=col.matrix, nquantiles=10)
plot(bivmap2,frame.plot=F,axes=F,box=F,add=F,legend=F,col=as.vector(col.matrix))
map(interior=F,add=T)


par(mfrow=c(2,1))
plot(bivmap2,frame.plot=F,axes=F,box=F,add=F,legend=F,col=as.vector(col.matrix))
map(interior=F,add=T)
plot(bivmap,frame.plot=F,axes=F,box=F,add=F,legend=F,col=as.vector(t(col.matrix)))
map(interior=F,add=T)

par(mfrow=c(1,1))

# The transposed colors:

colmat<-function(nquantiles=10, upperleft=rgb(255,255,0, maxColorValue=255), upperright=rgb(50,205,50, maxColorValue=255), bottomleft="red", bottomright="blue", xlab="x label", ylab="y label"){
  my.data<-seq(0,1,.01)
  my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
  my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
  my.pal.2<-findColours(my.class,c(upperright, bottomright))
  col.matrix<-matrix(nrow = 101, ncol = 101, NA)
  for(i in 1:101){
    my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
    col.matrix[102-i,]<-findColours(my.class,my.col)}
  col.matrix <- t(col.matrix)
  plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, xlab=xlab, ylab=ylab,cex.lab=1.3)
  for(i in 1:101){
    col.temp<-col.matrix[i-1,]
    points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
  seqs<-seq(0,100,(100/nquantiles))
  seqs[1]<-1
  col.matrix<-col.matrix[c(seqs), c(seqs)]}
plot(col.matrix)
col.matrix<-colmat(nquantiles=10,
                   xlab="PCoA1", ylab="PCoA2")


################## FINALLY...... bat evo distinctiveness

evo2 <- evol.distinct(tree2, type = 'fair.proportion') # 'equal.splits' # 'fair.proportion'
iucn.all <- iucn[iucn$binomial %in% gsub("_"," ", evo2$Species),]

iucn.2 <- left_join(iucn.all, evo2 %>% mutate(Species = gsub("_"," ", Species)), 
                    by = c('binomial' = 'Species'))

map.num <- fasterize(iucn.2, r, field = NULL, fun = 'count')
ed <- fasterize(iucn.2, r, field = 'w', fun = 'max')
writeRaster(ed, 'BatDistinctiveness.tif')
