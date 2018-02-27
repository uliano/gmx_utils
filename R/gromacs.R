require(tools)

unquote <- function(s){
  pos <- gregexpr(pattern = '"', s)
  substr(s, head(pos[[1]], n=1) +1 ,tail(pos[[1]], n=1) -1)
} 

uncomment <-function(s) {
  pos1 <- gregexpr(pattern = '/*', s)
  pos2 <- gregexpr(pattern = '*/', s)
  substr(s, head(pos1[[1]], n=1) +2 ,tail(pos2[[1]], n=1) -2)
}

col_decode <-function(c){
  color <- strsplit(c,'/*',fixed=TRUE)[[1]]
  value <- trimws(unquote(color[2]))
  tmp <- strsplit(unquote(color[1]), ' c ', fixed=TRUE)[[1]]
  key <- trimws(tmp[1])
  color <- trimws(tmp[2])
  list(key=key, color=color, value=value)
}

process_meta <- function(meta){
  uncommented <- character(0)
  for (i in meta) {
    j <- uncomment(i)
    if (length(j)== 0) {
      next
    }
    if (j[1] == ' ') {
      uncommented[lenght(uncommented)] <- paste(uncommented[lenght(uncommented)], j, sep='')
    } else if (grepl(":",j) ) {
      uncommented<-c(uncommented, j)
    }
  }
  result <- list()
  for (i in uncommented) {
    index <- regexpr(':',i)
    key <- substr(i, 1, index-1)
    value <- substr(i, index=1, nchar(i))
    if (substr(value, 1,1) == '"') {
      value <- unquote(value)
    }
    if (key %in% result) {
      result[key] <- paste(result[key],value, sep=' ')
    }
    result[key]<-value
  }
  result
}

read_xpm <- function(file_name) {
  f <- file(file_name,'r')
  meta <- character(0)
  line <- readLines(f, n=1, ok=FALSE)
  print(line)
  while (! startsWith(line, "static char *gromacs_xpm[]")){
    meta <- c(meta, line)
    line <- readLines(f, n=1, ok=FALSE)
    print(line)
  }
  dim <- as.integer(strsplit(unquote(readLines(f, n=1, ok=FALSE)),' +')[[1]])
  nx <- dim[1]
  ny <- dim[2]
  nc <- dim[3]
  nb <- dim[4]
  
  values<-list()
  colors<-list()
  for (i in 1:nc){
    line <- readLines(f, n=1, ok=FALSE)
    decoded <- col_decode(line)
    key <- decoded['key'][[1]]
    value <- decoded['value'][[1]]
    color <- decoded['color'][[1]]

    values[key] <- value
    colors[key] <- color
  }

  lines <- readLines(f)
  df <- data.frame(x=integer(), y=integer(), color=character(), value<-character())
  
  mat <-matrix(nrow=ny, ncol=nx)
  y<-ny
  for (line in lines){
    if (startsWith(line,'/*')){
      meta <-c(meta,line)
      next
    }
    j <- unquote(line)
    for (x in seq(from=1, to=nx*nb, by=nb)){
      code<-substr(j,1+(x-1)*nb,x*nb)
      color <- colors[code][[1]]
      value <- values[code][[1]]
      rgb <- col2rgb(color)
      df1 <- data.frame(x=x, y=y, color=color, value=value, rgb=rgb)
      df<-rbind(df,df1)
      mat[x,y]<-as.numeric(value)
      print(x)
      print(y)
    }
    y<-y-1
  }
  meta <- process_meta(meta)
  list(data=df, meta=meta, matrix=mat)
}


read_gro <-function(file_name) {
  f <- file(file_name)
  lines <- readLines(f, n = 3, warn = FALSE)
  line_len <- nchar(lines[3])
  close(f)
  n_atoms=as.integer(lines[2])
  df <- read.fwf(file_name, c(5,5,5,5,line_len-20), skip=2, n=n_atoms, col.names=c('rnum','rname','aname','anum_trunc','waste'))
  df$anum<-seq(nrow(df))
  if(any((df$anum %% 100000) != df$anum_trunc)) {
    stop("Error in ",file_name," atom numbering!")
  }
  df <- df[c('rnum','rname','anum','aname')]
}

read_many_files <- function(directories, file_name, category_name, categories, col_names,
                            directory=NULL, t0=0.0, dt=NULL ){
  if (length(directories) != length(categories)){
    stop("the number of directories:", directories, "should be equal to the number of categories:", categories)
  }
  if(is.null(directory)) {
    abs_directories <- file.path(directories)
  } else {
    abs_directories <- file.path(directories, directory)
  }
  existing <- file.exists(abs_directories)
  existing_directories <- abs_directories[existing]
  if(any(!existing)){
    missing <- abs_directories[!existing]
    stop("missing: ", missing)
  }
  files <- file.path(existing_directories,file_name)
  existing <- file.exists(files)
  files <- files[existing]
  categories <- categories[existing]
  
  result <- NULL
  
  for (i in 1:length(files)){
    extension <- file_ext(file_name)
    if (extension == 'dat'){
      if (is.null(dt)){
        stop("for .dat files dt is required!")
      }
      df <-read.table(files[i],skip=1,col.names=c('frame','value'))
      value <- df$value
      time=seq(from=t0, by=dt, along.with=value)
      const <- rep(categories[i],nrow(df))
      df <- data.frame(time, value, const)
      names(df) <- c(col_names, category_name)
    } else if (extension == 'xvg') {
      df <- read_xvg(files[i], column_names = col_names, const_column=c(category_name,categories[i]))
    } else if (extension == 'csv') {
      df <- read.table(files[i], sep=',', skip=1, col.names = col_names)
      df[[category_name]] <- rep(categories[i],nrow(df))
    } else {
      stop("unknown extension:", extension)
    }
    if (is.null(result)) {
      result <- df
    } else {
      result<- rbind(result, df)
    }
  }
  result
}

read_xvg <- function(file_name, column_names = NULL, const_column = NULL) {
  f <- file(file_name)
  lines <- readLines(f, n = 5000, warn = FALSE)
  close(f)
  meta_indexes <- grep("^[@|#]", lines)
  meta <- lines[meta_indexes]
  title_index <- grep(' title ', meta)
  title <- scan(text = meta[title_index], what = "", quiet = TRUE)[3]
  xaxis_index <- grep(' xaxis ', meta)
  xaxis <- scan(text = meta[xaxis_index], what = "", quiet = TRUE)[4]
  yaxis_index <- grep(' yaxis ', meta)
  col_names <- c(scan(text = xaxis, what = "", quiet = TRUE)[1])
  yaxis <- scan(text = meta[yaxis_index], what = "", quiet = TRUE)[4]
  type_index <- grep('TYPE ', meta)
  type <- scan(text = meta[type_index], what = "", quiet = TRUE)[2]
  series_indexes <-grep('@ s[0-9]',meta)
  n_series <- length(series_indexes)
  if (n_series == 0) {
    col_names <- c(col_names, scan(text = meta[yaxis_index], what = "", quiet = TRUE)[4])
  } else {
    for (index in series_indexes) {
      col_names <- c(col_names, scan(text = meta[index], what = "", quiet = TRUE)[4])
    }
  }
  if (! is.null(column_names)) {
    col_names <- column_names
  }
  df <- read.table(file_name, col.names = col_names, skip = length(meta_indexes))
  if (! is.null(const_column)) {
    df[[as.character(const_column[[1]])]] <- const_column[[2]]
  }
  attr(df, "title") <- title
  attr(df, "xlabel") <- xaxis
  attr(df, "ylabel") <- yaxis
  attr(df, "type") <- type
  df
}

read_2xvg <- function(file1, file2, column_names = NULL, var_col=NULL, const_col=NULL) {
  if(! is.null(var_col)) {
    name <- var_col[[1]]
    value1 <- var_col[[2]]
    value2 <- var_col[[3]]
    cc1 <- list(name,value1) 
    cc2 <- list(name,value2)
  } else {
    cc1 <- NULL
    cc2 <- NULL
  }
  df1 <- read_xvg(file1, column_names = column_names, const_column = cc1)
  df2 <- read_xvg(file2, column_names = column_names, const_column = cc2)
  df<-rbind(df1,df2)
  if(! is.null(const_col)) {
    df[[as.character(const_col[[1]])]] <- const_col[[2]]
  }
  df
}

make_labels <- function(df, extra_title=NULL) {
  title <- attr(df, "title")
  if (! is.null(extra_title)) {
    title <- paste(title, extra_title, sep=" ") 
  }
  labs(title=title, x=attr(df,"xlabel"),y=attr(df,"ylabel"))
}

pplot_hbonds <- function(df, category_name, delta=0.025) {
  df$cons <- ifelse(df$consenus == 'Present', 1.0, NA)
  categories <- levels(factor(df[[category_name]]))
  for (i in seq_along(categories)) {
    df$cons <- ifelse(df[[category_name]] == categories[i], df$cons*(0.2 - delta*(i-1)), df$cons)
  }
  p <- ggplot(df, aes(x=time, y=distance, color=Ligand)) +
    geom_point(alpha=0.2) +
    labs(title='H Bonds', x='Time (ps)', y='Distance (nm)') +
    geom_smooth(method="loess",span = span,se=F) +
    geom_point(aes(y=cons))+
    guides(colour = guide_legend(override.aes = list(alpha=1)))+
    facet_grid(interaction ~ .)
}


plot_hbonds <- function(directories, directory, bonds, category_name, categories, delta=-0.05, base=0.15, diam = 1, title=NULL) {
  the_title <- title
  df <- NULL
  file_names <- paste('hbond',bonds,sep='-')
  file_names <- paste(file_names, 'csv', sep='.')
  col_names <- c('time','distance','consenus')
  for (i in seq_along(file_names)) {
    dftemp <- read_many_files(directories, file_names[i], category_name, categories, col_names, directory=directory)
    dftemp$interaction=bonds[i]
    if (is.null(df)) {
      df <- dftemp
    } else {
      df <- rbind(df, dftemp)
    }
  }
  df$cons <- ifelse(df$consenus == 'Present', 1.0, NA)
  df[[category_name]]<-factor(df[[category_name]])
  df[[category_name]]<-factor(df[[category_name]], levels=sort(levels(df[[category_name]])))
  
  categories <- levels(factor(df[[category_name]]))
  for (i in seq_along(categories)) {
    df$cons <- ifelse(df[[category_name]] == categories[i], df$cons*(base + delta*(i-1)), df$cons)
  }
  p <- ggplot(df, aes(x=time, y=distance, color=Ligand)) +
    geom_point(alpha=0.2) +
    labs(x='Time (ps)', y='Distance (nm)') +
    geom_smooth(method="loess",span = span,se=F) +
    geom_point(aes(y=cons), alpha=0.05, size =diam, na.rm= TRUE)+
    guides(colour = guide_legend(override.aes = list(alpha=1)))+
    coord_cartesian(ylim=c(0,1))
  
  if (length(bonds) == 1) {
    if (is.null(title))
      the_title <- paste('H Bond') #, bonds[1])
  } else {
    the_title <- paste('H Bonds')
     
  }
  p <- p + facet_grid(interaction ~ .)
  p <- p + labs(title=the_title)
}

plot_saltbr <- function(directories, directory, bridges, category_name, categories, dt=NULL, title=NULL) {
  the_title <- title
  df <- NULL
  file_names <- paste('saltbr',bridges ,sep='-')
  file_names <- paste(file_names, 'dat', sep='.')
  col_names <- c('time','distance')
  for (i in seq_along(file_names)) {
    dftemp <- read_many_files(directories, file_names[i], category_name, categories, col_names, directory=directory, dt=dt)
    dftemp$interaction=bridges[i]
    if (is.null(df)) {
      df <- dftemp
    } else {
      df <- rbind(df, dftemp)
    }
  }
  df$distance <- df $ distance/10.0
  df[[category_name]]<-factor(df[[category_name]])
  df[[category_name]]<-factor(df[[category_name]], levels=sort(levels(df[[category_name]])))
  
  categories <- levels(factor(df[[category_name]]))
  p <- ggplot(df, aes(x=time, y=distance, color=Ligand)) +
    geom_point(alpha=0.2) +
    labs(x='Time (ps)', y='Distance (nm)') +
    geom_smooth(method="loess",span = span,se=F) +
    guides(colour = guide_legend(override.aes = list(alpha=1)))+
    coord_cartesian(ylim=c(0,1))
  
  if (length(bridges) == 1) {
    if (is.null(title))
      the_title <- paste('Salt Bridge') #, bridges[1])
  } else {
    the_title <- paste('Salt Bridges')
     
  }
  p <- p + facet_grid(interaction ~ .)
  p <- p + labs(title=the_title)
}


