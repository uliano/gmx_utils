library(ggplot2)
library(scales)

source('gromacs.R')

twocol<-scales::hue_pal()(2)
fourcols<-scales::hue_pal()(4)
#system('./generate_files.sh')
theme_set(theme_gray(base_size = 18))

span<-0.1 # decide quanto sono smoothed le linee dei grafici (+ grande = + smooth)

directories=c('S1PR1_Lys+S1P', 'S1PR1_Lys+FP')
file_name<-'saltbr-ASP279-LYS200.dat'
category_name='Ligand'
categories<-c('S1P','FP')
col_names <-c('time','distance')
df1 <-read_many_files(directories, file_name, category_name, categories, col_names, dt=100.0)
df1$interaction<-'279-200'
file_name<-'saltbr-GLU250-ARG78.dat'
df2 <-read_many_files(directories, file_name, category_name, categories, col_names, dt=100.0)
df2$interaction<-'250-78'
file_name<-'saltbr-GLU250-LYS75.dat'
df3 <-read_many_files(directories, file_name, category_name, categories, col_names, dt=100.0)
df3$interaction<-'250-75'
file_name<-'saltbr-GLU294-LYS41.dat'
df4 <-read_many_files(directories, file_name, category_name, categories, col_names, dt=100.0)
df4$interaction<-'294-41'
df<-rbind(df1,df2,df3,df4)
df$distance <- df$distance/10.0
p <- ggplot(df, aes(x=time, y=distance, color=Ligand)) +
  geom_point(alpha=0.2) +
  labs(title='Salt Bridges', x='Time (ps)', y='Distance (nm)') +
  geom_smooth(method="loess",span = span,se=F) +
  guides(colour = guide_legend(override.aes = list(alpha=1)))+
  facet_grid(interaction ~ .)
print(p)

file_name <-'mindist_i138_v258.xvg'
df1 <-read_many_files(directories, file_name, category_name, categories, col_names)
df1$interaction<-'3x46 - 6x37'
file_name <-'mindist_i138_y311.xvg'
df2 <-read_many_files(directories, file_name, category_name, categories, col_names)
df2$interaction<-'3x46 - 7x53'
df<-rbind(df1,df2)
p <- ggplot(df, aes(x=time, y=distance, color=Ligand)) +
  geom_point(alpha=0.2) +
  labs(title='Minimum Distance', x='Time (ps)', y='Distance (nm)') +
  geom_smooth(method="loess",span = span,se=F) +
  guides(colour = guide_legend(override.aes = list(alpha=1)))+
  facet_grid(interaction ~ .)
print(p)

file_name <-'ncontact_i138_v258.xvg'
df1 <-read_many_files(directories, file_name, category_name, categories, col_names)
df1$interaction<-'3x46 - 6x37'
file_name <-'ncontact_i138_y311.xvg'
df2 <-read_many_files(directories, file_name, category_name, categories, col_names)
df2$interaction<-'3x46 - 7x53'
df<-rbind(df1,df2)
p <- ggplot(df, aes(x=time, y=distance, color=Ligand)) +
  geom_point(alpha=0.2) +
  labs(title='Number of Contacts < 0.6 nm', x='Time (ps)', y='Number') +
  geom_smooth(method="loess",span = span,se=F) +
  guides(colour = guide_legend(override.aes = list(alpha=1)))+
  facet_grid(interaction ~ .)
print(p)

directory <- 'hbonds_tm'
file_name <- 'hbond-HIS28-CYS184.csv'
col_names <- c('time','distance','consenus')
df1 <-read_many_files(directories, file_name, category_name, categories, col_names, directory=directory)
df1$interaction='HIS28-CYS184'
file_name <- 'hbond-TYR29-TRP117.csv'
df2 <-read_many_files(directories, file_name, category_name, categories, col_names, directory=directory)
df2$interaction='TYR29-TRP117'

df<-rbind(df1,df2)

directories=c('S1PR1_Lys+S1P', 'S1PR1_Lys+FP')
directory <- 'hbonds_tm'
bonds = c('LYS200-SER190','LYS200-ASP279','LYS200-VAL280','LYS200-CYS282')
delta=-0.1
base=0.2
diam = 1
category_name='Ligand'
categories<-c('S1P','FP')

p<-plot_hbonds(directories, directory, bonds, category_name, categories)
print(p)

bonds<-c('HIS77-TRP71')
p<-plot_hbonds(directories, directory, bonds, category_name, categories)
print(p)

saltbr<-'saltbr'


bridges<-c('ASP18-LYS111')
p<-plot_saltbr(directories, saltbr, bridges, category_name, categories, dt=100)
print(p)

bridges<-c('ASP40-ARG292')
p<-plot_saltbr(directories, saltbr, bridges, category_name, categories, dt=100)
print(p)

bridges<-c('GLU294-LYS34')
p<-plot_saltbr(directories, saltbr, bridges, category_name, categories, dt=100)
print(p)

bridges<-c('GLU317-ARG78')
p<-plot_saltbr(directories, saltbr, bridges, category_name, categories, dt=100)
print(p)

bridges<-c('GLU317-ARG247')
p<-plot_saltbr(directories, saltbr, bridges, category_name, categories, dt=100)
print(p)
