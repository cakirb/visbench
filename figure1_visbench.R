setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

used.packages <- c("cowplot", "dplyr", "stringr", "tidyr", "ggplot2")
new.packages <- used.packages[!(used.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {install.packages(new.packages)}
lapply(used.packages, require, character.only = TRUE)


### cellxgene
data <- read.csv('data/cellxgene.csv')
names(data) <- sub("^X", "", names(data))

ram_cxg <- data %>% slice(1:10) %>% select(-1) %>%  
  "/"(1024) %>% slice(2,4,7,8,10) %>% stack() %>% 
  rename(memory = values, dataset = ind) %>% 
  mutate(memory = memory/1024)
tim_cxg  <- data %>% slice(11:20) %>% select(-1) %>% 
  slice(2,4,7,8,10) %>% stack() %>% 
  rename(time = values, dataset = ind)


### isee
# isee_sce
data <- read.csv('data/isee_sce.csv')
names(data) <- sub("^X", "", names(data))

ram_isee_sce <- data %>% slice(1:5) %>% select(-1) %>%
  stack() %>% rename(memory = values, dataset = ind) %>%
  mutate(memory = memory/1024)
tim_isee_sce  <- data %>% slice(6:10) %>% select(-1) %>%
  stack() %>% rename(time = values, dataset = ind)

# isee_loom
data <- read.csv('data/isee_loom.csv')
names(data) <- sub("^X", "", names(data))

ram_isee_loom <- data %>% slice(1:5) %>% select(-1) %>%
  stack() %>% rename(memory = values, dataset = ind) %>%
  mutate(memory = memory/1024)
tim_isee_loom  <- data %>% slice(6:10) %>% select(-1) %>%
  stack() %>% rename(time = values, dataset = ind)


### scsva
data <- read.csv('data/scsva.csv') %>% select(-X)
data$name <- str_remove(data$name, '_')
data$name <- factor(data$name, levels = unique(data$name))

ram_scsva <- data %>% select(name,memory) %>%
  rename(dataset = name) %>% mutate(memory = memory/1024)
rt_scsva <- data %>% select(name,time) %>%
  rename(dataset = name)


### loom-viewer
data <- read.csv('data/loo.csv', header = TRUE)
names(data) <- sub("^X", "", names(data))

ram_loom <- data %>% slice(1:5) %>% stack() %>%
  rename(memory = values, dataset = ind) %>% mutate(memory = memory/(1024*1024))
usr_loom <- data %>% slice(6:18) %>% slice(seq(1,20,3))
sys_loom <- data %>% slice(7:19) %>% slice(seq(1,20,3))
usys_loom<- usr_loom + sys_loom
usys_loom<- gather(as.data.frame(usys_loom), "dataset", "time")
rt_loom <- data %>% slice(8:20) %>% slice(seq(1,20,3)) %>% stack() %>% rename(time = values, dataset = ind)


### ucsc cell browser
data <- read.csv('data/ucsccb.csv', header = TRUE)
names(data) <- sub("^X", "", names(data))

ram_ucb <- data %>% slice(1:5) %>% stack() %>% 
  rename(memory = values, dataset = ind) %>% mutate(memory = memory/(1024*1024))
usr_ucb <- data %>% slice(6:10) %>% select(-1) %>% stack() %>% rename(time = values, dataset = ind) 
sys_ucb <- data %>% slice(11:15) %>% select(-1) %>% stack() %>% rename(time = values, dataset = ind) 
usys_ucb <- data %>% slice(16:20) %>% select(-1) %>% stack() %>% rename(time = values, dataset = ind) 


### single cell browser
data <- read.csv('data/scb.csv', header = TRUE)
names(data) <- sub("^X", "", names(data))

ram_scb <- data %>% slice(1:5) %>% stack() %>% 
  rename(memory = values, dataset = ind) %>% mutate(memory = memory/(1024*1024))
usr_scb <- data %>% slice(6:10) %>% select(-1) %>% stack() %>% rename(time = values, dataset = ind) 
sys_scb <- data %>% slice(11:15) %>% select(-1) %>% stack() %>% rename(time = values, dataset = ind) 
usys_scb <- data %>% slice(16:20) %>% select(-1) %>% stack() %>% rename(time = values, dataset = ind) 


### SCope
data <- read.csv('data/scope.csv', header = TRUE)
names(data) <- sub("^X", "", names(data))

ram_scope <- data %>% slice(1:5) %>% stack() %>% 
  rename(memory = values, dataset = ind) %>% mutate(memory = memory/(1024*1024))
usr_scope <- data %>% slice(6:10) %>% select(-1) %>% stack() %>% rename(time = values, dataset = ind)  
sys_scope <- data %>% slice(11:15) %>% select(-1) %>% stack() %>% rename(time = values, dataset = ind) 
usys_scope <- data %>% slice(16:20) %>% select(-1) %>% stack() %>% rename(time = values, dataset = ind) 




# merging results
tool <- c(rep('cellxgene',nrow(ram_cxg)), 
          rep('iSEE-SCE',nrow(ram_isee_sce)),
          rep('iSEE-loom',nrow(ram_isee_loom)),
          rep('scSVA',nrow(ram_scsva)),
          rep('loom-viewer',nrow(ram_loom)),
          rep('UCSC Cell Browser',nrow(ram_ucb)),
          rep('Single Cell Explorer',nrow(ram_scb)),
          rep('SCope',nrow(ram_scope))
          )
ram_all <- rbind(ram_cxg, ram_isee_sce, ram_isee_loom, ram_scsva, ram_loom, ram_ucb, ram_scb, ram_scope)
ram_all <- cbind(ram_all,tool)
ram_all <- ram_all %>% group_by(dataset, tool) %>% mutate(se = sd(memory)/sqrt(length(memory)), mean = mean(memory))

tim_all <- rbind(tim_cxg, tim_isee_sce, tim_isee_loom, rt_scsva, usys_loom, usys_ucb, usys_scb, usys_scope)
tim_all <- cbind(tim_all,tool)
tim_all <- tim_all %>% group_by(dataset, tool) %>% mutate(se = sd(time)/sqrt(length(time)), mean = mean(time)) 

dataset_names <- ram_all$dataset
dataset_names <- str_to_upper(dataset_names)
ram_all$dataset <- as.character(ram_all$dataset)
tim_all$dataset <- as.character(tim_all$dataset)

realval <-  c("5000","10000","25000","50000","100000","250000","500000","1000000","1500000","2000000")
k <- 1
for (i in unique(as.character(dataset_names))){
  ram_all$dataset[ram_all$dataset == i] <- realval[k]
  tim_all$dataset[tim_all$dataset == i] <- realval[k]
  k <- k+1
}
ram_all$dataset <- as.numeric(ram_all$dataset)
tim_all$dataset <- as.numeric(tim_all$dataset)

rm(list=setdiff(ls(), c("ram_all","tim_all","dataset_names","tool")))




# Fig.1a: Maximum RAM usage (x/y axes in log-scale)
p1 <- ggplot(ram_all, aes(x = dataset, y = mean, color = tool)) + geom_point(size = 2.5, position=position_dodge(w=0.02)) + 
  geom_path(size = 0.5) + labs(x = 'Number of cells', y = 'Preprocessing memory (GB)', color = 'Tools') + scale_y_log10() +
  scale_x_continuous(trans = "log10", labels = dataset_names, breaks = ram_all$dataset) + theme_bw() + 
  theme(text = element_text(size=18), axis.text.x=element_text(angle = -45, hjust = 0), axis.text.y = element_text(vjust = 0.1)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_brewer(palette = "Set2") + 
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se))

# Fig.1b: Start-up time (x/y axes in log-scale) 
p2 <- ggplot(tim_all, aes(x = dataset, y = mean, color = tool)) + geom_point(size = 2.5, position=position_dodge(w=0.02)) + 
  geom_path(size = 0.5) + labs(x = 'Number of cells', y = 'Preprocessing time (minutes)', color = 'Tools') + scale_y_log10() + 
  scale_x_continuous(trans = "log10", labels = dataset_names, breaks = tim_all$dataset) + theme_bw() + 
  theme(text = element_text(size=18), axis.text.x=element_text(angle = -45, hjust = 0), axis.text.y = element_text(vjust = 0.1)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_brewer(palette = "Set2") +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se))

# combining two plots
prow <- plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  align = 'vh',
  labels = c("a", "b"),
  label_size = 20,
  label_y = 1.03,
  nrow = 2)
legend_b <- get_legend(p1 + guides(color = guide_legend(nrow = 2)) + theme(legend.position = "bottom", text = element_text(size=18)))
plot_grid(prow, legend_b, ncol=1, rel_heights = c(3, .3)) # Give it 0.1/3 of the width of one plot (via rel_widths).

ggsave("figure1_visbench.pdf", h = 10)
ggsave("figure1_visbench.png", h = 10)