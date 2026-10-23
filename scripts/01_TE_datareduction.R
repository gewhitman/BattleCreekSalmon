#Trace Element data reduction

# Prepare the environment -------------------------------------------------
library(tidyverse)
library(here) #better working directory management
library(stringr) #working with strings
library(viridis)#Color scale
library(cowplot) #GGplot extensions
library(readxl)
library(openxlsx)
library(zoo)

#Clear the workspace
rm(list = ls())

# Read in data ------------------------------------------------------------
#read in exog scores
exog_readin <- read.xlsx(here('data','battle_creek_exog_scores.xlsx'), sheet="exog")%>%
  dplyr::select(lab_id, exog_score)

#Iolite timeseries data from UC Davis
paths_d<-list.files(path = here("data","TE_Davis","timeseries_data"),full.names = TRUE, recursive = TRUE)
filenames_d<- list.files(path = here("data","TE_Davis","timeseries_data"), recursive = TRUE)
iolite_readin_d <-data_frame(filename = filenames_d) %>% 
  mutate(file_contents = map(paths_d, ~read_csv (., skip = 1)))%>%
  unnest(cols = c(file_contents))

#Iolite timeseries data from Essex
paths_e<-list.files(path = here("data","TE_Essex","timeseries_data"),full.names = TRUE, recursive = TRUE)
filenames_e<- list.files(path = here("data","TE_Essex","timeseries_data"), recursive = TRUE)
iolite_readin_e <-data_frame(filename = filenames_e) %>% 
  mutate(file_contents = map(paths_e, ~read_csv (., skip = 1)))%>%
  unnest(cols = c(file_contents))

# #Output filenames for metadata
# filenames_out <- filenames_e %>%
#   data.frame()
# write.csv(filenames_out, "filenames_out.csv", row.names=F)

iolite_d <- iolite_readin_d %>%
dplyr::select(filename,'Elapsed Time', Mg24_ppm, Ca44_ppm, Mn55_ppm, Zn66_ppm, Sr88_ppm, Ba137_ppm)%>%
  mutate(Lab="UCDavis")%>%
  mutate(Join_ID=str_sub(filename,1,-5))

iolite_e <- iolite_readin_e %>%
  dplyr::select(filename,'Elapsed Time', Mg24_ppm, Ca44_ppm, Mn55_ppm, Zn66_ppm, Sr88_ppm, Ba137_ppm)%>%
  mutate(Lab="Essex")%>%
  mutate(Join_ID=str_sub(filename,1,-5))

iolite_data <-rbind(iolite_d,iolite_e)

#Laser runnotes
te_notes <- read.xlsx(here::here("data", "TSC_trace_elements_sample_notes_20240420.xlsx"),
                            sheet = "sample_notes")

# Clean data --------------------------------------------------------------
TE_data <- iolite_data%>%
  #Join runnotes
  left_join(te_notes, by="Join_ID")%>%
  #Calculate time interval and check 0.97
  mutate(timedif=`Elapsed Time` - lag(`Elapsed Time`))%>%
  drop_na(Fish_ID)%>%
  #Calculate distance
  group_by(run_id)%>%
  mutate(Distance=row_number()*0.97*run_speed_ums)%>%
  #Determine natal region
  mutate(Region=case_when(Distance >=Natal_start & Distance <=Natal_stop ~"Natal",
                          TRUE~"Profile"))%>%
  #Select columns to keep
  dplyr::select(Fish_ID,run_id, Distance, Natal_start, Natal_stop, Origin, Lifestage, Region, Lab,analysis_no, Comparison,
                Mg24_ppm, Ca44_ppm, Mn55_ppm, Zn66_ppm, Sr88_ppm, Ba137_ppm)%>%
 #Calculate ratios to Ca44 based on the ppm values
  mutate (Mg_Ca =Mg24_ppm/Ca44_ppm) %>%
  mutate (Mn_Ca =Mn55_ppm/Ca44_ppm) %>%
  mutate (Zn_Ca =Zn66_ppm/Ca44_ppm) %>%
  mutate (Sr_Ca =Sr88_ppm/Ca44_ppm) %>%
  mutate (Ba_Ca =Ba137_ppm/Ca44_ppm) %>%
  filter(Sr_Ca>0)%>%
  # Change to long format
  pivot_longer(names_to="Element", values_to="Value", 12:22)%>%
  #Remove fish with no data
  filter(!Fish_ID %in%c("BUT081612CS08",
                       "BUT092718CS03",
                       "BUT092518CS05",
                       "BUT100418CS01",
                       "BUT091818CS03",
                       "BUT091818CS02",
                       "AL040814-2"
                       
    
  ))%>%
  group_by(run_id, Element)%>%
  #Outlier rejection
  mutate(median = rollapply(Value, width = 20, FUN = median, partial = T)) %>%
  mutate(upper75prob = rollapply(Value, width = 20, FUN = quantile, partial = T, probs = 0.75, na.rm=TRUE),
         lower25prob = rollapply(Value, width = 20, FUN = quantile, partial = T, probs = 0.25, na.rm=TRUE),
         IQR = upper75prob - lower25prob,
         diff_flag = (Value - median)/IQR,
         despiked = ifelse(abs(diff_flag) > 2, median, Value)) %>%
  ungroup()

write.csv(TE_data, "outputs/TE_data.csv", row.names=F)

#Sample table
TE_samples_table <-TE_data %>%
  distinct(Fish_ID, run_id, analysis_no, Lab,Origin, Lifestage)%>%
  group_by(Lab,Origin, Lifestage)%>%
  mutate(n_fish=n())%>%
  distinct(Lab,Origin, Lifestage, n_fish)
write.csv(TE_samples_table, "outputs/TE_samples_table.csv", row.names=F)

#Plot with outliers for each measured element and ratio. Each fish on its own page. 
#Filter which elements to plot
element_plotter <-TE_data %>%
  filter(Element%in% c("Mg_Ca", "Mn_Ca", "Zn_Ca","Sr_Ca","Ba_Ca"))

fish_plotter <- unique(element_plotter$run_id)
fish_plotter #Show list of fish that will be plotted

#Create panel plot for each fish
pdf("outputs/BAT_TE.pdf",16, 8)
for(i in fish_plotter){
  plot_outlier<- element_plotter %>% filter (run_id==i)
  p <- ggplot(data=plot_outlier)+
    geom_point(aes(x=Distance,y=despiked, color=Element))+
    geom_smooth(aes(x=Distance,y=despiked, color=Element), method="loess", span = 0.1)+
    geom_vline(aes(xintercept=Natal_start))+
    geom_vline(aes(xintercept=Natal_stop))+
    facet_wrap(~Element, scales = "free")+
    theme_bw() + ggtitle(i)
  print(p)
}
dev.off()


#Lab comparison
TE_data_comb <- TE_data %>%
  filter(Element%in% c("Sr_Ca"))%>%
  filter(Comparison=="y")

#Plot comparison fish Sr/Ca
p_comp_allele <- ggplot(data=TE_data_comb)+
  geom_smooth(aes(x=Distance,y=despiked, group=run_id, color=Lab), span=0.2, se=F)+
  theme_bw() + 
  #Limit x axis to 1000um to make comparison more readable (Essex went to 1000 out)
  scale_x_continuous(limits=c(0,1000))+
  theme(legend.position = "top")+
  facet_wrap(Element~Fish_ID)
p_comp_allele
ggsave(plot=p_comp_allele, "outputs/comp_allele.png", width = 10, height = 10)


# Mean natal values -------------------------------------------------------
TE_natal <-element_plotter%>%
  filter(analysis_no=="1")%>%
  group_by(Fish_ID, Element, Region)%>%
  mutate(Value_mean=mean(despiked, na.rm=T))%>%
  filter(Region=="Natal")%>%
  ungroup()%>%
  distinct(Fish_ID, Origin, Lab, Lifestage,Element, Value_mean)


#Plot all juvenile fish all Elements by origin
p_allfish_org <- ggplot(data=TE_natal%>%filter(Lifestage=="juvenile"))+
  geom_boxplot(aes(x=Origin,y=Value_mean, fill=Origin), width=0.5)+
  geom_point(aes(x=Origin,y=Value_mean), alpha=0.5)+
  # scale_fill_manual(values=c("orange","steelblue", "darkolivegreen"))+
  theme_classic() + 
  theme (panel.background = element_rect(colour = "black"), legend.position="none")+
  ylab("Mean Natal Element/Ca")+
  facet_wrap(~Element, scales = "free")
p_allfish_org
ggsave(plot=p_allfish_org, "outputs/juvfish_org.png", width = 8, height = 4)


#Sr vs Ba plot
TE_natal_srba <-TE_natal %>%
  filter(Element %in% c("Ba_Ca","Sr_Ca"))%>%
  pivot_wider(names_from=Element, values_from=Value_mean)


#Plot only the Sr overlapping sources juvenile  fish by origin Sr vs Ba
p_assign_srba <- ggplot(data=TE_natal_srba%>% 
                          filter(Lifestage=="juvenile")
                        )+
  geom_point(aes(y=Ba_Ca,x=Sr_Ca, fill=Origin), shape=21, size=3,alpha=0.7) + 
  scale_fill_brewer(palette = "Set1")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position='right')+
  xlab("Sr/Ca")+
  ylab("Ba/Ca")
p_assign_srba
ggsave(plot=p_assign_srba, "outputs/p_assign_srba.png", width = 5, height = 3)

write.csv(TE_natal_srba, "outputs/TE_natal_srba.csv", row.names=F)


# Preparing juvenile data
TE_natal_srba_juv_check <- TE_natal%>%
  filter(Lifestage=="juvenile")%>%
  mutate(Site=factor(Origin))%>%
  mutate(lab_id=Fish_ID)%>%
  left_join(exog_readin, by="lab_id")

write.csv(TE_natal_srba_juv_check, "outputs/TE_natal_srba_juv_check.csv", row.names=F)


