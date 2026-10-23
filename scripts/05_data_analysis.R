#Data analysis
# Prepare the environment -------------------------------------------------
library(tidyverse)
library(here)
library(stringr)
library(viridis)
library(cowplot)
library(rstatix)
library(openxlsx)
library(nlstools) 
library(FSA)
library(wesanderson)
library(ggpmisc)
library(ggpubr)
library(boot)
library(ggrepel)

#Clear the workspace
rm(list = ls())


# Define parameters -------------------------------------------------------
colors <- c("BAT" = "#8DD3C7", "CNH" = "#B3DE69", 
            "DEE" = "#BEBADA", "MIL" = "#FB8072", 
            "SAC" = "#80B1D3", "BUT" = "#FDB462",
            "FEH" = "#FFFFB3", "NIH" = "#FCCDE5", 
            "MOK" = "#D9D9D9", "MOH" = "#BC80BD", 
            "STA" = "#CCEBC5", "TUO" = "#FFED6F",
            "OTH"="black")

#set theshold values for Natal and FW exit 
Sr8786_FWexit = 0.70785 # Mean Chipps value used in Phillis et al and Morais et al papers
Sr8786_Natalexit=0.70431 # Value of Battle creek

# Read in data ------------------------------------------------------------
#read Sr profiles
battle_sr <-read.csv(here('data','battle_data.csv'),  stringsAsFactors = F)

#read sr natal predictions
natal_data <- read.csv(here('outputs','battle_adult_assign_cwt_corrected.csv'),  stringsAsFactors = F) %>% 
  # filter out the unknown-origin sample (UNK) - adipose fin absent & CWT info absent & wild assigned by otolith
  filter(!.pred_class_final_cwt == "UNK") %>% 
  mutate(assign= case_when( .pred_class_final_cwt == "STA" ~ "OTH",
                           .pred_class_final_cwt == "TUO" ~ "OTH",
                           .pred_class_final_cwt == "MOK" ~ "OTH",
                           .pred_class_final_cwt == "NIH" ~ "OTH",
                           .pred_class_final_cwt == "AME" ~ "OTH",
                           .pred_class_final_cwt == "FEA" ~ "OTH",
                          # .pred_class_final_cwt == "FEH" ~ "OTH",
                           .pred_class_final_cwt == "MOH" ~ "OTH",
                           TRUE~ .pred_class_final_cwt),
         assignment=factor(assign, 
                           levels=c("BAT", "CNH", "DEE","MIL",
                                    "SAC","BUT","FEH","OTH"))) %>% 
  select(lab_id,assignment)

#prop bootstrap dataframe              
prop_boot <- read.csv(here('outputs','prop_boot.csv'),  stringsAsFactors = F) %>% 
  mutate(assign= case_when( NAT_LOC == "STA" ~ "OTH",
                            NAT_LOC == "TUO" ~ "OTH",
                            NAT_LOC == "MOK" ~ "OTH",
                            NAT_LOC == "NIH" ~ "OTH",
                            NAT_LOC == "AME" ~ "OTH",
                            NAT_LOC  == "FEA" ~ "OTH",
                           # NAT_LOC  == "FEH" ~ "OTH",
                            NAT_LOC  == "MOH" ~ "OTH",
                            TRUE~ NAT_LOC),
         year=factor(year),
         assignment=factor(assign,levels=c("BAT", "CNH", "DEE","MIL",
                                           "SAC","BUT","FEH","OTH")))


#Read in environmental data
SacIndexWY <- read.xlsx(here('data','SacIndexWY.xlsx')) %>% 
              mutate(WYT= as.factor(WYT))

#Read in additional data from USFW for Battle Creek
bc_usfw <-read.xlsx(here('data','Battle_Creek_otoliths_USFWS_sub.xlsx'), sheet="overview")

#RMIS
rmis_readin <-read.csv(here('data','CSV19422.csv'),  stringsAsFactors = F)
rmis <-rmis_readin %>%
  dplyr::select(tag_code=tag_code_or_release_id,   brood_year, hatchery_location_name)

#Carcass survey data
carcass_snorkel <- read.csv(here('data','Carcass_Snorkel.csv'),  stringsAsFactors = F)
carcass_snorkel_clipped <- read.csv(here('data','Carcass_Snorkel_clipped.csv'),  stringsAsFactors = F)
carcass_tagcode <- read.csv(here('data','snorkelcarcass_clipped_hatchcode.csv'),  stringsAsFactors = F)


# Work with the data ------------------------------------------------------
#Dataframe with full profiles
battle_sr_profile <- battle_sr %>%
  #Add natal information to the Sr profiles
  left_join(natal_data, by="lab_id")%>%
  #Remove fish with no natal assignments
  drop_na(assignment)%>%
  left_join(bc_usfw, by="lab_id")%>%
  mutate(Sample_type="Adult")%>%
  ungroup()%>%
  #Remove manually entered natal_exit
  dplyr::select(-natal_exit)%>%
  group_by(lab_id)%>%
  #Determine natal exit based on cutoff
   mutate(natal_exit_flag=case_when(
     assignment=="BAT" & Distance_um>=natal_stop & sr_spline >Sr8786_Natalexit ~ "natal_exit"))%>%
   group_by(lab_id,natal_exit_flag)%>%
   mutate(natal_exit=case_when(assignment=="BAT" ~min(Distance_um)))%>%
   ungroup()%>%
  group_by(lab_id)%>%
  #To fill all the rows with natal_exit value
  mutate(natal_exit=max(natal_exit))%>%
  #For the yearlings
  mutate(natal_exit=case_when(
    assignment=="BAT" & max(natal_exit)==0 ~ max(Distance_um),
    TRUE~natal_exit))%>%
  #Fix natal exit for a few profiles where the automatic detection failed
  #This is due to a small extra peak in noisy data
  mutate(natal_exit=case_when(lab_id=="SC_06_50026E" ~ 450,
                              lab_id=="SC_11_50013E" ~ 430,
                              lab_id=="SC_06_50042E" ~ 630,
                              lab_id=="SC_07_50034E" ~ 600,
                              lab_id=="SC_07_50044E" ~ 610,
                              lab_id=="SC_14_50041E" ~ 750,
                              lab_id=="SC_09_50006E" ~ 307, # changed manually natal exit from 320 to 307  
                              lab_id=="SC_13_50059E" ~ 300, # added by flora on 11/04 because value above 2000
                              lab_id== "SC_13_51007E" ~ 400,  # added by flora on 11/04  because value above 2000
                              lab_id== "SC_13_50021E_2" ~ 400, #to remove small bump earlier
    TRUE~natal_exit)) %>%
  #Assign migrant phenotype based on natal exit
  mutate(phenotype_mig = case_when (natal_exit <= 307 ~"Fry",
                                    natal_exit > 307 & natal_exit<= 699 ~"Subyearling",
                                    natal_exit >699  ~ "Yearling")) %>%
  mutate(phenotype_mig_rev=factor(phenotype_mig, levels=c("Yearling", "Subyearling", "Fry"))) %>% 
  ungroup() %>% 
  mutate(age_f=factor(age))%>%
  mutate(age_f_rev = factor(age, levels = c('4','3','2','1')))%>%
  #Determine freshwater exit
  mutate(fw_exit_flag=case_when(Distance_um>natal_exit & sr_spline >Sr8786_FWexit ~ "fw_exit"))%>%
  group_by(lab_id,fw_exit_flag)%>%
  mutate(fw_exit=min(Distance_um))%>%
  ungroup()%>%
  group_by(lab_id)%>%
  #To fill all the rows with fw_exit value
  mutate(fw_exit=max(fw_exit))%>%
  #Calculate amount of rearing
  mutate(natal_rearing_dis=natal_exit-natal_start,
         non_natal_rearing_dis=fw_exit-natal_exit)%>%
  mutate(non_natal_rearing_dis=case_when(non_natal_rearing_dis<0 ~0,
                                         TRUE ~non_natal_rearing_dis))%>%
  mutate(non_natal_rearing_prop=round(non_natal_rearing_dis/(natal_rearing_dis+non_natal_rearing_dis),2))%>%
  ungroup()%>%
  #Calculate outmigration years
  mutate(outmigration_year=year-age+1)%>%
  left_join(SacIndexWY, by = c("outmigration_year" = "WY")) %>%
  mutate(WY=factor(WYT, levels=c("W","AN","BN","D","C")))
write.csv(battle_sr_profile, "outputs/battle_sr_profile.csv", row.names=F)

# Natal origin analyses ---------------------------------------------------
#Dataframe for histograms (1 row per fish)
battle_sr_natal_including_adclip <-battle_sr_profile%>%
  distinct(lab_id,survey_year=year, assignment, age, age_f, age_f_rev, 
           FL_mm=fork_length, sex, adipose_fin, outmigration_year, WY,
           phenotype_mig, natal_exit, fw_exit, non_natal_rearing_prop)%>%
  drop_na(assignment) %>% 
  # Add info about carcass location
  left_join(bc_usfw, by = 'lab_id') %>% 
  mutate(Date = openxlsx::convertToDate(Date),
         Month=month(Date))

write.csv(battle_sr_natal_including_adclip, "outputs/battle_sr_natal_including_adclip.csv", row.names=F)

#Histogram plot of natal origins with ad clipped fish and carcass location
(p_histpredict_loc <- ggplot(data=battle_sr_natal_including_adclip)+
  geom_bar(aes(x=assignment, fill=Subwatershed), stat = "count", color="black")+
  geom_text(stat='count',aes(x=assignment,label=after_stat(count)),vjust=-1)+
  scale_fill_viridis(discrete = T, name="Natal origin", direction=-1)+
  theme_classic()+
  theme (panel.background = element_rect(colour = "black"), legend.position="right")+
  xlab("Predicted natal origin")+
  ylab("Number of Fish")+
    facet_wrap(~adipose_fin)
)
p_histpredict_loc
ggsave(plot=p_histpredict_loc, "outputs/Hist_predict_loc_withadclip.png", width = 15, height = 7)


# Age analyses ------------------------------------------------------------
####Scale and otolith age comparison
age_comparison <- battle_sr_profile %>%
  distinct(lab_id,oto_age=age, scale_age=Scale_Age, assignment,fork_length)%>%
  mutate(scale_age=as.numeric(scale_age))

## Age precision
ap.A <- agePrecision(~oto_age+scale_age,data=age_comparison)
summary(ap.A,what="difference")
summary(ap.A,what="precision")

ab.A <- ageBias(oto_age~scale_age,data=age_comparison)
summary(ab.A,what="symmetry")
plot(ab.A)

#Age size plot / Paper figure Appendix 9
(p_age_size_histo<-ggplot(data=battle_sr_natal_including_adclip  %>% drop_na(age_f) 
                          %>% drop_na(FL_mm)%>%
                            filter((assignment=="BAT"))) +
  geom_histogram(aes(x=FL_mm, fill=age_f_rev), color="black") +
  scale_x_continuous(name = "Fork Length (mm)", limits=c(0,NA)) +
  scale_y_continuous("Fish count")+
  scale_fill_manual(name="Age Class", values = c( "darkblue", "#0868AC","#7BCCC4","#F0F9E8"))+
  theme_bw()+
  theme(legend.position="right", axis.text = element_text(size = 8))
)
ggsave("outputs/age_size_histo.png",p_age_size_histo,  width = 8, height = 4)


#as a stacked plot with proportions
(p_ages_stacked <- ggplot(data= battle_sr_natal_including_adclip  %>% drop_na(age_f)) + 
  geom_bar(aes(x=as.factor(survey_year), fill = age_f_rev), position="fill",color = "black", width = 0.5) +
  scale_fill_manual(name="Age Class", values = c( "darkblue", "#0868AC","#7BCCC4","#F0F9E8"))+
  theme_classic() +
  scale_x_discrete(name="Survey year") + 
  scale_y_continuous(name="Proportion") 
)
ggsave(plot=p_ages_stacked, "outputs/p_ages_stacked.png", width = 8, height = 5)

age_prop <-battle_sr_natal_including_adclip  %>% drop_na(age_f) %>%
  filter(assignment=="BAT")%>%
  mutate(survey_year=as.factor(survey_year)) %>% 
  group_by(survey_year,age)%>%
  mutate(age_count=n())%>%
  distinct(age, age_count)%>%
  ungroup()%>%
  group_by(survey_year) %>% 
  mutate(age_sum=sum(age_count))%>%
  mutate(age_p=age_count/age_sum) %>% 
  ungroup()
write.csv(age_prop, "outputs/age_prop.csv", row.names=F)

age_prop_extended <- age_prop %>%
  complete(survey_year,age = c(2, 3, 4),fill = list(age_p = 0)) %>% 
  select(survey_year,age,age_p) %>% 
  pivot_wider( names_from = age, values_from = age_p, names_prefix = "age_p_" )

# Genotype vs Phenotypes --------------------------------------------------
battle_geno_pheno <-battle_sr_profile %>%
  distinct(lab_id, assignment, adipose_fin, run_pheno=run, run_gen= run_RBFWO_gen, CWT_Code) %>%
  filter(!run_gen=="No Call")%>%
  mutate(run_gen=case_when(run_gen=="SCS"~ "Spring",
                           run_gen=="FCS"~ "Fall"))%>%
  drop_na(run_gen)%>%
  mutate(gp_match=case_when(run_pheno==run_gen ~"match",
                            TRUE~"no_match"))

# CWT test ----------------------------------------------------------------
battle_cwt <-battle_sr_profile %>%
  distinct(lab_id, year, age, assignment, assignment, adipose_fin, tag_code=CWT_Code)%>%
  # mutate(tag_code=as.numeric(tag_code))%>%
  drop_na(tag_code)%>%
  left_join(rmis, by="tag_code")%>%
  filter(!tag_code=="NTD")%>%
  filter(!tag_code=="LOST TAG")%>%
  mutate(cwt_age=year-brood_year)%>%
  mutate(age_match=case_when(age==cwt_age ~"match",
                                TRUE~"no_match"))%>%
  mutate(cwt_origin=case_when(hatchery_location_name=="FEATHER R HATCHERY" ~"FEH",
                              TRUE~ hatchery_location_name
                              ))%>%
  mutate(origin_match=case_when(assignment==cwt_origin ~"match",
                            TRUE~"no_match"))

### CWT and otolith age comparison
## Age precision
ap.CWT <- agePrecision(~age+cwt_age,data=battle_cwt)
summary(ap.CWT,what="difference")
summary(ap.CWT,what="precision")

ab.CWT <- ageBias(age~cwt_age,data=battle_cwt)
summary(ab.CWT,what="symmetry")
plot(ab.CWT)

# USFWS carcass survey --------------------------------------------------------------------
carcass_summary <- carcass_snorkel_clipped %>% 
  group_by(SPAWN.YEAR) %>% 
  summarise(total_fish =n(),
            clipped_fish = sum(ADIPOSE == "Absent"),
            prop_clipped = clipped_fish/total_fish)

write.csv(carcass_summary, "outputs/snorkelcarcass_summary.csv", row.names=F)

carcass_clipped <- carcass_snorkel_clipped %>% 
filter(ADIPOSE =="Absent") %>%   #ADIPOSE %in% c("Absent",  "Unknown", "UNKNOWN")
  dplyr::select(ID.COMBINED,DATE,LOCATION,ADIPOSE,CWT,CWT_CODE)

TagCode <- carcass_tagcode %>% 
            rename(CWT_CODE = Tag.Code..Id, ID.COMBINED = Biosample.ID) %>% 
            mutate(CWT_CODE = as.character(CWT_CODE))

carcass_clipped_wtcode <- carcass_clipped %>% 
                          filter(CWT =="Y") %>% 
                          left_join(TagCode,by=c("ID.COMBINED","CWT_CODE")) %>% 
                          mutate(parsed_date = dmy(DATE),
                                 year= year(parsed_date))

write.csv(carcass_clipped_wtcode, "outputs/carcass_clippedCWTcode.csv", row.names=F)

carcass_wtcode_summary <-carcass_clipped_wtcode %>% 
  group_by(year) %>% 
  summarise(total_fish =n(),
            clipped_fish = sum(ADIPOSE == "Absent"),
            prop_clipped = clipped_fish/total_fish)



otoCNH_summary <- battle_sr_natal_including_adclip %>% 
              group_by(survey_year) %>% 
            summarise(total_fish =n(),
                      CNH_fish = sum(assignment == "CNH"),
            prop_CNH = CNH_fish/total_fish)

carcassoto_summary <- merge(carcass_wtcode_summary ,otoCNH_summary, by.x='year', by.y ='survey_year')

write.csv(carcassoto_summary, "outputs/carcassotoCNH_summary.csv", row.names=F)


# Figure of CNH proportions
CNH_data_long <- carcassoto_summary %>%
  dplyr::select(year, clipped_fish, CNH_fish) %>% #, prop_clipped, prop_CNH
  pivot_longer(cols = c( clipped_fish, CNH_fish),  #prop_clipped, prop_CNH
               names_to = "Type", 
               values_to = "Number")

(carcassoto_barplot <- ggplot(CNH_data_long, 
                              aes(x = factor(year), y = Number, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(x = "Spawn Year", y = "Number")+
    scale_fill_manual(values = c("clipped_fish" = "#1f77b4", "CNH_fish" = "#ff7f0e"),
                      labels = c("Clipped Feather Fish", "CNH_assigned Fish")) +
    theme_minimal()
) 
ggsave('outputs/p_carcassoto_CNH.jpg',plot=carcassoto_barplot, dpi = 350, height = 4, width = 6)

## CNH assignment through time
otoCNH_date <- battle_sr_natal_including_adclip %>% 
                filter(assignment == "CNH") %>% 
                mutate(month=month(Date)) %>% 
                count(month) 

(p_otoCNH_date <-  ggplot(otoCNH_date, aes(x = month, y = n)) +
    geom_col(fill = "steelblue") +
    labs(title = "Number of CNH origin Fish Sampled per Month",
         x = "Month",
         y = "Number of Fish") +
    geom_text(aes(label = n), vjust = -0.5, size = 3) +
    theme_minimal()
)
ggsave('outputs/p_otoCNH_date.jpg',plot=p_otoCNH_date, dpi = 350, height = 4, width = 6)


# Outmigration profiles ---------------------------------------------------
#Show which fish to plot
fish_ids <- unique(battle_sr_profile$lab_id)
fish_ids #Show list of fish that will be plotted
#Create pdf
pdf("outputs/sr_iso_trimmed_battle.pdf",16, 8)
for(i in fish_ids){
  plot_chem<- battle_sr_profile%>% filter (lab_id==i)
  p <- ggplot(data=plot_chem)+
    geom_ribbon(aes(x=Distance_um, ymin = sr_spline-(2*sr_spline_se), ymax =sr_spline+(2*sr_spline_se)),
                fill = "black", alpha=0.2)+
    geom_line(aes(x=Distance_um, y=sr_spline), color="black")+
    geom_hline(aes(yintercept=0.70918), color="steelblue")+ #Ocean
    geom_vline(aes(xintercept=core_region), color="orange", linetype="dashed")+
    geom_vline(aes(xintercept=natal_start), color="black", linetype="dashed")+
    geom_vline(aes(xintercept=natal_stop), color="black", linetype="dashed")+
    geom_vline(aes(xintercept=natal_exit), color="blue", linetype="dashed")+
    geom_vline(aes(xintercept=fw_exit), color="red", linetype="dashed")+
    geom_text(aes(x=100, y=0.710, label=assignment), check_overlap = T, size=8)+
    scale_x_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 20))+
    scale_y_continuous(expand=c(0,0),expression(""^87*Sr*"/"^86*Sr), limits= c(0.702, 0.711), 
                       breaks = scales::pretty_breaks(n = 8))+
    theme_bw() +
    theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    ggtitle(paste(i, plot_chem$phenotype_mig))
  print(p)
}
dev.off()

# Juvenile outmigration analysis -----------------------------------------------------
### Otolith Sr data
#Lineplots by phenotype for BAT fish
battle_org_pheno <-battle_sr_profile%>%
  filter(assignment=="BAT")%>%
  filter(!adipose_fin=="absent")%>%
  drop_na(phenotype_mig)%>%
  drop_na(age_f)%>%
  ungroup()

(p_sriso_lines <- ggplot(data=battle_org_pheno)+
    geom_line(aes(y=sr_spline, x=Distance_um, group=lab_id, color = phenotype_mig), alpha = .5, linewidth = 0.8)+
    geom_hline(yintercept = 0.70918, alpha = .8, linetype = "dashed")+ #Ocean
    scale_color_manual(name="Phenotype", values=c("#FF0000","#00A08A","#F2AD00"),
                       labels = c("Fry", "Subyearling", "Yearling"))+
    scale_x_continuous(expand=c(0,0), breaks = scales::pretty_breaks(n = 5),
                       limits=c(0, 1100), "Distance (µm) from Core")+
    scale_y_continuous(expand=c(0,0),expression(""^87*Sr*"/"^86*Sr), 
                       limits= c(0.702, 0.711), breaks = scales::pretty_breaks(n = 3))+
    theme_bw(base_size = 18)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          legend.position='none')+
    facet_wrap(~phenotype_mig, ncol=2)
)

# Phenotype hist by outmigration year
battle_org_year <-battle_org_pheno %>%
  mutate(year =as.factor(year)) %>% 
  distinct(lab_id,phenotype_mig, phenotype_mig_rev, year, outmigration_year,natal_exit)

(p_phenotype_year <-ggplot(data=battle_org_year) +
    geom_bar(aes(x=outmigration_year,fill=phenotype_mig_rev), color="black",  position="fill") + 
    scale_x_continuous("Outmigration Year")+ 
    scale_y_continuous("Proportion of Escapement")+
    scale_fill_manual(values=c("#F2AD00","#00A08A","#FF0000"), name="Phenotype")+
    theme_classic()+
    theme(legend.position="top", axis.text = element_text(size = 8))
)
ggsave('outputs/p_phenotype_year.jpg',plot=p_phenotype_year, dpi = 350, height = 4, width = 6)


#Combined plot / Paper Figure Appendix 8 
(phenoytpe_lines_and_bars<- plot_grid(p_sriso_lines, p_phenotype_year , labels = c('A', 'B'),
                                      nrow = 1)
)
ggsave(plot=phenoytpe_lines_and_bars , "outputs/phenoytpe_lines_and_bars .png", width = 10, height = 4)

# Phenotype overview
(p_phenotype_ov <-ggplot(data=battle_org_year) +
    geom_bar(aes(x=phenotype_mig,fill=phenotype_mig), color="black") + 
    geom_text(stat='count',aes(x=phenotype_mig,label=after_stat(count)),vjust=-1)+
    scale_x_discrete("Juvenile outmigration phenotype")+
    scale_y_continuous("Number of fish", limits=c(0,110))+
    scale_fill_manual(values=c("#FF0000","#00A08A","#F2AD00"), name="Phenotype")+
    theme_classic()+
    theme(legend.position="none", axis.text = element_text(size = 8))
)
ggsave('outputs/p_phenotype_ov.jpg',plot=p_phenotype_ov, dpi = 350, height = 3, width = 5)

# Summary table of phenotypes
battle_pheno_summary <- battle_org_year %>% 
  group_by(phenotype_mig) %>% 
  summarize(N=n(),
            Proportion = n() / nrow(battle_org_year))

#Phenotype proportions and hydrology
Battle_pheno_props_yearly <- battle_org_year %>%
                              group_by(outmigration_year, phenotype_mig) %>%
                             summarize(n_fish=n()) %>%
                             ungroup()%>%
                             group_by(outmigration_year)%>%
                            mutate(prop_pheno = n_fish/(sum(n_fish))) %>% 
                            left_join(SacIndexWY, by=c('outmigration_year'='WY')) 

smolt = filter(Battle_pheno_props_yearly , phenotype_mig=="Subyearling") %>%
  rename(percent_smolt =prop_pheno) %>%
  filter(n_fish>2)

fry = filter(Battle_pheno_props_yearly , phenotype_mig=="Fry") %>%
  rename(percent_fry =prop_pheno) 

yearling = filter(Battle_pheno_props_yearly , phenotype_mig=="Yearling") %>%
  rename(percent_yearl =prop_pheno) 

(p_flow <- ggplot(smolt, aes(x=sac_index, y=percent_smolt)) +
  geom_point(aes(size=n_fish))+
  geom_smooth(method='lm', se=F, color = 'black') +
  #stat_poly_eq( aes(label = paste(..rr.label..)), formula = y ~ x, parse = TRUE )+
  # geom_smooth(method = "glm",
  #             method.args = list(family = "binomial"),
  #             aes(weight = n_fish),
  #             se = FALSE,
  #             color="black",
  #             linewidth=1) +
 # xlim(4.9,15.2)+
#  ylim(0.1,0.7)+
  theme_classic(base_size=12)+
  geom_text_repel(aes(label=as.factor(outmigration_year)))+
  labs(x = "Sacramento Water Index", y = "Proportion of smolt migrants", size = "N samples")
)

### RST data from USFWS (Schraml et al. in prep)
RST_USFWSdata <-read.csv(here('data','UBCRST_BY2008-2022.csv'),  stringsAsFactors = F) |>
  filter(FWSRace == "S" & ForkLength >0) |> 
  mutate(date = as.Date(SampleDate, format = "%m/%d/%Y"),
         true_fl = ForkLength,
         week = week(date),
         week = if_else(week == 53, 52, week), # all December 30th or 31st
         month = month(date),
         day = ifelse(day(date) < 10, paste0("0", day(date)), day(date)),
         Outyear = BroodYear +1,
         fake_year = ifelse(month %in% 10:12, 1971, 1972),
         # create fake date for plotting all years on same x axis
         fake_date = ymd(paste0(fake_year, "-",month, "-", day))) 

### Download yearling cutoff based on FlowWest ruleset:
# https://srjpe.github.io/SRJPEdata/articles/yearling_ruleset.html#apply-cutoff-to-catch-data
Yearling_cutoff <- SRJPEdata::daily_yearling_ruleset %>% 
  filter(stream=="battle creek") %>% 
  mutate(fake_year = ifelse(month %in% 10:12, 1971, 1972),
         # create fake date for plotting all years on same x axis
         fake_date = ymd(paste0(fake_year, "-",month, "-", day))) %>% 
  dplyr::select(cutoff, fake_date) 

RST_USFWSdata <- RST_USFWSdata |> 
  left_join(Yearling_cutoff, by = c("fake_date")) |> 
  mutate(is_yearling = case_when(true_fl <= 45  ~ 'Fry',
                                 # true_fl > cutoff ~ "Yearling",
                                 true_fl > 110 ~ "Yearling",
                                 TRUE ~ 'Subyearling'))

### plot Battle RST
(pRST_Battle <- ggplot() + 
    geom_point(RST_USFWSdata, mapping = aes(x = fake_date, y = true_fl, colour = is_yearling))+
    scale_x_date(date_breaks = "1 month", date_labels = "%b",
                 limits = c(ymd("1971-10-01"), ymd("1972-09-30"))) +
    scale_y_continuous(limits=c(20,150),breaks= seq(20,150,20)) +
    labs(y = "fork length (mm)",x = "",color="Phenotype") +
    #ggtitle("Battle Creek") +
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle=45, hjust = 1), 
          text = element_text(size=20),
          legend.position = "none")+
    scale_color_manual(values = wes_palette("Darjeeling1", n = 3)) 
)

ggsave(plot=pRST_Battle , "outputs/RST_USFWS_Battle.png", width = 7, height = 5)

### Plot RST vs otolith juvenile phenotypes
rst_data_c <- RST_USFWSdata %>%
  #  filter(stream=="battle creek")%>%
  mutate(id=row_number())%>%
  mutate(phenotype=is_yearling,
         type="RST")%>%
  dplyr::select(id, phenotype, type)

oto_data_c <- battle_org_year %>%
  mutate(id=lab_id,
         phenotype=phenotype_mig,
         type="Otolith")%>%
  dplyr::select(id, phenotype, type)

rstoto=rbind(rst_data_c,oto_data_c)%>%
  mutate(type=factor(type, levels=c("RST","Otolith")))%>%
  mutate(phenotype=factor(phenotype, levels=c("Yearling","Subyearling","Fry")))

(p_rstoto <-ggplot(data=rstoto) +
    geom_bar(aes(x=type,fill=phenotype), color="black",  position="fill", width=0.8) + 
    scale_x_discrete("Type")+
    scale_y_continuous("Proportion of Phenotype")+
    scale_fill_manual(values=c("#F2AD00","#00A08A","#FF0000"), name="Phenotype")+
    theme_classic()+
    theme(legend.position="right", axis.text = element_text(size = 8))
)
ggsave('outputs/p_rstoto.jpg',plot=p_rstoto, dpi = 350, height = 4, width = 5)


## Histogram and density figures
(pRSTsizehist_pery <- ggplot(data=RST_USFWSdata, aes(x=true_fl)) +
   geom_histogram(binwidth = 10)+
   scale_x_continuous(limits=c(20,120), breaks=seq(20,120,40))+
   labs(x="fish size (mm)",y="Count")+
   ggtitle('Natal Exit - RST')+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"),
         plot.title = element_text(face = "bold"))+
   facet_grid(~Outyear)
)

ggsave('outputs/p_RSTsizehist_pery.jpg',plot=pRSTsizehist_pery, dpi = 350, height = 5, width = 13)

(pRSTsizehist_selecty <- ggplot(data=RST_USFWSdata %>% filter(Outyear %in% c(2009,2011:2013,2016,2018)), aes(x=true_fl)) +
    geom_histogram(binwidth = 10)+
    scale_x_continuous(limits=c(20,120), breaks=seq(20,120,20))+
    labs(x="fish size (mm)",y="Count")+
    ggtitle('Natal Exit - RST')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(face = "bold"))+
    facet_grid(~Outyear)
)
ggsave('outputs/p_RSTsizehist_selecty.jpg',plot=pRSTsizehist_selecty, dpi = 350, height = 5, width =9)

(pRSTsizehistdens <- ggplot(data=RST_USFWSdata, aes(x=true_fl)) +
    geom_histogram(aes(y = ..density..),colour = 1) +
    geom_density(adjust=2, size=0.8,colour="blue")+
    labs(x="fish size (mm)",y="Density")+
    ggtitle('Natal Exit - RST')+
    scale_x_continuous(limits=c(20,120), breaks=seq(20,120,20))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(face = "bold"))
)

ggsave('outputs/p_RSTsizehistdens.jpg',plot=pRSTsizehistdens, dpi = 350, height = 4, width = 5)

# Density plot colored by LH
(pRSTsizedens <- ggplot()+
    geom_density(data=RST_USFWSdata,aes(x=true_fl, fill = is_yearling,color=is_yearling),
                 alpha=0.3,adjust=1.3, size=1) +
    labs(x="fish size (mm)",y="Density")+
    ggtitle('Natal Exit - RST')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),
          #      axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
          legend.title = element_blank(),
          #     legend.position='none',
          plot.title = element_text(face = "bold"))+
    scale_x_continuous(limits=c(25,155), breaks=seq(25,155,25))+
    scale_color_manual( values=c("#FF0000","#00A08A","#F2AD00"))+
    scale_fill_manual(values=c("#FF0000","#00A08A","#F2AD00"))
)


## Calculation of FL at natal exit based on OR 
# Broken stick FL calibration model
calib =   read.csv( "data/or_fl_calibration_input_data_updated2025_with_predicted.csv") 
calib = subset(calib, select=c("Sample_ID","OR","FL")) %>% 
  drop_na(FL)
FL <- calib$FL
OR <- calib$OR
forced.intercept <- 30
res.lm <- lm(I(FL-forced.intercept) ~ 0 + OR)
res.bs <- segmented::segmented(res.lm, seg.Z = ~ 0 + OR)
FL.bs <- predict(res.bs) + forced.intercept
calib$Predict <- FL.bs 
calib$predicted_FL = case_when(OR<256 ~ (OR*0.021259)+31.08748, TRUE ~ (OR*0.165660)-5.8799)

# New back calculation method from Anna 
back_predict_OR <- function(pred_FL) {
  # invert both equations
  OR1 <- (pred_FL - 31.08748) / 0.021259   # candidate for OR < 256
  OR2 <- (pred_FL + 5.8799) / 0.165660     # candidate for OR >= 256
  
  # choose appropriate region
  if (OR1 < 256 & OR1 >= 0) {
    return(OR1)
  } else {
    return(OR2)
  }
}

OR_at_FL_45 <- back_predict_OR(45)
print(OR_at_FL_45)

OR_at_FL_110 <- back_predict_OR(110)
print(OR_at_FL_110)

# Natal exit FL estimate for Battle fish
battle_org_year <- battle_org_year  %>% 
  mutate(natal_exit_FL = case_when(natal_exit <256 ~ (natal_exit *0.021259)+31.08748, 
                                   TRUE ~ (natal_exit *0.165660)-5.8799) 
  )

## Histogram and density figures
(potosizehist <- ggplot(data=battle_org_year, aes(x=natal_exit_FL)) +
    geom_histogram(bins = 20, colour=1)+
    scale_x_continuous(limits=c(20,120), breaks=seq(20,120,10))+
    ggtitle('Natal Exit - Otoliths')+
    labs(x="fish size (mm)",y="Count")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(face = "bold"))
)

(potosizehist_pery <- ggplot(data=battle_org_year, aes(x=natal_exit_FL)) +
    geom_histogram(bins = 20, colour=1)+
    scale_x_continuous(limits=c(20,120), breaks=seq(20,120,40))+
    ggtitle('Natal Exit - Otoliths')+
    labs(x="fish size (mm)",y="Count")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(face = "bold"))+
    facet_grid(~outmigration_year)
)
ggsave('outputs/p_Otosizehist_pery.jpg',plot=potosizehist_pery, dpi = 350, height = 5, width = 13)

(potosizehist_selecty <- ggplot(data=battle_org_year %>% filter(outmigration_year %in% c(2009,2011:2013,2016,2018)), 
                                aes(x=natal_exit_FL)) +
    geom_histogram(bins = 20, colour=1)+
    scale_x_continuous(limits=c(20,120), breaks=seq(20,120,20))+
    ggtitle('Natal Exit - Otoliths')+
    labs(x="fish size (mm)",y="Count")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(face = "bold"))+
    facet_grid(~outmigration_year)
)

ggsave('outputs/p_Otosizehist_selecty.jpg',plot=potosizehist_selecty, dpi = 350, height = 5, width = 9)

(otosizehistdens <- ggplot(data=battle_org_year, aes(x=natal_exit_FL)) +
    geom_histogram(aes(y = ..density..),colour = 1) +
    geom_density(adjust=2,size=0.8,colour="blue")+
    labs(x="fish size (mm)",y="Density")+
    ggtitle('Natal Exit - Otoliths')+
    scale_x_continuous(limits=c(30,130), breaks=seq(30,130,20))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(face = "bold"))
)

(potosizedens <- ggplot(data=battle_org_year)+
    geom_density(aes(x=natal_exit_FL), 
                 #fill = phenotype_mig,color=phenotype_mig),
                 alpha=0.3,adjust=0.8, size=1) +
    labs(x="fish size (mm)",y="Density")+
    ggtitle('Natal Exit - Otoliths')+
    scale_x_continuous(limits=c(10,140), breaks=seq(10,140,20))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(face = "bold"))
)

ggsave('outputs/p_Otosizedens.jpg',plot=potosizedens, dpi = 350, height = 4, width = 5)

### Combined RST and oto data / Paper figure 9
(pRSTOtosizedens <- ggplot(data=battle_org_year)+
    geom_density(aes(x=natal_exit_FL,y=..scaled..),adjust=0.7, size=1, color='blue') +
    geom_density(data=RST_USFWSdata, aes(x=true_fl,y=..scaled..),alpha=0.5,adjust=1.5, size=1, color='red')+
    labs(x="fish size (mm)",y="Density")+
    ggtitle('Natal Exit ')+
    scale_x_continuous(limits=c(10,140), breaks=seq(10,140,20))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(face = "bold"))
  #    scale_color_manual(values = c("mediumvioletred","grey48","goldenrod2"))+
  #   scale_fill_manual(values = c("mediumvioletred","grey48","goldenrod2"))
)
ggsave('outputs/p_RST+Otosizedens.jpg',plot=pRSTOtosizedens, dpi = 350, height = 4, width = 5)

(p_rstoto <-ggplot(data=rstoto) +
   geom_bar(aes(x=type,fill=phenotype), color="black",  position="fill", width=0.6) + 
   scale_x_discrete("")+
   scale_y_continuous("")+
   scale_fill_manual(values=c("#F2AD00","#00A08A","#FF0000"), name="Phenotype")+
   theme_classic()+
   theme(legend.position="right", axis.text = element_text(size = 8))
)

(pRSTsizehist <- ggplot(data=RST_USFWSdata, aes(x=true_fl,fill = is_yearling,color=is_yearling)) +
    geom_histogram(binwidth = 10)+
    scale_x_continuous(limits=c(20,160), breaks=seq(20,160,20))+
    scale_y_continuous(expand=c(0,0))+
    labs(x="Fish size (mm)",y=" Fish count")+
    ggtitle('Natal Exit - RST')+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(face = "bold"), legend.position="none")+
    scale_color_manual( values=c("#FF0000","#00A08A","#F2AD00"))+
    scale_fill_manual(values=c("#FF0000","#00A08A","#F2AD00"))
)

(potosizehist <- ggplot(data=battle_org_year, aes(x=natal_exit_FL,fill =phenotype_mig,color=phenotype_mig)) +
    geom_histogram(binwidth = 10)+
    scale_x_continuous(limits=c(20,160), breaks=seq(20,160,20))+
    scale_y_continuous(expand=c(0,0))+
    ggtitle('Natal Exit - Otoliths')+
    labs(x="Fish size (mm)",y="Fish count")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(face = "bold"), legend.position="none") +
    scale_color_manual( values=c("#FF0000","#00A08A","#F2AD00"))+
    scale_fill_manual(values=c("#FF0000","#00A08A","#F2AD00"))
)


(pRSTOtoLHSize <- ggarrange(p_rstoto,pRSTsizehist,potosizehist,
                            # align="v",
                            # common.legend = TRUE,
                            nrow=3,
                            labels=c('a','b','c'))
)

ggsave('outputs/pRSTOtoLHSize.jpg',plot=pRSTOtoLHSize, dpi = 350, height =9, width =6)

# Stock-recruit analysis ---------------------------------------------------
### Create stock recruit table for paper
prop_boot_group <- prop_boot %>% 
  group_by(year, assignment) %>% 
  summarise(prop_as = sum(Mean_Proportion),
            prop_as_lowCI = sum(Lower_95_CI),
            prop_as_upperCI = sum(Upper_95_CI))

JuveEscap_data <-read.csv(here('data','RSTEscap_data.csv'),  stringsAsFactors = F) %>% 
  mutate(numyear=as.numeric(Year),
         Year=as.factor(Year)) %>%  
  left_join(SacIndexWY,by=c("Outyear"="WY")) %>% 
  left_join(age_prop_extended,by=c("Year"="survey_year")) %>% 
  left_join(prop_boot_group, by = c("Year"="year")) %>% 
  filter(assignment == "BAT"| Year %in% c(2002, 2003, 2004)) %>% 
  mutate(escap_year = as.numeric(as.character(Year)),
         Nadult_bat = Nadult*prop_as, 
         Nadult_bat_A2 = lead(Nadult * age_p_2 *prop_as, 2),
         Nadult_bat_A3 = lead(Nadult * age_p_3 *prop_as, 3),
         Nadult_bat_A4 = lead(Nadult * age_p_4 *prop_as, 4),
         cohort_lead = lead(Nadult * age_p_2 *prop_as, 2) + 
           lead(Nadult * age_p_3*prop_as, 3) + 
           lead(Nadult * age_p_4*prop_as, 4),
         cohort_lead_upperCI = lead(Nadult * age_p_2 *prop_as_upperCI, 2) + 
           lead(Nadult * age_p_3*prop_as_upperCI, 3) + 
           lead(Nadult * age_p_4*prop_as_upperCI, 4),
         cohort_lead_lowCI = lead(Nadult * age_p_2 *prop_as_lowCI, 2) + 
           lead(Nadult * age_p_3*prop_as_lowCI, 3) + 
           lead(Nadult * age_p_4*prop_as_lowCI, 4),
          RS_cohort_lead = cohort_lead/Nadult,
         RS_cohort_lead_upperCI = cohort_lead_upperCI/Nadult,
         RS_cohort_lead_lowCI = cohort_lead_lowCI/Nadult,
         se = (RS_cohort_lead_upperCI  - RS_cohort_lead_lowCI) / (2 * 1.96),
         weight = 1 / se^2,
         juveadult_ratio = Njuve/Nadult) 

write.csv(JuveEscap_data, "outputs/JuveEscap_data.csv", row.names=F)

# keep years with N samples >5
JuveEscap_data_over5 <- JuveEscap_data %>% 
  filter(Year %in% c(2006:2014)) 

write.csv(JuveEscap_data_over5, "outputs/JuveEscap_data_over5.csv", row.names=F)

### Make Nate extra figures for paper
JuveEscap_data_stock_recruit <- JuveEscap_data %>%
  filter(!is.na(Nadult), !is.na(Njuve), Nadult > 0, Njuve > 0)

(p_juvenilevsescape_BH <- ggplot(data=JuveEscap_data,aes(x=Nadult,y=Njuve))+
    geom_point(size=2)+
    stat_smooth(
      method = "nls",
      formula = y ~ (alpha * x) / (1 + beta * x),
      method.args = list(
        start = list(
          alpha = max(JuveEscap_data$Njuve / JuveEscap_data$Nadult, na.rm = TRUE),
          beta = 1 / max(JuveEscap_data$Nadult, na.rm = TRUE)
        )
      ),
      se = FALSE,
      color = "black" ) +
    theme_classic()+
    scale_y_continuous(name="Njuveniles")+
    scale_x_continuous(name="Nadults")+
    geom_text(data=JuveEscap_data,aes(label = Year),
              vjust=-0.5) 
    
)

(p_juvenilevsescape_Ricker <- ggplot(
  data = JuveEscap_data_stock_recruit,
  aes(x = Nadult, y = Njuve)
) +
    geom_point(size = 2) +
    stat_smooth(
      method = "nls",
      formula = y ~ alpha * x * exp(-beta * x),
      method.args = list(
        start = list(
          alpha = max(JuveEscap_data_stock_recruit$Njuve / JuveEscap_data_stock_recruit$Nadult, na.rm = TRUE),
          beta = 1 / max(JuveEscap_data_stock_recruit$Nadult, na.rm = TRUE)
        )
      ),
      se = FALSE,
      color = "black"
    ) +
    scale_y_continuous(name = "Njuveniles") +
    scale_x_continuous(name = "Nadults") +
    geom_text(
      aes(label = Year),
      vjust = -0.5
    ) +
    theme_classic()
)

ggsave('outputs/p_juvenilevsescape_BH.jpg',plot=p_juvenilevsescape_BH, dpi = 350, height = 5, width = 6)
ggsave('outputs/p_juvenilevsescape_Ricker.jpg',plot=p_juvenilevsescape_Ricker, dpi = 350, height = 5, width = 6)

(p_juveperspawn_WY <- ggplot(data=JuveEscap_data,
                             aes(x=sac_index,y=juveadult_ratio))+
    geom_point()+
    stat_poly_line(
      formula = y ~ poly(x, 2, raw = TRUE),
      se = FALSE,
      color = "black") +
    theme_classic()+
    scale_x_continuous(name="Water year index",
                       limits = c(3, 15),
                       breaks = 3:15,
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(name="Juvenile per spawner",
                       limits = c(0, 600))+ 
    geom_text(data=JuveEscap_data,aes(x=sac_index,y=juveadult_ratio,label = Year),
              vjust=-0.5) +
    theme_classic()
)

ggsave('outputs/p_juveperspawn_WY_v2.jpg',plot=p_juveperspawn_WY, dpi = 350, height = 5, width = 6)


(p_recruitperjuve <- ggplot(data=JuveEscap_data_over5,aes(y=cohort_lead,x=Njuve))+
    geom_point(size=2)+
    theme_classic()+
    scale_y_continuous(name="BAT-origin 2, 3 and 4 year adult returns",
                                          limits = c(40, 300),
                                          breaks = seq(40,300,20),
                                          expand = expansion(mult = c(0, 0)))+
    scale_x_continuous(name="Njuveniles",
                       limits = c(10000, 110000),
                       breaks = seq(10000, 110000,25000))+
    geom_text(data=JuveEscap_data_over5,aes(label = Year),
              vjust=-0.5) 
)

ggsave('outputs/p_recruitperjuve.jpg',plot=p_recruitperjuve, dpi = 350, height = 5, width = 6)

### Paper figure 8
(p_recperspawn_cohort_WY <- ggplot(data=JuveEscap_data_over5,
                                   aes(x=sac_index,y=RS_cohort_lead))+
    geom_point()+
    geom_errorbar(aes(
      ymin = RS_cohort_lead_lowCI,
      ymax = RS_cohort_lead_upperCI ),
      width = 0.15) +
    geom_smooth(aes(weight = weight),method = "lm",se = TRUE,color="black")+
    theme_classic()+
    scale_x_continuous(
      name = "Water year index",
      limits = c(4.5, 11),
      breaks = 4:11,
      expand = expansion(mult = c(0, 0)))+
    scale_y_continuous(limits = c(0, 3),
                       #breaks = seq(0,3,0.5),
                       expand = expansion(mult = c(0, 0)),
                       name="BAT-origin recruits per spawner 2, 3 and 4 years prior")+
    geom_text(data=JuveEscap_data_over5,aes(x=sac_index,y=RS_cohort_lead,label = Year),
              vjust=-0.5) +
    theme_classic()
)

ggsave('outputs/p_recruitperspawnervsWY_withCI.jpg',plot=p_recperspawn_cohort_WY, dpi = 350, height = 5, width = 6)

(p_recperspawn_spawn <- ggplot(data=JuveEscap_data,aes(x=Nadult,y=RS_cohort_lead))+
    geom_point()+
    theme_classic()+
    stat_poly_line(color="black") +
    stat_poly_eq() +
    scale_x_continuous(name="Spawners",limits = c(50, 300)) +
    scale_y_continuous(name="Recruits per spawner")+
    #geom_smooth(method='lm', formula= y~x)+
    geom_text(data=JuveEscap_data,aes(x=Nadult,y=RS_cohort_lead,label = Year),
              vjust=-0.5) +
    theme_classic()
)

(p_StockRecruit <- ggplot(data=JuveEscap_data,mapping=aes(x=Nadult,y=cohort_lead)) +
    geom_point(aes(color=sac_index),size=2.5) +
    geom_smooth(method = "gam", formula = y ~ s(x,k=3), se = TRUE,color="black")+
    #   geom_text(aes(x=Nspawn,y=Nrecruit_BAT,label = Year),vjust=-0.5) +
    theme_classic()+
    scale_x_continuous(name="Escapement",limits = c(50, 300)) +
    scale_y_continuous(name="Recruits")
)

# Spawner origin analysis -----------------------------------------------
# Plot total Battle escapement time series
Escap_data <-read.csv(here('data','BattleEscape.csv'),  stringsAsFactors = F)

(p_Battleescape <- ggplot(data = Escap_data %>% filter(Year > 1969), aes(x = Year, y = Count)) +
    geom_point(color = "darkorange", size = 2.5) +
    geom_line(color = "darkorange", linewidth = 1) +
    ggtitle('Battle Creek') +
    ylab("")+
    theme_bw()
)

ggsave('outputs/p_historicescape.jpg',plot=p_Battleescape ,dpi = 350, height = 4, width = 6)


battle_escape <- battle_sr_natal_including_adclip %>%
  mutate(year=factor(survey_year))%>%
  group_by(year, assignment, .drop=F)%>%
  summarise(n_year_as=n(), .groups = "drop_last") %>%
  group_by(year)%>%
  mutate(n_year=sum(n_year_as))%>%
  ungroup()%>%
  left_join(
    prop_boot_group %>%
      rename(
        prop_origin = prop_as,
        prop_origin_lowCI = prop_as_lowCI,
        prop_origin_upperCI = prop_as_upperCI
      ),
    by = c("year","assignment")
  ) %>%
  left_join(
    JuveEscap_data %>%
      distinct(Year, Nadult, numyear),
    by = c("year"="Year")
  ) %>%
  mutate(fish_expand=Nadult*prop_origin,
         fish_expand_lowCI=Nadult*prop_origin_lowCI,
         fish_expand_upperCI=Nadult*prop_origin_upperCI,
         numyear = as.numeric(as.character(year))) %>%
  ungroup()

write.csv(battle_escape, "outputs/battle_escape_perorigin.csv", row.names=F)

battle_escape_filtered <- battle_escape %>% 
                          filter(assignment == "BAT") 

battle_escape_over5 <- battle_escape %>% 
                       filter(n_year > 5) 
                                
# Bar plot Battle origin only and years with at least >5 otoliths analyzed
(p_natalbattle_exp <- ggplot(data=battle_escape_filtered %>% filter(year %in% c(2006:2014)))+
    geom_bar(aes(x = year, y = fish_expand), stat='identity', color="black")+
    scale_fill_viridis(discrete = T, name="Natal origin", direction=-1)+
    theme_classic()+
    theme (panel.background = element_rect(colour = "black"), legend.position="right")+
    ylab("Number of Adults")
)

# Bar plot per year
(p_nataladults_exp <- ggplot(data=battle_escape_over5)+
    geom_bar(aes(x = assignment, y = fish_expand, fill=assignment), stat='identity', color="black")+
    scale_fill_manual(values = colors, name="Natal origin")+
    theme_classic()+
    theme (panel.background = element_rect(colour = "black"), legend.position="right",
           axis.text.x=element_blank(),
           axis.ticks.x=element_blank())+
    xlab("Predicted natal origin")+
    ylab("Number of Adults")+
    facet_wrap(~year,scales='free_y')
  #facet_grid(~survey_year,scales='free_y') #Other per year option
)

ggsave('outputs/p_adultexp_peryear.jpg',plot=p_nataladults_exp, dpi = 350, height = 4, width = 8)

# Area plot / Paper figure 7
origin_levels <- c("OTH", "BUT", "MIL", "CNH", "FEH", "SAC", "DEE", "BAT")
origin_colors <- colors[origin_levels]

total_by_year <- JuveEscap_data %>%
  distinct(numyear, Nadult) %>%
  filter(!is.na(numyear), !is.na(Nadult)) %>%
  arrange(numyear)

data_all <- battle_escape %>%
  filter(n_year > 5) %>% 
  filter(!is.na(numyear), !is.na(fish_expand)) %>%
  mutate(assignment = as.character(assignment)) %>%
  filter(assignment %in% origin_levels) %>%
  mutate(
    assignment = factor(assignment, levels = origin_levels),
    assignment_stack = factor(assignment, levels = rev(origin_levels))
  )

battle_origin_area <- data_all %>%
  filter(assignment == "BAT")

(p_nataladults_area <-ggplot()+
    geom_area(data=total_by_year, aes(x=numyear,y=Nadult),
              stat = "identity", position = "identity",
              linewidth=0.5, fill="grey88", colour = NA)+
    geom_area(data=data_all,
              aes(x = numyear, y = fish_expand,
                  group = assignment_stack, colour = assignment, fill = assignment),
              position = "stack",
              linewidth = 1, alpha = 0.25) +
    scale_fill_manual(values=origin_colors, breaks = origin_levels, name="Natal origin", drop = FALSE)+
    scale_colour_manual(values=origin_colors, breaks = origin_levels, name="Natal origin", drop = FALSE)+
    guides(fill = guide_legend(nrow = 2, byrow = TRUE),
           colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_classic()+  
    theme(legend.position = "bottom")+
    scale_x_continuous(limits=c(2005,2021), breaks=seq(2005,2021,1), "Survey year")+
    scale_y_continuous(limits = c(0, 850), expand = expansion(mult = c(0, 0)))+
    ylab("Number of Adult spawners")
)

(p_natalBAT_area <-ggplot()+
    geom_area(data=total_by_year,
              aes(x=numyear,y=Nadult),
              stat = "identity", position = "identity",
              linewidth=0.5, fill="grey88", colour = NA)+
    geom_ribbon(data=battle_origin_area,
                aes(x = numyear, ymin = fish_expand_lowCI, ymax = fish_expand_upperCI),
                fill = origin_colors["BAT"], alpha = 0.25) +
    geom_line(data=battle_origin_area,
              aes(x = numyear, y = fish_expand, colour = assignment),
              linewidth=1) +
    geom_line(data=battle_origin_area,
              aes(x = numyear, y = fish_expand_lowCI, colour = assignment),
              linetype = "dashed", linewidth=1) +
    geom_line(data=battle_origin_area,
              aes(x = numyear, y = fish_expand_upperCI, colour = assignment),
              linetype = "dashed", linewidth=1) +
    scale_colour_manual(values=origin_colors, breaks = origin_levels, name="Natal origin", drop = FALSE)+
    theme_classic()+
    theme(legend.position = "none")+
    scale_x_continuous(limits=c(2005,2021), breaks=seq(2005,2021,1), "Survey year")+
    scale_y_continuous(limits = c(0, 850), expand = expansion(mult = c(0, 0)))+
    ylab("Number of Adult spawners")
)


(p_natalarea <- ggarrange(p_nataladults_area,p_natalBAT_area,
                          nrow=1,ncol=2,
                          legend="bottom",common.legend = TRUE,
                          font.label = list(size = 18),
                          labels=c('a','b'))
)

ggsave('outputs/p_natalarea_withFEH.jpg',plot=p_natalarea, dpi = 350, height = 5, width = 13)

# Predict lambda and escapement numbers for future years ------------------------------------------------------------------
#### Using Ricker stock-recruit model
# Fit Ricker stock-recruit model with Water Year Index (WYindex) as covariate
ricker_wy_model <- nls(
  Nrecruit ~ Nadult * exp(r * (1 - Nadult / K) + beta * WYindex),
  data = data,
  start = list(r = 0.1, K = max(data$Nadult), beta = 0)
)

# Model summary
summary(ricker_wy_model)

# Bootstrap to estimate confidence intervals
ricker_boot <- function(data, indices) {
  d <- data[indices, ]
  fit <- tryCatch({
    nls(Nrecruit ~ Nadult * exp(r * (1 - Nadult / K) + beta * WYindex),
        data = d, start = list(r = 0.1, K = max(d$Nadult), beta = 0))
  }, error = function(e) return(rep(NA, 3)))
  
  if (any(is.na(fit))) return(c(NA, NA, NA))
  coef(fit)
}

set.seed(123)
boot_results <- boot(data = data, statistic = ricker_boot, R = 1000)

# Confidence intervals
r_ci <- boot.ci(boot_results, type = "perc", index = 1)
K_ci <- boot.ci(boot_results, type = "perc", index = 2)
beta_ci <- boot.ci(boot_results, type = "perc", index = 3)

print(r_ci)
print(K_ci)
print(beta_ci)

# Future prediction function
predict_future <- function(r, K, beta, WYmean, start_Nadult, years = 20) {
  result <- data.frame(Year = 1:years, Nadult = NA, Predicted_Recruits = NA, Lambda = NA)
  result$Nadult[1] <- start_Nadult
  
  for (i in 1:years) {
    S <- result$Nadult[i]
    W <- WYmean
    R <- S * exp(r * (1 - S / K) + beta * W)
    result$Predicted_Recruits[i] <- R
    if (i < years) {
      result$Nadult[i + 1] <- R
    }
  }
  result$Lambda <- c(result$Nadult[-1] / result$Nadult[-years], NA)
  return(result)
}

# Define scenarios
scenarios <- data.frame(
  Scenario = c("Average", "Wet", "Dry"),
  WYmean = c(mean(data$WYindex), max(data$WYindex) + 1, min(data$WYindex) - 1)
)

# Predict for each scenario
params <- coef(ricker_wy_model)

all_scenarios <- lapply(1:nrow(scenarios), function(i) {
  predict_future(
    r = params["r"],
    K = params["K"],
    beta = params["beta"],
    WYmean = scenarios$WYmean[i],
    start_Nadult = tail(data$Nadult, 1),
    years = 20
  ) %>% mutate(Scenario = scenarios$Scenario[i])
})

results <- bind_rows(all_scenarios)

# Plot Adult Abundance by Scenario
ggplot(results, aes(x = Year, y = Nadult, color = Scenario)) +
  geom_line(size = 1.2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Predicted Adult Salmon Abundance under Water Year Scenarios",
    y = "Predicted Adults"
  )

# Plot Lambda by Scenario
ggplot(results, aes(x = Year, y = Lambda, color = Scenario)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Predicted Lambda under Water Year Scenarios",
    y = expression(lambda)
  )

# ### Testing stock-recruit curves
# ricker_sr_dat <- JuveEscap_data %>%
#   filter(Year %in% c(2009:2015)) %>%
#   filter(!is.na(Nadult), !is.na(cohort_lead), Nadult > 0, cohort_lead > 0) %>%
#   mutate(retperesc = cohort_lead / Nadult)
# 
# svR <- srStarts(cohort_lead ~ Nadult, data = ricker_sr_dat, type = "Ricker", plot = TRUE)
# #fixed=list(a=3,b=0.008))
# #svBH <- srStarts(return~escapement,data=StockRecruit_data,type="BevertonHolt",plot=TRUE)
# 
# rckr <- srFuns("Ricker")
# srR <- nls(
#   log(cohort_lead) ~ log(rckr(Nadult, a, b)),
#   data = ricker_sr_dat,
#   start = svR,
#   control = nls.control(maxiter = 200, warnOnly = TRUE)
# )
# 
# ricker_boot <- function(dat, i) {
#   boot_dat <- dat[i, ]
#   boot_fit <- try(
#     nls(
#       log(cohort_lead) ~ log(rckr(Nadult, a, b)),
#       data = boot_dat,
#       start = as.list(coef(srR)),
#       control = nls.control(maxiter = 200, warnOnly = TRUE)
#     ),
#     silent = TRUE
#   )
#   
#   if (inherits(boot_fit, "try-error")) {
#     return(c(a = NA_real_, b = NA_real_))
#   }
#   
#   coef(boot_fit)
# }
# 
# set.seed(123)
# bootR <- boot(data = ricker_sr_dat, statistic = ricker_boot, R = 1000)
# boot_coef <- as.data.frame(bootR$t) %>%
#   setNames(names(coef(srR))) %>%
#   drop_na()
# 
# if (nrow(boot_coef) < 20) {
#   stop("Ricker bootstrap produced fewer than 20 converged fits. Try fewer years, different starting values, or plot the fitted curve without bootstrap confidence intervals.")
# }
# 
# cbind(
#   estimates = coef(srR),
#   t(apply(boot_coef, 2, quantile, probs = c(0.025, 0.975)))
# )
# 
# x <- seq(min(ricker_sr_dat$Nadult, na.rm = TRUE), max(ricker_sr_dat$Nadult, na.rm = TRUE), length.out = 199)        # many S for prediction
# pR <- rckr(x, a = coef(srR)["a"], b = coef(srR)["b"])           # predicted mean R
# LCI <- UCI <- numeric(length(x))
# 
# for(i in 1:length(x)) {             # CIs for mean R @ each S
#   tmp <- apply(boot_coef, MARGIN = 1, function(z) rckr(x[i], a = z["a"], b = z["b"]))
#   LCI[i] <- quantile(tmp,0.025)
#   UCI[i] <- quantile(tmp,0.975)
# }
