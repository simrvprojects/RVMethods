# create a set of sample pedigrees for examples

library(SimRVPedigree)
load("C:/Users/cnieuwoudt/Google Drive/sfuthesis-2.2.1/R code/Part 3 Simulation/RV_peds_com.rda")

ped_sum = summary(RV_peds)
length(ped_sum$family_info$FamID[ped_sum$family_info$numAffected > 3])
keep_fams = ped_sum$family_info$FamID[ped_sum$family_info$numAffected > 3]


#cherry pick some informative pedigrees
i = 14
plot(reassign_gen(RV_peds[RV_peds$FamID == keep_fams[i], ]))


keep_fams = keep_fams[c(10, 17, 33, 37, 38)]

RV_peds = RV_peds[RV_peds$FamID %in% keep_fams, ]


study_pedigrees = RV_peds[1, ]
study_pedigrees = study_pedigrees[-1, ]


#NOTE: must export affect_onlyPed and assign_gen to run next chunk

for (i in unique(RV_peds$FamID)) {
 study_pedigrees = rbind(study_pedigrees,
                         affected_onlyPed(RV_peds[RV_peds$FamID == i, ]))
}

plot(study_pedigrees[study_pedigrees$FamID == 58, ])

study_pedigrees = study_pedigrees[order(study_pedigrees$FamID,study_pedigrees$Gen,study_pedigrees$ID), ]
head(study_pedigrees)
row.names(study_pedigrees) = NULL
head(study_pedigrees)

save(study_pedigrees, file="data/study_pedigrees.rdata", compress='xz')
