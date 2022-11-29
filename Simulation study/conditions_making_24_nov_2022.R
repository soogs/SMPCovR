# conditions for simulation (discovR) #
# initiated: 24-nov-2022 #

# 6 Sep 2022 simulation is the same as the simulation before (16 September 2022). 
# I only do "notall3" and "notall2"

# but the model selection is different. Refer to the word document at 9-nov-2022

# 24-nov-2022 is the same as the 9-nov-2022 simulation, just with 50 replicates instead of 20

# And I have only 20 replicates to make it faster

# 1. conditions ####
conditions <- list(dimension = c("low", "high"),
                   signal_level =  c(0.5, 0.9),
                   Jy = c(5, 20),
                   py_pattern = c("notall3", "notall2"),
                   reps = c(1:50))

condition_df <- data.frame(dimension = NA,
                           signal_level = NA,
                           Jy = NA,
                           py_pattern = NA,
                           reps = NA)

counts <- 0
for (dimensionz in 1:2){
  for (signal_levelz in 1:2){
    for (jyz in 1:2){
      for (patternz in 1:2){
        for (repsz in 1:50){
          
          counts <- counts + 1
          
          dimension_now <- conditions$dimension[dimensionz]
          signal_level_now <- conditions$signal_level[signal_levelz]
          Jy_now <- conditions$Jy[jyz]
          py_pattern_now <- conditions$py_pattern[patternz]
          reps_now <- conditions$reps[repsz]
          
          condition_df[counts,] <- c(dimension_now,
                                     signal_level_now,
                                     Jy_now,
                                     py_pattern_now,
                                     reps_now)
          
          print(counts) 
        }
      }
    }
  }
}

