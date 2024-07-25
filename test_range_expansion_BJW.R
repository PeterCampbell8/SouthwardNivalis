#Test for range expansion

setwd("FiveDecadesCounted")

library(readr)
library(tidyr)


# read in data
counts50 <- read_csv("All50sCounts/All50sCounts.csv")
counts60 <- read_csv("All60sCounts/All60sCounts.csv")
counts70 <- read_csv("All70sCounts/All70sCounts.csv")
counts80 <- read_csv("All80sCounts/All80sCounts.csv")
counts90 <- read_csv("All90sCounts/All90sCounts.csv")

# combine into one dataframe
identical(counts50$fid, counts60$fid)
identical(counts50$fid, counts70$fid)
identical(counts50$fid, counts80$fid)
identical(counts50$fid, counts90$fid)
counts <- cbind(counts80[,1:9], counts50[,7:10], counts60[,7:10], counts70[,7:10], counts80[,10:13], counts90[,7:10])

# average number of nivalis per background
sum(counts$NivalisCou)/sum(counts$Background)

# average number of nivalis per background, only hexes where >= nivalis
success.rate <- sum(counts[counts$NivalisCou > 0,]$NivalisCou)/sum(counts[counts$NivalisCou > 0,]$Background)

# number of trials (BackPre80) needed to get 1 nivalis
trials.needed <- 1 / success.rate






# number of trials needed to have X% confidence in getting a weasel if the true proportion is Y
sampleSize <- function(p, sig.level) {
  # p: proportion of occurrences in the population
  # sig.level: significance level
  require(stats)
  pval <- 1
  trials <- 0
  while (pval > sig.level) {
    trials <- trials + 1
    pval <- binom.test(0, trials, p, alternative = "less")$p.value
  }
  return(trials)
}

for (x in c(0.01, 0.05, 0.1, 0.2)) {
  for (y in c(0.005, 0.01, 0.015, 0.02)) {
    n <- sampleSize(p = x, sig.level = y)
    print(paste(x, y, n, sep = " "))
  }
}

sampleSize(p = success.rate, 0.05)






# number of hexes where BackPre80 is greater than trials.needed
nrow(counts[counts$BackPre80 >= trials.needed,])

# histograms of number of nivalis per hex
hist(counts$NivalisCou)
hist(counts[counts$NivalisCou > 0,]$NivalisCou)

# histograms of number of background per hex
hist(counts$Background)
hist(counts[counts$Background > 0,]$Background)


# binomial test for every hex where there is background sampling both pre-XXs and post-XXs
for (i in 1:nrow(counts)) {
  if (i == 1) {
    counts$p.value50 <- rep(NA, nrow(counts))
    counts$p.value60 <- rep(NA, nrow(counts))
    counts$p.value70 <- rep(NA, nrow(counts))
    counts$p.value80 <- rep(NA, nrow(counts))
    counts$p.value90 <- rep(NA, nrow(counts))
  }
  if (counts$BackPre50[i] != 0 & counts$BackPost50[i] !=0) {
    counts[i, "p.value50"] <- binom.test(counts$NivalPre50[i], counts$BackPre50[i], counts$NivalPost50[i]/counts$BackPost50[i], alternative = "less" )$p.value
  }
  if (counts$BackPre60[i] != 0 & counts$BackPost60[i] !=0) {
    counts[i, "p.value60"] <- binom.test(counts$NivalPre60[i], counts$BackPre60[i], counts$NivalPost60[i]/counts$BackPost60[i], alternative = "less" )$p.value
  }
  if (counts$BackPre70[i] != 0 & counts$BackPost70[i] !=0) {
    counts[i, "p.value70"] <- binom.test(counts$NivalPre70[i], counts$BackPre70[i], counts$NivalPost70[i]/counts$BackPost70[i], alternative = "less" )$p.value
  }
  if (counts$BackPre80[i] != 0 & counts$BackPost80[i] !=0) {
    counts[i, "p.value80"] <- binom.test(counts$NivalPre80[i], counts$BackPre80[i], counts$NivalPost80[i]/counts$BackPost80[i], alternative = "less" )$p.value
  }
  if (counts$BackPre90[i] != 0 & counts$BackPost90[i] !=0) {
    counts[i, "p.value90"] <- binom.test(counts$NivalPre90[i], counts$BackPre90[i], counts$NivalPost90[i]/counts$BackPost90[i], alternative = "less" )$p.value
  }
  
}


# add column for status of each hex
# Version using trials.needed as minimum sample size required, but only for pre-decade
for (i in 1:nrow(counts)) {
  if (i == 1) {
    counts$status50 <- rep(NA, nrow(counts))
    counts$status60 <- rep(NA, nrow(counts))
    counts$status70 <- rep(NA, nrow(counts))
    counts$status80 <- rep(NA, nrow(counts))
    counts$status90 <- rep(NA, nrow(counts))
    
  }
  if (counts$BackPre50[i] < trials.needed & counts$BackPost50[i] < trials.needed) {
    counts[i, "status50"] <- "not enough sampling pre50s or post50s"
  } else if (counts$BackPre50[i] < trials.needed) {
    counts[i, "status50"] <- "not enough sampling pre50s"
  } else if (counts$BackPost50[i] < trials.needed & counts$NivalPost50[i] < 1) {
    counts[i, "status50"] <- "not enough sampling post50s"
  } else if (replace_na(counts$p.value50[i], 1) <= 0.05) {
    counts[i, "status50"] <- "range expansion"
  } else if (counts$NivalPre50[i] == 0 & counts$NivalPost50[i] > 0) {
    counts[i, "status50"] <- "range expansion but non-significant"
  }
  if (counts$NivalPre50[i] > 0) {
    counts[i, "status50"] <- "exist pre50"
  } else if (counts$NivalPost50[i] < 1 & counts$BackPost50[i] > trials.needed) {
    counts[i, "status50"] <- "enough sampling post 50s but not documented"
  }
  if (counts$BackPre60[i] < trials.needed & counts$BackPost60[i] < trials.needed) {
    counts[i, "status60"] <- "not enough sampling pre60s or post60s"
  } else if (counts$BackPre60[i] < trials.needed) {
    counts[i, "status60"] <- "not enough sampling pre60s"
  } else if (counts$BackPost60[i] < trials.needed & counts$NivalPost60[i] < 1) {
    counts[i, "status60"] <- "not enough sampling post60s"
  } else if (replace_na(counts$p.value60[i], 1) <= 0.05) {
    counts[i, "status60"] <- "range expansion"
  } else if (counts$NivalPre60[i] == 0 & counts$NivalPost60[i] > 0) {
    counts[i, "status60"] <- "range expansion but non-significant"
  }
  if (counts$NivalPre60[i] > 0) {
    counts[i, "status60"] <- "exist pre60"
  } else if (counts$NivalPost60[i] < 1 & counts$BackPost60[i] > trials.needed) {
    counts[i, "status60"] <- "enough sampling post 60s but not documented"
  }
  if (counts$BackPre70[i] < trials.needed & counts$BackPost70[i] < trials.needed) {
    counts[i, "status70"] <- "not enough sampling pre70s or post70s"
  } else if (counts$BackPre70[i] < trials.needed) {
    counts[i, "status70"] <- "not enough sampling pre70s"
  } else if (counts$BackPost70[i] < trials.needed & counts$NivalPost70[i] < 1) {
    counts[i, "status70"] <- "not enough sampling post70s"
  } else if (replace_na(counts$p.value70[i], 1) <= 0.05) {
    counts[i, "status70"] <- "range expansion"
  } else if (counts$NivalPre70[i] == 0 & counts$NivalPost70[i] > 0) {
    counts[i, "status70"] <- "range expansion but non-significant"
  }
  if (counts$NivalPre70[i] > 0) {
    counts[i, "status70"] <- "exist pre70"
  } else if (counts$NivalPost70[i] < 1 & counts$BackPost70[i] > trials.needed) {
    counts[i, "status70"] <- "enough sampling post 70s but not documented"
  }
  if (counts$BackPre80[i] < trials.needed & counts$BackPost80[i] < trials.needed) {
    counts[i, "status80"] <- "not enough sampling pre80s or post80s"
  } else if (counts$BackPre80[i] < trials.needed) {
    counts[i, "status80"] <- "not enough sampling pre80s"
  } else if (counts$BackPost80[i] < trials.needed & counts$NivalPost80[i] < 1) {
    counts[i, "status80"] <- "not enough sampling post80s"
  } else if (replace_na(counts$p.value80[i], 1) <= 0.05) {
    counts[i, "status80"] <- "range expansion"
  } else if (counts$NivalPre80[i] == 0 & counts$NivalPost80[i] > 0) {
    counts[i, "status80"] <- "range expansion but non-significant"
  }
  if (counts$NivalPre80[i] > 0) {
    counts[i, "status80"] <- "exist pre80"
  } else if (counts$NivalPost80[i] < 1 & counts$BackPost80[i] > trials.needed) {
    counts[i, "status80"] <- "enough sampling post 80s but not documented"
  }
  if (counts$BackPre90[i] < trials.needed & counts$BackPost90[i] < trials.needed) {
    counts[i, "status90"] <- "not enough sampling pre90s or post90s"
  } else if (counts$BackPre90[i] < trials.needed) {
    counts[i, "status90"] <- "not enough sampling pre90s"
  } else if (counts$BackPost90[i] < trials.needed & counts$NivalPost90[i] < 1) {
    counts[i, "status90"] <- "not enough sampling post90s"
  } else if (replace_na(counts$p.value90[i], 1) <= 0.05) {
    counts[i, "status90"] <- "range expansion"
  } else if (counts$NivalPre90[i] == 0 & counts$NivalPost90[i] > 0) {
    counts[i, "status90"] <- "range expansion but non-significant"
  }
  if (counts$NivalPre90[i] > 0) {
    counts[i, "status90"] <- "exist pre90"
  } else if (counts$NivalPost90[i] < 1 & counts$BackPost90[i] > trials.needed) {
    counts[i, "status90"] <- "enough sampling post 90s but not documented"
  }
}




# Version using 50, 75, and 100 as minimum trials needed both pre- and post-decade
for (tn in c(50, 75, 100)) {
  for (i in 1:nrow(counts)) {
    if (i == 1) {
      a <- rep(NA, nrow(counts))
      b <- rep(NA, nrow(counts))
      c <- rep(NA, nrow(counts))
      d <- rep(NA, nrow(counts))
      e <- rep(NA, nrow(counts))
      counts <- cbind(counts, a, b, c, d, e)
    }
    if (counts$BackPre50[i] < trials.needed & counts$BackPost50[i] < trials.needed) {
      counts[i, "a"] <- "not enough sampling pre50s or post50s"
    } else if (counts$BackPre50[i] < trials.needed) {
      counts[i, "a"] <- "not enough sampling pre50s"
    } else if (counts$BackPost50[i] < trials.needed) {
      counts[i, "a"] <- "not enough sampling post50s"
    } else if (replace_na(counts$p.value50[i], 1) <= 0.05) {
      counts[i, "a"] <- "range expansion"
    } else if (counts$NivalPre50[i] == 0 & counts$NivalPost50[i] > 0) {
      counts[i, "a"] <- "range expansion but non-significant"
    }
    if (counts$NivalPre50[i] > 0) {
      counts[i, "a"] <- "exist pre50"
    } else if (counts$NivalPost50[i] < 1 & counts$BackPost50[i] > trials.needed) {
      counts[i, "a"] <- "enough sampling post 50s but not documented"
    }
    if (counts$BackPre60[i] < trials.needed & counts$BackPost60[i] < trials.needed) {
      counts[i, "b"] <- "not enough sampling pre60s or post60s"
    } else if (counts$BackPre60[i] < trials.needed) {
      counts[i, "b"] <- "not enough sampling pre60s"
    } else if (counts$BackPost60[i] < trials.needed) {
      counts[i, "b"] <- "not enough sampling post60s"
    } else if (replace_na(counts$p.value60[i], 1) <= 0.05) {
      counts[i, "b"] <- "range expansion"
    } else if (counts$NivalPre60[i] == 0 & counts$NivalPost60[i] > 0) {
      counts[i, "b"] <- "range expansion but non-significant"
    }
    if (counts$NivalPre60[i] > 0) {
      counts[i, "b"] <- "exist pre60"
    } else if (counts$NivalPost60[i] < 1 & counts$BackPost60[i] > trials.needed) {
      counts[i, "b"] <- "enough sampling post 60s but not documented"
    }
    if (counts$BackPre70[i] < trials.needed & counts$BackPost70[i] < trials.needed) {
      counts[i, "c"] <- "not enough sampling pre70s or post70s"
    } else if (counts$BackPre70[i] < trials.needed) {
      counts[i, "c"] <- "not enough sampling pre70s"
    } else if (counts$BackPost70[i] < trials.needed) {
      counts[i, "c"] <- "not enough sampling post70s"
    } else if (replace_na(counts$p.value70[i], 1) <= 0.05) {
      counts[i, "c"] <- "range expansion"
    } else if (counts$NivalPre70[i] == 0 & counts$NivalPost70[i] > 0) {
      counts[i, "c"] <- "range expansion but non-significant"
    }
    if (counts$NivalPre70[i] > 0) {
      counts[i, "c"] <- "exist pre70"
    } else if (counts$NivalPost70[i] < 1 & counts$BackPost70[i] > trials.needed) {
      counts[i, "c"] <- "enough sampling post 70s but not documented"
    }
    if (counts$BackPre80[i] < trials.needed & counts$BackPost80[i] < trials.needed) {
      counts[i, "d"] <- "not enough sampling pre80s or post80s"
    } else if (counts$BackPre80[i] < trials.needed) {
      counts[i, "d"] <- "not enough sampling pre80s"
    } else if (counts$BackPost80[i] < trials.needed) {
      counts[i, "d"] <- "not enough sampling post80s"
    } else if (replace_na(counts$p.value80[i], 1) <= 0.05) {
      counts[i, "d"] <- "range expansion"
    } else if (counts$NivalPre80[i] == 0 & counts$NivalPost80[i] > 0) {
      counts[i, "d"] <- "range expansion but non-significant"
    }
    if (counts$NivalPre80[i] > 0) {
      counts[i, "d"] <- "exist pre80"
    } else if (counts$NivalPost80[i] < 1 & counts$BackPost80[i] > trials.needed) {
      counts[i, "d"] <- "enough sampling post 80s but not documented"
    }
    if (counts$BackPre90[i] < trials.needed & counts$BackPost90[i] < trials.needed) {
      counts[i, "e"] <- "not enough sampling pre90s or post90s"
    } else if (counts$BackPre90[i] < trials.needed) {
      counts[i, "e"] <- "not enough sampling pre90s"
    } else if (counts$BackPost90[i] < trials.needed) {
      counts[i, "e"] <- "not enough sampling post90s"
    } else if (replace_na(counts$p.value90[i], 1) <= 0.05) {
      counts[i, "e"] <- "range expansion"
    } else if (counts$NivalPre90[i] == 0 & counts$NivalPost90[i] > 0) {
      counts[i, "e"] <- "range expansion but non-significant"
    }
    if (counts$NivalPre90[i] > 0) {
      counts[i, "e"] <- "exist pre90"
    } else if (counts$NivalPost90[i] < 1 & counts$BackPost90[i] > trials.needed) {
      counts[i, "e"] <- "enough sampling post 90s but not documented"
    }
  }
  names(counts)[names(counts) %in% c("a","b","c","d","e")] <- c(paste0("status50min", tn),
                                                                paste0("status60min", tn),
                                                                paste0("status70min", tn),
                                                                paste0("status80min", tn),
                                                                paste0("status90min", tn))
}





# Version using trials.needed as minimum sample size required, and if background samples is lower than 
# this, than use success.rate instead of counts$NivalPostXX[i]/counts$BackPostXX[i]


# binomial test for every hex where there is background sampling both pre-XXs and post-XXs
for (i in 1:nrow(counts)) {
  if (i == 1) {
    counts$p.value50.min76 <- rep(NA, nrow(counts))
    counts$p.value60.min76 <- rep(NA, nrow(counts))
    counts$p.value70.min76 <- rep(NA, nrow(counts))
    counts$p.value80.min76 <- rep(NA, nrow(counts))
    counts$p.value90.min76 <- rep(NA, nrow(counts))
  }
  if (counts$BackPre50[i] != 0 & counts$BackPost50[i] !=0) {
    if (counts$BackPost50[i] > trials.needed) {
      counts[i, "p.value50.min76"] <- binom.test(counts$NivalPre50[i], counts$BackPre50[i], counts$NivalPost50[i]/counts$BackPost50[i], alternative = "less" )$p.value
    } else {
      counts[i, "p.value50.min76"] <- binom.test(counts$NivalPre50[i], counts$BackPre50[i], success.rate, alternative = "less" )$p.value
    }
  }
  if (counts$BackPre60[i] != 0 & counts$BackPost60[i] !=0) {
    if (counts$BackPost60[i] > trials.needed) {
      counts[i, "p.value60.min76"] <- binom.test(counts$NivalPre60[i], counts$BackPre60[i], counts$NivalPost60[i]/counts$BackPost60[i], alternative = "less" )$p.value
    } else {
      counts[i, "p.value60.min76"] <- binom.test(counts$NivalPre60[i], counts$BackPre60[i], success.rate, alternative = "less" )$p.value
    }
  }
  if (counts$BackPre70[i] != 0 & counts$BackPost70[i] !=0) {
    if (counts$BackPost70[i] > trials.needed) {
      counts[i, "p.value70.min76"] <- binom.test(counts$NivalPre70[i], counts$BackPre70[i], counts$NivalPost70[i]/counts$BackPost70[i], alternative = "less" )$p.value
    } else {
      counts[i, "p.value70.min76"] <- binom.test(counts$NivalPre70[i], counts$BackPre70[i], success.rate, alternative = "less" )$p.value
    }
  }
  if (counts$BackPre80[i] != 0 & counts$BackPost80[i] !=0) {
    if (counts$BackPost80[i] > trials.needed) {
      counts[i, "p.value80.min76"] <- binom.test(counts$NivalPre80[i], counts$BackPre80[i], counts$NivalPost80[i]/counts$BackPost80[i], alternative = "less" )$p.value
    } else {
      counts[i, "p.value80.min76"] <- binom.test(counts$NivalPre80[i], counts$BackPre80[i], success.rate, alternative = "less" )$p.value
    }
  }
  if (counts$BackPre90[i] != 0 & counts$BackPost90[i] !=0) {
    if (counts$BackPost90[i] > trials.needed) {
      counts[i, "p.value90.min76"] <- binom.test(counts$NivalPre90[i], counts$BackPre90[i], counts$NivalPost90[i]/counts$BackPost90[i], alternative = "less" )$p.value
    } else {
      counts[i, "p.value90.min76"] <- binom.test(counts$NivalPre90[i], counts$BackPre90[i], success.rate, alternative = "less" )$p.value
    }
  }
}

# Assign status
for (i in 1:nrow(counts)) {
  if (i == 1) {
    counts$status50.min76 <- rep(NA, nrow(counts))
    counts$status60.min76 <- rep(NA, nrow(counts))
    counts$status70.min76 <- rep(NA, nrow(counts))
    counts$status80.min76 <- rep(NA, nrow(counts))
    counts$status90.min76 <- rep(NA, nrow(counts))
    
  }
  if (counts$BackPre50[i] < trials.needed & counts$BackPost50[i] < trials.needed) {
    counts[i, "status50.min76"] <- "not enough sampling pre50s or post50s"
  } else if (counts$BackPre50[i] < trials.needed) {
    counts[i, "status50.min76"] <- "not enough sampling pre50s"
  } else if (counts$BackPost50[i] < trials.needed & counts$NivalPost50[i] < 1) {
    counts[i, "status50.min76"] <- "not enough sampling post50s"
  } else if (replace_na(counts$p.value50.min76[i], 1) <= 0.05) {
    counts[i, "status50.min76"] <- "range expansion"
  } else if (counts$NivalPre50[i] == 0 & counts$NivalPost50[i] > 0) {
    counts[i, "status50.min76"] <- "range expansion but non-significant"
  }
  if (counts$NivalPre50[i] > 0) {
    counts[i, "status50.min76"] <- "exist pre50"
  } else if (counts$NivalPost50[i] < 1 & counts$BackPost50[i] > trials.needed) {
    counts[i, "status50.min76"] <- "enough sampling post 50s but not documented"
  }
  if (counts$BackPre60[i] < trials.needed & counts$BackPost60[i] < trials.needed) {
    counts[i, "status60.min76"] <- "not enough sampling pre60s or post60s"
  } else if (counts$BackPre60[i] < trials.needed) {
    counts[i, "status60.min76"] <- "not enough sampling pre60s"
  } else if (counts$BackPost60[i] < trials.needed & counts$NivalPost60[i] < 1) {
    counts[i, "status60.min76"] <- "not enough sampling post60s"
  } else if (replace_na(counts$p.value60.min76[i], 1) <= 0.05) {
    counts[i, "status60.min76"] <- "range expansion"
  } else if (counts$NivalPre60[i] == 0 & counts$NivalPost60[i] > 0) {
    counts[i, "status60.min76"] <- "range expansion but non-significant"
  }
  if (counts$NivalPre60[i] > 0) {
    counts[i, "status60.min76"] <- "exist pre60"
  } else if (counts$NivalPost60[i] < 1 & counts$BackPost60[i] > trials.needed) {
    counts[i, "status60.min76"] <- "enough sampling post 60s but not documented"
  }
  if (counts$BackPre70[i] < trials.needed & counts$BackPost70[i] < trials.needed) {
    counts[i, "status70.min76"] <- "not enough sampling pre70s or post70s"
  } else if (counts$BackPre70[i] < trials.needed) {
    counts[i, "status70.min76"] <- "not enough sampling pre70s"
  } else if (counts$BackPost70[i] < trials.needed & counts$NivalPost70[i] < 1) {
    counts[i, "status70.min76"] <- "not enough sampling post70s"
  } else if (replace_na(counts$p.value70.min76[i], 1) <= 0.05) {
    counts[i, "status70.min76"] <- "range expansion"
  } else if (counts$NivalPre70[i] == 0 & counts$NivalPost70[i] > 0) {
    counts[i, "status70.min76"] <- "range expansion but non-significant"
  }
  if (counts$NivalPre70[i] > 0) {
    counts[i, "status70.min76"] <- "exist pre70"
  } else if (counts$NivalPost70[i] < 1 & counts$BackPost70[i] > trials.needed) {
    counts[i, "status70.min76"] <- "enough sampling post 70s but not documented"
  }
  if (counts$BackPre80[i] < trials.needed & counts$BackPost80[i] < trials.needed) {
    counts[i, "status80.min76"] <- "not enough sampling pre80s or post80s"
  } else if (counts$BackPre80[i] < trials.needed) {
    counts[i, "status80.min76"] <- "not enough sampling pre80s"
  } else if (counts$BackPost80[i] < trials.needed & counts$NivalPost80[i] < 1) {
    counts[i, "status80.min76"] <- "not enough sampling post80s"
  } else if (replace_na(counts$p.value80.min76[i], 1) <= 0.05) {
    counts[i, "status80.min76"] <- "range expansion"
  } else if (counts$NivalPre80[i] == 0 & counts$NivalPost80[i] > 0) {
    counts[i, "status80.min76"] <- "range expansion but non-significant"
  }
  if (counts$NivalPre80[i] > 0) {
    counts[i, "status80.min76"] <- "exist pre80"
  } else if (counts$NivalPost80[i] < 1 & counts$BackPost80[i] > trials.needed) {
    counts[i, "status80.min76"] <- "enough sampling post 80s but not documented"
  }
  if (counts$BackPre90[i] < trials.needed & counts$BackPost90[i] < trials.needed) {
    counts[i, "status90.min76"] <- "not enough sampling pre90s or post90s"
  } else if (counts$BackPre90[i] < trials.needed) {
    counts[i, "status90.min76"] <- "not enough sampling pre90s"
  } else if (counts$BackPost90[i] < trials.needed & counts$NivalPost90[i] < 1) {
    counts[i, "status90.min76"] <- "not enough sampling post90s"
  } else if (replace_na(counts$p.value90.min76[i], 1) <= 0.05) {
    counts[i, "status90.min76"] <- "range expansion"
  } else if (counts$NivalPre90[i] == 0 & counts$NivalPost90[i] > 0) {
    counts[i, "status90.min76"] <- "range expansion but non-significant"
  }
  if (counts$NivalPre90[i] > 0) {
    counts[i, "status90.min76"] <- "exist pre90"
  } else if (counts$NivalPost90[i] < 1 & counts$BackPost90[i] > trials.needed) {
    counts[i, "status90.min76"] <- "enough sampling post 90s but not documented"
  }
}


# number of hexes that fall into each category:
table(counts$status50)
table(counts$status60)
table(counts$status70)
table(counts$status80)
table(counts$status90)

colnames(counts[,c(35:54, 60:64)])
any(is.na(counts[,c(35:54, 60:64)]))


# write out csv file
write.csv(counts, file = "AllCountsCategories.csv", quote = F, row.names = F)

#--------------------------------------------------------------------
#Test for clumping and south-ness--------------------------------
setwd("FiveDecadesCounted")

library(readr)
library(geosphere)
library(ggplot2)

set.seed(234890572)

# Use only hexes whose centroids are within the Great Plains
# (excludes some significant hexes just on the border)
#GP_centroids <- read_csv("All60sCountsClippedReprojectedCentroidCoordinatesGreatPlainsIntersection.csv")

# Use hexes who overlap at all with Great Plains
# (more inclusive of hexes right on the edge of the Great Plains)
GP_centroids <- read_csv("All60sCountsClippedReprojectedHexGreatPlainsIntersection.csv")


hist(GP_centroids$Latitude, breaks = 20)

# Average latitude of centroids with significant range expansion
sig_lat <- GP_centroids[GP_centroids$AllCountsCategories_status60.min76 == "range expansion",]$Latitude
num_sig <- length(sig_lat)
ave_sig_lat <- mean(sig_lat)


######################################
# Bootstrap sampling to test whether
# average latitude of significant 
# cells is farther south than expected 
# by chance
######################################

# Make function to do bootstrapping
weasel_boot <- function (x, length, replicates) {
  ave_boot_lats <- c()
  for (i in 1:replicates) {
    boot_samp <- sample(x$Latitude, length, replace = F)
    ave_boot_samp <- mean(boot_samp)
    ave_boot_lats <- c(ave_boot_lats, ave_boot_samp)
  }
  ave_boot_lats
}

# Do the bootstrapping and look at it
ave_boot_lats <- weasel_boot(GP_centroids, num_sig, 100000)
hist(ave_boot_lats, main = "", xlab = "Average latitude (degrees) of random resamples")
abline(v=ave_sig_lat, col = "black", lty = 2)
abline(v=quantile(ave_boot_lats, probs = 0.05), col = "blue")


# Look at quantiles for various significance levels
quantile(ave_boot_lats, probs = c(0, 0.01, 0.025, 0.05))
# Get p-value for observed average latitude of significant range expansion cells
percentile <- ecdf(ave_boot_lats)
percentile(ave_sig_lat)





######################################
# Bootstrap sampling to test whether
# average distance between significant
# centroids is smaller than expected
# by chance
######################################



# Average pairwise distance of centroids with significant range expansion

# Dataframe of only significant range expansion centroids
GP_centroids_sig <- GP_centroids[GP_centroids$AllCountsCategories_status60.min76 == "range expansion",]


# Make function that calculates average pairwise distance between any set of coordinates
mean_pairwise_dist <- function(x, bootstrap = F, length = NULL, replicates = 1) {
  ave_dist <- c()
  for (r in 1:replicates) {
    if (r%%1000==0) {
      print(r)
    }
    dists <- c()
    y <- x
    if (bootstrap) {
      y <- x[sample(nrow(x), size = length, replace = F),]
    }
    
    for (i in 1:nrow(y)) {
      j <- i+1
      while (j <= nrow(y)) {
        # centroid points are WGS84, distGeo calculates shortest distance between two points on a WGS84 ellipsoid
        # outputs distance in meters
        dist <- distm(c(y$Longitude[i], y$Latitude[i]), 
                      c(y$Longitude[j], y$Latitude[j]),
                      fun = distGeo)   
        dists <- c(dists, dist)
        j <- j+1
      }
    }
    ave_dist <- c(ave_dist, mean(dists))
    #print(dists)
    #hist(dists)
  }
  return(ave_dist)
}

# To run the function on the real data, plug in GP_centroids_sig, no bootstraps, replicates = 1
ave_sig_dist <- mean_pairwise_dist(GP_centroids_sig, bootstrap = F, replicates = 1)

# To run the funciton on bootstraps, plug in GP_centroids, num_sig, bootstraps = T, replicates)
ave_boot_dists <- mean_pairwise_dist(GP_centroids, bootstrap = T, length = num_sig, replicates = 100000)



#Look at bootstrap distribution
hist(ave_boot_dists, main = "", xlab = "Average pairwise distance (meters) of random resamples")
abline(v=ave_sig_dist, col = "black", lty = 2)
abline(v=quantile(ave_boot_dists, probs = 0.05), col = "blue")


# Look at quantiles for various significance levels
quantile(ave_boot_dists, probs = c(0, 0.01, 0.025, 0.05))
# Get p-value for observed average latitude of significant range expansion cells
percentile_dists <- ecdf(ave_boot_dists)
percentile_dists(ave_sig_dist)

